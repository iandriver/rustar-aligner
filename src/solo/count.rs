//! UMI deduplication and raw count-matrix output (Phase 14.4).
//!
//! Collates the per-read `(cell, UMI, gene)` records produced during alignment
//! into a sparse per-cell, per-gene count matrix:
//!   1. resolve deferred 1MM_multi cell barcodes via the count+quality posterior
//!      (STAR `SoloReadFeature_inputRecords.cpp`: weight = exactCount·10^(−q/10));
//!   2. group reads by `(cell, gene)` and collapse UMIs per `--soloUMIdedup`
//!      (STAR `SoloFeature_collapseUMIall.cpp`);
//!   3. write `Solo.out/Gene/raw/{matrix.mtx, barcodes.tsv, features.tsv}` in
//!      CellRanger-compatible MatrixMarket layout (features × barcodes, 1-based).

use crate::error::Error;
use crate::solo::whitelist::CbWhitelist;
use crate::solo::{SoloContext, SoloCountRecord};
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;
use std::str::FromStr;

// ---------------------------------------------------------------------------
// UMI deduplication
// ---------------------------------------------------------------------------

/// `--soloUMIdedup` method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiDedup {
    /// Count distinct UMI sequences (no error correction).
    Exact,
    /// No collapsing — count every read.
    NoDedup,
    /// Collapse all UMIs within Hamming-1 transitively (connected components).
    OneMmAll,
    /// UMI-tools directional, `count_hub >= 2*count_leaf + 0`.
    OneMmDirectional,
    /// UMI-tools directional original, `count_hub >= 2*count_leaf - 1`.
    OneMmDirectionalUmiTools,
    /// CellRanger 2–4 1MM collapse: each UMI is corrected to a higher-count
    /// 1MM neighbor (non-transitive); count = distinct corrected UMIs.
    OneMmCr,
}

impl FromStr for UmiDedup {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Exact" => Ok(Self::Exact),
            "NoDedup" => Ok(Self::NoDedup),
            "1MM_All" => Ok(Self::OneMmAll),
            "1MM_Directional" => Ok(Self::OneMmDirectional),
            "1MM_Directional_UMItools" => Ok(Self::OneMmDirectionalUmiTools),
            "1MM_CR" => Ok(Self::OneMmCr),
            _ => Err(format!(
                "unknown soloUMIdedup '{s}'; expected Exact, NoDedup, 1MM_All, 1MM_Directional, 1MM_Directional_UMItools, or 1MM_CR"
            )),
        }
    }
}

/// `--soloUMIfiltering`: removal of UMIs that map to multiple genes within a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiFiltering {
    /// No multi-gene UMI filtering.
    None,
    /// Remove lower-count gene assignments of a multi-gene UMI; if every gene
    /// has a single read, drop the UMI entirely (STAR `MultiGeneUMI`).
    MultiGeneUmi,
    /// CellRanger > 3.0 variant: keep only the highest-read-count gene for a
    /// multi-gene UMI (ties retained), without the all-singletons drop.
    MultiGeneUmiCr,
}

impl FromStr for UmiFiltering {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "-" | "None" => Ok(Self::None),
            // MultiGeneUMI_All behaves like MultiGeneUMI for the count matrix.
            "MultiGeneUMI" | "MultiGeneUMI_All" => Ok(Self::MultiGeneUmi),
            "MultiGeneUMI_CR" => Ok(Self::MultiGeneUmiCr),
            _ => Err(format!(
                "unknown soloUMIfiltering '{s}'; expected -, None, MultiGeneUMI, MultiGeneUMI_CR, or MultiGeneUMI_All"
            )),
        }
    }
}

/// True if packed UMIs `a` and `b` (length `len`) differ at exactly one base.
fn hamming1(a: u64, b: u64, len: usize) -> bool {
    let x = a ^ b;
    let mut diff = 0u32;
    for i in 0..len {
        if (x >> (2 * i)) & 0b11 != 0 {
            diff += 1;
            if diff > 1 {
                return false;
            }
        }
    }
    diff == 1
}

/// Deduplicate the UMIs observed for one `(cell, gene)` pair into a molecule
/// count. `umis` maps each packed UMI to its read multiplicity.
#[allow(clippy::implicit_hasher)] // always called with the default hasher
pub fn dedup_count(umis: &HashMap<u64, u32>, method: UmiDedup, umi_len: usize) -> u64 {
    match method {
        UmiDedup::Exact => umis.len() as u64,
        UmiDedup::NoDedup => umis.values().map(|&c| u64::from(c)).sum(),
        UmiDedup::OneMmAll => connected_components(umis, umi_len),
        UmiDedup::OneMmDirectional => directional(umis, umi_len, 0),
        UmiDedup::OneMmDirectionalUmiTools => directional(umis, umi_len, -1),
        UmiDedup::OneMmCr => cellranger_1mm(umis, umi_len),
    }
}

/// 1MM_CR: CellRanger's 1-mismatch UMI collapse (STAR `umiArrayCorrect_CR`).
/// UMIs are sorted ascending by `(count, umi)`; each UMI is corrected to the
/// LAST (highest-count) 1MM neighbor with a strictly later sort position — i.e.
/// its highest-count 1MM neighbor. Correction is non-transitive (it points to
/// the neighbor's raw UMI, not its corrected value); the molecule count is the
/// number of distinct corrected UMIs.
fn cellranger_1mm(umis: &HashMap<u64, u32>, umi_len: usize) -> u64 {
    let mut items: Vec<(u64, u32)> = umis.iter().map(|(&u, &c)| (u, c)).collect();
    // Ascending by count, then by UMI value (mirrors funCompareSolo1 ordering,
    // so the inner scan from the end meets higher-count neighbors first).
    items.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut corrected: Vec<u64> = Vec::with_capacity(n);
    for iu in 0..n {
        let mut corr = items[iu].0;
        let mut iuu = n;
        while iuu > iu + 1 {
            iuu -= 1;
            if hamming1(items[iu].0, items[iuu].0, umi_len) {
                corr = items[iuu].0;
                break;
            }
        }
        corrected.push(corr);
    }
    let distinct: std::collections::HashSet<u64> = corrected.into_iter().collect();
    distinct.len() as u64
}

/// 1MM_All: number of connected components when UMIs within Hamming-1 are
/// merged transitively (union-find).
fn connected_components(umis: &HashMap<u64, u32>, umi_len: usize) -> u64 {
    let keys: Vec<u64> = umis.keys().copied().collect();
    let n = keys.len();
    if n <= 1 {
        return n as u64;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if hamming1(keys[i], keys[j], umi_len) {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut roots = std::collections::HashSet::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        roots.insert(r);
    }
    roots.len() as u64
}

/// 1MM_Directional: a lower-count UMI within Hamming-1 of a hub whose count
/// satisfies `count_hub >= 2*count_leaf + dir_count_add` is absorbed; the
/// molecule count is the number of surviving (non-absorbed) UMIs.
fn directional(umis: &HashMap<u64, u32>, umi_len: usize, dir_count_add: i64) -> u64 {
    // Sort by count desc, then by UMI value for determinism.
    let mut items: Vec<(u64, u32)> = umis.iter().map(|(&u, &c)| (u, c)).collect();
    items.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut absorbed = vec![false; n];
    for i in 0..n {
        if absorbed[i] {
            continue;
        }
        let hub_count = i64::from(items[i].1);
        for j in 0..n {
            if i == j || absorbed[j] {
                continue;
            }
            let leaf_count = i64::from(items[j].1);
            if leaf_count <= hub_count
                && hub_count >= 2 * leaf_count + dir_count_add
                && hamming1(items[i].0, items[j].0, umi_len)
            {
                absorbed[j] = true;
            }
        }
    }
    (n - absorbed.iter().filter(|&&a| a).count()) as u64
}

// ---------------------------------------------------------------------------
// Cell-barcode multi-match resolution (deferred 1MM_multi)
// ---------------------------------------------------------------------------

/// Resolve a 1MM_multi cell barcode to a single whitelist index using the
/// count+quality posterior: weight = `(exactCount[cand] + pseudocount) · 10^(−q/10)`
/// where `q` is the mismatch-position Phred score. `pseudocount` is 1 for the
/// `*_pseudocounts` match types (CellRanger ≥ 3.0). Returns the argmax, or
/// `None` if no candidate has positive weight.
fn resolve_multi_cb(
    candidates: &[crate::solo::whitelist::CbCandidate],
    exact_counts: &[u64],
    pseudocount: f64,
) -> Option<u32> {
    let mut best: Option<(u32, f64)> = None;
    let mut total = 0.0f64;
    for c in candidates {
        let prior = *exact_counts.get(c.wl_index as usize).unwrap_or(&0) as f64 + pseudocount;
        let q = f64::from(c.mismatch_qual.saturating_sub(33)); // Phred+33 → Phred
        let weight = prior * 10f64.powf(-q / 10.0);
        total += weight;
        match best {
            Some((_, w)) if w >= weight => {}
            _ => best = Some((c.wl_index, weight)),
        }
    }
    match best {
        Some((idx, w)) if total > 0.0 && w > 0.0 => Some(idx),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Matrix assembly + output
// ---------------------------------------------------------------------------

/// Build and stream the raw count matrix to `matrix_path` in one per-cell pass,
/// returning the number of non-zero entries written.
///
/// Mirrors STAR's `SoloFeature_collapseUMIall.cpp`: the flat record list is
/// sorted by cell barcode so each cell's reads are contiguous, then **one cell
/// is processed at a time** (Step 1 — peak build memory is a single cell's
/// `umi → gene` maps, not a global `cell → umi → gene` nest over all records).
///
/// Step 2 (streaming output): each cell's `gene → count` entries are written
/// straight to a temporary MatrixMarket body as they are produced — the global
/// `cell → (gene → count)` map is never materialized. `nnz` is counted on the
/// fly; the final `matrix.mtx` is the header (`rows cols nnz`) followed by the
/// temp body (the BySJout temp-file pattern). So matrix-output memory is bounded
/// by one cell regardless of how many cells the raw whitelist matrix spans.
///
/// Records are sorted by cb (ascending column), and each cell's genes are
/// emitted ascending, so entries come out in the same order as before.
#[allow(clippy::too_many_arguments)]
/// Per-cell summary collected while streaming the matrix: reads (records before
/// UMI dedup), UMIs (deduped column sum), and genes detected (nonzero entries).
#[derive(Clone, Copy)]
pub struct CellStat {
    pub n_reads: u64,
    pub n_umis: u64,
    pub n_genes: u32,
}

/// What `stream_matrix` returns alongside the written matrix.
pub struct MatrixStats {
    pub nnz: usize,
    /// One entry per barcode that received ≥1 UMI (the raw, unfiltered set).
    pub cells: Vec<CellStat>,
    /// Distinct genes with a nonzero count anywhere in the raw matrix.
    pub genes_detected: u32,
}

#[allow(clippy::too_many_arguments)]
fn stream_matrix(
    ctx: &SoloContext,
    recorder: &crate::solo::SoloRecorder,
    method: UmiDedup,
    filtering: UmiFiltering,
    umi_len: usize,
    pseudocount: f64,
    matrix_path: &Path,
    n_features: usize,
    n_barcodes: usize,
) -> Result<MatrixStats, Error> {
    let dir = matrix_path.parent().unwrap_or_else(|| Path::new("."));
    let mut body_tmp = tempfile::Builder::new()
        .prefix(".matrix_body")
        .tempfile_in(dir)
        .map_err(|e| Error::io(e, dir))?;
    let mut nnz = 0usize;
    let mut cell_stats: Vec<CellStat> = Vec::new();
    let mut gene_seen = vec![false; n_features];

    {
        let mut body = std::io::BufWriter::new(body_tmp.as_file_mut());

        // Move records out of the recorder; fold in resolved 1MM_multi cells.
        let mut records = std::mem::take(&mut *recorder.records.lock().unwrap());
        let exact_counts = ctx.whitelist.exact_count_snapshot();
        let multi = std::mem::take(&mut *recorder.multi_records.lock().unwrap());
        for m in &multi {
            if let Some(cb) = resolve_multi_cb(&m.candidates, &exact_counts, pseudocount) {
                records.push(SoloCountRecord {
                    cb,
                    umi: m.umi,
                    gene: m.gene,
                });
            }
        }
        drop(multi);

        // Group each cell's reads together so we can process + free one at a time.
        records.sort_unstable_by_key(|r| r.cb);

        let mut i = 0;
        while i < records.len() {
            let cb = records[i].cb;

            // umi → gene → read multiplicity, for this cell only.
            let mut umi_genes: HashMap<u64, HashMap<u32, u32>> = HashMap::new();
            let mut j = i;
            while j < records.len() && records[j].cb == cb {
                let r = &records[j];
                *umi_genes
                    .entry(r.umi)
                    .or_default()
                    .entry(r.gene)
                    .or_insert(0) += 1;
                j += 1;
            }

            // (gene → (umi → read_count)) after multi-gene UMI filtering.
            let mut gene_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::new();
            for (&umi, genes) in &umi_genes {
                for (&gene, &rc) in filter_multi_gene_umi(genes, filtering) {
                    *gene_umis.entry(gene).or_default().entry(umi).or_insert(0) += rc;
                }
            }

            // Collapse UMIs per gene, then emit this cell's entries gene-ascending.
            let mut cell_entries: Vec<(u32, u64)> = Vec::with_capacity(gene_umis.len());
            for (&gene, umis) in &gene_umis {
                let count = dedup_count(umis, method, umi_len);
                if count > 0 {
                    cell_entries.push((gene, count));
                }
            }
            cell_entries.sort_unstable_by_key(|&(g, _)| g);
            // Per-cell summary: reads = records (j-i), genes = nonzero entries,
            // UMIs = sum of deduped counts.
            let n_reads = (j - i) as u64;
            let n_genes = cell_entries.len() as u32;
            let mut n_umis = 0u64;
            for (g, c) in cell_entries {
                n_umis += c;
                gene_seen[g as usize] = true;
                writeln!(body, "{} {} {}", g + 1, cb + 1, c)
                    .map_err(|e| Error::io(e, matrix_path))?;
                nnz += 1;
            }
            if n_umis > 0 {
                cell_stats.push(CellStat {
                    n_reads,
                    n_umis,
                    n_genes,
                });
            }

            i = j;
        }
        body.flush().map_err(|e| Error::io(e, matrix_path))?;
    }

    // Final matrix.mtx = MatrixMarket header (now that nnz is known) + temp body.
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(matrix_path).map_err(|e| Error::io(e, matrix_path))?,
    );
    writeln!(out, "%%MatrixMarket matrix coordinate integer general")
        .map_err(|e| Error::io(e, matrix_path))?;
    writeln!(out, "%").map_err(|e| Error::io(e, matrix_path))?;
    writeln!(out, "{n_features} {n_barcodes} {nnz}").map_err(|e| Error::io(e, matrix_path))?;
    let mut body_read = body_tmp.reopen().map_err(|e| Error::io(e, matrix_path))?;
    std::io::copy(&mut body_read, &mut out).map_err(|e| Error::io(e, matrix_path))?;
    out.flush().map_err(|e| Error::io(e, matrix_path))?;
    let genes_detected = gene_seen.iter().filter(|&&s| s).count() as u32;
    Ok(MatrixStats {
        nnz,
        cells: cell_stats,
        genes_detected,
    })
}

/// Apply `--soloUMIfiltering` to the gene→read_count map of a single UMI,
/// returning the surviving (gene, read_count) entries.
fn filter_multi_gene_umi(genes: &HashMap<u32, u32>, filtering: UmiFiltering) -> Vec<(&u32, &u32)> {
    if filtering == UmiFiltering::None || genes.len() <= 1 {
        return genes.iter().collect();
    }
    let max = genes.values().copied().max().unwrap_or(0);
    match filtering {
        // STAR MultiGeneUMI: threshold = max (or 2 if max==1, dropping all
        // single-read multi-gene UMIs); keep genes with read_count >= threshold.
        UmiFiltering::MultiGeneUmi => {
            let thresh = if max == 1 { 2 } else { max };
            genes.iter().filter(|&(_, &rc)| rc >= thresh).collect()
        }
        // CellRanger > 3.0: keep the highest-read-count gene(s); no singleton drop.
        UmiFiltering::MultiGeneUmiCr => genes.iter().filter(|&(_, &rc)| rc >= max).collect(),
        UmiFiltering::None => unreachable!(),
    }
}

/// CellRanger-2.2 knee threshold on per-barcode UMI totals (STARsolo's default
/// `--soloCellFilter CellRanger2.2 3000 0.99 10`). Returns the minimum UMI count
/// for a barcode to be called a cell.
fn knee_cr22(umis_desc: &[u64], n_expected: usize, max_pct: f64, max_min_ratio: f64) -> u64 {
    if umis_desc.is_empty() {
        return 0;
    }
    let idx = ((n_expected as f64 * (1.0 - max_pct)).round() as usize).min(umis_desc.len() - 1);
    let robust_max = umis_desc[idx] as f64;
    (robust_max / max_min_ratio).ceil() as u64
}

/// Median of an ascending-sorted slice (0 if empty).
fn median_sorted(sorted: &[u64]) -> u64 {
    let n = sorted.len();
    if n == 0 {
        0
    } else if n % 2 == 1 {
        sorted[n / 2]
    } else {
        u64::midpoint(sorted[n / 2 - 1], sorted[n / 2])
    }
}

/// Write the raw gene-count matrix + `Summary.csv` for a finished solo run.
/// No-op (with a warning) when there is no explicit whitelist.
pub fn write_gene_matrix(
    ctx: &SoloContext,
    params: &crate::params::Parameters,
    align_stats: &crate::stats::AlignmentStats,
) -> Result<(), Error> {
    let CbWhitelist::List { sorted, .. } = &ctx.whitelist else {
        log::warn!(
            "STARsolo: --soloCBwhitelist None matrix output is not yet supported (Phase 14.4); skipping matrix"
        );
        return Ok(());
    };

    let method: UmiDedup = params
        .solo_umi_dedup
        .first()
        .map_or("1MM_All", String::as_str)
        .parse()
        .unwrap_or(UmiDedup::OneMmAll);
    let filtering: UmiFiltering = params
        .solo_umi_filtering
        .first()
        .map_or("-", String::as_str)
        .parse()
        .unwrap_or(UmiFiltering::None);
    // `*_pseudocounts` CB-match types add 1 to the posterior prior.
    let pseudocount = if params.solo_cb_match_wl_type.contains("pseudocounts") {
        1.0
    } else {
        0.0
    };
    let umi_len = params.solo_umi_len as usize;

    let solo_dir = params
        .solo_out_file_names
        .first()
        .cloned()
        .unwrap_or_else(|| "Solo.out/".to_string());
    let features_name = params
        .solo_out_file_names
        .get(1)
        .cloned()
        .unwrap_or_else(|| "features.tsv".to_string());
    let barcodes_name = params
        .solo_out_file_names
        .get(2)
        .cloned()
        .unwrap_or_else(|| "barcodes.tsv".to_string());
    let matrix_name = params
        .solo_out_file_names
        .get(3)
        .cloned()
        .unwrap_or_else(|| "matrix.mtx".to_string());

    // Global mapping funnel (shared across features). The region tallies are
    // CellRanger-style positional bins over uniquely-mapped reads, populated only
    // when both Gene and GeneFull run (otherwise the split is unavailable).
    use std::sync::atomic::Ordering;
    let total_reads = align_stats.total_reads.load(Ordering::Relaxed);
    let mapped_unique = align_stats.uniquely_mapped.load(Ordering::Relaxed);
    let mapped_multi = align_stats.multi_mapped.load(Ordering::Relaxed);
    let valid_barcodes = ctx.stats.yes_exact.load(Ordering::Relaxed)
        + ctx.stats.yes_one_mm.load(Ordering::Relaxed)
        + ctx.stats.yes_mult_mm.load(Ordering::Relaxed);
    let reads_of = |f: crate::solo::SoloFeature| -> u64 {
        ctx.features
            .iter()
            .position(|&x| x == f)
            .map_or(0, |i| ctx.feature_reads[i].load(Ordering::Relaxed))
    };
    let have_funnel = ctx.features.contains(&crate::solo::SoloFeature::Gene)
        && ctx.features.contains(&crate::solo::SoloFeature::GeneFull);
    let region = have_funnel.then(|| RegionFunnel {
        exonic: ctx.region_stats.exonic.load(Ordering::Relaxed),
        intronic: ctx.region_stats.intronic.load(Ordering::Relaxed),
        intergenic: ctx.region_stats.intergenic.load(Ordering::Relaxed),
        antisense: ctx.region_stats.antisense.load(Ordering::Relaxed),
    });

    // One {prefix}{soloOutFileNames[0]}<feature>/raw/ directory per feature
    // (Gene, GeneFull, …), each fed from its own recorder.
    for (feature, recorder) in ctx.features.iter().zip(&ctx.recorders) {
        let feature_dir = params.output_path(&format!("{solo_dir}{}/", feature.dir_name()));
        let raw_dir = feature_dir.join("raw");
        std::fs::create_dir_all(&raw_dir).map_err(|e| Error::io(e, &raw_dir))?;

        write_features(&raw_dir.join(&features_name), &ctx.gene_ann.gene_ids)?;
        write_barcodes(&raw_dir.join(&barcodes_name), &ctx.whitelist, sorted.len())?;
        let mstats = stream_matrix(
            ctx,
            recorder,
            method,
            filtering,
            umi_len,
            pseudocount,
            &raw_dir.join(&matrix_name),
            ctx.gene_ann.gene_ids.len(),
            sorted.len(),
        )?;

        log::info!(
            "STARsolo: wrote {}/raw matrix to {} ({} genes × {} barcodes, {} entries)",
            feature.dir_name(),
            raw_dir.display(),
            ctx.gene_ann.gene_ids.len(),
            sorted.len(),
            mstats.nnz,
        );

        write_summary(
            &feature_dir.join("Summary.csv"),
            feature.dir_name(),
            &mstats,
            total_reads,
            valid_barcodes,
            mapped_unique,
            mapped_multi,
            reads_of(*feature),
            region,
        )?;
        log::info!("STARsolo: wrote {}/Summary.csv", feature.dir_name());
    }
    Ok(())
}

/// CellRanger-style positional mapping bins over uniquely-mapped reads.
#[derive(Clone, Copy)]
struct RegionFunnel {
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    antisense: u64,
}

/// Write a CellRanger/STARsolo-style `Summary.csv` for one feature: the
/// sequencing/mapping funnel (genome → exonic → intronic → intergenic, antisense)
/// plus per-cell UMI/gene statistics over the CR2.2-knee-called cells.
#[allow(clippy::too_many_arguments)]
fn write_summary(
    path: &Path,
    feature_name: &str,
    mstats: &MatrixStats,
    total_reads: u64,
    valid_barcodes: u64,
    mapped_unique: u64,
    mapped_multi: u64,
    feature_mapped: u64,
    region: Option<RegionFunnel>,
) -> Result<(), Error> {
    let frac = |num: u64| -> f64 {
        if total_reads == 0 {
            0.0
        } else {
            num as f64 / total_reads as f64
        }
    };

    // Cell calling: CR2.2 knee on per-barcode UMI totals.
    let mut umis_desc: Vec<u64> = mstats.cells.iter().map(|c| c.n_umis).collect();
    umis_desc.sort_unstable_by(|a, b| b.cmp(a));
    let thr = knee_cr22(&umis_desc, 3000, 0.99, 10.0);
    let cells: Vec<&CellStat> = mstats.cells.iter().filter(|c| c.n_umis >= thr).collect();
    let n_cells = cells.len();

    // Totals across all barcodes (for sequencing saturation + fraction-in-cells).
    let total_reads_counted: u64 = mstats.cells.iter().map(|c| c.n_reads).sum();
    let total_umis_all: u64 = mstats.cells.iter().map(|c| c.n_umis).sum();
    let saturation = if total_reads_counted > 0 {
        1.0 - total_umis_all as f64 / total_reads_counted as f64
    } else {
        0.0
    };

    // Per-cell aggregates over called cells.
    let reads_in_cells: u64 = cells.iter().map(|c| c.n_reads).sum();
    let umis_in_cells: u64 = cells.iter().map(|c| c.n_umis).sum();
    let mut reads_sorted: Vec<u64> = cells.iter().map(|c| c.n_reads).collect();
    let mut umis_sorted: Vec<u64> = cells.iter().map(|c| c.n_umis).collect();
    let mut genes_sorted: Vec<u64> = cells.iter().map(|c| c.n_genes as u64).collect();
    reads_sorted.sort_unstable();
    umis_sorted.sort_unstable();
    genes_sorted.sort_unstable();
    let mean = |sum: u64| -> u64 {
        if n_cells == 0 {
            0
        } else {
            sum / n_cells as u64
        }
    };

    use std::fmt::Write as _;
    let mut out = String::new();
    let mut row = |k: &str, v: String| {
        let _ = writeln!(out, "{k},{v}");
    };
    row("Number of Reads", total_reads.to_string());
    row(
        "Reads With Valid Barcodes",
        format!("{:.6}", frac(valid_barcodes)),
    );
    row("Sequencing Saturation", format!("{saturation:.6}"));
    row(
        "Reads Mapped to Genome: Unique+Multiple",
        format!("{:.6}", frac(mapped_unique + mapped_multi)),
    );
    row(
        "Reads Mapped to Genome: Unique",
        format!("{:.6}", frac(mapped_unique)),
    );
    row(
        &format!("Reads Mapped to {feature_name}: Unique {feature_name}"),
        format!("{:.6}", frac(feature_mapped)),
    );
    // CellRanger-style positional funnel over uniquely-mapped reads (each region
    // counted by where the read falls, independent of strand; antisense is a
    // separate orientation metric). Available only with Gene + GeneFull.
    if let Some(r) = region {
        row(
            "Reads Mapped Confidently to Exonic Regions",
            format!("{:.6}", frac(r.exonic)),
        );
        row(
            "Reads Mapped Confidently to Intronic Regions",
            format!("{:.6}", frac(r.intronic)),
        );
        row(
            "Reads Mapped Confidently to Intergenic Regions",
            format!("{:.6}", frac(r.intergenic)),
        );
        row(
            "Reads Mapped Antisense to Gene",
            format!("{:.6}", frac(r.antisense)),
        );
    }
    row("Estimated Number of Cells", n_cells.to_string());
    row(
        &format!("Unique Reads in Cells Mapped to {feature_name}"),
        reads_in_cells.to_string(),
    );
    row(
        "Fraction of Unique Reads in Cells",
        format!(
            "{:.6}",
            if total_reads_counted > 0 {
                reads_in_cells as f64 / total_reads_counted as f64
            } else {
                0.0
            }
        ),
    );
    row("Mean Reads per Cell", mean(reads_in_cells).to_string());
    row(
        "Median Reads per Cell",
        median_sorted(&reads_sorted).to_string(),
    );
    row("UMIs in Cells", umis_in_cells.to_string());
    row("Mean UMI per Cell", mean(umis_in_cells).to_string());
    row(
        "Median UMI per Cell",
        median_sorted(&umis_sorted).to_string(),
    );
    row(
        &format!("Mean {feature_name} per Cell"),
        mean(genes_sorted.iter().sum()).to_string(),
    );
    row(
        &format!("Median {feature_name} per Cell"),
        median_sorted(&genes_sorted).to_string(),
    );
    row(
        &format!("Total {feature_name} Detected"),
        mstats.genes_detected.to_string(),
    );

    std::fs::write(path, out).map_err(|e| Error::io(e, path))?;
    Ok(())
}

/// `features.tsv`: `gene_id <TAB> gene_name <TAB> "Gene Expression"` (CellRanger
/// v3 layout). We have no gene names, so the id is repeated.
fn write_features(path: &Path, gene_ids: &[String]) -> Result<(), Error> {
    let mut f =
        std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| Error::io(e, path))?);
    for id in gene_ids {
        writeln!(f, "{id}\t{id}\tGene Expression").map_err(|e| Error::io(e, path))?;
    }
    f.flush().map_err(|e| Error::io(e, path))
}

/// `barcodes.tsv`: one barcode per line in sorted whitelist order (the same
/// order the matrix columns are indexed by).
///
/// This lists the full whitelist (millions of lines), so it MUST be buffered —
/// an unbuffered writer issues one syscall per line and dominates runtime,
/// especially over a virtiofs mount. Barcodes are unpacked into a reused scratch
/// buffer to avoid a `String` allocation per line.
fn write_barcodes(path: &Path, whitelist: &CbWhitelist, n: usize) -> Result<(), Error> {
    use std::io::Write as _;
    let mut f =
        std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| Error::io(e, path))?);
    let len = whitelist.barcode_len();
    let mut line: Vec<u8> = Vec::with_capacity(len + 1);
    for i in 0..n {
        line.clear();
        whitelist.unpack_barcode_into(i as u32, &mut line);
        line.push(b'\n');
        f.write_all(&line).map_err(|e| Error::io(e, path))?;
    }
    f.flush().map_err(|e| Error::io(e, path))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::encode_base;
    use crate::solo::whitelist::pack_barcode;

    #[test]
    fn median_sorted_odd_even_empty() {
        assert_eq!(median_sorted(&[]), 0);
        assert_eq!(median_sorted(&[5]), 5);
        assert_eq!(median_sorted(&[1, 2, 3]), 2);
        assert_eq!(median_sorted(&[10, 20, 30, 40]), 25); // midpoint(20,30)
    }

    #[test]
    fn knee_cr22_threshold() {
        // 100 cells at 1000 UMI, then a long ambient tail at 10.
        let mut umis: Vec<u64> = vec![1000; 100];
        umis.extend(std::iter::repeat_n(10u64, 5000));
        umis.sort_unstable_by(|a, b| b.cmp(a));
        // robust max = umis[round(3000*0.01)] = umis[30] = 1000; thr = 1000/10 = 100.
        let thr = knee_cr22(&umis, 3000, 0.99, 10.0);
        assert_eq!(thr, 100);
        let cells = umis.iter().filter(|&&u| u >= thr).count();
        assert_eq!(cells, 100); // the 100 real cells, none of the ambient tail
    }

    fn umi(s: &str) -> u64 {
        match pack_barcode(&s.bytes().map(encode_base).collect::<Vec<_>>()) {
            crate::solo::whitelist::PackResult::NoN(p) => p,
            _ => panic!("N in test UMI"),
        }
    }

    fn counts(pairs: &[(&str, u32)]) -> HashMap<u64, u32> {
        pairs.iter().map(|&(s, c)| (umi(s), c)).collect()
    }

    #[test]
    fn dedup_method_parsing() {
        assert_eq!("1MM_All".parse::<UmiDedup>().unwrap(), UmiDedup::OneMmAll);
        assert_eq!("Exact".parse::<UmiDedup>().unwrap(), UmiDedup::Exact);
        assert_eq!("NoDedup".parse::<UmiDedup>().unwrap(), UmiDedup::NoDedup);
        assert!("bogus".parse::<UmiDedup>().is_err());
    }

    #[test]
    fn exact_counts_distinct_umis() {
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::Exact, 4), 3);
    }

    #[test]
    fn nodedup_sums_reads() {
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::NoDedup, 4), 9);
    }

    #[test]
    fn one_mm_all_merges_neighbors() {
        // AAAA–AAAC are Hamming-1 (one component); TTTT separate → 2 molecules.
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 2);
    }

    #[test]
    fn one_mm_all_transitive_chain() {
        // AAAA–AAAC–AACC chain: all one component even though AAAA/AACC are 2 apart.
        let c = counts(&[("AAAA", 1), ("AAAC", 1), ("AACC", 1)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 1);
    }

    #[test]
    fn directional_absorbs_low_count_neighbor() {
        // hub AAAA count 5 absorbs AAAC count 1 (5 >= 2*1+0); TTTT survives.
        let c = counts(&[("AAAA", 5), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectional, 4), 2);
        // Equal counts are NOT absorbed (5 >= 2*5 is false).
        let c2 = counts(&[("AAAA", 5), ("AAAC", 5)]);
        assert_eq!(dedup_count(&c2, UmiDedup::OneMmDirectional, 4), 2);
    }

    #[test]
    fn directional_umitools_threshold() {
        // count_hub >= 2*leaf - 1: hub 3 absorbs leaf 2 (3 >= 3). Directional(0)
        // would not (3 >= 4 false).
        let c = counts(&[("AAAA", 3), ("AAAC", 2)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectionalUmiTools, 4), 1);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectional, 4), 2);
    }

    #[test]
    fn cellranger_1mm_collapses_neighbor() {
        // AAAA (5) and AAAC (1) are 1MM → low-count corrected to high-count →
        // 1 molecule. TTTT separate → 2 total.
        let c = counts(&[("AAAA", 5), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmCr, 4), 2);
        assert_eq!("1MM_CR".parse::<UmiDedup>().unwrap(), UmiDedup::OneMmCr);
    }

    #[test]
    fn cellranger_1mm_non_transitive() {
        // Chain AAAA(1)–AAAC(2)–AACC(4): each corrects to its highest-count 1MM
        // neighbor. AAAA→AAAC (only neighbor), AAAC→AACC, AACC→self. Corrected
        // set {AAAC, AACC, AACC} → 2 molecules (NOT 1 like the transitive All).
        let c = counts(&[("AAAA", 1), ("AAAC", 2), ("AACC", 4)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmCr, 4), 2);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 1);
    }

    #[test]
    fn umi_filtering_parsing() {
        assert_eq!("-".parse::<UmiFiltering>().unwrap(), UmiFiltering::None);
        assert_eq!(
            "MultiGeneUMI_CR".parse::<UmiFiltering>().unwrap(),
            UmiFiltering::MultiGeneUmiCr
        );
        assert!("bogus".parse::<UmiFiltering>().is_err());
    }

    #[test]
    fn multi_gene_umi_cr_keeps_top_gene() {
        // UMI maps to gene 0 (3 reads) and gene 1 (1 read). CR keeps only gene 0.
        let mut genes = HashMap::new();
        genes.insert(0u32, 3u32);
        genes.insert(1u32, 1u32);
        let kept = filter_multi_gene_umi(&genes, UmiFiltering::MultiGeneUmiCr);
        assert_eq!(kept.len(), 1);
        assert_eq!(*kept[0].0, 0);
        // Plain MultiGeneUMI with all-singletons drops the UMI entirely.
        let mut single = HashMap::new();
        single.insert(0u32, 1u32);
        single.insert(1u32, 1u32);
        assert_eq!(
            filter_multi_gene_umi(&single, UmiFiltering::MultiGeneUmi).len(),
            0
        );
    }

    #[test]
    fn resolve_multi_prefers_higher_prior() {
        use crate::solo::whitelist::CbCandidate;
        let cands = vec![
            CbCandidate {
                wl_index: 0,
                mismatch_pos: 1,
                mismatch_qual: b'I',
            },
            CbCandidate {
                wl_index: 1,
                mismatch_pos: 2,
                mismatch_qual: b'I',
            },
        ];
        // Same quality → higher exact-count prior wins.
        assert_eq!(resolve_multi_cb(&cands, &[10, 3], 0.0), Some(0));
        assert_eq!(resolve_multi_cb(&cands, &[3, 10], 0.0), Some(1));
        // No prior signal and no pseudocount → rejected.
        assert_eq!(resolve_multi_cb(&cands, &[0, 0], 0.0), None);
        // Pseudocount gives every candidate positive weight → argmax accepted.
        assert!(resolve_multi_cb(&cands, &[0, 0], 1.0).is_some());
    }
}
