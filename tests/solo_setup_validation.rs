//! CLI behaviour of the solo barcode-setup checks.
//!
//! A whitelist for the wrong chemistry, or a CB/UMI geometry that does not match
//! the barcode read, produces no error of its own: the matrices come out nearly
//! empty and the run exits 0. These tests pin the two guards against that, at the
//! level a user meets them — the process exit status and what it printed — rather
//! than the arithmetic behind them, which is unit-tested in `src/lib.rs`.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const CB: &str = "AAAACCCCGGGGTTTT";
const UMI: &str = "ACGTACGTAC";
const GENE: (usize, usize) = (2_000, 2_600);
const READ_LEN: usize = 60;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut seq = Vec::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

/// A 20 kb single-contig genome with one gene, and its GTF.
fn write_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let genome = lcg_seq(88888, 20_000);
    let fasta = dir.join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(&genome).unwrap();
    writeln!(f).unwrap();

    let gtf = dir.join("genes.gtf");
    let mut g = fs::File::create(&gtf).unwrap();
    writeln!(
        g,
        "chr1\tsyn\texon\t{}\t{}\t.\t+\t.\tgene_id \"G1\"; transcript_id \"G1_T1\";",
        GENE.0 + 1,
        GENE.1
    )
    .unwrap();
    (fasta, gtf)
}

fn build_index(fasta: &Path, gtf: &Path, genome_dir: &Path) {
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--genomeSAindexNbases",
            "7",
            "--runThreadN",
            "1",
        ])
        .assert()
        .success();
}

/// `n` reads from the gene, every one carrying `cb` + [`UMI`].
fn write_reads(dir: &Path, cb: &str, n: usize) {
    let genome = lcg_seq(88888, 20_000);
    let mut cf = fs::File::create(dir.join("cdna.fq")).unwrap();
    let mut bf = fs::File::create(dir.join("barcode.fq")).unwrap();
    for i in 0..n {
        let start = GENE.0 + 20 + (i % 5) * 7;
        writeln!(cf, "@read{i}").unwrap();
        cf.write_all(&genome[start..start + READ_LEN]).unwrap();
        writeln!(cf, "\n+\n{}", "I".repeat(READ_LEN)).unwrap();
        writeln!(
            bf,
            "@read{i}\n{cb}{UMI}\n+\n{}",
            "I".repeat(cb.len() + UMI.len())
        )
        .unwrap();
    }
}

fn write_whitelist(path: &Path, barcodes: &[&str]) {
    let mut w = fs::File::create(path).unwrap();
    for b in barcodes {
        writeln!(w, "{b}").unwrap();
    }
}

/// Run solo against `whitelist` with an explicit `--soloCBlen`, returning the
/// finished command assertion for the caller to judge.
fn run(
    dir: &Path,
    genome_dir: &Path,
    gtf: &Path,
    tag: &str,
    whitelist: &Path,
    cb_len: &str,
) -> std::process::Output {
    let out = dir.join(format!("out_{tag}"));
    fs::create_dir_all(&out).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            dir.join("cdna.fq").to_str().unwrap(),
            dir.join("barcode.fq").to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            whitelist.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            cb_len,
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "10",
            "--soloFeatures",
            "Gene",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outSAMtype",
            "None",
            "--runThreadN",
            "1",
            "--outFileNamePrefix",
            &format!("{}/", out.display()),
        ])
        // Pin the level: the warning is `warn`, which survives an inherited
        // RUST_LOG=warn, but an inherited `error` would hide it and the test
        // would pass for the wrong reason.
        .env("RUST_LOG", "info")
        .output()
        .unwrap()
}

/// A whitelist whose barcodes cannot be the configured length is refused before
/// any alignment work, because no read could ever match it.
#[test]
fn a_whitelist_of_the_wrong_barcode_length_is_rejected() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);
    write_reads(root, CB, 40);

    // 16-base whitelist, --soloCBlen 12.
    let wl = root.join("wl16.txt");
    write_whitelist(&wl, &[CB]);
    let out = run(root, &genome_dir, &gtf, "wronglen", &wl, "12");

    assert!(
        !out.status.success(),
        "a whitelist that can never match must fail the run"
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("16-base barcodes") && err.contains("--soloCBlen is 12"),
        "the error must name both lengths: {err}"
    );
}

/// A whitelist of the right length that simply holds different barcodes is a
/// valid run — it must not fail — but it earns the diagnostic.
#[test]
fn a_whitelist_that_matches_nothing_warns_without_failing() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);
    write_reads(root, CB, 40);

    // Same length as the reads' barcode, but none of them is CB.
    let wl = root.join("wl_other.txt");
    write_whitelist(
        &wl,
        &["TTTTGGGGCCCCAAAA", "GGGGTTTTAAAACCCC", "CCCCAAAATTTTGGGG"],
    );
    let out = run(root, &genome_dir, &gtf, "nomatch", &wl, "16");

    assert!(
        out.status.success(),
        "a low match rate is a diagnostic, not a failure: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("carry a whitelisted barcode"),
        "the diagnostic must be emitted: {err}"
    );
    // The advice has to name the layout this run selected, not assume one.
    assert!(
        err.contains("--soloBarcodeMate 0"),
        "the warning must name the barcode layout: {err}"
    );
    // Claims the review found overreaching must not come back.
    assert!(
        !err.contains("will be close to empty"),
        "must not assert the matrix is necessarily empty: {err}"
    );
    assert!(
        !err.contains("does not match a 5'"),
        "must not claim 3'/5' whitelists are categorically incompatible: {err}"
    );
}

/// The matching whitelist on the same reads stays silent, so the diagnostic is
/// not simply always on.
#[test]
fn the_right_whitelist_produces_no_diagnostic() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);
    write_reads(root, CB, 40);

    let wl = root.join("wl_right.txt");
    write_whitelist(&wl, &[CB, "TTTTGGGGCCCCAAAA"]);
    let out = run(root, &genome_dir, &gtf, "match", &wl, "16");

    assert!(out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        !err.contains("carry a whitelisted barcode"),
        "a matching whitelist must not warn: {err}"
    );
}
