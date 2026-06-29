use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

/// Global allocator override — see the comment in `Cargo.toml`
/// next to the `mimalloc` dependency for the rationale. The
/// `#[global_allocator]` attribute applies to the whole binary; no
/// further work needed at allocation sites.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    // mimalloc returns freed pages to the OS aggressively by default (~10 ms purge
    // delay). On this high-throughput aligner that means constant page re-faulting +
    // kernel zeroing — ~12% of alignment time (clear_page_erms) in profiling. We
    // default to never-purge (-1) for speed (~13% faster). The tradeoff: mimalloc then
    // retains a resident free-page reserve that inflates Max RSS (reusable scratch, not
    // data — the actual heap footprint is ~3 GB, vs a reclaimable ~25 GB mmap'd index).
    // Set RUSTAR_PURGE_DELAY_MS to a finite value (e.g. 1000) to cap that reserve at a
    // small wall-time cost. `mi_option_set` overrides the cached default, so it takes
    // effect even though mimalloc initialized before main(). `mi_option_purge_delay` is
    // option 15 in mimalloc v2.x (not a named const in libmimalloc-sys 0.1.49).
    // SAFETY: mi_option_t is a plain c_int; trivial FFI call.
    const MI_OPTION_PURGE_DELAY: std::os::raw::c_int = 15;
    let purge_ms: std::os::raw::c_long = std::env::var("RUSTAR_PURGE_DELAY_MS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(-1);
    unsafe { libmimalloc_sys::mi_option_set(MI_OPTION_PURGE_DELAY, purge_ms) };

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    rustar_aligner::run(&Parameters::parse())
}
