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
    // kernel zeroing — ~12% of alignment time (clear_page_erms) in profiling. Disable
    // purging so freed pages are retained and reused; the extra heap footprint is
    // negligible next to the mmap'd index (measured: Max RSS unchanged, wall ~8%
    // faster). `mi_option_set` overrides the cached default, so this takes effect even
    // though mimalloc initialized before main(). `mi_option_purge_delay` is option 15
    // in mimalloc v2.x (not exposed as a named const by libmimalloc-sys 0.1.49);
    // value -1 = never purge. SAFETY: mi_option_t is a plain c_int; trivial FFI call.
    const MI_OPTION_PURGE_DELAY: std::os::raw::c_int = 15;
    unsafe { libmimalloc_sys::mi_option_set(MI_OPTION_PURGE_DELAY, -1) };

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    rustar_aligner::run(&Parameters::parse())
}
