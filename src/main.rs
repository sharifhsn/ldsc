mod bed;
mod cli;
mod cts_annot;
#[cfg(feature = "gpu")]
mod gpu;
mod h2;
mod irwls;
mod jackknife;
mod l2;
mod la;
mod make_annot;
mod munge;
mod parse;
mod regressions;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;
use tracing_subscriber::EnvFilter;

use cli::{Cli, Command};

fn init_tracing(level: &str) {
    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(format!("ldsc={}", level)));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .init();
}

fn is_subcommand(arg: &str) -> bool {
    matches!(
        arg,
        "munge-sumstats" | "l2" | "h2" | "rg" | "make-annot" | "cts-annot"
    )
}

fn inject_h2_subcommand(mut args: Vec<String>) -> Vec<String> {
    if args.len() < 2 {
        return args;
    }
    let first = args[1].as_str();
    if !is_subcommand(first) && args.iter().any(|a| a == "--h2") {
        args.insert(1, "h2".to_string());
    } else if !is_subcommand(first) && args.iter().any(|a| a == "--l2") {
        args.insert(1, "l2".to_string());
    }
    args
}

fn main() -> Result<()> {
    let args = inject_h2_subcommand(std::env::args().collect());
    let cli = Cli::parse_from(args);
    init_tracing(&cli.log_level);

    if let Some(n) = cli.rayon_threads
        && let Err(err) = ThreadPoolBuilder::new().num_threads(n).build_global()
    {
        eprintln!("warning: failed to set Rayon thread pool size: {}", err);
    }

    if let Some(n) = cli.polars_threads {
        // Setting environment variables is unsafe in Rust 2024 due to potential
        // data races if other threads read env concurrently.
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", n.to_string());
        }
    }

    match cli.command {
        Command::MungeSumstats(args) => munge::run(args),
        Command::L2(args) => l2::run(args),
        Command::H2(args) => regressions::run_h2(args),
        Command::Rg(args) => regressions::run_rg(args),
        Command::MakeAnnot(args) => make_annot::run(args),
        Command::CtsAnnot(args) => cts_annot::run(args),
    }
}
