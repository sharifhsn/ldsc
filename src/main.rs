mod cli;
mod blas;
mod irwls;
mod jackknife;
mod ldscore;
mod make_annot;
mod munge;
mod parse;
mod regressions;

use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;

use cli::{Cli, Command};

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(n) = cli.rayon_threads {
        if let Err(err) = ThreadPoolBuilder::new().num_threads(n).build_global() {
            eprintln!("warning: failed to set Rayon thread pool size: {}", err);
        }
    }

    if let Some(n) = cli.polars_threads {
        std::env::set_var("POLARS_MAX_THREADS", n.to_string());
    }

    blas::set_openblas_threads(cli.blas_threads);

    match cli.command {
        Command::MungeSumstats(args) => munge::run(args),
        Command::Ldscore(args) => ldscore::run(args),
        Command::H2(args) => regressions::run_h2(args),
        Command::Rg(args) => regressions::run_rg(args),
        Command::MakeAnnot(args) => make_annot::run(args),
    }
}
