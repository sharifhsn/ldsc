mod cli;
mod irwls;
mod jackknife;
mod ldscore;
mod make_annot;
mod munge;
mod parse;
mod regressions;

use anyhow::Result;
use clap::Parser;

use cli::{Cli, Command};

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::MungeSumstats(args) => munge::run(args),
        Command::Ldscore(args) => ldscore::run(args),
        Command::H2(args) => regressions::run_h2(args),
        Command::Rg(args) => regressions::run_rg(args),
        Command::MakeAnnot(args) => make_annot::run(args),
    }
}
