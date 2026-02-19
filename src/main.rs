mod cli;
mod parse;
mod munge;
mod ldscore;
mod irwls;
mod jackknife;
mod regressions;
mod make_annot;

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
