mod cli;
mod commands;
mod report;

use std::io::{self, Write};
use std::process::ExitCode;

fn run() -> Result<(), String> {
    let command = cli::parse(std::env::args_os().skip(1))?;
    let mut output = io::BufWriter::new(io::stdout().lock());
    commands::execute(command, &mut output)?;
    output.flush().map_err(|error| error.to_string())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("aiips: {error}");
            ExitCode::FAILURE
        }
    }
}
