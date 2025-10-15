use clap::{Parser, Subcommand};
use log::info;

use scnado::barcodes;

#[derive(Parser, Debug)]
struct FindBarcodesArgs {
    #[arg(short, long, help = "Input FASTQ file")]
    input: String,
    #[arg(
        short,
        long,
        help = "Barcode file. CSV file with columns 'barcode_type' and 'barcode_sequence'"
    )]
    barcodes: String,
    #[arg(short, long, help = "Output file")]
    output: String,
    #[arg(long, help = "Slack left", default_value = "0")]
    slack_left: Option<usize>,
    #[arg(long, help = "Slack right", default_value = "0")]
    slack_right: Option<usize>,
    #[arg(
        short,
        long,
        help = "Number of allowed mismatches",
        default_value = "0"
    )]
    n_missmatches: Option<usize>,
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Verbosity level
    #[arg(short, long, required = false, default_value = "2")]
    verbose: u8,
}

#[derive(Subcommand)]
enum Commands {
    /// Find barcodes in FASTQ file
    FindBarcodes(FindBarcodesArgs),
}

fn main() {
    colog::init();
    log::set_max_level(log::LevelFilter::Info);
    log::info!("Starting barcode detection");
    let cli = Cli::parse();

    match &cli.command {
        Commands::FindBarcodes(args) => {
            let slack =
                barcodes::Slack::new(args.slack_left.unwrap_or(0), args.slack_right.unwrap_or(0));
            info!("Finding barcodes in {}", args.input);
            if let Err(e) = barcodes::run(
                &args.input,
                &args.barcodes,
                &args.output,
                Some(slack),
                args.n_missmatches,
            ) {
                eprintln!("Error: {:?}", e);
                std::process::exit(1);
            }
        }
    }
}
