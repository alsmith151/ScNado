use clap::{Parser, Subcommand};
// use log::info;

use scnado::barcodes;

#[derive(Parser, Debug)]
struct FindBarcodesArgs {
    #[arg(long, help = "Input FASTQ file")]
    input_r1: String,
    #[arg(long, help = "Input FASTQ file (R2)")]
    input_r2: String,
    #[arg(
        long,
        help = "Barcode file. CSV file with columns 'barcode_type' and 'barcode_sequence'"
    )]
    barcodes: String,
    #[arg(short, long, help = "Output FASTQ file prefix")]
    output: String,
    #[arg(
        short,
        long,
        help = "Number of allowed mismatches",
        default_value = "0"
    )]
    n_missmatches: Option<usize>,

    #[arg(
        long,
        help = "Allow 'N' characters in barcode matching",
        default_value_t = false
    )]
    enable_n_matching: bool,

    #[arg(long, help = "Allow missing barcodes", default_value_t = false)]
    allow_missing_barcodes: bool,
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
            if let Err(e) = barcodes::run(
                &args.input_r1,
                &args.input_r2,
                &args.barcodes,
                args.output.clone(),
                args.n_missmatches.unwrap_or(0),
                args.enable_n_matching,
                args.allow_missing_barcodes,
            ) {
                eprintln!("Error: {:?}", e);
                std::process::exit(1);
            }
        }
    }
}
