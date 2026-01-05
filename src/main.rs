use clap::{Parser, Subcommand};

mod barcodes;
mod fragments;

#[derive(Parser, Debug)]
struct FindBarcodesArgs {
    #[arg(long, help = "Input FASTQ file")]
    read1: String,
    #[arg(long, help = "Input FASTQ file (R2)")]
    read2: String,
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
    enable_n_to_match: bool,

    #[arg(long, help = "Allow missing barcodes", default_value_t = false)]
    allow_missing_barcodes: bool,
}

#[derive(Debug, Parser)]
struct Fragments {
    #[arg(long, help = "Input BAM file")]
    bam: String,
    #[arg(long, help = "Output fragments file")]
    output: String,
    #[arg(
        long,
        help = "Regex to extract barcode from read name. Capture group 1 is used."
    )]
    barcode_regex: Option<String>,
    #[arg(long, help = "Shift plus", default_value_t = 0)]
    shift_plus: i64,
    #[arg(long, help = "Shift minus", default_value_t = 0)]
    shift_minus: i64,
    #[arg(
        long,
        help = "If set extend each read to a fragment by adding this length"
    )]
    fragment_length_extension: Option<usize>,
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
    /// Extract fragments from BAM file
    Fragments(Fragments),
}

fn main() {
    colog::init();
    log::set_max_level(log::LevelFilter::Info);
    log::info!("Starting barcode detection");
    let cli = Cli::parse();

    match &cli.command {
        Commands::FindBarcodes(args) => {
            if let Err(e) = barcodes::run(
                &args.read1,
                &args.read2,
                &args.barcodes,
                args.output.clone(),
                args.n_missmatches.unwrap_or(0),
                args.enable_n_to_match,
                args.allow_missing_barcodes,
            ) {
                eprintln!("Error: {:?}", e);
                std::process::exit(1);
            }
        }
        Commands::Fragments(args) => {
            if let Err(e) = fragments::run(
                &args.bam,
                args.barcode_regex
                    .as_ref()
                    .map(|s| regex::Regex::new(s).unwrap()),
                &args.output,
                args.shift_plus,
                args.shift_minus,
                args.fragment_length_extension,
            ) {
                eprintln!("Error: {:?}", e);
                std::process::exit(1);
            }
        }
    }
}
