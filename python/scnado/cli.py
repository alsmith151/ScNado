import argparse
import sys
from .cat import process_cat, analyze_cat_dataset
from .rna import process_rna
from .multiome import integrate_cat_rna

def main():
    parser = argparse.ArgumentParser(description="scnado Python CLI")
    subparsers = parser.add_subparsers(dest="command")
    
    # Add subcommands as needed for the pipeline
    # For now, we'll let Snakemake call the functions directly via 'script:' or 'run:'
    # but having a CLI is good practice.
    
    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
