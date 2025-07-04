"""
Protein Orthology Comparison Pipeline

This script serves as a command-line interface for analyzing orthologous protein groups
across species using precomputed FASTA and orthology inference data.

Runs the `run_table()` function which performs:
   - Merging and comparison of orthologous groups from multiple inference tools.
   - Dynamic selection of the most compact reference group for species comparison.
   - Evaluation of orthology consistency across tools.
   - Optional integration of tree-based analyses.

Usage example:
    python3 main.py -fd /path/to/fasta -o /path/to/output -i /path/to/proteins.csv -s /path/to/species.csv
    - t /path/to/tree.nwk -pf /path/to/fasta_protein

Arguments:
    -fd / --fasta_dir        : Directory containing translated protein FASTA files.
    -o  / --output_dir       : Directory to save analysis results.
    -i  / --interest_file    : CSV file with proteins of interest.
    -s  / --species_file     : CSV file with species/sample metadata.
    -t  / --tree_file        : Optional path to a phylogenetic tree in Newick format.
    -pf / --protein_fasta    : Optional path to a combined protein FASTA file.

Dependencies:
    - Python libraries: argparse, shutil, pandas.
    - External tools for ortholog inference must be run beforehand.
    - `run_table()` function (imported from `table.py`) performs core analysis.

Author(s): Romain DAGUERRE
"""

from table import run_table

import argparse
from argparse import RawTextHelpFormatter

def _set_args():
    """
    This function sets the parameters provided by the user.
    It handles paths for FASTA input, BUSCO database, outputs, and processing options.
    """
    parser = argparse.ArgumentParser(
        description="Run the \n\n"
                    "Example usage:\n"
                    "python3 main.py -fd /path/to/fasta -db /path/to/busco -o /path/to/output -csv /path/to/results.csv",
        formatter_class=RawTextHelpFormatter
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-fd", "--fasta_dir",
        metavar="PATH",
        required=True,
        help="Directory containing FASTA files."
    )
    required.add_argument(
        "-o", "--output_dir",
        metavar="PATH",
        required=True,
        help="Directory to save all outputs."
    )
    required.add_argument(
        "-i", "--interest_file",
        metavar="PATH",
        required=True,
        help="CSV path with the information from interest proteins."
    )
    required.add_argument(
        "-s", "--species_file",
        metavar="PATH",
        required=True,
        help="Species file."
    )

    optional.add_argument(
        "-t", "--tree_file",
        metavar="PATH",
        help="Path to the phylogenetic tree (active this option if needed)."
    )

    optional.add_argument(
        "-pf", "--protein_fasta",
        metavar="PATH",
        help="Optional path to all protein FASTA file."
    )

    args = vars(parser.parse_args())

    global fasta_dir
    fasta_dir = args["fasta_dir"]
    global output_dir
    output_dir = args["output_dir"]
    global tree_file
    tree_file = args["tree_file"]
    global interest_file
    interest_file = args["interest_file"]
    global species_file
    species_file = args["species_file"]
    global protein_fasta
    protein_fasta = args["protein_fasta"]


if __name__ == "__main__":

    _set_args()

    print(protein_fasta)
    if protein_fasta:
        run_table(fasta_dir, output_dir, tree_file, interest_file, species_file, protein_fasta)
    else:
        run_table(fasta_dir, output_dir, tree_file, interest_file, species_file, None)