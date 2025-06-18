"""
FASTA Processing Pipeline

This script provides a comprehensive pipeline for managing a collection of FASTA files,
with optional integration of BUSCO analysis and ortholog inference.

Pipeline steps:
1. Detect and translate nucleotide sequences into protein sequences using EMBOSS Transeq.
   - Clean and format FASTA headers.
   - Rename FASTA files based on a reference CSV file.
   - Add organism-specific prefixes to FASTA headers.
   - Update the CSV file accordingly.

2. (Optional) Assess genome completeness with BUSCO.
   - Runs BUSCO for each translated FASTA using a specified BUSCO database.
   - Collects and renames summary figures to a common output directory.

3. (Optional) Infer ortholog groups using three tools: SonicParanoid, Orthologer, and OrthoFinder.
   - Parses the orthogroups from each tool's output.
   - Formats orthogroup data into standardized CSV files.
   - Only orthogroups containing at least two proteins are retained.

Usage example:
    python3 main.py -fd /path/to/fasta -db /path/to/busco -o /path/to/output -s /path/to/species.csv

Arguments:
    -fd / --fasta_dir      : Input directory containing nucleotide FASTA files.
    -db / --busco_db       : Path to the BUSCO lineage dataset (unless --no_busco is specified).
    -o  / --output_dir     : Directory to save all outputs (translated files, figures, results).
    -s  / --species_file   : CSV file with species/sample metadata.
    -t  / --force_type     : Optional flag to force the type of sequence (e.g. for BUSCO input).
    -b  / --no_busco       : Optional flag to skip BUSCO analysis.
    -ortho / --no_orthology: Optional flag to skip orthology inference and downstream analysis.

Dependencies:
    - Python libraries: argparse, os, glob, shutil, pandas.
    - EMBOSS Transeq (must be available in PATH).
    - External tools for ortholog detection and BUSCO must be installed.

Author(s): Romain DAGUERRE
"""

from translate_fasta import run_fasta_translation
from Busco import run_busco_analysis
from orthologs import run_orthologs
from output_orthologs import run_analysis

import argparse
import os
import glob
import shutil
from argparse import RawTextHelpFormatter


def find_busco_figures(root_dir):
    """
    Search for all 'busco_figure.png' files located in 'my_summaries' subdirectories
    under a given root directory.

    Parameters
    ----------
    root_dir : str
        Root directory containing 'my_summaries' folders.

    Returns
    -------
    list
        List of paths to the located 'busco_figure.png' files.
    """
    search_pattern = os.path.join(root_dir, "my_summaries*/busco_figure.png")
    return glob.glob(search_pattern)


def copy_busco_figures(busco_files, target_dir):
    """
    Copy all 'busco_figure.png' files to the target directory. If duplicates exist,
    the files are renamed with incremental suffixes.

    Parameters
    ----------
    busco_files : list
        List of file paths to copy.
    target_dir : str
        Destination directory where the files will be copied.
    """
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for file in busco_files:
        base_name = os.path.basename(file)
        target_file = os.path.join(target_dir, base_name)

        counter = 2
        while os.path.exists(target_file):
            new_name = f"{os.path.splitext(base_name)[0]}_{counter}{os.path.splitext(base_name)[1]}"
            target_file = os.path.join(target_dir, new_name)
            counter += 1

        shutil.copy(file, target_file)


def _set_args():
    """
    Set and parse command-line arguments provided by the user.

    This function manages parameters for FASTA inputs, BUSCO database,
    analysis outputs, and user-defined options.
    """
    parser = argparse.ArgumentParser(
        description="Run FASTA translation and optionally BUSCO and orthology analyses.\n\n"
        "Example usage:\n"
        "python3 main.py -fd /path/to/fasta -db /path/to/busco -o /path/to/output -csv /path/to/results.csv",
        formatter_class=RawTextHelpFormatter,
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-fd",
        "--fasta_dir",
        metavar="PATH",
        required=True,
        help="Directory containing FASTA files.",
    )
    required.add_argument(
        "-o",
        "--output_dir",
        metavar="PATH",
        required=True,
        help="Directory to save all outputs.",
    )
    required.add_argument(
        "-s", "--species_file", metavar="PATH", required=True, help="species file."
    )

    optional.add_argument(
        "-db",
        "--busco_db",
        metavar="PATH",
        help="Path to the BUSCO database (required unless --no_busco is set).",
    )
    optional.add_argument(
        "-t",
        "--force_type",
        action="store_true",
        help="Force the type (activate this option if needed).",
    )
    optional.add_argument(
        "-b", "--no_busco", action="store_true", help="Skip BUSCO analysis."
    )
    optional.add_argument(
        "-ortho",
        "--no_orthology",
        action="store_true",
        help="Skip orthology inference and analysis.",
    )

    args = vars(parser.parse_args())

    global fasta_dir
    fasta_dir = args["fasta_dir"]
    global busco_db
    busco_db = args["busco_db"]
    global output_dir
    output_dir = args["output_dir"]
    global species_file
    species_file = args["species_file"]

    global force_type
    force_type = False
    if args["force_type"]:
        force_type = True
    global no_busco
    no_busco = args["no_busco"]

    global no_orthology
    no_orthology = False
    if args["no_orthology"]:
        no_orthology = True

    if not no_busco and not busco_db:
        parser.error(
            "The --busco_db (-db) argument is required unless --no_busco is specified."
        )


if __name__ == "__main__":

    _set_args()

    fasta_output = os.path.join(output_dir, "Fasta")

    # Step 1: Translate FASTA sequences using species metadata
    run_fasta_translation(fasta_dir, fasta_output, species_file)

    # Step 2: Run BUSCO analysis unless explicitly skipped
    if not no_busco:
        output_Busco = os.path.join(output_dir, "Busco")
        run_busco_analysis(busco_db, output_Busco, species_file, fasta_output, force_type)

        busco_files = find_busco_figures(output_Busco)
        target_directory = os.path.join(output_dir, "Busco_results")
        copy_busco_figures(busco_files, target_directory)
    else:
        print("BUSCO analysis skipped")

    # Step 3: Run orthology inference unless explicitly skipped
    if not no_orthology:
        input_sonic = os.path.join(output_dir, "sonicparanoid")
        input_orthologer = os.path.join(output_dir, "orthologer")
        input_orthofinder = os.path.join(output_dir, "orthofinder")

        run_orthologs(fasta_output, input_sonic, input_orthologer, input_orthofinder)
        run_analysis(input_sonic, input_orthologer, input_orthofinder, output_dir)
    else:
        print("Orthology analysis skipped")
