"""
BUSCO Completeness Assessment Pipeline for Proteomes

This script provides an automated pipeline to assess the completeness of proteomes using the BUSCO tool.

Pipeline steps:
- Parse an input CSV file to extract organism metadata and proteome file paths.
- Run BUSCO on each proteome file in parallel using multithreading.
- Organize results into batches and generate summary plots with `gen_plot.py`.
- Extract and append BUSCO completeness scores to the original CSV file.

Dependencies:
- BUSCO must be installed and available in the system PATH.
- A compatible BUSCO database must be downloaded and accessible.

Author(s): Romain DAGUERRE

"""

import subprocess
import os
import csv
import sys
from concurrent.futures import ThreadPoolExecutor
import math
import pandas as pd
import re


def parse_csv_to_dict(csv_file):
    """Parse the input CSV file to extract organisms and associated files.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing the proteomes to analyze.

    Returns
    -------
    dict
        Dictionary with organism names as keys and metadata (excluding 'Organism') as values.
    """
    organism_dict = {}

    with open(csv_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        headers = reader.fieldnames

        if not headers or "Organism" not in headers or "File" not in headers:
            raise ValueError(
                "CSV file must contain at least 'Organism' and 'File' columns."
            )

        for row in reader:
            organism_name = row["Organism"]
            organism_dict[organism_name] = {
                key: value for key, value in row.items() if key != "Organism"
            }

    return organism_dict


def run_busco(proteome_file, busco_db, output_dir, organism):
    """Run BUSCO on a proteome file.

    Parameters
    ----------
    proteome_file : str
        Path to the FASTA proteome file.
    busco_db : str
        BUSCO database to use.
    output_dir : str
        Directory to save BUSCO results.
    organism : str
        Name of the organism being analyzed.
    """
    if not os.path.exists(proteome_file):
        raise FileNotFoundError(f"The file {proteome_file} was not found.")

    os.makedirs(output_dir, exist_ok=True)

    output_name = organism.replace(" ", "_")

    busco_command = [
        "busco",
        "-f",
        "-i", proteome_file,
        "-o", output_name,  # Juste le nom
        "-l", busco_db,
        "-m", "proteins",
        "--cpu", "8",
        "--out_path", output_dir  # Le chemin séparé
    ]

    try:
        with open(os.devnull, "w") as FNULL:
            subprocess.run(busco_command, stdout=FNULL, stderr=FNULL)
        print(f"BUSCO completed for {organism}")
    except subprocess.CalledProcessError as e:
        print(f"BUSCO error for {organism}: {e}")


def busco_analysis(
    spec_names, busco_db, output_dir, batch_index, type_csv, force_type=False
):
    """Aggregate BUSCO results for a batch of organisms and generate a summary plot.

    Parameters
    ----------
    spec_names : list
        List of organism names in the current batch.
    busco_db : str
        BUSCO database used.
    output_dir : str
        Directory containing BUSCO result folders.
    batch_index : int
        Index of the current batch.
    type_csv : str
        Path to the original input CSV.
    force_type : bool, optional
        Optional flag for forcing the type in the visualization tool (default: False).
    """
    os.makedirs(f"{output_dir}/my_summaries{batch_index}", exist_ok=True)

    for spec_name in spec_names:
        normalized_name = spec_name.replace(" ", "_")
        summary_file = f"{output_dir}/{normalized_name}/short_summary.specific.{busco_db}.{normalized_name}.txt"
        subprocess.run(
            ["cp", summary_file, f"{output_dir}/my_summaries{batch_index}/"], check=True
        )

    args = [
        "python3",
        "gen_plot.py",
        "-wd",
        f"{output_dir}/my_summaries{batch_index}",
        "-csv",
        f"{type_csv}",
    ]

    if force_type:
        args.append("-t")

    subprocess.run(args, check=True)


def multiprocess_busco(
    batch_size,
    organisms,
    proteome_dict,
    busco_db,
    output_dir,
    proteome_dir,
    type_csv,
    force_type=False,
    max_workers=5,
):
    """Run BUSCO in parallel using threads on multiple proteomes.

    Parameters
    ----------
    batch_size : int
        Number of organisms per batch.
    organisms : list
        List of organism names to analyze.
    proteome_dict : dict
        Dictionary linking each organism to its metadata.
    busco_db : str
        BUSCO database to use.
    output_dir : str
        Output directory for BUSCO results.
    proteome_dir : str
        Directory containing the proteome files.
    type_csv : str
        Path to the original CSV file for visualization purposes.
    force_type : bool, optional
        Whether to force a specific visualization type (default: False).
    max_workers : int, optional
        Maximum number of threads to use (default: 5).
    """
    with ThreadPoolExecutor(max_workers) as executor:
        for batch_num, i in enumerate(range(0, len(organisms), batch_size), start=1):
            batch = organisms[i : i + batch_size]

            futures = []
            for organism in batch:
                if "File" not in proteome_dict[organism]:
                    print(f"Error: 'File' key missing for {organism}. Skipping.")
                    continue

                file_name = proteome_dict[organism]["File"]

                if "Database" in proteome_dict[organism]:
                    database_name = proteome_dict[organism]["Database"]
                    proteome_file = os.path.join(proteome_dir, database_name, file_name)
                else:
                    proteome_file = os.path.join(proteome_dir, file_name)

                if not os.path.exists(proteome_file):
                    print(f"File not found: {proteome_file} for {organism}. Skipping.")
                    continue

                futures.append(
                    executor.submit(
                        run_busco, proteome_file, busco_db, output_dir, organism
                    )
                )

            for future in futures:
                future.result()

            busco_analysis(batch, busco_db, output_dir, batch_num, type_csv, force_type)


def best_batch_size(n, max_batch=70):
    """Determine the optimal batch size for BUSCO processing.

    Parameters
    ----------
    n : int
        Total number of items to process.
    max_batch : int, optional
        Maximum number of batches (default: 70).

    Returns
    -------
    int
        Optimal batch size.
    """
    min_b = math.ceil(n / max_batch)
    return math.ceil(n / min_b)


def count_files_in_folder(folder_path):
    """Count the number of files in a given folder.

    Parameters
    ----------
    folder_path : str
        Path to the folder.

    Returns
    -------
    int
        Number of files in the folder.
    """
    return sum(
        1
        for entry in os.listdir(folder_path)
        if os.path.isfile(os.path.join(folder_path, entry))
    )


def add_busco(csv_file, busco_dir):
    """Add BUSCO completeness scores to the input CSV file.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file to update.
    busco_dir : str
        Path to the directory containing BUSCO result files.
    """
    df = pd.read_csv(csv_file, sep=",")
    busco_data = {}

    for root, _, files in os.walk(busco_dir):
        for file in files:
            if file.endswith(".txt"):
                file_path = os.path.join(root, file)
                fasta_file = None
                busco_score = "NA"

                with open(file_path, "r") as f:
                    for line in f:
                        if "for file" in line:
                            match_fasta = re.search(r"for file (.+\.fasta)", line)
                            if match_fasta:
                                fasta_file = os.path.basename(match_fasta.group(1))

                        # Extract the BUSCO score
                        match_score = re.search(r"C:(\d+\.\d+)%", line)
                        if match_score:
                            busco_score = float(match_score.group(1))

                if fasta_file:
                    busco_data[fasta_file] = busco_score

    df["BUSCO_score"] = df["File"].apply(lambda x: busco_data.get(x, "NA"))
    df.to_csv(csv_file, sep=",", index=False)


def run_busco_analysis(
    busco_db, output_dir, proteome_csv_path, proteome_dir, force_type
):
    """Launch the full BUSCO analysis pipeline on a set of proteomes.

    Parameters
    ----------
    busco_db : str
        Name of the BUSCO database to use.
    output_dir : str
        Directory where BUSCO results will be saved.
    proteome_csv_path : str
        Path to the CSV file containing proteome metadata.
    proteome_dir : str
        Directory containing the proteome files.
    force_type : bool
        Whether to force a specific visualization type.
    """
    if not os.path.exists(proteome_csv_path):
        print(f"Error: file '{proteome_csv_path}' does not exist.")
        sys.exit(1)

    if not os.path.exists(proteome_dir):
        print(f"Error: directory '{proteome_dir}' does not exist.")
        sys.exit(1)

    proteome_dict = parse_csv_to_dict(proteome_csv_path)
    organisms = sorted(proteome_dict.keys())
    n_files = count_files_in_folder(proteome_dir)
    batch_size = best_batch_size(n_files)

    multiprocess_busco(
        batch_size,
        organisms,
        proteome_dict,
        busco_db,
        output_dir,
        proteome_dir,
        proteome_csv_path,
        force_type,
    )

    add_busco(proteome_csv_path, output_dir)
