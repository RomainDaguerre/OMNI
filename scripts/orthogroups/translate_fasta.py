"""
FASTA Processing Pipeline

This script provides a comprehensive pipeline for managing a collection of FASTA files.

Pipeline steps :
- Detect and translate nucleotide sequences into protein sequences using EMBOSS Transeq.
- Clean and format FASTA headers.
- Rename FASTA files based on a reference CSV file.
- Add organism-specific prefixes (particles) to FASTA headers.
- Update the CSV file accordingly with renamed files.

Dependencies:
- EMBOSS Transeq must be installed and available in the system PATH.


Author(s): Romain DAGUERRE

"""

import os
import sys
import subprocess
import pandas as pd
import shutil


def is_nucleotide_sequence(sequence, threshold=0.9):
    """
    Detects whether a sequence is DNA or protein.

    Parameters
    ----------
    sequence : str
        String containing a biological sequence.
    threshold : float
        Proportion threshold to consider a sequence as DNA (default 90%).

    Returns
    -------
    boolean
        True if DNA sequence, False if protein.
    """
    sequence = sequence.upper()
    dna_chars = set("ATCGN")
    valid_count = sum(1 for char in sequence if char in dna_chars)
    return (valid_count / len(sequence)) >= threshold


def should_translate_fasta(fasta_path):
    """
    Checks whether a FASTA file is nucleotide and needs translation.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    boolean
        True if DNA and should be translated, False if protein.
    """
    with open(fasta_path, "r") as file:
        for line in file:
            if not line.startswith(">"):
                return is_nucleotide_sequence(line.strip())
    return False

def reformat_fasta_one_line(fasta_file):
    """
    Reformat a FASTA file so that each sequence is on a single line.

    This function:
    - Reads a FASTA file.
    - Concatenates multiline sequences into a single line per sequence.
    - Keeps the header lines (starting with '>') intact.
    - Overwrites the original file with the reformatted content.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file to reformat.
    """
    with open(fasta_file, "r") as f:
        lines = f.readlines()

    reformatted_lines = []
    current_seq = []

    for line in lines:
        if line.startswith(">"):
            if current_seq:
                reformatted_lines.append("".join(current_seq) + "\n")
                current_seq = []
            reformatted_lines.append(line.strip() + "\n")
        else:
            current_seq.append(line.strip())

    if current_seq:
        reformatted_lines.append("".join(current_seq) + "\n")

    with open(fasta_file, "w") as f:
        f.writelines(reformatted_lines)


def translate_fasta(fasta_path):
    """
    Translate a DNA FASTA file into protein using EMBOSS Transeq and overwrite the original file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"File {fasta_path} does not exist.")

    temp_output = fasta_path + ".tmp"

    command = [
        "transeq",
        "-sequence",
        fasta_path,
        "-outseq",
        temp_output,
        "-frame",
        "6",
    ]

    try:
        subprocess.run(command, check=True)
        os.replace(temp_output, fasta_path)
        reformat_fasta_one_line(fasta_path) 
        print(f"Translation completed: {fasta_path}")
    except subprocess.CalledProcessError as e:
        print(f"Transeq error on {fasta_path}: {e}")
        raise


def process_fasta_files(fasta_directory):
    """
    Loop through all FASTA files in a folder and translate those that contain nucleotide sequences.

    Parameters
    ----------
    fasta_directory : str
        Directory containing the FASTA files.
    """
    for filename in os.listdir(fasta_directory):
        if filename.endswith((".fasta", ".fa", ".fna")):
            fasta_path = os.path.join(fasta_directory, filename)
            if should_translate_fasta(fasta_path):
                print(f"Translation in progress: {filename}")
                translate_fasta(fasta_path)
            else:
                print(f"Already translated: {filename}")


def clean_fasta_file(file_path):
    """
    Ensure FASTA headers start with a single '>'.

    Parameters
    ----------
    file_path : str
        Path to the FASTA file.
    """
    with open(file_path, "r") as infile:
        lines = infile.readlines()

    cleaned_lines = []
    for line in lines:
        if line.startswith(">"):
            line = ">" + line.lstrip(">").strip() + "\n"
        cleaned_lines.append(line)

    with open(file_path, "w") as outfile:
        outfile.writelines(cleaned_lines)


def rename_fasta_files(source_folder, output_folder, csv_file):
    """
    Rename and organize FASTA files based on a reference CSV.

    This function performs the following:
    - Reads a CSV containing two columns: 'Previous_file' (original FASTA filename) 
      and 'Organism' (species name).
    - Renames each FASTA file in the source folder using the 'Organism' name, formatted 
      in lowercase with underscores and '.fasta' extension.
    - Copies the renamed files into the output folder.
    - Adds or updates the 'File' column in the CSV with the new filenames.
    - Skips and does not copy files that are not listed in the CSV.
    - If the 'File' column already exists in the CSV and all expected files are present 
      in the output folder, the function does nothing.

    Parameters
    ----------
    source_folder : str
        Path to the folder containing the original FASTA files.
    output_folder : str
        Path where the renamed FASTA files will be saved.
    csv_file : str
        Path to the CSV file listing the original filenames and corresponding organisms.
    """
    df = pd.read_csv(csv_file, sep=",")
    os.makedirs(output_folder, exist_ok=True)

    existing_files = set(os.listdir(source_folder))

    if "File" in df.columns:
        expected_files = set(df["File"].dropna())
        if expected_files.issubset(set(os.listdir(output_folder))):
            print("All files from 'File' column are already in output folder.")
            return

    file_mapping = dict(zip(df["Previous_file"], df["Organism"]))
    renamed_files = {}

    for filename in existing_files:
        old_path = os.path.join(source_folder, filename)

        if filename in file_mapping:
            new_name = file_mapping[filename].lower().replace(" ", "_") + ".fasta"
            new_path = os.path.join(output_folder, new_name)
            shutil.copyfile(old_path, new_path)
            renamed_files[filename] = new_name
            print(f"Copied and renamed: {filename} → {new_name}")
        else:
            print(f"Skipped (not in CSV): {filename}")

    df["File"] = df["Previous_file"].map(renamed_files).fillna(df.get("File", ""))
    df.to_csv(csv_file, index=False)


def add_particule_to_fasta(fasta_file, particle):
    """
    Add a prefix (particle) to each FASTA header if not already present.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.
    particle : str
        Prefix to add to each sequence identifier.
    """
    with open(fasta_file, "r") as f:
        lines = f.readlines()

    with open(fasta_file, "w") as f:
        for line in lines:
            if line.startswith(">"):
                header = line[1:].strip()
                if not header.startswith(f"{particle}_"):
                    line = f">{particle}_{header}\n"
            f.write(line)


def process_fasta_folder(folder_path, csv_file):
    """
    Apply cleaning and prefixing modifications to all FASTA files.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing FASTA files.
    csv_file : str
        Path to the CSV file listing proteomes.
    """
    df = pd.read_csv(csv_file, sep=",")
    fasta_particle_map = dict(zip(df["File"], df["Particule"]))

    for filename in os.listdir(folder_path):
        fasta_path = os.path.join(folder_path, filename)
        if not os.path.isfile(fasta_path):
            continue

        if filename in fasta_particle_map:
            clean_fasta_file(fasta_path)
            add_particule_to_fasta(fasta_path, fasta_particle_map[filename])


def run_fasta_translation(source_folder, output_folder, csv_file):
    """
    Run the full processing pipeline on a directory of FASTA files.

    Parameters
    ----------
    source_folder : str
        Path to the directory containing FASTA files.
    ouput_foler : str
        Path to thr directory output.
    csv_file : str
        Path to the CSV file listing proteomes and their mappings.
    """
    if not os.path.exists(source_folder):
        print(f"Error: source folder '{source_folder}' does not exist.")
        sys.exit(1)

    rename_fasta_files(source_folder, output_folder, csv_file)
    process_fasta_folder(output_folder, csv_file)
    process_fasta_files(output_folder)