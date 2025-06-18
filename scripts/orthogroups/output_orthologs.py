"""
Orthologs groups Analyze Pipeline

This script analyze results from orthologs groups methods.

Pipeline steps :
- Recover orthologs groups results from the 3 methods (SonicParanoid, Orthologer, OrthoFinder).
- Write orthologs groups on an output file.


Author(s): Romain DAGUERRE

"""

import os
import pandas as pd
import glob


def get_ortholog_groups_directory(file_path):
    """
    Retrieves the path to the 'ortholog_groups.tsv' file from the SonicParanoid log.

    Parameters
    ----------
    file_path : str
        Path to the SonicParanoid output directory.

    Returns
    -------
    str
        Full path to the 'ortholog_groups.tsv' file, or None if not found.
    """
    file_info = os.path.join(file_path, "last_run_info.txt")

    with open(file_info, "r") as f:
        for line in f:
            if line.startswith("Directory with ortholog groups:"):
                directory = line.strip().split(":", 1)[1].strip()
                return os.path.join(directory, "ortholog_groups.tsv")
    return None


def analyze_sonic_groups(file_path, output_dir, sep="\t"):
    """
    Parses SonicParanoid ortholog groups and outputs a formatted CSV file.

    Parameters
    ----------
    file_path : str
        Path to the SonicParanoid results directory.
    output_dir : str
        Directory to save the output CSV.
    sep : str, optional
        Field separator used in the input file (default: tab).
    """
    output_sonic = os.path.join(output_dir, "sonicParanoid.csv")
    orthogroups_file = get_ortholog_groups_directory(file_path)

    df = pd.read_csv(orthogroups_file, sep=sep)
    species_columns = df.columns[4:]
    formatted_data = []

    for _, row in df.iterrows():
        group_id = row["group_id"]
        group_data = {"Busco_Group": group_id}
        protein_count = 0

        for species in species_columns:
            value = row[species]
            if pd.isna(value) or value.strip() in ["", "*"]:
                group_data[species] = "Na"
            else:
                group_data[species] = value
                protein_count += 1

        if protein_count >= 2:
            formatted_data.append(group_data)

    df_output = pd.DataFrame(formatted_data)
    df_output.to_csv(output_sonic, sep=",", index=False)


def analyze_orthologer_groups(file_path, output_dir, sep=" "):
    """
    Parses Orthologer ortholog group output and formats it into a CSV file.

    Parameters
    ----------
    file_path : str
        Path to the Orthologer output file.
    output_dir : str
        Directory to save the output CSV.
    sep : str, optional
        Field separator used in the input file (default: space).
    """
    output_orthologer = os.path.join(output_dir, "orthologer.csv")
    df = pd.read_csv(file_path, sep=sep, header=None, comment="#")
    df.columns = [
        "group_id",
        "species",
        "e_value",
        "start",
        "end",
        "score",
        "value",
        "val",
        "p_value",
    ]

    groups = {}

    for _, row in df.iterrows():
        group_id = row["group_id"]
        species_code = row["species"].split("_")[0]

        if group_id not in groups:
            groups[group_id] = {"Busco_Group": group_id, "total_proteins": 0}

        if species_code in groups[group_id]:
            groups[group_id][species_code] += f",{row['species']}"
        else:
            groups[group_id][species_code] = row["species"]

        groups[group_id]["total_proteins"] += 1

    df_output = pd.DataFrame.from_dict(groups, orient="index")
    df_output = df_output[df_output["total_proteins"] >= 2]
    df_output = df_output.drop(columns=["total_proteins"])
    df_output = df_output.fillna("Na")
    df_output.to_csv(output_orthologer, sep=",", index=False)


def analyze_orthofinder_groups(file_path, output_dir, sep="\t"):
    """
    Parses OrthoFinder orthogroup TSV and creates a formatted CSV file.

    Parameters
    ----------
    file_path : str
        Path to the OrthoFinder Orthogroups.tsv file.
    output_dir : str
        Directory to save the output CSV.
    sep : str, optional
        Field separator used in the input file (default: tab).
    """
    output_orthofinder = os.path.join(output_dir, "orthofinder.csv")
    df = pd.read_csv(file_path, delimiter=sep)
    df = df.loc[:, ~df.columns.str.contains("OG", case=False, na=False)]
    species_columns = df.columns[1:]

    def count_proteins(row):
        return sum(
            len(proteins.split(",")) if isinstance(proteins, str) else 0
            for proteins in row[species_columns]
        )

    df["protein_count"] = df.apply(count_proteins, axis=1)
    df.insert(0, "Busco_Group", range(1, len(df) + 1))
    df.drop(columns=["protein_count"], inplace=True)
    df.to_csv(output_orthofinder, sep=",", index=False)


def run_analysis(output_sonic, output_orthologer, output_orthofinder, output_dir):
    """
    Runs all orthogroup analyses for SonicParanoid, Orthologer, and OrthoFinder.

    Parameters
    ----------
    output_sonic : str
        Path to the SonicParanoid results directory.
    output_orthologer : str
        Path to the Orthologer results directory.
    output_orthofinder : str
        Path to the OrthoFinder results directory.
    output_dir : str
        Output directory where all formatted CSVs will be stored.
    """
    output_orthogroups = os.path.join(output_dir, "orthogroups")
    os.makedirs(output_orthogroups, exist_ok=True)

    orthologer_file = os.path.join(output_orthologer, "Results/project_orthogroups.txt")
    orthofinder_matches = glob.glob(
        os.path.join(output_orthofinder, "*/Orthogroups/Orthogroups.tsv")
    )

    if orthofinder_matches:
        orthofinder_file = orthofinder_matches[0]
    else:
        raise FileNotFoundError(
            "Orthogroups.tsv not found in output_orthofinder subdirectories."
        )

    analyze_sonic_groups(output_sonic, output_orthogroups)
    analyze_orthologer_groups(orthologer_file, output_orthogroups)
    analyze_orthofinder_groups(orthofinder_file, output_orthogroups)
