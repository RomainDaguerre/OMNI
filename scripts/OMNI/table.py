"""
Protein Orthology Comparison Pipeline

This script serves as a command-line interface for analyzing orthologous protein groups
across species using precomputed FASTA and orthology inference data.

Runs the `run_table()` function which performs:
   - Merging and comparison of orthologous groups from multiple inference tools.
   - Dynamic selection of the most compact reference group for species comparison.
   - Evaluation of orthology consistency across tools.
   - Optional integration of tree-based analyses.

Dependencies:
    - Python libraries.

Author(s): Romain DAGUERRE
"""


import pandas as pd
import os
import csv
from Bio.Blast.Applications import NcbiblastpCommandline
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
from ete3 import Tree, TreeStyle, faces
from Bio import Entrez, SeqIO, pairwise2
from matplotlib.colors import ListedColormap
from collections import defaultdict
import re
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from matplotlib.patches import Wedge, Circle
from PIL import Image
from matplotlib.patches import Patch



def fetch_sequence_ncbi_single(protein_id, protein_name, output_fasta, email="romain.daguerre1119@gmail.com"):
    """
    Fetch a protein sequence from NCBI and save it in FASTA format.

    This function:
    - Uses the NCBI Entrez API to search for a protein using its ID.
    - Retrieves the corresponding protein sequence in FASTA format.
    - Cleans the FASTA to have a single-line sequence with a custom header.

    Parameters
    ----------
    protein_id : str
        Identifier of the protein to search in NCBI (e.g., UniProt or RefSeq ID).
    
    protein_name : str
        Descriptive name of the protein (used for print messages).
    
    output_fasta : str
        Path to the output FASTA file where the sequence will be saved.
    
    email : str, optional
        Email address required by NCBI to use Entrez services (default is provided).
    """
    Entrez.email = email  # Mandatory for using ENTREZ

    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

    try:
        print(f"Recherche de {protein_id} ({protein_name})...")
        search_handle = Entrez.esearch(db="protein", term=protein_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if search_results["IdList"]:
            protein_ncbi_id = search_results["IdList"][0]

            fetch_handle = Entrez.efetch(
                db="protein", id=protein_ncbi_id, rettype="fasta", retmode="text"
            )
            fasta = fetch_handle.read()
            fetch_handle.close()

            # Cleaning
            lines = fasta.strip().split('\n')
            header = f">{protein_id}"
            sequence = "".join(lines[1:])

            with open(output_fasta, "w") as f_out:
                f_out.write(header + "\n" + sequence + "\n")

        else:
            print(f"Aucune séquence trouvée pour {protein_id}")

    except Exception as e:
        print(f"Erreur lors de la récupération de {protein_id} : {e}")


def read_protein_interest(file_path):
    """
    Read a CSV file containing proteins of interest and return a nested dictionary.

    This function:
    - Reads a CSV file where each row contains a species name, a protein code, and a protein name.
    - Organizes the data into a nested dictionary of the form {species: {protein_code: protein_name}}.

    Parameters
    ----------
    file_path : str
        Path to the CSV file containing protein information. The file must contain the columns:
        'espece', 'code_proteine', and 'nom_proteine'.

    Returns
    -------
    dict
        A nested dictionary mapping species to protein codes and corresponding names.
    """
    proteins = {}
    with open(file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            espece = row['espece']
            code_proteine = row['code_proteine']
            nom_proteine = row['nom_proteine']
            
            if espece not in proteins:
                proteins[espece] = {}
            proteins[espece][code_proteine] = nom_proteine

    return proteins


def load_groups(orthogroup_file):
    """
    Load orthogroup data and return a dictionary mapping group IDs to sets of gene identifiers.

    This function:
    - Reads a CSV file where the first column contains BUSCO group IDs and the remaining columns contain gene lists per species.
    - Returns a dictionary with group IDs as keys and sets of associated genes as values.

    Parameters
    ----------
    orthogroup_file : str
        Path to the CSV file containing orthogroup data.

    Returns
    -------
    dict
        A dictionary mapping group IDs to sets of gene identifiers (set of str).
    """
    df = pd.read_csv(orthogroup_file, sep=',')
    orthogroup_dict = {}
    for _, row in df.iterrows():
        busco_group = row['Busco_Group']
        genes = set()        
        for col in df.columns[1:]:
            if pd.notna(row[col]) and row[col] != 'Missing':
                gene_list = row[col].replace(',', ' ').split()
                genes.update(gene_list)
        if genes:
            orthogroup_dict[busco_group] = genes
    return orthogroup_dict


def run_blastp(query_fasta, db_fasta, output_file):
    """
    Run BLASTP to compare protein sequences against a protein database.

    This function:
    - Creates a BLAST protein database from the provided FASTA file if it doesn't already exist.
    - Executes a BLASTP search using the query sequences.
    - Outputs the BLAST results in tabular format (outfmt 6).

    Parameters
    ----------
    query_fasta : str
        Path to the FASTA file containing the query protein sequences.

    db_fasta : str
        Path to the FASTA file to be used as the BLAST protein database.

    output_file : str
        Path to the output file where BLAST results will be written.
    """
    db_name = db_fasta.replace(".fasta", "")
    
    if not os.path.exists(db_name + ".pin"):
        os.system(f"makeblastdb -in {db_fasta} -dbtype prot -out {db_name}")
    
    blastp_cline = NcbiblastpCommandline(query=query_fasta, db=db_name, outfmt=6, out=output_file)
    stdout, stderr = blastp_cline()


def parse_blast_results(blast_output):
    """
    Parse a BLAST output file and extract the best hit per query sequence.

    This function:
    - Reads a BLAST tabular output file (outfmt 6).
    - For each query sequence, identifies the hit with the highest bit score.
    - Returns a dictionary mapping each query ID to a tuple of (subject ID, bit score).

    Parameters
    ----------
    blast_output : str
        Path to the BLAST output file in tabular format.

    Returns
    -------
    dict
        A dictionary where keys are query sequence IDs and values are tuples
        (best subject ID, best bit score).
    """
    best_hits = {}

    with open(blast_output, 'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            query_id, subject_id = cols[0], cols[1]
            evalue, bit_score = map(float, cols[10:12])

            # Store the best hit (highest bit score)
            if query_id not in best_hits or bit_score > best_hits[query_id][1]:
                best_hits[query_id] = (subject_id, bit_score)

    return best_hits


def parse_hits_in_orthogroups(ortholog_dict, results):
    """
    Search for the best hit of each protein within ortholog groups.

    Parameters
    ----------
    ortholog_dict : dict
        Dictionary mapping group IDs to lists of protein IDs.
    results : dict
        BLAST results in the form {protein_name: {"best_hits": (hit, score)}}.

    Returns
    -------
    dict
        Dictionary mapping each protein name to a dictionary containing:
        - "best_hit": best hit protein ID,
        - "score": bit score of the best hit,
        - "ortholog_group": ID of the ortholog group where the best hit was found,
          or None if not found.
    """
    ortholog_groups = {}

    for protein_name, data in results.items():
        best_hit, score = data["best_hits"]
        found_group = False

        for group_id, members in ortholog_dict.items():
            if best_hit in members:
                ortholog_groups[protein_name] = {
                    "best_hit": best_hit,
                    "score": score,
                    "ortholog_group": group_id
                }
                found_group = True
                break

        if not found_group:
            ortholog_groups[protein_name] = {
                "best_hit": best_hit,
                "score": score,
                "ortholog_group": None
            }

    return ortholog_groups


def retrieve_orthologous_proteins(ortholog_dict, ortholog_groups):
    """
    Retrieve the proteins belonging to each ortholog group in ortholog_dict.

    Parameters
    ----------
    ortholog_dict : dict
        Dictionary containing ortholog groups with their associated proteins.
    ortholog_groups : dict
        Dictionary from previous results containing the ortholog group assigned to each best hit.

    Returns
    -------
    dict
        A dictionary mapping each ortholog group to the proteins that belong to it.
    """
    group_proteins = {}

    for query_id, data in ortholog_groups.items():
        group_id = data["ortholog_group"]

        if group_id is not None:
            if group_id in ortholog_dict:
                group_proteins[group_id] = ortholog_dict[group_id]
            else:
                # If the group does not exist in ortholog_dict, assign an empty list
                group_proteins[group_id] = []

    return group_proteins


def merge_orthologous_groups_by_protein(method1_orthologs, method1_groups,
                                        method2_orthologs, method2_groups):
    """
    Merges orthologous groups for each query protein by combining results from two orthology methods.

    Parameters:
    - method1_orthologs : dict : Orthologous group assignments per query protein from method 1
    - method1_groups : dict : Proteins associated with each orthologous group from method 1
    - method2_orthologs : dict : Orthologous group assignments per query protein from method 2
    - method2_groups : dict : Proteins associated with each orthologous group from method 2

    Returns:
    - dict : Merged dictionary where each query protein is associated with the unique set of orthologs from both methods
    """
    merged_proteins = {}

    all_query_ids = set(method1_orthologs.keys()).union(set(method2_orthologs.keys()))

    for query_id in all_query_ids:
        merged_proteins[query_id] = set()

        # Add proteins from method 1
        if query_id in method1_orthologs:
            group_id = method1_orthologs[query_id]["ortholog_group"]
            if group_id in method1_groups:
                merged_proteins[query_id].update(method1_groups[group_id])

        # Add proteins from method 2
        if query_id in method2_orthologs:
            group_id = method2_orthologs[query_id]["ortholog_group"]
            if group_id in method2_groups:
                merged_proteins[query_id].update(method2_groups[group_id])

    return {query_id: list(proteins) for query_id, proteins in merged_proteins.items()}


def remove_duplicate_groups(merged_proteins_by_query):
    """
    Removes redundant orthologous groups that contain exactly the same proteins.

    Parameters:
    - merged_proteins_by_query : dict : Dictionary where each query protein is associated with a set of orthologous proteins

    Returns:
    - dict : Dictionary without duplicated groups
    """
    unique_groups = {}
    query_to_group = {}

    for query_id, proteins in merged_proteins_by_query.items():
        proteins_set = frozenset(proteins)

        if proteins_set not in unique_groups:
            unique_groups[proteins_set] = query_id

        query_to_group[query_id] = unique_groups[proteins_set]

    final_groups = {query_id: list(proteins) for proteins, query_id in unique_groups.items()}

    return final_groups


def load_species_mapping(species_file):
    """
    Loads a CSV file containing species information and creates a dictionary mapping particles to organism names.

    The CSV file must contain the following columns: 'Organism', 'File', 'Particule', and 'BUSCO_score'.

    Parameters
    ----------
    species_file : str
        Path to the CSV file with species information.

    Returns
    -------
    dict
        Dictionary mapping 'Particule' values to corresponding 'Organism' names.
    """
    df = pd.read_csv(species_file)
    species_dict = pd.Series(df['Organism'].values, index=df['Particule']).to_dict()

    return species_dict


def analyze_orthologous_groups(final_merged_groups, method_dict):
    """
    Analyzes orthologous groups from a specified method by listing the proteins 
    associated with each query protein and showing their composition, 
    including those not found in any group.

    Parameters
    ----------
    final_merged_groups : dict
        Dictionary of merged groups with associated proteins for each query.

    method_dict : dict
        Dictionary of orthologous groups from a method (e.g., SonicParanoid, Orthologer, OrthoFinder),
        where keys are group IDs and values are lists of proteins.

    Returns
    -------
    tuple
        - dict: Detailed information for each group found, including:
            - 'queries': list of associated query proteins
            - 'count': number of matching proteins found in the group
            - 'proteins': full list of proteins in the group
        - dict: Proteins that were not found in any orthologous group, grouped by query protein.
    """
    group_details = {}
    proteins_in_no_group = defaultdict(set)

    for query_id, proteins in final_merged_groups.items():
        for protein in proteins:
            found_in_group = False
            for group_id, method_proteins in method_dict.items():
                if protein in method_proteins:
                    if group_id not in group_details:
                        group_details[group_id] = {
                            "queries": set(),
                            "count": 0,
                            "proteins": set()
                        }

                    group_details[group_id]["queries"].add(query_id)
                    group_details[group_id]["proteins"].update(method_proteins)
                    group_details[group_id]["count"] += 1
                    found_in_group = True
                    break

            # Add proteins that aren't find in a group
            if not found_in_group:
                proteins_in_no_group[query_id].add(protein)

    for group_id in group_details:
        group_details[group_id]["queries"] = list(group_details[group_id]["queries"])
        group_details[group_id]["proteins"] = list(group_details[group_id]["proteins"])

    proteins_in_no_group = {k: list(v) for k, v in proteins_in_no_group.items()}

    return group_details, proteins_in_no_group


def extract_species(protein_id):
    """
    Extracts the species tag from a protein ID (e.g., Pfalc_, TGME49_, etc.).

    Parameters
    ----------
    protein_id : str
        The protein identifier from which to extract the species prefix.

    Returns
    -------
    str or None
        The extracted species tag (prefix before the underscore) if found, otherwise None.
    """
    match = re.match(r"^([A-Za-z0-9]+)_", protein_id)
    return match.group(1) if match else None


def get_sequence_from_fasta(fasta_path, protein_id):
    """
    Retrieves the amino acid sequence of a given protein ID from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file containing protein sequences.
    protein_id : str
        Identifier of the protein whose sequence is to be retrieved.

    Returns
    -------
    str or None
        The amino acid sequence of the protein if found, otherwise None.
    """
    for record in SeqIO.parse(fasta_path, "fasta"):
        if protein_id in record.id:
            return str(record.seq)
    return None


def calculate_similarity(seq1, seq2):
    """
    Calculates a simple similarity score between two protein sequences using global alignment.

    Parameters
    ----------
    seq1 : str
        First protein sequence.
    seq2 : str
        Second protein sequence.

    Returns
    -------
    float
        Normalized similarity score between 0 and 1, based on the number of identical matches
        in the best global alignment. Returns 0 if either sequence is empty or None.
    """
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    return alignments / max(len(seq1), len(seq2)) if seq1 and seq2 else 0


def compute_group_similarity_score(group_proteins, query_seq, fasta_dir_species, species_df):
    """
    Computes the average similarity score between a query sequence and a set of orthologous proteins.

    Parameters
    ----------
    group_proteins : list
        List of protein IDs from an orthologous group.
    query_seq : str
        Amino acid sequence of the query protein.
    fasta_dir_species : str
        Path to the directory containing species-specific FASTA files.
    species_df : pandas.DataFrame
        DataFrame containing species metadata with at least 'Particule' and 'File' columns.

    Returns
    -------
    float
        Average similarity score between the query sequence and all valid group proteins.
        Returns 0 if no sequences could be compared.
    """
    scores = []
    for prot in group_proteins:
        species_code = extract_species(prot)
        species_row = species_df[species_df["Particule"] == species_code]
        if not species_row.empty:
            fasta_file = os.path.join(fasta_dir_species, species_row.iloc[0]["File"])
            seq = get_sequence_from_fasta(fasta_file, prot)
            if seq:
                sim = calculate_similarity(query_seq, seq)
                scores.append(sim)
    return sum(scores) / len(scores) if scores else 0


def get_seq_fasta(fasta_dir, protein_id):
    """
    Searches for a protein sequence in all FASTA files within the specified directory.

    Parameters
    ----------
    fasta_dir : str
        Path to the directory containing FASTA files.
    protein_id : str
        ID of the protein to search for.

    Returns
    -------
    str or None
        The protein sequence if found, otherwise None.
    """
    for file_name in os.listdir(fasta_dir):
        if not file_name.endswith(".fasta"):
            continue
        fasta_path = os.path.join(fasta_dir, file_name)
        for record in SeqIO.parse(fasta_path, "fasta"):
            if protein_id in record.id:
                return str(record.seq)
    return None


def build_tree(sequences, tree_dir):
    """
    Builds a phylogenetic tree from a dictionary of protein sequences using MAFFT for alignment
    and FastTree for tree inference.

    Parameters
    ----------
    sequences : dict
        Dictionary where keys are protein IDs and values are sequences (strings).
    tree_dir : str
        Path to the directory where intermediate and output files will be saved.

    Returns
    -------
    str
        Path to the resulting Newick tree file.
    """
    os.makedirs(tree_dir, exist_ok=True)
    raw_fasta = os.path.join(tree_dir, "all_seqs_raw.fasta")
    aligned_fasta = os.path.join(tree_dir, "aligned_seqs.fasta")
    tree_path = os.path.join(tree_dir, "tree.nwk")

    seq_records = [
        SeqRecord(Seq(seq), id=f"{id}", description="")
        for id, seq in sequences.items()
        if seq is not None and seq != ""
    ]
    SeqIO.write(seq_records, raw_fasta, "fasta")

    # Step 1: Align sequences using MAFFT
    mafft_cmd = f"mafft --auto {raw_fasta} > {aligned_fasta}"
    subprocess.run(mafft_cmd, shell=True, check=True)

    # Step 2: Build phylogenetic tree using FastTree
    fasttree_cmd = f"fasttree {aligned_fasta} > {tree_path}"
    subprocess.run(fasttree_cmd, shell=True, check=True)

    return tree_path


def compute_z_scores_and_ratios(long_branch_scores, orphans, df):
    """
    Computes Z-scores for protein branch lengths and identifies the lowest Z-score
    per species.

    Parameters
    ----------
    long_branch_scores : dict
        Dictionary mapping protein IDs to their long branch scores.
    orphans : set
        Set of orphan protein IDs to exclude from Z-score reference computation.
    df : pandas.DataFrame
        DataFrame containing species metadata with at least the columns 'Organism' and 'Particule'.

    Returns
    -------
    z_scores : dict
        Dictionary mapping each protein ID to its Z-score.
    min_z_by_species : dict
        Dictionary mapping each species name to the minimum Z-score observed among its proteins.
    """
    # Compute mean and standard deviation using only non-orphan proteins
    non_orphan_scores = [v for k, v in long_branch_scores.items() if k not in orphans]
    mean_non_orphans = np.mean(non_orphan_scores)
    std_non_orphans = np.std(non_orphan_scores)

    # Compute Z-scores for all proteins
    z_scores = {
        k: (v - mean_non_orphans) / std_non_orphans
        for k, v in long_branch_scores.items()
    }

    # Determine minimum Z-score per species
    min_z_by_species = {}
    for _, row in df.iterrows():
        organism = row["Organism"]
        prefix = row["Particule"]
        matched = {k: z for k, z in z_scores.items() if k.startswith(prefix)}

        if matched:
            min_prot = min(matched, key=matched.get)
            min_z_by_species[organism] = matched[min_prot]

    return z_scores, min_z_by_species


def load_species_and_protein_codes(species_file, prot_interest_file):
    """
    Load species data and map protein names to protein codes.

    Parameters
    ----------
    species_file : str
        Path to the CSV file containing species information.
    prot_interest_file : str
        Path to the CSV file mapping protein names to protein codes.

    Returns
    -------
    tuple
        species_df : pandas.DataFrame
            DataFrame with species information.
        query_to_code : dict
            Dictionary mapping protein names (str) to protein codes (str).
    """
    species_df = pd.read_csv(species_file)
    query_to_code = pd.read_csv(prot_interest_file, sep=",").set_index("nom_proteine")["code_proteine"].to_dict()
    return species_df, query_to_code


def gather_common_species(reference_groups, analyzed_groups):
    """
    Extract common species present in analyzed groups keyed by query protein.

    Parameters
    ----------
    reference_groups : iterable
        Iterable of reference group IDs to analyze.
    analyzed_groups : dict
        Dictionary with group_id keys and group information values.

    Returns
    -------
    defaultdict
        Nested defaultdict mapping query protein -> group_id -> set of species.
    """
    common_species = defaultdict(lambda: defaultdict(set))
    for group_id in reference_groups:
        group_info = analyzed_groups[group_id]
        query = group_info["queries"][0]
        group_species = set(extract_species(p) for p in group_info["proteins"] if extract_species(p))
        common_species[query][group_id] = group_species
    return common_species


def sort_groups(matching_group_ids, analyzed_groups, query_seq, fasta_dir_species, species_df, sort_method):
    """
    Sort matching groups by either species count or sequence similarity.

    Parameters
    ----------
    matching_group_ids : list
        List of group IDs to sort.
    analyzed_groups : dict
        Dictionary with group data including proteins.
    query_seq : str
        The sequence of the query protein.
    fasta_dir_species : str
        Directory containing species fasta sequences.
    species_df : pandas.DataFrame
        DataFrame containing species info.
    sort_method : str
        Sorting method, either "species_count" or "similarity".

    Returns
    -------
    list
        Sorted list of group IDs according to the chosen method.
    """
    if sort_method == "species_count":
        print("Tri par nombre d'espèces")
        sorted_groups = sorted(
            matching_group_ids,
            key=lambda gid: len(analyzed_groups[gid]["proteins"]),
            reverse=True
        )
    elif sort_method == "similarity":
        print("Tri par similarité de séquence")
        similarity_scores = {
            gid: compute_group_similarity_score(analyzed_groups[gid]["proteins"], query_seq, fasta_dir_species, species_df)
            for gid in matching_group_ids
        }
        sorted_groups = sorted(matching_group_ids, key=lambda gid: similarity_scores[gid], reverse=True)
    else:
        sorted_groups = matching_group_ids
    return sorted_groups


def merge_groups_for_query(query, ref_groups, analyzed_groups, sorted_matching_groups, fusion_count_start=0):
    """
    Merge groups that do not share any common species for each query protein.

    Parameters
    ----------
    common_species : defaultdict
        Nested dictionary query -> group_id -> species set.
    analyzed_groups : dict
        Dictionary with group info including proteins and queries.
    sort_method : str
        Method to sort matching groups ('species_count' or 'similarity').
    fasta_dir_query : str
        Directory containing query protein fasta sequences.
    fasta_dir_species : str
        Directory containing species fasta sequences.
    species_df : pandas.DataFrame
        DataFrame with species info.
    query_to_code : dict
        Mapping from query protein names to protein codes.

    Returns
    -------
    tuple
        merged_no_common : dict
            Dictionary with merged groups keyed by query and merged group IDs.
        fusion_count : int
            Number of fusions performed.
    """
    merged_groups = defaultdict(dict)
    fusion_count = fusion_count_start
    already_merged = set()

    for ref_group_id, ref_species_set in ref_groups.items():
        print(f"Groupe de référence {ref_group_id} ({len(ref_species_set)} espèces) :")

        current_group_id = ref_group_id
        current_proteins_by_group = defaultdict(list)
        current_species = ref_species_set.copy()

        for prot in analyzed_groups[ref_group_id]["proteins"]:
            current_proteins_by_group[ref_group_id].append(prot)

        for test_group_id in sorted_matching_groups:
            if test_group_id == ref_group_id or test_group_id in already_merged:
                continue

            test_info = analyzed_groups[test_group_id]
            test_species_set = set(extract_species(p) for p in test_info["proteins"] if extract_species(p))
            shared_species = current_species & test_species_set

            print(f"    ↪ Comparé avec groupe {test_group_id} ({len(test_species_set)} espèces) :")
            if shared_species:
                print(f"Espèces communes : {shared_species} → Pas de fusion")
                continue
            else:
                print(f"Aucune espèce commune → Fusion")
                for prot in test_info["proteins"]:
                    current_proteins_by_group[test_group_id].append(prot)
                current_species |= test_species_set
                already_merged.add(test_group_id)
                current_group_id = f"fusion_{fusion_count}"

        if isinstance(current_group_id, str) and current_group_id.startswith("fusion_"):
            fusion_count += 1

        merged_groups[query][current_group_id] = dict(current_proteins_by_group)

    return merged_groups, fusion_count


def filter_orphans(proteins_in_no_group, merged_no_common):
    """
    Identify orphans (proteins not in any group) whose species are not already represented.

    Parameters
    ----------
    merged_no_common : dict
        Dictionary of merged groups keyed by query and group IDs.
    proteins_in_no_group : dict
        Dictionary mapping query proteins to lists of orphan proteins.

    Returns
    -------
    list
        List of tuples (query, orphan_protein) for orphans with unique species.
    """
    filtered_orphans = []
    for query in proteins_in_no_group:
        if query not in merged_no_common:
            continue
        for orphan in proteins_in_no_group[query]:
            orphan_species = extract_species(orphan)
            existing_species = {
                extract_species(p)
                for group in merged_no_common[query].values()
                for plist in group.values()
                for p in plist
            }
            if orphan_species not in existing_species:
                filtered_orphans.append((query, orphan))
    print(f"Orphelins : {filtered_orphans}")
    return filtered_orphans


def process_orphans_and_build_tree(merged_no_common, proteins_in_no_group, fasta_dir_species, species_df, output_dir):
    """
    For each query, combine group protein sequences and orphans sequences,
    build a phylogenetic tree, calculate long branch scores and z-scores,
    and update merged_no_common with orphan info.

    Parameters
    ----------
    merged_no_common : dict
        Merged groups keyed by query and group IDs.
    proteins_in_no_group : dict
        Dictionary mapping queries to orphan proteins.
    fasta_dir_species : str
        Directory containing species fasta sequences.
    output_dir : str
        Directory to store output trees.
    species_df : pandas.DataFrame
        DataFrame with species information.

    Returns
    -------
    tuple
        merged_no_common : dict
            Updated dictionary including orphan information.
        min_by_species : dict
            Computed minimum z-scores or ratios by species.
    """
    min_by_species = None
    for query in merged_no_common:
        if query not in proteins_in_no_group:
            continue

        group_seqs_dict = {
            p: get_seq_fasta(fasta_dir_species, p)
            for group in merged_no_common[query].values()
            for plist in group.values()
            for p in plist
        }

        orphans = [
            orphan for orphan in proteins_in_no_group[query]
            if extract_species(orphan) not in {
                extract_species(p)
                for group in merged_no_common[query].values()
                for plist in group.values()
                for p in plist
            }
        ]

        orphan_seqs_dict = {
            orphan: get_seq_fasta(fasta_dir_species, orphan)
            for orphan in orphans
        }

        combined_seqs = {**group_seqs_dict, **orphan_seqs_dict}

        output_tree = os.path.join(output_dir, "tree")
        tree = build_tree(combined_seqs, output_tree)

        phykit_cmd = f"phykit long_branch_score {tree} -v"
        result = subprocess.run(phykit_cmd, shell=True, capture_output=True, text=True, check=True)

        long_branch_scores = {}
        for line in result.stdout.strip().split('\n'):
            parts = line.split()
            if len(parts) == 2:
                name, score = parts
                long_branch_scores[name] = float(score)

        z_scores, min_by_species = compute_z_scores_and_ratios(long_branch_scores, orphans, species_df)

        # Adding orphans
        all_orphan_ids = list(orphan_seqs_dict.keys())
        for orphan in all_orphan_ids:
            last_group_id = max(merged_no_common[query].keys(), key=lambda x: x.startswith("fusion_"))
            merged_no_common[query][last_group_id].setdefault("orphans", []).append(orphan)

    return merged_no_common, min_by_species


def compare_species_between_groups(reference_groups, analyzed_groups, proteins_in_no_group, species_file, fasta_dir_query, fasta_dir_species, prot_interest_file, output_dir, sort_method="species_count"):
    """
    Main function to merge groups with no common species and handle orphans.

    This function orchestrates the process by:
    - Loading species and protein code data.
    - Extracting common species per group.
    - Sorting and merging groups without shared species.
    - Filtering orphans that represent unique species.
    - Building phylogenetic trees with combined sequences.
    - Computing long branch scores and z-scores for orphan evaluation.

    Parameters
    ----------
    reference_groups : iterable
        Reference group IDs to consider.
    analyzed_groups : dict
        Dictionary of group information.
    proteins_in_no_group : dict
        Dictionary mapping query proteins to orphan proteins.
    species_file : str
        Path to species CSV file.
    fasta_dir_query : str
        Directory of query protein fasta files.
    fasta_dir_species : str
        Directory of species fasta files.
    prot_interest_file : str
        CSV file mapping protein names to codes.
    output_dir : str
        Directory where output files and trees will be saved.
    sort_method : str, optional
        Sorting method for group merging ('species_count' or 'similarity'), by default 'species_count'.

    Returns
    -------
    tuple
        merged_no_common : dict
            Dictionary of merged groups including orphans.
        min_by_species : dict
            Minimum z-scores or ratios computed per species for orphans.
    """
    species_df, query_to_code = load_species_and_protein_codes(species_file, prot_interest_file)
    common_species = gather_common_species(reference_groups, analyzed_groups)
    merged_no_common = defaultdict(dict)
    fusion_count = 0

    for query, ref_groups in common_species.items():
        print(f"\n▶ Protéine requête : {query}")

        protein_code = query_to_code.get(query)
        print(protein_code)
        if not protein_code:
            print(f"  ⚠️  Aucun code_proteine trouvé pour {query}")
            continue

        query_seq = get_sequence_from_fasta(fasta_dir_query, protein_code)

        matching_group_ids = [
            gid for gid, info in analyzed_groups.items()
            if "queries" in info and query in info["queries"]
        ]

        sorted_matching_groups = sort_groups(matching_group_ids, analyzed_groups, query_seq, fasta_dir_species, species_df, sort_method)

        merged_query_groups, fusion_count = merge_groups_for_query(query, ref_groups, analyzed_groups, sorted_matching_groups, fusion_count)
        merged_no_common.update(merged_query_groups)

    filter_orphans(proteins_in_no_group, merged_no_common)
    merged_no_common, min_by_species = process_orphans_and_build_tree(merged_no_common, proteins_in_no_group, fasta_dir_species, species_df, output_dir)

    return merged_no_common, min_by_species


def write_final_groups_with_sequences(nom_prot, new_groups, output_results, fasta_dir):
    """
    Writes the final ortholog groups along with their FASTA sequences to a file.

    Parameters
    ----------
    nom_prot : str
        Name of the reference protein or query.
    new_groups : dict
        Dictionary of grouped proteins organized by group, cluster, and subcluster.
    output_results : str
        Path to the main output directory.
    fasta_dir : str
        Path to the directory containing FASTA files for each species.

    Output
    ------
    A FASTA file is created in the 'final_groups' subdirectory containing all sequences 
    for the given groups, formatted with subcluster and protein identifiers.
    """
    dir_fgroups = os.path.join(output_results, "final_groups")
    os.makedirs(dir_fgroups, exist_ok=True)

    output_file = os.path.join(dir_fgroups, f"{nom_prot}_orthogroups.fasta")

    with open(output_file, "w") as f:
        for group, id_dict in new_groups.items():
            f.write(f">{group}\n")
            for cluster_id, subcluster in id_dict.items():
                for sub_id, seq_list in subcluster.items():
                    for protein_id in seq_list:
                        sequence = get_seq_fasta(fasta_dir, protein_id)
                        if sequence:
                            f.write(f">{sub_id}|{protein_id}\n{sequence}\n\n")
                        else:
                            print(f"Sequence not found for {protein_id}")


def load_protein_mapping(csv_file):
    """
    Loads a CSV file containing the mapping between protein codes and their names.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file. The file must contain at least two columns:
        'code_proteine' and 'nom_proteine'.

    Returns
    -------
    dict
        A dictionary mapping each protein code to its corresponding protein name.
    """
    mapping = {}
    with open(csv_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            mapping[row["code_proteine"]] = row["nom_proteine"]
    return mapping


def get_busco_scores(csv_file):
    """
    Loads a CSV file containing BUSCO scores and returns a dictionary
    mapping each organism to its corresponding BUSCO score.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file. The file must contain the columns:
        'Organism' and 'BUSCO_score'.

    Returns
    -------
    dict
        A dictionary where keys are organism names and values are BUSCO scores.
    """
    df = pd.read_csv(csv_file)
    return dict(zip(df['Organism'], df['BUSCO_score']))


def generate_presence_absence_table_single_group(orthologs_group, species_dict):
    """
    Generates a presence/absence matrix with ortholog group identifiers.

    Each cell contains:
    - 0 if the species is not present in any group
    - 1, 2, 3... indicating the group index the species belongs to
    - n+1 for the 'no_group' (orphans)

    Parameters
    ----------
    orthologs_group : dict
        Nested dictionary of ortholog groups structured as:
        {query_id: {fusion_id: {group_id: list of proteins}}}
    species_dict : dict
        Mapping from species tag (e.g., 'Pfalc') to species full name 
        (e.g., 'Plasmodium falciparum').

    Returns
    -------
    pd.DataFrame
        A DataFrame where rows correspond to query proteins and columns
        indicate the group assignment per species.
    """
    all_species = sorted(set(species_dict.values()))
    data = []

    for query_id, fusions in orthologs_group.items():
        merged_groups = {}
        for fusion_id, group_dict in fusions.items():
            for group_id, proteins in group_dict.items():
                merged_groups[group_id] = proteins

        species_to_group = {}
        group_index = 1

        for group_id, proteins in merged_groups.items():
            for prot in proteins:
                species = species_dict.get(prot[:5])
                if species and species not in species_to_group:
                    if group_id == "orphans":
                        species_to_group[species] = -1  # Temporarily mark as orphan
                    else:
                        species_to_group[species] = group_index
            if group_id != "orphans":
                group_index += 1

        # Assign orphan group to the next group index
        for species, value in species_to_group.items():
            if value == -1:
                species_to_group[species] = group_index

        row = [f"{query_id}_fusion"]
        for species in all_species:
            row.append(species_to_group.get(species, 0))

        data.append(row)

    columns = ['Query Protein'] + all_species
    df = pd.DataFrame(data, columns=columns)
    return df


def merge_ortholog_groups(sonic_ortholog, orthofinder_ortholog, orthologer_ortholog):
    """
    Merges ortholog group data from three different sources by associating 
    each protein ID with its ortholog group from each source.

    Parameters
    ----------
    sonic_ortholog : dict
        Dictionary of SonicParanoid ortholog groupings: {protein_id: {"ortholog_group": group_data}}.
    orthofinder_ortholog : dict
        Dictionary of OrthoFinder ortholog groupings.
    orthologer_ortholog : dict
        Dictionary of Orthologer ortholog groupings.

    Returns
    -------
    dict
        A merged dictionary where each protein ID maps to a sub-dictionary of 
        source-to-group associations:
        {
            protein_id: {
                "sonic": group_data,
                "orthofinder": group_data,
                "orthologer": group_data
            }
        }
    """
    merged_dict = {}

    for ortho_type, ortholog_data in {
        "sonic": sonic_ortholog,
        "orthofinder": orthofinder_ortholog,
        "orthologer": orthologer_ortholog
    }.items():
        for prot_id, data in ortholog_data.items():
            if prot_id not in merged_dict:
                merged_dict[prot_id] = {}
            merged_dict[prot_id][ortho_type] = data["ortholog_group"]

    return merged_dict


def rename_keys(orthogroups, mapping):
    """
    Rename the keys of a dictionary using a provided mapping.

    Parameters
    ----------
    orthogroups : dict
        Dictionary with original keys to rename.
    mapping : dict
        Dictionary mapping old keys to new keys.

    Returns
    -------
    dict
        New dictionary with keys renamed according to the mapping.
        Keys not found in the mapping remain unchanged.
    """
    return {mapping.get(k, k): v for k, v in orthogroups.items()}


def generate_presence_absence_table(orthologs_sonic_group, orthologs_orthologer_group, orthologs_orthofinder_group, species_dict, merged_ortholog_groups):
    """
    Generate a presence/absence table with original query_ids first, then renamed at the end.

    Parameters
    ----------
    orthologs_sonic_group : dict
        Dictionary of Sonic ortholog groups by query_id.
    orthologs_orthologer_group : dict
        Dictionary of Orthologer groups by query_id.
    orthologs_orthofinder_group : dict
        Dictionary of OrthoFinder groups by query_id.
    species_dict : dict
        Dictionary mapping species prefixes to full species names.
    merged_ortholog_groups : dict
        Dictionary mapping query_id to its group IDs for each method.

    Returns
    -------
    pd.DataFrame
        DataFrame showing presence (1) or absence (0) of species per query_id,
        with query_ids renamed by method suffix at the end.
    """
    all_species = set(species_dict.values())
    query_sequences = set(orthologs_sonic_group.keys()).union(set(orthologs_orthologer_group.keys()), set(orthologs_orthofinder_group.keys()))

    data = []
    
    for query_id in query_sequences:
        species_in_sonic = {species_dict.get(prot[:5], None) for prot in orthologs_sonic_group.get(query_id, [])}
        species_in_orthologer = {species_dict.get(prot[:5], None) for prot in orthologs_orthologer_group.get(query_id, [])}
        species_in_orthofinder = {species_dict.get(prot[:5], None) for prot in orthologs_orthofinder_group.get(query_id, [])}
        
        species_in_sonic.discard(None)
        species_in_orthologer.discard(None)
        species_in_orthofinder.discard(None)
        
        row = [query_id]
        
        for species in all_species:
            in_combined = 1 if (species in species_in_sonic or species in species_in_orthologer or species in species_in_orthofinder) else 0
            row.append(in_combined)

        data.append(row)

    columns = ['Query Protein'] + list(all_species)
    df = pd.DataFrame(data, columns=columns)

    renamed_queries = {}
    for query_id, methods in merged_ortholog_groups.items():
        for method, group_id in methods.items():
            new_query_name = f"{query_id}_{method}"
            renamed_queries[group_id] = new_query_name

    df['Query Protein'] = df['Query Protein'].map(renamed_queries).fillna(df['Query Protein'])
    return df

def sort_key(query_protein):
    """
    Generate a sorting key for query protein identifiers.

    This function:
    - Splits a query protein string by underscores.
    - Uses the first two segments as the protein identifier.
    - Uses the third segment (method suffix) if present.
    - Returns a tuple (protein_id, method) for consistent sorting.

    Parameters
    ----------
    query_protein : str
        The query protein identifier in the format 'ProteinID_Method' or 'ProteinID_SubID_Method'.

    Returns
    -------
    tuple
        A tuple (protein_id, method) where:
        - protein_id is the combined first two segments (e.g., 'Prot_XYZ').
        - method is the third segment if it exists, otherwise an empty string.
    """
    parts = query_protein.split("_")
    protein_id = parts[0] + "_" + parts[1]
    method = parts[2] if len(parts) > 2 else ""
    return (protein_id, method)

'''
def run_table(fasta_dir, output_dir, tree_file, prot_interest_file, species_file, fasta_protein):

    output_orthogroup_file = os.path.join(output_dir, "orthogroups")
    sonic_file_name = "sonicParanoid.csv"
    output_sonic_file = os.path.join(output_orthogroup_file, sonic_file_name)
    orthologer_file_name = "orthologer.csv"
    output_orthologer_file = os.path.join(output_orthogroup_file, orthologer_file_name)
    orthofinder_file_name = "orthofinder.csv"
    output_orthofinder_file = os.path.join(output_orthogroup_file, orthofinder_file_name)

    tree = Tree(tree_file, format=1) 

    proteins_dict = read_protein_interest(prot_interest_file)

    sonic_dict = load_groups(output_sonic_file)
    orthologer_dict = load_groups(output_orthologer_file)
    orthofinder_dict = load_groups(output_orthofinder_file)

    species_dict = load_species_mapping(species_file)
        
    dict_esp_fusion = {}

    for espece, proteins in proteins_dict.items():
        db_fasta = os.path.join(fasta_dir, f"{espece}.fasta")
        
        for protein, nom_prot in proteins.items():

            output_results = os.path.join(output_dir, nom_prot)
            blast_output_dir = os.path.join(output_results, "blast_results")
            os.makedirs(blast_output_dir, exist_ok=True)
            
            if fasta_protein == None :
                fasta_interest_dir = os.path.join(output_results, f"{nom_prot}.fasta")
                fetch_sequence_ncbi_single(protein, nom_prot, fasta_interest_dir)
            else :
                fasta_interest_dir = os.path.join(fasta_protein, f"{nom_prot}.fasta")

            if os.path.exists(fasta_interest_dir) and os.path.exists(db_fasta):
                blast_output = os.path.join(blast_output_dir, f"{espece}_blast.txt")
                print(f"Lancement de BLASTP pour {espece}...")
            
                print(fasta_interest_dir)
                
                results = {}
                run_blastp(fasta_interest_dir, db_fasta, blast_output)

                best_hits = parse_blast_results(blast_output)

                results[nom_prot] = {
                    "best_hits": best_hits.get(protein, [])
                }
                ortholog_sonic = parse_hits_in_orthogroups(sonic_dict, results)
                orthologs_sonic_group = retrieve_orthologous_proteins(sonic_dict, ortholog_sonic)

                ortholog_orthologer = parse_hits_in_orthogroups(orthologer_dict, results)
                orthologs_orthologer_group = retrieve_orthologous_proteins(orthologer_dict, ortholog_orthologer)

                ortholog_orthofinder = parse_hits_in_orthogroups(orthofinder_dict, results)
                orthologs_orthofinder_group = retrieve_orthologous_proteins(orthofinder_dict, ortholog_orthofinder)

                # Sélectionner dynamiquement la méthode avec le plus petit groupe total
                group_sizes = {
                    "sonic": sum(len(g) for g in orthologs_sonic_group.values()),
                    "orthologer": sum(len(g) for g in orthologs_orthologer_group.values()),
                    "orthofinder": sum(len(g) for g in orthologs_orthofinder_group.values())
                }

                smallest_method = min(group_sizes, key=group_sizes.get)

                if smallest_method == "sonic":
                    merged_proteins_by_query = merge_orthologous_groups_by_protein(
                    ortholog_orthologer, orthologs_orthologer_group, ortholog_orthofinder, orthologs_orthofinder_group)
                    final_merged_groups = remove_duplicate_groups(merged_proteins_by_query)

                    ref_groups = orthologs_sonic_group
                    analysis_groups, no_group_proteins = analyze_orthologous_groups(final_merged_groups, sonic_dict)
                
                elif smallest_method == "orthologer":
                    merged_proteins_by_query = merge_orthologous_groups_by_protein(
                    ortholog_sonic, orthologs_sonic_group, ortholog_orthofinder, orthologs_orthofinder_group)
                    final_merged_groups = remove_duplicate_groups(merged_proteins_by_query)

                    ref_groups = orthologs_orthologer_group
                    analysis_groups, no_group_proteins = analyze_orthologous_groups(final_merged_groups, orthologer_dict)

                else:  # orthofinder
                    merged_proteins_by_query = merge_orthologous_groups_by_protein(
                    ortholog_orthologer, orthologs_orthologer_group, ortholog_sonic, orthologs_sonic_group)
                    final_merged_groups = remove_duplicate_groups(merged_proteins_by_query)

                    ref_groups = orthologs_orthofinder_group
                    analysis_groups, no_group_proteins = analyze_orthologous_groups(final_merged_groups, orthofinder_dict)


                # Comparaison avec la méthode à plus petits groupes
                new_groups, max_z_by_species = compare_species_between_groups(
                    ref_groups, analysis_groups, no_group_proteins,
                    species_file, fasta_interest_dir, fasta_dir,
                    prot_interest_file, output_results, "similarity"
                )

                orthogroup output

                def load_fasta_mapping(particule_file):
                    """Charge le fichier CSV et associe chaque particule à un fichier FASTA."""
                    fasta_mapping = {}
                    with open(particule_file, mode='r') as file:
                        reader = csv.DictReader(file, delimiter=',')
                        for row in reader:
                            organism = row['Organism']
                            file_name = row['File']
                            particle = row['Particule']
                            fasta_mapping[particle] = file_name
                    return fasta_mapping

                def extract_sequences_by_group(new_groups, fasta_dir, csv_file):
                    """
                    Récupère les séquences des protéines dans les groupes fusionnés, en parcourant les fichiers FASTA.
                    
                    :param new_groups: dictionnaire {query_protein: {group_id: {'from_groups': [...], 'proteins': [...]}}}
                    :param fasta_dir: répertoire contenant les fichiers FASTA
                    :param csv_file: fichier CSV mappant les "particles" aux noms de fichiers FASTA
                    :return: un dictionnaire des séquences par protéine requête puis par groupe
                    """
                    group_sequences = {}

                    fasta_mapping = load_fasta_mapping(csv_file)

                    for query_protein, groups in new_groups.items():
                        group_sequences[query_protein] = {}

                        for group_id, group_data in groups.items():
                            group_sequences[query_protein][group_id] = []

                            for identifiant in group_data:
                                for protein in group_data[identifiant]:
                                    particle = protein.split('_')[0]

                                    if particle in fasta_mapping:
                                        fasta_file = fasta_mapping[particle]
                                        fasta_path = os.path.join(fasta_dir, fasta_file)

                                        if os.path.exists(fasta_path):
                                            found = False
                                            for record in SeqIO.parse(fasta_path, "fasta"):
                                                if record.id == protein:
                                                    group_sequences[query_protein][group_id].append((record.id, str(record.seq)))
                                                    found = True
                                                    break
                                            if not found:
                                                print(f"Identifiant {protein} non trouvé dans {fasta_file}")
                                        else:
                                            print(f"Fichier FASTA inexistant : {fasta_path}")
                                    else:
                                        print(f"Particule '{particle}' non trouvée dans le fichier CSV.")

                    return group_sequences

                def write_fasta_for_groups(group_sequences, fasta_output_dir):
                    """
                    Écrit les séquences des groupes orthologues dans des fichiers FASTA distincts.

                    Parameters:
                    - group_sequences : dict : {group_id: [(seq_id, sequence), ...]}
                    - fasta_output_dir : str : Chemin du dossier de sortie pour les fichiers FASTA
                    """
                    # Création du dossier de sortie s'il n'existe pas
                    os.makedirs(fasta_output_dir, exist_ok=True)

                    for requete, group in group_sequences.items():
                        for group_id, sequences in group.items():
                            fasta_path = os.path.join(fasta_output_dir, f"{group_id}.fasta")

                            with open(fasta_path, "w") as fasta_file:
                                for seq_id, sequence in sequences:
                                    fasta_file.write(f">{seq_id}\n{sequence}\n")

                def extract_sequences_by_group(new_groups, fasta_dir, csv_file):

                    group_sequences = defaultdict(list)

                    fasta_mapping = load_fasta_mapping(csv_file)
                    for groups in new_groups.items():
                        for prot in groups[1]:
                            particle = prot.split('_')[0]

                            if particle in fasta_mapping:
                                fasta_file = fasta_mapping[particle]
                                fasta_path = os.path.join(fasta_dir, fasta_file)

                                if os.path.exists(fasta_path):
                                    found = False
                                    for record in SeqIO.parse(fasta_path, "fasta"):
                                        if record.id == prot:
                                            group_sequences[groups[0]].append((record.id, str(record.seq)))
                                            found = True
                                            break
                                    if not found:
                                        print(f"Identifiant {prot} non trouvé dans {fasta_file}")
                                else:
                                    print(f"Fichier FASTA inexistant : {fasta_path}")
                            else:
                                print(f"Particule '{particle}' non trouvée dans le fichier CSV.")

                    return group_sequences


                def write_fasta_for_groups(group_sequences, fasta_output_dir):
                    """
                    Écrit les séquences des groupes orthologues dans des fichiers FASTA distincts.

                    Parameters:
                    - group_sequences : dict : {group_id: [(seq_id, sequence), ...]}
                    - fasta_output_dir : str : Chemin du dossier de sortie pour les fichiers FASTA
                    """
                    # Création du dossier de sortie s'il n'existe pas
                    os.makedirs(fasta_output_dir, exist_ok=True)

                    for group_id, sequences in group_sequences.items():
                        fasta_path = os.path.join(fasta_output_dir, f"{group_id}.fasta")
                        
                        with open(fasta_path, "w") as fasta_file:
                            for seq_id, sequence in sequences:
                                fasta_file.write(f">{seq_id}\n{sequence}\n")

                fasta_output_dir = "/media/isabelle_florent_linux/hard_disk/donnees/test/Ouput/Fasta_orthologer_PCR7"

                # Extraire les séquences des orthologues
                group_sequences = extract_sequences_by_group(orthologs_orthofinder_group, fasta_dir, species_file)
                write_fasta_for_groups(group_sequences, fasta_output_dir)

                fin test

                output_file = os.path.join(output_results, f"final_groups/{nom_prot}_orthogroups.txt")
                dir_fgroups = os.path.join(output_results, f"final_groups")
                if not os.path.exists(dir_fgroups) : 
                    os.makedirs(dir_fgroups)
                with open(output_file, "w") as f:
                    for group, id_dict in new_groups.items():
                        f.write(f">{group}\n")
                        for cluster_id, subcluster in id_dict.items():
                            for sub_id, seq_list in subcluster.items():
                                f.write(f">{sub_id}\n")
                                f.write("\n".join(map(str, seq_list)) + "\n\n")
                
                #write_final_groups_with_sequences(nom_prot, new_groups, output_results, fasta_dir)
                protein_mapping = load_protein_mapping(prot_interest_file)
                presence_absence_single_df = generate_presence_absence_table_single_group(new_groups, species_dict)

                merged_orthogroups = merge_ortholog_groups(ortholog_sonic, ortholog_orthofinder, ortholog_orthologer)
                renamed_orthogroups = rename_keys(merged_orthogroups, protein_mapping)

                presence_absence_df = generate_presence_absence_table(orthologs_sonic_group, orthologs_orthologer_group, orthologs_orthofinder_group, species_dict, renamed_orthogroups)
                presence_combined_df = pd.concat([presence_absence_df, presence_absence_single_df], axis=0)

                PB CSUI
                target_column = "Cystoisospora suis"
                target_prefix = "Csui"
                group_map = {}  # {protein_id: group_index}

                # Création d'une liste unique des groupes fusion_0 pour tous les new_groups
                group_names_ordered = []  # Pour donner un index cohérent
                seen_group_names = set()

                # Étape 1 : Indexation cohérente de tous les groupes
                for group_data in new_groups.values():
                    if "fusion_0" in group_data:
                        for group_name in group_data["fusion_0"]:
                            if group_name not in seen_group_names:
                                group_names_ordered.append(group_name)
                                seen_group_names.add(group_name)

                # Étape 2 : Exploration de tous les groupes pour chercher les protéines Csui
                for outer_key, group_data in new_groups.items():
                    fusion_groups = group_data.get("fusion_0", {})
                    for group_name, group_content in fusion_groups.items():
                        group_index = group_names_ordered.index(group_name) + 1  # décalage dynamique

                        if isinstance(group_content, dict):
                            # Sous-groupes (ex: 32674: [...])
                            for prot_list in group_content.values():
                                for prot in prot_list:
                                    if prot.startswith(target_prefix):
                                        group_map[prot] = group_index
                        elif isinstance(group_content, list):
                            # Groupe plat (ex: no_group)
                            for prot in group_content:
                                if prot.startswith(target_prefix):
                                    group_map[prot] = group_index


                if group_map:
                    group_value = list(group_map.values())[0]
                    presence_combined_df[target_column] = group_value
                    print(presence_combined_df[target_column])
                

                fusion_row = presence_combined_df[presence_combined_df["Query Protein"] == f"{nom_prot}_fusion"].iloc[0]

                fusion_dict = {}

                for species in presence_combined_df.columns[1:]:
                    fusion_value = fusion_row[species]
                    z_score = max_z_by_species.get(species, None)

                    fusion_dict[species] = {
                        "fusion": fusion_value,
                        "z_score": z_score
                    }

                dict_esp_fusion[nom_prot] = fusion_dict

                busco_scores = get_busco_scores(species_file)
                busco_scores_line = [busco_scores.get(leaf, None) for leaf in presence_combined_df.columns[1:]]
                busco_df = pd.DataFrame({'Species': presence_combined_df.columns[1:], 'BUSCO_score': busco_scores_line})
                busco_df.set_index('Species', inplace=True)
                busco_df = busco_df.astype(int)

                sorted_df = presence_combined_df.copy()
                sorted_df["Query Protein"] = sorted_df["Query Protein"].astype(str)
                sorted_df = sorted_df.sort_values(by="Query Protein", key=lambda x: x.map(sort_key))

                sorted_columns = sorted(presence_combined_df.columns[1:])
                sorted_df = sorted_df[["Query Protein"] + sorted_columns]

                df_numeric = sorted_df.set_index("Query Protein").astype(int)

                leaf_names = [leaf.name.replace("'", "").strip() for leaf in tree.iter_leaves()]
                for leaf in tree.iter_leaves():
                    leaf.name = leaf.name.replace("'", "").strip()
                df = df_numeric[leaf_names]

                busco_df = busco_df.loc[leaf_names]

                methods = ["orthofinder", "orthologer", "sonic", "fusion"]
                all_proteins = sorted({
                    "_".join(prot.split("_")[:-1])
                    for prot in df.index
                    if prot.count("_") >= 2 and prot.split("_")[-1] in methods
                })

                output_images_dir = os.path.join(output_results, "plot")
                os.makedirs(output_images_dir, exist_ok=True)

                for prot in all_proteins:
                    selected_rows = [f"{prot}_{method}" for method in methods]
                    df_subset = df.loc[df.index.intersection(selected_rows)]
                    df_subset = df_subset.reindex(selected_rows)

                    fig, axes = plt.subplots(1, 3, figsize=(30, 20), gridspec_kw={'width_ratios': [1.3, 0.1, 1.8], 'wspace': 0.0001})

                    def master_ly(node):
                        if node.is_leaf():
                            F = faces.TextFace(node.name, fsize=150)
                            F.hz_align = 2
                            faces.add_face_to_node(F, node, column=1, position="aligned")

                        node.img_style["vt_line_width"] = 50  # Épaisseur de la branche verticale
                        node.img_style["hz_line_width"] = 50  # Épaisseur de la branche horizontale
                        node.img_style["vt_line_color"] = "black"
                        node.img_style["hz_line_color"] = "black"

                    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
                    tree_style = TreeStyle()
                    tree_style.show_leaf_name = False
                    tree_style.draw_guiding_lines = True
                    tree_style.guiding_lines_color = "black"
                    tree_style.guiding_lines_type = 0
                    tree_style.scale = 1000
                    tree_style.layout_fn = master_ly
                    tree.render(temp_file.name, w=2000, tree_style=tree_style)
                    img = plt.imread(temp_file.name)
                    axes[0].imshow(img)
                    axes[0].axis("off")

                    sns.heatmap(busco_df, cmap="RdYlGn", annot=True, fmt="d", cbar=False, ax=axes[1], linewidths=0.3, linecolor="black", annot_kws={'size': 8})
                    axes[1].set_yticklabels([])
                    axes[1].set_ylabel("")
                    axes[1].set_xticklabels(["BUSCO Score"], rotation=45)
                    axes[1].tick_params(axis='both', length=0)

                    max_group_id = df_subset.to_numpy().max()
                    if pd.isna(max_group_id) or max_group_id < 1:
                        max_group_id = 1

                    fixed_colors = {
                        0: "#FFFFFF", # Blanc pour "Absence"
                        1: "#228B22", # Vert foncé pour "Groupe de référence"
                        'intermediate': sns.color_palette("husl", max(0, max_group_id - 2)),
                        max_group_id: "#e31a1c" # Rouge foncé pour "No_group"
                    }
                    

                    color_list = [fixed_colors[0], fixed_colors[1]] + fixed_colors['intermediate'] + [fixed_colors[max_group_id]]
                    custom_cmap = ListedColormap(color_list)

                    ax = sns.heatmap(df_subset.T, cmap=custom_cmap, linewidth=0.3,
                                    linecolor="black", cbar=True, xticklabels=True,
                                    yticklabels=True, annot=False, ax=axes[2], vmin=0, vmax=max_group_id, annot_kws={'size': 20})

                    # Ajout des z-scores sur la ligne "fusion"
                    for col_idx, leaf in enumerate(leaf_names):
                        if leaf not in max_z_by_species:
                            continue  # On saute si pas de z-score pour cette espèce

                        z_score = max_z_by_species[leaf]

                        text_color = "white"
                        
                        if z_score <= 1:
                            n_stars = 0
                        elif 1 < z_score <= 2:
                            n_stars = 2
                        elif 2 < z_score <= 3:
                            n_stars = 5
                        else:  # z > 3
                            n_stars = 10

                        if n_stars > 0:
                            stars = "★" * n_stars
                            s = f"{stars} {z_score:.2f} {stars}"
                        else:
                            s = f"{z_score:.2f}"

                        axes[2].text(
                            y=col_idx + 0.5,
                            x=df_subset.index.get_loc(f"{prot}_fusion") + 0.5,
                            s=s,
                            ha='center',
                            va='center',
                            fontsize=6,
                            color=text_color,
                            fontweight='bold'
                        )

                    axes[2].set_xlabel("Request protein")
                    axes[2].tick_params(axis="x", labelrotation=45)
                    axes[2].set_yticklabels([])
                    axes[2].set_yticks([])

                    colorbar = ax.collections[0].colorbar
                    ticks = list(range(0, max_group_id + 1))
                    ticklabels = []

                    for i in ticks:
                        if i == 0:
                            ticklabels.append("Absence")
                        elif i == 1:
                            ticklabels.append("Reference group")
                        elif i == max_group_id:
                            ticklabels.append("Orphans")
                        else:
                            ticklabels.append(f"Group {i - 1}")

                    colorbar.set_ticks(ticks)
                    colorbar.set_ticklabels(ticklabels)

                    output_path = os.path.join(output_images_dir, f"{prot}.png")
                    plt.savefig(output_path, dpi=300, bbox_inches='tight')

                    fig.subplots_adjust(left=0.01, right=0.97, top=0.99, bottom=0.15)
'''

def load_all_orthogroups(output_orthogroup_file):
    """
    Load all orthogroup mappings from SonicParanoid, Orthologer, and OrthoFinder.

    This function reads the orthogroup CSV files from each method and loads them 
    as dictionaries mapping orthogroup IDs to sets of protein IDs.

    Parameters
    ----------
    output_orthogroup_file : str
        Path to the directory containing the three orthogroup result files.

    Returns
    -------
    tuple of dict
        Tuple containing dictionaries for SonicParanoid, Orthologer, and OrthoFinder.
    """
    return (
        load_groups(os.path.join(output_orthogroup_file, "sonicParanoid.csv")),
        load_groups(os.path.join(output_orthogroup_file, "orthologer.csv")),
        load_groups(os.path.join(output_orthogroup_file, "orthofinder.csv"))
    )


def get_interest_fasta(protein, nom_prot, output_results, fasta_protein):
    """
    Retrieve or locate the FASTA file for the protein of interest.

    This function fetches the protein sequence from NCBI if no local FASTA directory
    is provided, or retrieves the path to an existing file.

    Parameters
    ----------
    protein : str
        Protein identifier (e.g., NCBI or UniProt ID).
    nom_prot : str
        Local name assigned to the protein.
    output_results : str
        Path to the output directory for results.
    fasta_protein : str or None
        Optional path to a directory containing pre-downloaded FASTA files.

    Returns
    -------
    str
        Path to the FASTA file of the protein of interest.
    """
    if fasta_protein is None:
        fasta_interest_dir = os.path.join(output_results, f"{nom_prot}.fasta")
        fetch_sequence_ncbi_single(protein, nom_prot, fasta_interest_dir)
    else:
        fasta_interest_dir = os.path.join(fasta_protein, f"{nom_prot}.fasta")
    return fasta_interest_dir


def run_and_parse_blast(espece, protein, nom_prot, fasta_interest_dir, db_fasta, output_results):
    """
    Run BLASTP for a protein against a database and parse the best hits.

    This function executes BLASTP with the input FASTA and database, saves the result,
    and parses the best matching proteins.

    Parameters
    ----------
    espece : str
        Species name or identifier used for output naming.
    protein : str
        Protein ID used to filter best hits.
    nom_prot : str
        Local name assigned to the protein.
    fasta_interest_dir : str
        Path to the protein FASTA file.
    db_fasta : str
        Path to the BLAST-formatted FASTA database.
    output_results : str
        Output directory to store BLAST results.

    Returns
    -------
    dict
        Dictionary with the protein name as key and best hits as values.
    """
    blast_output_dir = os.path.join(output_results, "blast_results")
    os.makedirs(blast_output_dir, exist_ok=True)

    blast_output = os.path.join(blast_output_dir, f"{espece}_blast.txt")
    run_blastp(fasta_interest_dir, db_fasta, blast_output)
    best_hits = parse_blast_results(blast_output)
    return {nom_prot: {"best_hits": best_hits.get(protein, [])}}


def choose_best_orthogroup_method(best_hits, sonic_dict, orthologer_dict, orthofinder_dict,
                                  species_file, fasta_interest_dir, fasta_dir,
                                  prot_interest_file, output_results):
    """
    Select the best orthogroup inference method and compute merged ortholog groups.

    This function:
    - Compares ortholog group sizes from each method (SonicParanoid, Orthologer, OrthoFinder).
    - Uses the method with the smallest total group size as the reference.
    - Merges ortholog groups from the other two methods.
    - Removes duplicate groups and compares species representation.
    - Computes final ortholog groups with interspecies comparisons.

    Parameters
    ----------
    best_hits : dict
        Dictionary of best protein hits from BLASTP.
    sonic_dict : dict
        Orthogroup mapping from SonicParanoid.
    orthologer_dict : dict
        Orthogroup mapping from Orthologer.
    orthofinder_dict : dict
        Orthogroup mapping from OrthoFinder.
    species_file : str
        Path to the file listing species information.
    fasta_interest_dir : str
        Path to the FASTA file of the protein of interest.
    fasta_dir : str
        Path to the directory containing all species' protein FASTA files.
    prot_interest_file : str
        Path to the protein interest metadata file.
    output_results : str
        Directory to store the output results.

    Returns
    -------
    tuple
        A tuple containing:
        - new_groups: Final merged ortholog groups after filtering.
        - ref_groups: Ortholog groups from the smallest method.
        - max_z_by_species: Dictionary with maximum z-score similarity by species.
    """
    ortholog_sets = {}
    sizes = {}

    for method, method_dict in zip(["sonic", "orthologer", "orthofinder"],
                                   [sonic_dict, orthologer_dict, orthofinder_dict]):
        parsed = parse_hits_in_orthogroups(method_dict, best_hits)
        group = retrieve_orthologous_proteins(method_dict, parsed)
        ortholog_sets[method] = (parsed, group)
        sizes[method] = sum(len(g) for g in group.values())

    smallest_method = min(sizes, key=sizes.get)
    other_methods = [m for m in ["sonic", "orthologer", "orthofinder"] if m != smallest_method]

    merged_proteins_by_query = merge_orthologous_groups_by_protein(
        ortholog_sets[other_methods[0]][0], ortholog_sets[other_methods[0]][1],
        ortholog_sets[other_methods[1]][0], ortholog_sets[other_methods[1]][1]
    )

    final_merged_groups = remove_duplicate_groups(merged_proteins_by_query)

    analysis_groups, no_group_proteins = analyze_orthologous_groups(
        final_merged_groups,
        eval(f"{smallest_method}_dict")
    )

    ref_groups = ortholog_sets[smallest_method][1]

    new_groups, max_z_by_species = compare_species_between_groups(
        ref_groups, analysis_groups, no_group_proteins,
        species_file, fasta_interest_dir, fasta_dir,
        prot_interest_file, output_results, "similarity"
    )

    return new_groups, ref_groups, max_z_by_species


def save_final_groups(nom_prot, new_groups, output_results):
    """
    Save the final orthologous groups into a formatted text file.

    This function:
    - Creates a directory to store the final orthogroups if it doesn't already exist.
    - Writes each orthogroup with its associated subclusters and sequences into a file.

    Parameters
    ----------
    nom_prot : str
        Name of the protein of interest.
    new_groups : dict
        Dictionary containing final orthologous group information, typically structured by group IDs, cluster IDs, and sequence lists.
    output_results : str
        Path to the main output directory.

    Returns
    -------
    str
        Path to the saved orthogroup file.
    """
    dir_fgroups = os.path.join(output_results, "final_groups")
    os.makedirs(dir_fgroups, exist_ok=True)

    output_file = os.path.join(dir_fgroups, f"{nom_prot}_orthogroups.txt")
    with open(output_file, "w") as f:
        for group, id_dict in new_groups.items():
            f.write(f">{group}\n")
            for cluster_id, subcluster in id_dict.items():
                for sub_id, seq_list in subcluster.items():
                    f.write(f">{sub_id}\n")
                    f.write("\n".join(map(str, seq_list)) + "\n\n")
    return output_file


def build_presence_absence_tables(best_hits, new_groups, species_dict,
                                  orthologs_ref, sonic_dict, orthologer_dict, orthofinder_dict,
                                  prot_interest_file):
    """
    Construct presence/absence tables from orthologous group data.

    This function:
    - Loads protein mapping data.
    - Generates a presence/absence table for the final (custom) orthogroups.
    - Reconstructs orthologous groups using SonicParanoid, Orthologer, and OrthoFinder.
    - Merges and renames orthogroups for unified presence/absence tracking.
    - Combines all information into a final DataFrame.

    Parameters
    ----------
    best_hits : dict
        Dictionary mapping query proteins to their best BLAST hits.
    new_groups : dict
        Final orthologous groups constructed through custom logic.
    species_dict : dict
        Dictionary mapping species codes to full names or other identifiers.
    orthologs_ref : dict
        Reference ortholog group (smallest set) used for comparison.
    sonic_dict : dict
        Dictionary of orthogroups from SonicParanoid.
    orthologer_dict : dict
        Dictionary of orthogroups from Orthologer.
    orthofinder_dict : dict
        Dictionary of orthogroups from OrthoFinder.
    prot_interest_file : str
        Path to the file mapping protein identifiers to query names.

    Returns
    -------
    pd.DataFrame
        Combined presence/absence table for custom and tool-based orthogroups.
    """
    protein_mapping = load_protein_mapping(prot_interest_file)

    presence_absence_single_df = generate_presence_absence_table_single_group(new_groups, species_dict)

    ortholog_sonic = parse_hits_in_orthogroups(sonic_dict, best_hits)
    orthologs_sonic_group = retrieve_orthologous_proteins(sonic_dict, ortholog_sonic)

    ortholog_orthologer = parse_hits_in_orthogroups(orthologer_dict, best_hits)
    orthologs_orthologer_group = retrieve_orthologous_proteins(orthologer_dict, ortholog_orthologer)

    ortholog_orthofinder = parse_hits_in_orthogroups(orthofinder_dict, best_hits)
    orthologs_orthofinder_group = retrieve_orthologous_proteins(orthofinder_dict, ortholog_orthofinder)

    merged_orthogroups = merge_ortholog_groups(
        ortholog_sonic, ortholog_orthofinder, ortholog_orthologer)

    renamed_orthogroups = rename_keys(merged_orthogroups, protein_mapping)

    presence_absence_df = generate_presence_absence_table(
        orthologs_sonic_group, orthologs_orthologer_group, orthologs_orthofinder_group, species_dict, renamed_orthogroups
    )

    return pd.concat([presence_absence_df, presence_absence_single_df], axis=0)


def analyze_fusion_groups(nom_prot, new_groups, presence_df, max_z_by_species):
    """
    Analyze fusion groups to detect species with putative gene fusion events.

    This function:
    - Identifies proteins belonging to fusion groups specific to the target species.
    - Assigns a unique group index to each fusion group.
    - Annotates the presence/absence table with fusion group values for the target species.
    - Constructs a dictionary summarizing fusion group and z-score values by species.

    Parameters
    ----------
    nom_prot : str
        Name of the query protein used to identify fusion events.
    new_groups : dict
        Dictionary of orthologous groups including potential fusion group information.
    presence_df : pd.DataFrame
        DataFrame containing presence/absence data across species.
    max_z_by_species : dict
        Dictionary mapping species to their maximum Z-score similarity values.

    Returns
    -------
    dict
        Dictionary summarizing fusion group membership and similarity Z-scores per species.
    """
    target_prefix = "Csui"
    target_column = "Cystoisospora suis"

    group_names_ordered = []
    group_map = {}
    seen_group_names = set()

    for group_data in new_groups.values():
        if "fusion_0" in group_data:
            for group_name in group_data["fusion_0"]:
                if group_name not in seen_group_names:
                    group_names_ordered.append(group_name)
                    seen_group_names.add(group_name)

    for outer_key, group_data in new_groups.items():
        fusion_groups = group_data.get("fusion_0", {})
        for group_name, group_content in fusion_groups.items():
            group_index = group_names_ordered.index(group_name) + 1
            if isinstance(group_content, dict):
                for prot_list in group_content.values():
                    for prot in prot_list:
                        if prot.startswith(target_prefix):
                            group_map[prot] = group_index
            elif isinstance(group_content, list):
                for prot in group_content:
                    if prot.startswith(target_prefix):
                        group_map[prot] = group_index

    fusion_dict = {}
    if group_map:
        group_value = list(group_map.values())[0]
        presence_df[target_column] = group_value

        fusion_row = presence_df[presence_df["Query Protein"] == f"{nom_prot}_fusion"].iloc[0]
        for species in presence_df.columns[1:]:
            fusion_value = fusion_row[species]
            z_score = max_z_by_species.get(species, None)
            fusion_dict[species] = {
                "fusion": fusion_value,
                "z_score": z_score
            }

    return fusion_dict


def plot_presence_absence_matrices(presence_df, busco_scores, tree, output_results, max_z_by_species):
    """
    Visualize presence/absence matrices of orthologous groups with phylogenetic tree and BUSCO scores.

    This function:
    - Sorts and structures a presence/absence matrix of proteins across species.
    - Displays a phylogenetic tree, BUSCO completeness scores, and a heatmap of orthogroup memberships.
    - Highlights fusion groups with custom color codes and annotates significant Z-score similarities.
    - Saves per-protein visualizations in the output directory.

    Parameters
    ----------
    presence_df : pd.DataFrame
        DataFrame with presence/absence information per query protein and species.
    busco_scores : dict
        Dictionary mapping species names to their BUSCO completeness scores (int).
    tree : ete3.Tree
        Phylogenetic tree object representing species relationships.
    output_results : str
        Path to the directory where the resulting plots will be saved.
    max_z_by_species : dict
        Dictionary mapping species names to their maximum Z-score for fusion group similarity.

    Returns
    -------
    tuple
        A tuple containing:
        - pd.DataFrame: DataFrame of BUSCO scores indexed by species.
        - list: Ordered list of species names (leaf names from the tree).
    """
    # Sort presence/absence dataframe based on custom logic
    sorted_df = presence_df.copy()
    sorted_df["Query Protein"] = sorted_df["Query Protein"].astype(str)
    sorted_df = sorted_df.sort_values(by="Query Protein", key=lambda x: x.map(sort_key))

    sorted_columns = sorted(presence_df.columns[1:])
    sorted_df = sorted_df[["Query Protein"] + sorted_columns]
    df_numeric = sorted_df.set_index("Query Protein").astype(int)

    # Clean leaf names from the tree
    leaf_names = [leaf.name.replace("'", "").strip() for leaf in tree.iter_leaves()]
    for leaf in tree.iter_leaves():
        leaf.name = leaf.name.replace("'", "").strip()

    df = df_numeric[leaf_names]
    busco_df = pd.DataFrame({'Species': df.columns, 'BUSCO_score': [busco_scores.get(sp, None) for sp in df.columns]})
    busco_df.set_index('Species', inplace=True)
    busco_df = busco_df.astype(int)

    methods = ["orthofinder", "orthologer", "sonic", "fusion"]
    all_proteins = sorted({
        "_".join(prot.split("_")[:-1])
        for prot in df.index
        if prot.count("_") >= 2 and prot.split("_")[-1] in methods
    })

    output_images_dir = os.path.join(output_results, "plot")
    os.makedirs(output_images_dir, exist_ok=True)

    for prot in all_proteins:
        selected_rows = [f"{prot}_{method}" for method in methods]
        df_subset = df.loc[df.index.intersection(selected_rows)]
        df_subset = df_subset.reindex(selected_rows)

        # Create the figure layout: tree + BUSCO + heatmap
        fig, axes = plt.subplots(1, 3, figsize=(30, 20), gridspec_kw={'width_ratios': [1.3, 0.1, 1.8], 'wspace': 0.0001})

        # Tree rendering style with thick black branches and aligned names
        def master_ly(node):
            if node.is_leaf():
                F = faces.TextFace(node.name, fsize=150)
                F.hz_align = 2
                faces.add_face_to_node(F, node, column=1, position="aligned")

            node.img_style["vt_line_width"] = 50
            node.img_style["hz_line_width"] = 50
            node.img_style["vt_line_color"] = "black"
            node.img_style["hz_line_color"] = "black"

        # Render tree to temporary file and load image into plot
        temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        tree_style = TreeStyle()
        tree_style.show_leaf_name = False
        tree_style.draw_guiding_lines = True
        tree_style.guiding_lines_color = "black"
        tree_style.guiding_lines_type = 0
        tree_style.scale = 1000
        tree_style.layout_fn = master_ly
        tree.render(temp_file.name, w=2000, tree_style=tree_style)
        img = plt.imread(temp_file.name)
        axes[0].imshow(img)
        axes[0].axis("off")

        # Plot BUSCO scores
        sns.heatmap(busco_df, cmap="RdYlGn", annot=True, fmt="d", cbar=False, ax=axes[1], linewidths=0.3, linecolor="black", annot_kws={'size': 8})
        axes[1].set_yticklabels([])
        axes[1].set_ylabel("")
        axes[1].set_xticklabels(["BUSCO Score"], rotation=45)
        axes[1].tick_params(axis='both', length=0)

        max_group_id = df_subset.to_numpy().max()
        if pd.isna(max_group_id) or max_group_id < 1:
            max_group_id = 1

        fixed_colors = {
            0: "#FFFFFF", #Absence
            1: "#228B22", # Reference group
            'intermediate': sns.color_palette("husl", max(0, max_group_id - 2)),
            max_group_id: "#e31a1c" # Orphans
        }
        
        color_list = [fixed_colors[0], fixed_colors[1]] + fixed_colors['intermediate'] + [fixed_colors[max_group_id]]
        custom_cmap = ListedColormap(color_list)

        # Plot presence/absence heatmap
        ax = sns.heatmap(df_subset.T, cmap=custom_cmap, linewidth=0.3,
                        linecolor="black", cbar=True, xticklabels=True,
                        yticklabels=True, annot=False, ax=axes[2], vmin=0, vmax=max_group_id, annot_kws={'size': 20})

        # Annotate z-scores for fusion rows
        for col_idx, leaf in enumerate(leaf_names):
            if leaf not in max_z_by_species:
                continue

            z_score = max_z_by_species[leaf]
            text_color = "white"
            
            if z_score <= 1:
                n_stars = 0
            elif 1 < z_score <= 2:
                n_stars = 2
            elif 2 < z_score <= 3:
                n_stars = 5
            else:  # z > 3
                n_stars = 10

            if n_stars > 0:
                stars = "★" * n_stars
                s = f"{stars} {z_score:.2f} {stars}"
            else:
                s = f"{z_score:.2f}"

            axes[2].text(
                y=col_idx + 0.5,
                x=df_subset.index.get_loc(f"{prot}_fusion") + 0.5,
                s=s,
                ha='center',
                va='center',
                fontsize=6,
                color=text_color,
                fontweight='bold'
            )

        axes[2].set_xlabel("Request protein")
        axes[2].tick_params(axis="x", labelrotation=45)
        axes[2].set_yticklabels([])
        axes[2].set_yticks([])

        colorbar = ax.collections[0].colorbar
        ticks = list(range(0, max_group_id + 1))
        ticklabels = []

        for i in ticks:
            if i == 0:
                ticklabels.append("Absence")
            elif i == 1:
                ticklabels.append("Reference group")
            elif i == max_group_id:
                ticklabels.append("Orphans")
            else:
                ticklabels.append(f"Group {i - 1}")

        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels(ticklabels)

        # Save figure
        output_path = os.path.join(output_images_dir, f"{prot}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        fig.subplots_adjust(left=0.01, right=0.97, top=0.99, bottom=0.15)

        return busco_df, leaf_names

    
def run_table(fasta_dir, output_dir, tree_file, prot_interest_file, species_file, fasta_protein):
    """
        Execute a comprehensive analysis pipeline for orthogroups and gene fusion detection
        across multiple species, integrating phylogenetic context and BUSCO completeness scores.

        This function performs the following major steps:
        - Loads input data: species tree, protein of interest, species mappings, BUSCO scores,
        and orthogroup assignments from various methods.
        - Iterates over species and proteins of interest to:
            * Extract relevant FASTA sequences.
            * Run BLAST searches and parse results.
            * Select the best orthogroup assignment by integrating multiple orthology methods.
            * Build presence/absence tables combining all orthogroup data.
            * Analyze putative gene fusion events.
            * Plot presence/absence matrices alongside BUSCO completeness scores on the species tree.
        - Aggregates fusion and similarity (z-score) data across all proteins and species.
        - Generates a multi-panel figure showing:
            * The species phylogenetic tree.
            * BUSCO scores per species.
            * Heatmaps and symbolic representations of fusion and similarity patterns,
            including detailed annotations for fusion groups and confidence (z-scores).
        - Saves the final figure in the specified output directory.

        Parameters
        ----------
        fasta_dir : str
            Directory containing FASTA files for each species.
        output_dir : str
            Directory where output files and plots will be saved.
        tree_file : str
            File path to the species phylogenetic tree in Newick format.
        prot_interest_file : str
            File containing the proteins of interest organized by species.
        species_file : str
            File containing species mapping and BUSCO score information.
        fasta_protein : str
            Path to the master FASTA file containing protein sequences of interest.

        Returns
        -------
        None
            Outputs results and plots directly to the output directory.
    """
    output_orthogroup_file = os.path.join(output_dir, "orthogroups")
    tree = Tree(tree_file, format=1)

    # Load proteins of interest per species
    proteins_dict = read_protein_interest(prot_interest_file)
    # Load species name mapping
    species_dict = load_species_mapping(species_file)
    # Load BUSCO completeness scores for species
    busco_scores = get_busco_scores(species_file)

    # Load orthogroup assignments from different methods
    sonic_dict, orthologer_dict, orthofinder_dict = load_all_orthogroups(output_orthogroup_file)

    dict_esp_fusion = {}

    # Iterate over proteins of interest
    for espece, proteins in proteins_dict.items():
        db_fasta = os.path.join(fasta_dir, f"{espece}.fasta")
        for protein, nom_prot in proteins.items():
            output_results = os.path.join(output_dir, nom_prot)
            os.makedirs(output_results, exist_ok=True)

            # Extract fasta sequences related to protein of interest
            fasta_interest_dir = get_interest_fasta(protein, nom_prot, output_results, fasta_protein)

            if not os.path.exists(fasta_interest_dir) or not os.path.exists(db_fasta):
                continue

            # Run BLAST and parse best hits results
            best_hits = run_and_parse_blast(espece, protein, nom_prot, fasta_interest_dir, db_fasta, output_results)

            # Choose best orthogroup assignment combining multiple orthology methods
            final_merged_groups, ref_groups, max_z_by_species = choose_best_orthogroup_method(
                best_hits, sonic_dict, orthologer_dict, orthofinder_dict,
                species_file, fasta_interest_dir, fasta_dir,
                prot_interest_file, output_results
            )

            # Save final orthogroup assignments to file
            output_file = save_final_groups(nom_prot, final_merged_groups, output_results)

            # Build presence/absence matrices combining data from orthogroups and species info
            presence_combined_df = build_presence_absence_tables(
                best_hits, final_merged_groups, species_dict,
                ref_groups, sonic_dict, orthologer_dict, orthofinder_dict,
                prot_interest_file
            )

            # Analyze gene fusion groups and detect fusion events
            dict_esp_fusion[nom_prot] = analyze_fusion_groups(
                nom_prot, final_merged_groups, presence_combined_df, max_z_by_species
            )

            # Plot presence/absence matrices with BUSCO completeness on the species tree
            busco_df, leaf_names = plot_presence_absence_matrices(
                presence_combined_df, busco_scores, tree, output_results, max_z_by_species
            )

    zscore_df = pd.DataFrame({
        protein: {
            species: values['z_score']
            for species, values in species_dict.items()
        }
        for protein, species_dict in dict_esp_fusion.items()
    })

    fusion_df = pd.DataFrame({
        protein: {
            species: values['fusion']
            for species, values in species_dict.items()
        }
        for protein, species_dict in dict_esp_fusion.items()
    })

    zscore_df = zscore_df.loc[leaf_names]
    fusion_df = fusion_df.loc[leaf_names]

    n_cols = len(zscore_df.columns)
    n_rows = len(zscore_df.index)

    # Setup matplotlib figure with 3 panels: tree, BUSCO heatmap, fusion/zscore heatmap
    fig_final, axes_final = plt.subplots(
        1, 3,
        figsize=(22, 30),
        gridspec_kw={'width_ratios': [0.88, 0.05, 0.14], 'wspace': 0.0}
    )

    def master_ly(node):
        if node.is_leaf():
            F = faces.TextFace(node.name, fsize=150)
            F.hz_align = 2
            faces.add_face_to_node(F, node, column=1, position="aligned")

        node.img_style["vt_line_width"] = 30
        node.img_style["hz_line_width"] = 30
        node.img_style["vt_line_color"] = "black"
        node.img_style["hz_line_color"] = "black"

    # Render tree to temporary file and load image into plot
    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.draw_guiding_lines = True
    tree_style.guiding_lines_color = "black"
    tree_style.guiding_lines_type = 0
    tree_style.scale = 1000
    tree_style.layout_fn = master_ly
    tree.render(temp_file.name, w=2000, dpi=600, tree_style=tree_style)
    img = Image.open(temp_file.name)
    img = img.crop(img.getbbox())
    axes_final[0].imshow(img)
    axes_final[0].axis("off")

    # plot BUSCO scores
    sns.heatmap(busco_df, cmap="RdYlGn", annot=True, fmt="d", cbar=False, ax=axes_final[1], linewidths=0.3, linecolor="black", annot_kws={'size': 10})
    axes_final[1].set_xticks([0])
    axes_final[1].set_yticklabels([])
    axes_final[1].set_ylabel("")
    axes_final[1].set_xticklabels(["BUSCO Score"], rotation=45)
    axes_final[1].tick_params(axis='both', length=0)

    # Configure the third subplot for fusion and z-score visualization
    n_rows, n_cols = zscore_df.shape
    axes_final[2].set_xlim(0, n_cols)
    axes_final[2].set_ylim(0, n_rows)

    axes_final[2].set_xticks([x + 0.5 for x in range(n_cols)])
    axes_final[2].set_xticklabels(zscore_df.columns, fontsize=12, rotation=90)
    axes_final[2].xaxis.tick_top()

    box = axes_final[2].get_position()
    axes_final[2].set_position([box.x0, box.y0, box.width, box.height])

    axes_final[2].set_facecolor('white')
    list_max_group = []

    # Loop over species (columns) to draw fusion and z-score markers
    for j, species in enumerate(zscore_df.columns):
        fusion_values = [
            dict_esp_fusion.get(species, {}).get(protein, {}).get("fusion", 0)
            for protein in zscore_df.index
        ]
        max_group_id = max(fusion_values)
        list_max_group.append(max_group_id)

        intermediate_ids = list(range(2, max_group_id))
        fixed_colors = {
            0: "#FFFFFF",            # Absence: white
            1: "#228B22",            # Reference group: green
            "intermediate": sns.color_palette("husl", len(intermediate_ids)),  # Intermediate groups
            max_group_id: "#e31a1c"  # Orphans: red
        }

        # Loop over proteins (rows) to add circles and wedges indicating fusion and similarity
        for i, protein in enumerate(zscore_df.index):
            info = dict_esp_fusion.get(species, {}).get(protein, {})
            fusion = info.get("fusion", 0)
            z = info.get("z_score", None)
            center = (j + 0.5, n_rows - i - 0.5)

            if fusion == 0:
                edge_color = fixed_colors[0]
            elif fusion == 1:
                edge_color = fixed_colors[1]
            elif fusion == max_group_id:
                edge_color = fixed_colors[max_group_id]
            else:
                try:
                    idx = intermediate_ids.index(fusion)
                    edge_color = fixed_colors["intermediate"][idx]
                except ValueError:
                    edge_color = "#000000"  # fallback in black if error

            # Draw circle/wedge depending on z-score thresholds
            if pd.isna(z):
                # White circle
                circle = Circle(center, 0.3, facecolor='none', edgecolor='black')
                axes_final[2].add_patch(circle)

            elif z < 2:
                # Full circle
                circle = Circle(center, 0.3, facecolor=edge_color, edgecolor=edge_color)
                axes_final[2].add_patch(circle)

            elif 2 <= z < 3:
                # Half circle
                circle = Circle(center, 0.3, edgecolor=edge_color, facecolor='none', linewidth=1)
                axes_final[2].add_patch(circle)

                wedge = Wedge(center, 0.3, 270, 90, facecolor=edge_color, edgecolor='none')
                axes_final[2].add_patch(wedge)

            elif z >= 3:
                # Quarter circle
                circle = Circle(center, 0.3, edgecolor=edge_color, facecolor='none', linewidth=1)
                axes_final[2].add_patch(circle)

                wedge = Wedge(center, 0.3, 0, 90, facecolor=edge_color, edgecolor='none')
                axes_final[2].add_patch(wedge)

    for x in range(n_cols + 1):
        axes_final[2].axvline(x, color='black', lw=0.5)
    for y in range(n_rows + 1):
        axes_final[2].axhline(y, color='black', lw=0.5)

    axes_final[2].set_aspect('equal')  # square cases
    
    axes_final[2].set_yticks([])
    axes_final[2].tick_params(which="minor", length=0)
    axes_final[2].set_xlabel("")

    global_max_group_id = max(list_max_group)
    intermediate_ids = list(range(2, global_max_group_id))

    palette_inter = sns.color_palette("husl", len(intermediate_ids))

    legend_elements = [
        Patch(facecolor="#FFFFFF", edgecolor='black', label="Absence"),
        Patch(facecolor="#228B22", edgecolor='black', label="Reference group"),
    ]

    for idx, color in enumerate(palette_inter):
        legend_elements.append(Patch(facecolor=color, edgecolor='black', label=f"Group {idx + 1}"))

    # Orphans
    legend_elements.append(Patch(facecolor="#e31a1c", edgecolor='black', label="Orphans"))

    # Add legend
    axes_final[2].legend(
        handles=legend_elements,
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        borderaxespad=0.,
        frameon=False,
        fontsize=10
    )

    plt.tight_layout()

    output_fig = os.path.join(output_dir, f"graph_final.png")
    plt.savefig(output_fig, dpi=400, bbox_inches='tight')