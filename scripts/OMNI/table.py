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
from matplotlib.patches import Wedge, Circle, Rectangle
from matplotlib.collections import PatchCollection
from PIL import Image
from matplotlib.patches import Patch


def fetch_sequence_ncbi_single(protein_id, protein_name, output_fasta, email="romain.daguerre1119@gmail.com"):
    Entrez.email = email  # Obligatoire pour utiliser Entrez

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

            # Nettoyage : une seule ligne de séquence
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
    """Lit le fichier CSV contenant les protéines d'intérêt et retourne un dictionnaire {espece: {code_proteine: nom_proteine}}."""
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
    """Charge le fichier orthogroup et crée un dictionnaire des groupes."""
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
    """Exécute BLASTP pour comparer les séquences d'intérêt avec la base de données."""
    db_name = db_fasta.replace(".fasta", "")
    
    if not os.path.exists(db_name + ".pin"):
        os.system(f"makeblastdb -in {db_fasta} -dbtype prot -out {db_name}")
    
    blastp_cline = NcbiblastpCommandline(query=query_fasta, db=db_name, outfmt=6, out=output_file)
    stdout, stderr = blastp_cline()

def parse_blast_results(blast_output):
    """Lit le fichier BLAST et extrait le meilleur hit et celui correspondant à l'ID d'intérêt."""
    best_hits = {}  # Stocke le meilleur hit {query_id: (best_hit_id, score)}

    with open(blast_output, 'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue  # Éviter les lignes incomplètes

            query_id, subject_id = cols[0], cols[1]
            identity, alignment_length, mismatches, gap_opens = map(float, cols[2:6])
            q_start, q_end, s_start, s_end = map(int, cols[6:10])
            evalue, bit_score = map(float, cols[10:12])

            # Stocker le meilleur hit (le plus grand bit score)
            if query_id not in best_hits or bit_score > best_hits[query_id][1]:
                best_hits[query_id] = (subject_id, bit_score)

    return best_hits

def parse_hits_in_orthogroups(ortholog_dict, results):
    """
    Recherche le meilleur hit de chaque protéine dans les groupes orthologues.

    Parameters:
    - ortholog_dict : dict : Dictionnaire {group_id: [list of protein ids]}
    - results : dict : Résultats du BLAST {nom_prot: {"best_hits": (hit, score)}}

    Returns:
    - dict : {nom_prot: {"best_hit": ..., "score": ..., "ortholog_group": ...}}
    """
    ortholog_groups = {}

    for nom_prot, data in results.items():
        best_hit, score = data["best_hits"]
        found_group = False

        for group_id, members in ortholog_dict.items():
            if best_hit in members:
                ortholog_groups[nom_prot] = {
                    "best_hit": best_hit,
                    "score": score,
                    "ortholog_group": group_id
                }
                found_group = True
                break

        if not found_group:
            ortholog_groups[nom_prot] = {
                "best_hit": best_hit,
                "score": score,
                "ortholog_group": None
            }

    return ortholog_groups

def retrieve_orthologous_proteins(ortholog_dict, ortholog_groups):
    """
    Fonction qui récupère les protéines de chaque groupe orthologue dans ortholog_dict.

    Parameters:
    - ortholog_dict : dict : Dictionnaire contenant les groupes orthologues avec leurs protéines
    - ortholog_groups : dict : Dictionnaire des résultats précédents contenant les groupes orthologues pour chaque best_hit

    Returns:
    - dict : Un dictionnaire avec les groupes orthologues et les protéines qui y appartiennent
    """
    # Initialisation d'un dictionnaire pour stocker les résultats
    group_proteins = {}

    # Parcours des groupes orthologues dans les résultats précédents
    for query_id, data in ortholog_groups.items():
        group_id = data["ortholog_group"]

        if group_id is not None:
            # Si un groupe est trouvé, récupère toutes les protéines de ce groupe dans ortholog_dict
            if group_id in ortholog_dict:
                group_proteins[group_id] = ortholog_dict[group_id]
            else:
                # Si le groupe n'existe pas dans ortholog_dict, on assigne une valeur vide
                group_proteins[group_id] = []

    return group_proteins

def merge_orthologous_groups_by_protein(ortholog_orthologer, orthologs_orthologer_group, 
                                        ortholog_orthofinder, orthologs_orthofinder_group):
    """
    Fusionne les groupes d'orthologues pour chaque protéine requête en combinant les résultats d'Orthologer et d'OrthoFinder.

    Parameters:
    - ortholog_orthologer : dict : Groupes orthologues pour chaque protéine requête selon Orthologer
    - orthologs_orthologer_group : dict : Protéines associées aux groupes d'Orthologer
    - ortholog_orthofinder : dict : Groupes orthologues pour chaque protéine requête selon OrthoFinder
    - orthologs_orthofinder_group : dict : Protéines associées aux groupes d'OrthoFinder

    Returns:
    - dict : Dictionnaire fusionné où chaque protéine requête est associée aux protéines orthologues uniques des deux méthodes
    """
    merged_proteins = {}

    for query_id in set(ortholog_orthologer.keys()).union(set(ortholog_orthofinder.keys())):
        merged_proteins[query_id] = set()

        # Ajout des protéines trouvées par Orthologer
        if query_id in ortholog_orthologer:
            group_id = ortholog_orthologer[query_id]["ortholog_group"]
            if group_id in orthologs_orthologer_group:
                merged_proteins[query_id].update(orthologs_orthologer_group[group_id])

        # Ajout des protéines trouvées par OrthoFinder
        if query_id in ortholog_orthofinder:
            group_id = ortholog_orthofinder[query_id]["ortholog_group"]
            if group_id in orthologs_orthofinder_group:
                merged_proteins[query_id].update(orthologs_orthofinder_group[group_id])

    # Convertir les sets en listes pour la sortie finale
    return {query_id: list(proteins) for query_id, proteins in merged_proteins.items()}

def remove_duplicate_groups(merged_proteins_by_query):
    """
    Supprime les groupes d'orthologues redondants ayant exactement les mêmes protéines.

    Parameters:
    - merged_proteins_by_query : dict : Dictionnaire où chaque protéine requête est associée à un ensemble de protéines orthologues

    Returns:
    - dict : Dictionnaire sans groupes dupliqués (groupes identiques fusionnés)
    """
    unique_groups = {}
    query_to_group = {}

    for query_id, proteins in merged_proteins_by_query.items():
        # Convertir en set immuable (clé hashable) pour identifier les groupes identiques
        proteins_set = frozenset(proteins)

        if proteins_set not in unique_groups:
            unique_groups[proteins_set] = query_id  # Associer cet ensemble à un ID de requête

        # Associer la requête à l'ID du groupe unique conservé
        query_to_group[query_id] = unique_groups[proteins_set]

    # Recréer un dictionnaire final
    final_groups = {query_id: list(proteins) for proteins, query_id in unique_groups.items()}

    return final_groups

def load_species_mapping(species_file):
    """
    Charge le fichier CSV contenant les informations sur les espèces et crée un dictionnaire.
    Le fichier CSV doit contenir les colonnes : Organism, File, Particule, BUSCO_score.
    """
    # Charger le fichier CSV
    df = pd.read_csv(species_file)
    
    # Créer un dictionnaire avec le nom de l'organisme comme clé et le code 'Particule' comme valeur
    species_dict = pd.Series(df['Organism'].values, index=df['Particule']).to_dict()

    return species_dict

def sonic_analysis(final_merged_groups, sonic_dict):
    """
    Analyse les groupes Sonicparanoid en listant les protéines associées 
    à chaque protéine requête et en affichant leur composition, y compris 
    celles qui ne sont dans aucun groupe Sonicparanoid.

    Parameters:
    - final_merged_groups : dict : Dictionnaire des groupes fusionnés avec leurs protéines
    - sonic_dict : dict : Dictionnaire des groupes Sonicparanoid et leurs protéines

    Returns:
    - dict : Dictionnaire détaillant chaque groupe Sonicparanoid avec :
        - Les protéines requêtes associées
        - Le nombre de protéines correspondantes
        - La liste complète des protéines de Sonicparanoid pour ce groupe
    """
    group_details = {}
    proteins_in_no_group = defaultdict(set)

    for query_id, proteins in final_merged_groups.items():
        for protein in proteins:
            found_in_group = False
            for group_id, sonic_proteins in sonic_dict.items():
                if protein in sonic_proteins:
                    if group_id not in group_details:
                        group_details[group_id] = {"queries": set(), "count": 0, "proteins": set()}

                    # Ajouter la protéine requête et les protéines associées
                    group_details[group_id]["queries"].add(query_id)
                    group_details[group_id]["proteins"].update(sonic_proteins)
                    group_details[group_id]["count"] += 1
                    found_in_group = True
                    break  # Sortir dès que le groupe est trouvé

            if not found_in_group:
                # Ajouter les protéines qui ne sont dans aucun groupe
                proteins_in_no_group[query_id].add(protein)

    # Convertir les ensembles en listes pour affichage
    for group_id in group_details:
        group_details[group_id]["queries"] = list(group_details[group_id]["queries"])
        group_details[group_id]["proteins"] = list(group_details[group_id]["proteins"])

    proteins_in_no_group = {k: list(v) for k, v in proteins_in_no_group.items()}

    return group_details, proteins_in_no_group

def extract_species(protein_id):
    """Extrait le tag d'espèce depuis un ID de protéine (ex: Pfalc_, TGME49_, etc.)."""
    match = re.match(r"^([A-Za-z0-9]+)_", protein_id)
    return match.group(1) if match else None

def get_sequence_from_fasta(fasta_path, protein_id):
    for record in SeqIO.parse(fasta_path, "fasta"):
        if protein_id in record.id:
            return str(record.seq)
    return None

def calculate_similarity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    return alignments / max(len(seq1), len(seq2)) if seq1 and seq2 else 0

def compute_group_similarity_score(group_proteins, query_seq, fasta_dir_species, species_df):
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
    return sum(scores)/len(scores) if scores else 0

def get_seq_fasta(fasta_dir, protein_id):
    """
    Recherche une séquence dans tous les fichiers fasta du dossier donné.
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
    os.makedirs(tree_dir, exist_ok=True)
    raw_fasta = os.path.join(tree_dir, "all_seqs_raw.fasta")
    aligned_fasta = os.path.join(tree_dir, "aligned_seqs.fasta")
    tree_path = os.path.join(tree_dir, "tree.nwk")

    # Écrire les séquences originales dans un FASTA
    seq_records = [
        SeqRecord(Seq(seq), id=f"{id}", description="")
        for id, seq in sequences.items()
        if seq is not None and seq != ""
    ]
    SeqIO.write(seq_records, raw_fasta, "fasta")

    # Étape 1 : alignement avec MAFFT
    mafft_cmd = f"mafft --auto {raw_fasta} > {aligned_fasta}"
    subprocess.run(mafft_cmd, shell=True, check=True)

    # Étape 2 : arbre avec FastTree
    fasttree_cmd = f"fasttree {aligned_fasta} > {tree_path}"
    subprocess.run(fasttree_cmd, shell=True, check=True)

    return tree_path


def compute_z_scores_and_ratios(long_branch_scores, orphans, df):
    # Séparer les scores
    non_orphan_scores = [v for k, v in long_branch_scores.items() if k not in orphans]

    # Moyenne et écart-type des non-orphelins
    mean_non_orphans = np.mean(non_orphan_scores)
    std_non_orphans = np.std(non_orphan_scores)

    # Calcul des z-scores pour toutes les protéines
    z_scores = {
        k: (v - mean_non_orphans) / std_non_orphans
        for k, v in long_branch_scores.items()
    }

    # Dictionnaire des max z-scores par espèce
    min_z_by_species = {}

    for _, row in df.iterrows():
        organism = row["Organism"]
        prefix = row["Particule"]

        # Sélectionner les protéines avec le bon préfixe
        matched = {k: z for k, z in z_scores.items() if k.startswith(prefix)}

        if matched:
            min_prot = min(matched, key=matched.get)
            min_z_by_species[organism] = matched[min_prot]

    return z_scores, min_z_by_species


def compare_species_between_groups(reference_groups, analyzed_groups, proteins_in_no_group, species_file, fasta_dir_query, fasta_dir_species, prot_interest_file, output_dir, sort_method="species_count"):
    """
    Fusionne les groupes sans espèces communes.
    sort_method : 'species_count' (par défaut) ou 'similarity'
    """
    species_df = pd.read_csv(species_file)
    query_to_code = pd.read_csv(prot_interest_file, sep=",").set_index("nom_proteine")["code_proteine"].to_dict()

    common_species = defaultdict(lambda: defaultdict(set))
    merged_no_common = defaultdict(dict)
    fusion_count = 0

    for group_id in reference_groups:
        group_info = analyzed_groups[group_id]
        query = group_info["queries"][0]
        group_species = set(extract_species(p) for p in group_info["proteins"] if extract_species(p))
        common_species[query][group_id] = group_species

    for query, ref_groups in common_species.items():
        print(f"\n▶ Protéine requête : {query}")

        # Récupérer le code_proteine à partir du nom de la protéine (query)
        protein_code = query_to_code.get(query)
        print(protein_code)
        if not protein_code:
            print(f"  ⚠️  Aucun code_proteine trouvé pour {query}")
            continue

        # Charger la séquence de la protéine requête à partir de son code_proteine
        query_seq = get_sequence_from_fasta(fasta_dir_query, protein_code)

        matching_group_ids = [
            gid for gid, info in analyzed_groups.items()
            if "queries" in info and query in info["queries"]
        ]
        already_merged = set()

        # Choisir méthode de tri
        if sort_method == "species_count":
            print("Tri par nombre d'espèces")
            sorted_matching_groups = sorted(
                matching_group_ids,
                key=lambda gid: len(analyzed_groups[gid]["proteins"]),
                reverse=True
            )
        elif sort_method == "similarity":
            print("Tri par similarité de séquence")
            similarity_scores = {
                gid: compute_group_similarity_score(analyzed_groups[gid]["proteins"], query_seq, fasta_dir_species, species_df)
                for gid in matching_group_ids  # Tri des groupes comparés (matching_group_ids)
            }
            # Tri
            sorted_matching_groups = sorted(matching_group_ids, key=lambda gid: similarity_scores[gid], reverse=True)
        else:
            sorted_matching_groups = matching_group_ids

        for ref_group_id, ref_species_set in ref_groups.items():
            print(f"  🔬 Groupe de référence {ref_group_id} ({len(ref_species_set)} espèces) :")

            current_group_id = ref_group_id
            current_proteins_by_group = defaultdict(list)
            current_species = ref_species_set.copy()

            for prot in analyzed_groups[ref_group_id]["proteins"]:
                current_proteins_by_group[ref_group_id].append(prot)

            for test_group_id in sorted_matching_groups:  # Parcours des groupes triés
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

            merged_no_common[query][current_group_id] = dict(current_proteins_by_group)

        filtered_orphans = []
        for query in proteins_in_no_group:
            if query not in merged_no_common:
                continue
            for orphan in proteins_in_no_group[query]:
                orphan_species = extract_species(orphan)

                # Récupérer toutes les espèces déjà présentes dans les groupes de ce query
                existing_species = {
                    extract_species(p)
                    for group in merged_no_common[query].values()
                    for plist in group.values()
                    for p in plist
                }

                if orphan_species not in existing_species:
                    filtered_orphans.append((query, orphan))

        print(f"Orphelins : {filtered_orphans}")

        no_groups = []

        for query in merged_no_common:
            if query not in proteins_in_no_group:
                continue

            group_seqs_dict = {
                p: get_seq_fasta(fasta_dir_species, p)
                for group in merged_no_common[query].values()
                for plist in group.values()
                for p in plist
            }

            # Sélectionner les orphelins dont l'espèce n'est pas déjà représentée
            orphans = [
                orphan for orphan in proteins_in_no_group[query]
                if extract_species(orphan) not in {
                    extract_species(p)
                    for group in merged_no_common[query].values()
                    for plist in group.values()
                    for p in plist
                }
            ]

            # Séquences des orphelins sous forme de dict {id: seq}
            orphan_seqs_dict = {
                orphan: get_seq_fasta(fasta_dir_species, orphan)
                for orphan in orphans
            }

            # Fusionner les deux dictionnaires
            combined_seqs = {**group_seqs_dict, **orphan_seqs_dict}

            output_tree = os.path.join(output_dir, "tree")
            tree = build_tree(combined_seqs, output_tree)

            '''
            t = Tree(tree)
            ts = TreeStyle()
            ts.show_leaf_name = True
            t.show(tree_style=ts)
            '''
            
            phykit_cmd = f"phykit long_branch_score {tree} -v"
            result = subprocess.run(phykit_cmd, shell=True, capture_output=True, text=True, check=True)

            #spurious_output = subprocess.check_output(["phykit", "spurious_seq", tree]).decode().strip()

            # Parser les résultats dans un dictionnaire {leaf_name: score}
            long_branch_scores = {}
            for line in result.stdout.strip().split('\n'):
                parts = line.split()
                if len(parts) == 2:
                    name, score = parts
                    long_branch_scores[name] = float(score)
            
            z_scores, min_by_species = compute_z_scores_and_ratios(long_branch_scores, orphans, species_df)

            all_orphan_ids = list(orphan_seqs_dict.keys())
            for orphan in all_orphan_ids:
                # if is_long_branch(tree, orphan, all_orphan_ids):
                #     no_groups_rejected.append(orphan)
                # else:
                #     no_groups.append(orphan)
                no_groups.append(orphan)

            merged_no_common[query][current_group_id]["orphans"] = no_groups

    return merged_no_common, min_by_species

def write_final_groups_with_sequences(nom_prot, new_groups, output_results, fasta_dir):
    """
    Écrit les groupes finaux avec leurs séquences FASTA dans un fichier.
    """
    # Créer le dossier si nécessaire
    dir_fgroups = os.path.join(output_results, "final_groups")
    os.makedirs(dir_fgroups, exist_ok=True)

    # Définir le fichier de sortie
    output_file = os.path.join(dir_fgroups, f"{nom_prot}_orthogroups.fasta")

    with open(output_file, "w") as f:
        for group, id_dict in new_groups.items():
            f.write(f">{group}\n")  # Nom du groupe (commentaire)
            for cluster_id, subcluster in id_dict.items():
                for sub_id, seq_list in subcluster.items():
                    for protein_id in seq_list:
                        sequence = get_seq_fasta(fasta_dir, protein_id)
                        if sequence:
                            f.write(f">{sub_id}|{protein_id}\n{sequence}\n\n")
                        else:
                            print(f"Séquence non trouvée pour {protein_id}")


def load_protein_mapping(csv_file):
    mapping = {}
    with open(csv_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            mapping[row["code_proteine"]] = row["nom_proteine"]
    return mapping

def get_busco_scores(csv_file):
    df = pd.read_csv(csv_file)
    return dict(zip(df['Organism'], df['BUSCO_score']))

def generate_presence_absence_table_single_group(orthologs_group, species_dict):
    """
    Génère un tableau de présence/absence avec identifiants de groupes orthologues.

    Chaque cellule contient :
    - 0 si l'espèce n'est présente dans aucun groupe
    - 1, 2, 3... selon le groupe auquel appartient l'espèce
    - n+1 pour le groupe 'no_group'

    Parameters:
    - orthologs_group : dict : {query_id: dict(fusion_id: dict(group_id: list(proteins)))}
    - species_dict : dict : {'Pfalc': 'Plasmodium falciparum', ...}

    Returns:
    - pd.DataFrame : DataFrame avec identifiants de groupe pour chaque espèce
    """
    all_species = sorted(set(species_dict.values()))
    data = []

    for query_id, fusions in orthologs_group.items():
        # Fusionner tous les groupes dans un seul dictionnaire
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
                        # on le marque à part, temporairement
                        species_to_group[species] = -1
                    else:
                        species_to_group[species] = group_index
            if group_id != "orphans":
                group_index += 1

        # Après avoir terminé, on assigne no_group à group_index (le suivant)
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
    """Fusionne les données des trois sources en associant chaque ID à ses groupes orthologues."""
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
    return {mapping.get(k, k): v for k, v in orthogroups.items()}

def generate_presence_absence_table(orthologs_sonic_group, orthologs_orthologer_group, orthologs_orthofinder_group, species_dict, merged_ortholog_groups):
    """
    Génère un tableau de présence/absence avec les query_id d'abord inchangés, puis renommés à la fin.
    """
    all_species = set(species_dict.values())
    query_sequences = set(orthologs_sonic_group.keys()).union(set(orthologs_orthologer_group.keys()), set(orthologs_orthofinder_group.keys()))

    data = []
    
    for query_id in query_sequences:
        # Récupérer les espèces présentes pour chaque groupe
        species_in_sonic = {species_dict.get(prot[:5], None) for prot in orthologs_sonic_group.get(query_id, [])}
        species_in_orthologer = {species_dict.get(prot[:5], None) for prot in orthologs_orthologer_group.get(query_id, [])}
        species_in_orthofinder = {species_dict.get(prot[:5], None) for prot in orthologs_orthofinder_group.get(query_id, [])}
        
        species_in_sonic.discard(None)
        species_in_orthologer.discard(None)
        species_in_orthofinder.discard(None)
        
        row = [query_id]  # Garde d'abord le query_id original
        
        for species in all_species:
            in_combined = 1 if (species in species_in_sonic or species in species_in_orthologer or species in species_in_orthofinder) else 0
            row.append(in_combined)

        data.append(row)

    # Création du DataFrame
    columns = ['Query Protein'] + list(all_species)
    df = pd.DataFrame(data, columns=columns)

    # Renommage des query_id à la fin
    renamed_queries = {}
    for query_id, methods in merged_ortholog_groups.items():
        for method, group_id in methods.items():
            new_query_name = f"{query_id}_{method}"
            renamed_queries[group_id] = new_query_name

    # Appliquer le renommage sur la colonne "Query Protein"
    df['Query Protein'] = df['Query Protein'].map(renamed_queries).fillna(df['Query Protein'])
    return df

def sort_key(query_protein):
    parts = query_protein.split("_")
    protein_id = parts[0] + "_" + parts[1]
    method = parts[2] if len(parts) > 2 else ""
    return (protein_id, method)


def run_table(fasta_dir, output_dir, tree_file, prot_interest_file, species_file):

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

            fasta_interest_dir = os.path.join(output_results, f"{nom_prot}.fasta")
            fetch_sequence_ncbi_single(protein, nom_prot, fasta_interest_dir)

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

                merged_proteins_by_query = merge_orthologous_groups_by_protein(
                    ortholog_orthologer, orthologs_orthologer_group, ortholog_orthofinder, orthologs_orthofinder_group
                )
                final_merged_groups = remove_duplicate_groups(merged_proteins_by_query)

                sonic_group_analysis, no_group_proteins = sonic_analysis(final_merged_groups, sonic_dict)
                new_groups, max_z_by_species = compare_species_between_groups(orthologs_sonic_group, sonic_group_analysis, no_group_proteins, species_file, fasta_interest_dir, fasta_dir, prot_interest_file, output_results, "species_count")

                '''orthogroup output'''
                '''

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

                '''
                '''fin test'''

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

                '''PB CSUI'''
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
                    #plt.show()
                    #plt.close(fig)

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

        node.img_style["vt_line_width"] = 30  # Épaisseur de la branche verticale
        node.img_style["hz_line_width"] = 30  # Épaisseur de la branche horizontale
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
    #tree.render(temp_file.name, w=400, h =800, tree_style=tree_style)
    tree.render(temp_file.name, w=2000, dpi=600, tree_style=tree_style)
    img = Image.open(temp_file.name)
    img = img.crop(img.getbbox())
    axes_final[0].imshow(img)
    axes_final[0].axis("off")

    sns.heatmap(busco_df, cmap="RdYlGn", annot=True, fmt="d", cbar=False, ax=axes_final[1], linewidths=0.3, linecolor="black", annot_kws={'size': 10})
    axes_final[1].set_xticks([0])
    axes_final[1].set_yticklabels([])
    axes_final[1].set_ylabel("")
    axes_final[1].set_xticklabels(["BUSCO Score"], rotation=45)
    axes_final[1].tick_params(axis='both', length=0)

    n_rows, n_cols = zscore_df.shape
    axes_final[2].set_xlim(0, n_cols)
    axes_final[2].set_ylim(0, n_rows)

    axes_final[2].set_xticks([x + 0.5 for x in range(n_cols)])
    axes_final[2].set_xticklabels(zscore_df.columns, fontsize=12, rotation=90)
#    axes_final[2].tick_params(axis='x', length=0, pad=8)
    axes_final[2].xaxis.tick_top()

    # for label in axes_final[2].get_xticklabels():
    #     label.set_bbox(dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.4'))

    box = axes_final[2].get_position()
    axes_final[2].set_position([box.x0, box.y0, box.width, box.height])

    # Couleurs de fond
    axes_final[2].set_facecolor('white')

    list_max_group = []

    # Ajouter des formes selon les z-scores
    for j, species in enumerate(zscore_df.columns):

        fusion_values = [
            dict_esp_fusion.get(species, {}).get(protein, {}).get("fusion", 0)
            for protein in zscore_df.index
        ]
        max_group_id = max(fusion_values)
        list_max_group.append(max_group_id)

        intermediate_ids = list(range(2, max_group_id))
        fixed_colors = {
            0: "#FFFFFF",
            1: "#228B22",
            "intermediate": sns.color_palette("husl", len(intermediate_ids)),
            max_group_id: "#e31a1c"
        }

        for i, protein in enumerate(zscore_df.index):
            info = dict_esp_fusion.get(species, {}).get(protein, {})
            fusion = info.get("fusion", 0)
            z = info.get("z_score", None)
            center = (j + 0.5, n_rows - i - 0.5)

            # 3. Couleur de bordure selon fusion
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
                    edge_color = "#000000"  # fallback noir si problème

            if pd.isna(z):
                # Rond blanc vide
                circle = Circle(center, 0.3, facecolor='none', edgecolor='black')
                axes_final[2].add_patch(circle)

            elif z < 2:
                # Rond noir plein
                circle = Circle(center, 0.3, facecolor=edge_color, edgecolor=edge_color)
                axes_final[2].add_patch(circle)

            elif 2 <= z < 3:
                # Cercle vide avec bordure
                circle = Circle(center, 0.3, edgecolor=edge_color, facecolor='none', linewidth=1)
                axes_final[2].add_patch(circle)

                # Remplissage supérieur avec un wedge noir
                wedge = Wedge(center, 0.3, 270, 90, facecolor=edge_color, edgecolor='none')
                axes_final[2].add_patch(wedge)

            elif z >= 3:
                # Quart de cercle (haut gauche rempli)
                circle = Circle(center, 0.3, edgecolor=edge_color, facecolor='none', linewidth=1)
                axes_final[2].add_patch(circle)

                wedge = Wedge(center, 0.3, 0, 90, facecolor=edge_color, edgecolor='none')
                axes_final[2].add_patch(wedge)

    # Quadrillage
    for x in range(n_cols + 1):
        axes_final[2].axvline(x, color='black', lw=0.5)
    for y in range(n_rows + 1):
        axes_final[2].axhline(y, color='black', lw=0.5)

    axes_final[2].set_aspect('equal')  # Assure des cases carrées
    
    axes_final[2].set_yticks([])
    axes_final[2].tick_params(which="minor", length=0)
    axes_final[2].set_xlabel("")

    global_max_group_id = max(list_max_group)
    intermediate_ids = list(range(2, global_max_group_id))

    # Palette cohérente avec le nombre de groupes intermédiaires
    palette_inter = sns.color_palette("husl", len(intermediate_ids))

    # Créer la légende
    legend_elements = [
        Patch(facecolor="#FFFFFF", edgecolor='black', label="Absence"),
        Patch(facecolor="#228B22", edgecolor='black', label="Reference group"),
    ]

    # Ajouter les groupes intermédiaires
    for idx, color in enumerate(palette_inter):
        legend_elements.append(Patch(facecolor=color, edgecolor='black', label=f"Group {idx + 1}"))

    # Groupe maximal (No_group)
    legend_elements.append(Patch(facecolor="#e31a1c", edgecolor='black', label="Orphans"))

    # Ajouter la légende à axes_final[2]
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
