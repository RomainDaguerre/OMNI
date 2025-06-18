"""
Orthologs groups Processing Pipeline

This script run 3 orthologs groups methods :
    - SonicParanoid
    - Orthologer
    - Orthofinder

Dependencies:
- SonicParanoid must be installed and available in the system PATH (https://gitlab.com/salvo981/sonicparanoid2).
- Orthologer must be installed and available in the system PATH (https://github.com/drostlab/orthologr).
- OrthoFinder must be installed and available in the system PATH (https://github.com/davidemms/OrthoFinder).


Author(s): Romain DAGUERRE

"""

import subprocess
import time
import os


def run_sonicparanoid(input_dir, output_dir, threads=8):
    """
    Runs SonicParanoid on a directory containing FASTA files.

    Parameters
    ----------
    input_dir : str
        Path to the directory containing input FASTA files.
    output_dir : str
        Path to store SonicParanoid output.
    threads : int, optional
        Number of threads to use (default: 8).
    """

    command = [
        "sonicparanoid",
        "-i",
        input_dir,
        "-o",
        output_dir,
        "-t",
        str(threads),
        "--no-compress",
        "-dmnd",
        "fast",
        "--min-bitscore",
        "50",
        "--overwrite",
    ]

    start_time = time.perf_counter()

    try:
        subprocess.run(command, check=True)
        elapsed_time = time.perf_counter() - start_time
        print(f"SonicParanoid terminé avec succès en {elapsed_time:.2f} secondes !")
    except subprocess.CalledProcessError as e:
        elapsed_time = time.perf_counter() - start_time
        print(
            f"Erreur lors de l'exécution de SonicParanoid après {elapsed_time:.2f} secondes : {e}"
        )


def run_orthologer(work_dir, input_dir, project_name):
    """
    Runs Orthologer on a directory containing FASTA files.

    Parameters
    ----------
    work_dir : str
        Working directory where Orthologer will store intermediate and result files.
    input_dir : str
        Directory containing input FASTA files.
    project_name : str
        Name of the Orthologer project.
    """

    if not os.path.exists(work_dir):
        print(f"Création du dossier : {work_dir}")
        os.makedirs(work_dir, exist_ok=True)

    create_cmd = ["orthologer", "-c", "create", "-w", work_dir, "-p", project_name]
    subprocess.run(create_cmd, check=True)

    import_cmd = [
        "orthologer",
        "-c",
        "import",
        "-w",
        work_dir,
        "-d",
        input_dir,
        "-p",
        project_name,
    ]
    subprocess.run(import_cmd, check=True)

    run_cmd = ["orthologer", "-c", "run", "-w", work_dir, "-p", project_name]
    subprocess.run(run_cmd, check=True)


def run_orthofinder(input_dir, output_dir):
    """
    Runs OrthoFinder on a directory containing FASTA files.

    Parameters
    ----------
    input_dir : str
        Path to the directory with input FASTA files.
    output_dir : str
        Path to store OrthoFinder output.
    """

    start_time = time.time()

    command = ["orthofinder", "-f", input_dir, "-o", output_dir, "-t", str(4), "-og"]
    subprocess.run(command)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Le temps de calcul est : {elapsed_time:.2f} secondes")


def run_orthologs(input_dir, output_sonic, output_orthologer, output_orthofinder):
    """
    Executes the full orthology inference pipeline using SonicParanoid, Orthologer, and OrthoFinder.

    Parameters
    ----------
    input_dir : str
        Directory containing FASTA files.
    output_sonic : str
        Directory to store SonicParanoid output.
    output_orthologer : str
        Directory to store Orthologer output.
    output_orthofinder : str
        Directory to store OrthoFinder output.
    """

    run_sonicparanoid(input_dir, output_sonic)
    run_orthologer(output_orthologer, input_dir, "project")
    run_orthofinder(input_dir, output_orthofinder)
