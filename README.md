# OMNI (Orthologous group Merger for Numerous Inference tools)

<p align="center">
  <img alt="Made with Python" src="https://img.shields.io/badge/Made%20with-Python-1f425f.svg?color=%23539fc9">
</p>

## Repository Structure

```
.
├── Fasta/                           # Raw FASTA files of proteomes from different species
├── Output/                          # Output files from orthogroup
├── fasta_protein/                   # Request proteins sequences in FASTA format
├── particule_species_Busco.csv      # Species data
├── phylogenie_apicomplexa.nwk       # Phylogenetic tree in Newick format
├── protein.csv                      # Request proteins data
```

## Requirements

- Python 3.8+ 
- [EMBOSS](https://emboss.sourceforge.net/)
- [BUSCO](https://busco.ezlab.org/)
- [SonicParanoid](https://github.com/fenderglass/SonicParanoid)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [Orthologer](https://github.com/drostlab/orthologr)
- `pandas`, `argparse`, `os`, `glob`, `shutil`, etc.

---

## Usage

### 1. Orthogroups Processing Pipeline

This pipeline translates sequences, cleans headers, runs BUSCO, and infers orthogroups using multiple tools.

```bash
python3 main.py \
    -fd ./Fasta \
    -db /path/to/busco_db \
    -o ./Output \
    -s particule_species_Busco.csv
```

**Additional options :**
| Option | Description |
|--------|-------------|
| `-t`   | If proteome type not defined |
| `-b`   | Skip BUSCO analysis |
| `-ortho` | Skip ortholog inference and analysis |

---

## 2. 🧠 OMNI Tool (Ontology-based Metabolic Network Integration)

This tool serves as a command-line interface for analyzing orthologous protein groups across species using precomputed FASTA files and orthology inference results :

- Merging and comparing orthologous groups across multiple inference tools.
- Dynamically selecting the most compact reference group for species comparison.
- Evaluating the consistency of orthology predictions between tools.
- Optionally integrating a phylogenetic tree for enhanced analysis.

```bash
python3 main.py \
    -fd ./Fasta \
    -o ./Output \
    -i ./interest_proteins.csv \
    -s ./particule_species_Busco.csv \
    -t ./tree.nwk \
    -pf ./combined_protein.fasta
```

**Additional options :**
| Option | Description |
|--------|-------------|
| `-t`   | Path to a phylogenetic tree in Newick format |
| `-pf`   | Path to a combined FASTA file with all protein sequences |

---

## Author

**Romain Daguerre**
UMR 7245 - CNRS MNHM - Parasites et Protistes Libres
Master Bioinformatique — Université Paris Cité


---

## Contact

For any questions, suggestions, or contributions :  
📧 romain.daguerre1119 [at] gmail.com
