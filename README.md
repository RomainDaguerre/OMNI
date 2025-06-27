# OMNI (Orthologous group Merger for Numerous Inference tools)


## Repository Structure

```
.
├── Fasta/                           # Raw FASTA files used in the analysis
├── Output/                          # Output files from orthogroup and phylogeny analyses
├── fasta_protein/                   # Protein sequences in FASTA format
├── particule_species_Busco.csv      # BUSCO results for different species
├── phylogenie_apicomplexa.nwk       # Phylogenetic tree in Newick format
├── prot_test.csv                    # Protein test data
```

## Requirements

- Python 3.8+<p align="center"><img alt="Made with Python" src="https://img.shields.io/badge/Made%20with-Python-1f425f.svg?color=%23539fc9"></p>
- EMBOSS (transeq dans le `$PATH`)
- [BUSCO](https://busco.ezlab.org/)
- [SonicParanoid](https://github.com/fenderglass/SonicParanoid)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- Orthologer
- `pandas`, `argparse`, `os`, `glob`, `shutil`, etc.

---

## 👤 Author

**Romain Daguerre**
UMR 7245 - CNRS MNHM - Parasites et Protistes Libres
Master Bioinformatique — Université Paris Cité


---
