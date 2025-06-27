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

## Author

**Romain Daguerre**
UMR 7245 - CNRS MNHM - Parasites et Protistes Libres
Master Bioinformatique — Université Paris Cité


---

## Contact

For any questions, suggestions, or contributions :  
📧 romain.daguerre1119 [at] gmail.com
