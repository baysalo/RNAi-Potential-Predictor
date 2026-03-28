# RNAi-Potential-Predictor
RNAi Potential Predictor: Bacterial-Fungal Cross-Kingdom Silencing
# Bacterial-Fungal RNAi Potential Predictor 🧬

This repository contains an R-based bioinformatics pipeline to identify **cross-kingdom RNAi triggers** from bacterial genomes targeting fungal pathogens. It identifies 21-nt sequences with high fungal complementarity and zero host off-target risk.

## 🔬 How it Works
The pipeline simulates the natural process of RNA interference (RNAi) where small RNAs (siRNAs) guide the RISC complex to degrade specific mRNA targets in a pathogen.

1.  **In Silico Dicing:** Fragments bacterial genes into 21-mer virtual siRNAs.
2.  **Pathogen Mapping:** Identifies siRNAs with 100% complementarity to the target fungal transcriptome.
3.  **Host Filtering:** Discards any sequences with homology (up to 2 mismatches) to the host genome (e.g., the plant) to ensure biosecurity.
4.  **Target Annotation:** Maps the final safe siRNAs back to specific fungal transcript IDs.

## 🛠 Prerequisites
Install the required Bioconductor packages in R:
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
