# RNA-Seq Analysis of *Staphylococcus aureus* Persisters

**Goal:** Reproduce the RNA-Seq differential expression analysis from
 [Nature Communications (2020)](https://www.nature.com/articles/s41467-020-15966-7)
*“Intracellular Staphylococcus aureus persisters upon antibiotic exposure.”*

---

## Project Overview
This project aims to build a **reproducible workflow (Nextflow)** using **containers (Docker/Apptainer)**
to analyze RNA-Seq data and reproduce the paper’s results.

### Biological context
- Pathogen: *Staphylococcus aureus* NCTC 8325
- Conditions: intracellular persisters vs control (3 replicates each)
- Dataset: GEO [GSE139659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139659)

### Technologies
- Workflow manager: **Nextflow DSL2**
- Containers: **Docker** or **Apptainer**
- Tools: FastQC, Cutadapt, Bowtie, FeatureCounts, DESeq2

---

## Repository Structure

rna-seq-s-aureus-persisters/
├── main.nf
├── nextflow.config
├── run.sh
├── containers/
├── conf/
├── resources/
│ ├── genome/
│ └── samplesheets/
├── modules/
├── scripts/
├── results/
└── reports/


---

## 👥 Team

- Turbot Alan
- Beaumatin Yohan
- Faramus François
- Bouteville Tristan

