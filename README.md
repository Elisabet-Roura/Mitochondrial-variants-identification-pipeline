# Mitochondrial-variants-identification-pipeline
--- 
This repository contains the code developed as part of my Master’s Thesis (TFM):  
**“Design of a Pipeline for the Identification of Mitochondrial DNA Variants in a Sudden Cardiac Death Cohort”**  
(*“Diseño de una Pipeline para la Identificación de Variantes en ADN Mitocondrial de una Cohorte de Muerte Súbita Cardíaca”*).  

Author: Elisabet Roura 
Master’s Thesis – Universidad Internacional de la Rioja  
Academic Year: 2025 


## Overview

This repository contains:

- **`mtdna_pipeline.py`** – Main Python script to detect mtDNA variants.
- **`download_cohort_files.py`** – Utility script to download cohort sample files form database for downstream analysis.
- **`Simulations/`** – Directory containing scripts to generate synthetic datasets used to evaluate pipeline performance across different sequencing coverages and heteroplasmy proportions.
- **`Modules/`** – Supporting modules for core pipeline functionality.
- **`Modules/Annotation_files`** –  Annotation files adapted from the [mitoHPC](https://github.com/mitoNGS/mitoHPC) repository, used to annotate the final VCF files.
- **`LICENSE`** – GPL-3.0 license governing usage and distribution.


## Acknowledgements

Some resources and annotation files included in this project were adapted from the [mitoHPC repository](https://github.com/mitoNGS/mitoHPC). I would like to acknowledge the authors of mitoHPC for making their work publicly available and reusable for academic purposes. 
