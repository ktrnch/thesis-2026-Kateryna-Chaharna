# Genotype-phenotype simulator "Gen2PhenSim" 

Genotype-to-phenotype simulator from PLINK `.raw` input.

## Overview

- Supports quantitative and binary phenotypes
- Supports additive/pure interaction modes
- Optional LD-aware SNP inclusion
- Configurable heritability/noise

## Repository Structure
- `main.py` — loads config and runs simulation
- `converter.py` — generates config JSON files
- `simulator/`
  - `simulate.py` — core simulation pipeline
  - `ld.py`, `sampling.py`, `io_utils.py`, `transformer.py` — helpers
- `data/configs/` — archived configs (`<prefix>_config.json`)
- `data/runs/<prefix>/` — simulation outputs per run

## Requirements
- Python 3
- `numpy`
- `pandas`

## Configuration
Main fields include:
- run settings (`output_prefix`, `seed`, `phenotype_number`, `n_indiv_max`)
- causal SNP pool settings
- LD settings (`ld_mode`, `ld_threshold`, `chromosome_file`)
- interaction settings (`second_inter_num`, `third_inter_num`, etc.)
- heritability/noise settings (`target_h2`, `noise_factor`)
- output mapping settings

## How to Run
1. Create/update config:
   - `python3 converter.py`
2. Run simulation:
   - `python3 main.py`

## Outputs
Each run is written to `data/runs/<output_prefix>/`:
- `phenotype_without_noise.txt`
- `phenotype_with_noise.txt`
- `final_phenotype.txt`
- `phenotype_snps.csv`
- `<prefix>_SNP_ASSIGNMENTS.txt`
- `<prefix>_weights*.tsv`
- `<prefix>_h2_report.tsv` (quantitative)
- `<prefix>_meta.json`

Config snapshot:
- `data/configs/<prefix>_config.json`

# Reproducibility

Simulated data can be reprodused using the same seed number. Besides there is configs archive. It stores data provided by user for the simulation.
`data/configs/` — archived configs (`<prefix>_config.json`) can be identified based on prefix. 

# Limitations
Simulator only accepts PLINK.raw files and relies on its structure. 

# How to cite this simulator
