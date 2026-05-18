# Thesis Analysis and Reproducibility Repository

This repository is a complementary resource to [Gen2PhenSim](https://github.com/genodeco/Phen-sim), containing simulation configurations, run outputs, and downstream analysis results.

## Overview

This repository stores:
- Configuration files and run outputs from Gen2PhenSim simulator used in the thesis
- Results from Regenie v3.4 genome-wide association study (GWAS) analysis
- Data extraction and processing scripts for analyzing simulation results

## Repository Structure

- `extract_data.py`          Main script for extracting and processing results
- `README.md`                This file
- `data/`
   - `configs/`              Gen2PhenSim configuration files for all runs
      - `*.json`             Individual configuration files (various scenarios and seeds)
   - `runs/`                 Output directories from Gen2PhenSim simulations
- `data_for_regenie/`         GWAS data processed for Regenie analysis
   - `data_filtered.*`         Binary PLINK format files (filtered dataset)
   - `data_step2_filtered.* `   Binary BGEN format files (step 2 analysis)
   - `data.fam`              Family file
- `results/`                 Analysis outputs
   - `extracted_data.csv`    Processed data from extract_data.py
   - `manhattan_plots/`      Visualization outputs
   - `regenie/`              Regenie GWAS results


## Requirements

- Python 3.8+
- `numpy` - Numerical computing
- `pandas` - Data manipulation and analysis
- `matplotlib` (optional) - For visualization
- Regenie v3.4+ - For GWAS analysis (if running GWAS)

## Usage

### Extracting Data

Run the data extraction script(`extract_data.py `) to process simulation results.

The script will:
- Read simulation outputs from `data/runs/`
- Process configuration files from `data/configs/`
- Generate `results/extracted_data.csv` with processed results

### GWAS Analysis

Regenie analysis results are stored in `results/regenie/`. Processed data files for Regenie are located in `data_for_regenie/`.
Phenotype files are stored in `data/runs/`

## Configuration Files

Configuration files in `data/configs/` follow the naming convention:
- `{model}_config.json` - Base configuration (e.g., `addit2_config.json`, `b_lin_config.json`)
- `{model}_ld_config.json` - Configuration with linkage disequilibrium
- `{model}_seed{N}_config.json` - Configuration with specific random seed


## Outputs

Key output files:
- `results/extracted_data.csv` - Compiled results from all simulation runs
- `results/manhattan_plots/` - Manhattan plots and other visualizations
- `results/regenie/` - GWAS results including association statistics and p-values


