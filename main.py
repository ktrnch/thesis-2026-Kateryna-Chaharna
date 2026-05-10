#!/usr/bin/env python3
import json
import sys
from pathlib import Path

from simulator.simulate import simulate

RAW = "data/10K_SNP_1000G_real.raw"
#upgrade to easier config handling
CONFIG_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else Exception("Please provide a config file path as an argument")
try:
    with open(CONFIG_PATH, "r") as config_file:
        CONFIG = json.load(config_file)
except FileNotFoundError:
    print(f"Config file not found: {CONFIG_PATH}")
    sys.exit(1)
except Exception as e:
    print(f"Error loading config file: {e}")
    sys.exit(1)

if __name__ == "__main__":
    prefix = CONFIG["output_prefix"]
    output_dir = Path("data") / "runs" / prefix
    config_archive_dir = Path("data") / "configs"

    # Handle different config formats
    n_indiv_max = CONFIG.get("n_indiv_max", CONFIG.get("n_individuals", None))
    if "quant_range" in CONFIG:
        quant_range_min = CONFIG["quant_range"][0]
        quant_range_max = CONFIG["quant_range"][1]
    else:
        quant_range_min = CONFIG.get("quant_range_min", 0.0)
        quant_range_max = CONFIG.get("quant_range_max", 1.0)

    simulate(
        raw_path=RAW,
        output_prefix=prefix,
        output_dir=str(output_dir),
        config_source_path=str(CONFIG_PATH),
        config_archive_dir=str(config_archive_dir),

        seed=CONFIG["seed"],
        phenotype_type=CONFIG["phenotype_type"],
        phenotype_number=CONFIG["phenotype_number"],
        n_indiv_max=n_indiv_max,
        causal_var_pool_num=CONFIG["causal_var_pool_num"],
        causal_var_min=CONFIG["causal_var_min"],
        causal_var_max=CONFIG["causal_var_max"],
        causal_snps=CONFIG["causal_snps"],
        ld_mode=CONFIG["ld_mode"],
        ld_threshold=CONFIG["ld_threshold"],
        chromosome_file=CONFIG["chromosome_file"],
        second_inter_num=CONFIG["second_inter_num"],
        third_inter_num=CONFIG["third_inter_num"],
        dom_num=CONFIG["dom_num"],
        rec_num=CONFIG["rec_num"],
        int_var_max=CONFIG["int_var_max"],
        interactions_mode=CONFIG["interactions_mode"],
        separate_interaction_weights=CONFIG["separate_interaction_weights"],
        interaction_weight_scale=CONFIG["interaction_weight_scale"],
        force_defined_interactions=CONFIG["force_defined_interactions"],
        target_h2=CONFIG["target_h2"],
        noise_factor=CONFIG["noise_factor"],
        quant_range_min=quant_range_min,
        quant_range_max=quant_range_max,
        percentile_or_threshold=CONFIG["percentile_or_threshold"],
    )