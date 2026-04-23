#!/usr/bin/env python3
#load data
import json
from simulator.simulate import simulate

RAW = "10K_SNP_1000G_real.raw"
with open("configs.json", "r") as linear_only:
    CONFIG = json.load(linear_only)

#simulate
if __name__ == "__main__":
    simulate(
    raw_path                         = RAW,
    output_prefix                    = CONFIG["output_prefix"],
    seed                             = CONFIG["seed"],
    phenotype_type                   = CONFIG["phenotype_type"],
    phenotype_number                 = CONFIG["phenotype_number"],
    n_indiv_max                      = CONFIG["n_indiv_max"],
    causal_var_pool_num              = CONFIG["causal_var_pool_num"],
    causal_var_min                   = CONFIG["causal_var_min"],
    causal_var_max                   = CONFIG["causal_var_max"],
    causal_snps                      = CONFIG["causal_snps"],
    ld_mode                          = CONFIG["ld_mode"],
    ld_threshold                     = CONFIG["ld_threshold"],
    chromosome_file                  = CONFIG["chromosome_file"],   
    second_inter_num                 = CONFIG["second_inter_num"],
    third_inter_num                  = CONFIG["third_inter_num"],
    dom_num                          = CONFIG["dom_num"],
    rec_num                          = CONFIG["rec_num"],
    int_var_max                      = CONFIG["int_var_max"],
    interactions_mode                = CONFIG["interactions_mode"],
    separate_interaction_weights     = CONFIG["separate_interaction_weights"],
    interaction_weight_scale         = CONFIG["interaction_weight_scale"],
    force_defined_interactions       = CONFIG["force_defined_interactions"],
    target_h2                        = CONFIG["target_h2"],
    noise_factor                     = CONFIG["noise_factor"],
    quant_range_min                  = CONFIG["quant_range_min"],
    quant_range_max                  = CONFIG["quant_range_max"],
    percentile_or_threshold          = CONFIG["percentile_or_threshold"],
)
