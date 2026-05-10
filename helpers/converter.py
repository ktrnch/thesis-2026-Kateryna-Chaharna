import json
from pathlib import Path

CONFIG = dict(
    output_prefix="b_lin",
    seed=7,  # same seed = same weights for the same param.
    phenotype_type="binary",  # "quantitative" or "binary"
    phenotype_number=1,
    n_indiv_max=None,  # keep first N individuals for memory friendliness (None = use all)

    causal_var_pool_num=100,  # pre selected SNPs
    causal_var_min=100,  # min num of SNPs forming phenotype
    causal_var_max=100,
    causal_snps=None,  # OR comma-separated snp names string e. g. SNP166, SNPn

    # LD mode
    ld_mode= False,
    ld_threshold=1000,
    chromosome_file="data/10K_SNP_1000G_real.snpinfo.tsv",

    second_inter_num=0,
    third_inter_num=0,
    dom_num=0,
    rec_num=0,
    int_var_max=0,

    # interaction knobs (ignored if no interactions)
    interactions_mode="pure",  # "additive" or "pure"
    separate_interaction_weights=False,
    interaction_weight_scale=1.0,
    force_defined_interactions=False,

    # heritability / noise (for quantitative)
    target_h2=None,
    noise_factor=0.0,  # ignored if target_h2 is set

    # output mapping phenotype
    quant_range_min=0.0,
    quant_range_max=1.0,
    percentile_or_threshold=80.0,  # for binary only
)

prefix = CONFIG["output_prefix"]
config_path = Path(f"configs_{prefix}.json")

with open(config_path, "w") as configs:
    json.dump(CONFIG, configs, indent=2)

print(f"[done] wrote {config_path}")