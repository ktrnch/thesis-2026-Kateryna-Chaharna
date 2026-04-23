import json

CONFIG = dict(
    output_prefix              = "lin3_h2p20_q01",
    seed                       = 7, #same seed = same weights for the same param.
    phenotype_type             = "quantitative",  #"quantitative" or "binary"
    phenotype_number           = 3,
    n_indiv_max                = 3,  #keep first N individuals for memory friendliness (None = use all)


    causal_var_pool_num        = 100, #pre selected SNPs
    causal_var_min             = 100, # min num of SNPs forming phenotype
    causal_var_max             = 100,
    causal_snps                = None,  # OR comma-separated snp names string e. g. SNP166, SNPn
    ##how to make it more controlable for pleiotropy

    #LD mode
    ld_mode = True,
    ld_threshold = 1000,
    chromosome_file = "10K_SNP_1000g_real.snpinfo.tsv",

    second_inter_num           = 0,    # <-- no interactions
    third_inter_num            = 0,
    dom_num                    = 0,
    rec_num                    = 0,
    int_var_max                = 0,

    # interaction knobs (ignored if no interactions)
    interactions_mode          = "additive",     # "additive" or "pure" #check
    separate_interaction_weights = False, # check for how it is derived
    interaction_weight_scale   = 1.0,
    force_defined_interactions = False,

    # heritability / noise (for quantitative) #check how it influences
    target_h2                  = 0.2,  # exact on observed scale
    noise_factor               = 0.0,  # ignored if target_h2 is set

    # output mapping phenotype # give eqaitions
    quant_range_min            = 0.0,
    quant_range_max            = 1.0,
    percentile_or_threshold    = 80.0, # for binary only
)

with open("configs.json", "w") as linear_only:
    json.dump(CONFIG, linear_only, indent=2)

