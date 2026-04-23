# Genotype-phenotype simulator "GenPhenSim" 

The code takes PLINK.raw file as indut and generates phenotypes based on it

## Overwiew

This is genotype to phenotype simulator, writen in prthon.
User can set
    - Interaction mode (aditive or pure)
    - LD mode
    - targetheritability
    - output type (binary/quantitative)

## Repository structure

- simulator.py - main simulation script.
- converter.py - builds configs.json.
- configs.json - runtime configuration.

- input data files:
    * .raw - esential
    Structure: first 6 PLINK metadata columns + SNP columns.
    * .snpinfo.tsv - for LD mode
    Structure: SNP_id chromosome, position

## Reqirements

### Python 3
### Packages
    - numpy
    - pandas

## Configurations
User can create config file using converter. User can set following:
    - core run setings (output, prefics, seed, num of individuals)
    - casual SNP pool setings
    - LD threshold, if LD mode is set
    - interactions
    - heritability
    - Output mapping

    ### Common pitfalls
    - Causal_var_max is more or eqal causal_var_pool_num. 
    - Causal_var_min s more or eqal causal_var_max.
    - The force_defined_interactions are not set, but total number of interactions is less than maximal number of interactions
    - heritability is set, but phenotype is not "qantitative"
    - The noice factor is not zero, but heritability is set

## How to run?

## Outputs
        phenotype_without_noise.txt
        phenotype_with_noise.txt
        final_phenotype.txt
        <prefix>_h2_report.tsv
        <prefix>_weights*.tsv
        <prefix>_SNP_ASSIGNMENTS.txt
        <prefix>_meta.json

#Reproducibility

#Limitations

#How to cite this simulator






