
# this file conteins functions to andlyse GWAS performance
# This function calculates which SNPs are causal for chosen phenotype 
# based on step2 output and the threshold for p-value, and then saves them in a file.
# This function compares them with the ones that are in the file with true causal SNPs

import pandas as pd

def extract_p_values(step2_output):
    # Read the step2 output file
    df = pd.read_csv(step2_output, sep=' ')
    #print(df.head())  # Print the first few rows to check the structure of the data
    
    # Extract the SNPs and their corresponding log10p
    snp_log10p_values = df[['ID', 'LOG10P']].copy()
    # Convert log10p to p-values
    snp_log10p_values['P_VALUE'] = 10 ** (-snp_log10p_values['LOG10P'])
    snp_p_values = snp_log10p_values[['ID', 'P_VALUE']]

    return snp_p_values

def identify_causal_snps(snp_p_values, p_value_threshold):
    # Identify causal SNPs based on the p-value threshold
    causal_snps = snp_p_values[snp_p_values['P_VALUE'] < p_value_threshold]['ID'].tolist()
    return causal_snps

def compare_causal_snps(identified_causal_snps, true_causal_snps_file):
    # Read the true causal SNPs from the file
    true_causal_snps_df = pd.read_csv(true_causal_snps_file, sep='\t', header=None)
    # First row is list of causal SNPs without the first column, so we take the first row and convert it to a list
    true_causal_snps = true_causal_snps_df.iloc[0, 1:].dropna().tolist()  # Skip the first column which is not a SNP, drop NaN for empty
    #print(f"True Causal SNPs: {true_causal_snps}")  # Print the true causal SNPs for verification 
    # Compare identified causal SNPs with true causal SNPs
    true_positives = set(identified_causal_snps) & set(true_causal_snps)
    false_positives = set(identified_causal_snps) - set(true_causal_snps)
    false_negatives = set(true_causal_snps) - set(identified_causal_snps)

    return true_positives, false_positives, false_negatives

def extract_prefix_from_filename(filename):
    prefix = filename.split('.')[0] # Remove the file extension
    prefix = prefix.split('/')[-1] # Get the last part of the path, which is the filename without extension
    return prefix


# Lists of paths to the step2 output files
step2_linear_outputs = ['results/regenie/b_lin_ld_step2_Phenotype_1.regenie', 'results/regenie/b_lin_step2_Phenotype_1.regenie',
                 'results/regenie/q_lin_ld_step2_Phenotype_1.regenie', 'results/regenie/q_lin_step2_Phenotype_1.regenie' ]  

step2_additive_outputs = ['results/regenie/addit2/addit2_step2_Phenotype_1.regenie',
                          'results/regenie/addit2_ld/addit2_ld_step2_Phenotype_1.regenie', 
                          'results/regenie/addit2nd3/addit2nd3_step2_Phenotype_1.regenie',
                          'results/regenie/addit2nd3_ld/addit2nd3_ld_step2_Phenotype_1.regenie']

p_value_threshold = 0.05  # Threshold for identifying causal SNPs
#list of paths to the files containing true causal SNPs for each phenotype, in the same order as the step2 output files

true_causal_snps_files = ['data/runs/b_lin_ld/b_lin_ld_SNP_ASSIGNMENTS.txt', 'data/runs/b_lin/b_lin_SNP_ASSIGNMENTS.txt',
                          'data/runs/q_lin_ld/q_lin_ld_SNP_ASSIGNMENTS.txt', 'data/runs/q_lin/q_lin_SNP_ASSIGNMENTS.txt',
                          'data/runs/addit2/addit2_SNP_ASSIGNMENTS.txt', 'data/runs/addit2_ld/addit2_ld_SNP_ASSIGNMENTS.txt',
                          'data/runs/addit2nd3/addit2nd3_SNP_ASSIGNMENTS.txt', 'data/runs/addit2nd3_ld/addit2nd3_ld_SNP_ASSIGNMENTS.txt']


def main():
    results_list  = []
    for step2_output, true_causal_snps_file in zip(step2_linear_outputs + step2_additive_outputs, true_causal_snps_files):
        print(f"Processing file: {step2_output}")
        try:
            snp_p_values = extract_p_values(step2_output)
            identified_causal_snps = identify_causal_snps(snp_p_values, p_value_threshold)
            true_positives, false_positives, false_negatives = compare_causal_snps(identified_causal_snps, true_causal_snps_file)
            prefix = extract_prefix_from_filename(step2_output)
            results_list.append({ 
                'Conditions': prefix,
                "True Positives" :len(true_positives),
                "False Positives": len(false_positives), 
                "False Negatives": len(false_negatives),
                "True Causal SNPs": len(true_positives) + len(false_negatives)   
            })

        except FileNotFoundError as e:
            print(f"File not found: {e}")
            print("Skipping this file.\n")
        except Exception as e:
            print(f"Error processing {step2_output}: {e}")
            print("Skipping this file.\n")
    # Create a DataFrame from the results list and save it to a CSV file
    results_df = pd.DataFrame(results_list)
    results_df.to_csv('results/gwas_performance.csv', index=False)


if __name__ == "__main__":
    main()