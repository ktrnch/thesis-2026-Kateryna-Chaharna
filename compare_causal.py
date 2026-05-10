
# this file conteins functions to andlyse GWAS performance
# This function calculates which SNPs are causal for chosen phenotype 
# based on step2 output and the threshold for p-value, and then saves them in a file.
# This function compares them with the ones that are in the file with true causal SNPs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def extract_data_from_files(file_names: list, p_value_threshold=5e-8):
    step2_linear_outputs = []
    true_causal_snps = []
    for file_name in file_names:
        if file_name.startswith('results/regenie'):
            df = pd.read_csv(file_name, sep=r'\s+', engine='python')
            step2_linear_outputs.append(df)
        elif file_name.startswith('data/runs'):
            df = pd.read_csv(file_name)
            if 'Selected_SNPs' in df.columns:
                selected = df['Selected_SNPs'].dropna().astype(str).tolist()
                true_causal_snps.extend([snp.strip() for row in selected for snp in str(row).split(',') if snp.strip()])
            else:
                true_causal_snps.extend(df.iloc[:, 0].dropna().astype(str).tolist())
    snp_log10p_values = step2_linear_outputs[0][['ID', 'LOG10P']].copy()
    snp_log10p_values['P_VALUE'] = 10 ** (-snp_log10p_values['LOG10P'])
    snp_p_values = snp_log10p_values[['ID', 'P_VALUE']]
    identified_causal_snps = snp_p_values[snp_p_values['P_VALUE'] < p_value_threshold]['ID'].tolist()
    total_num_snps = len(snp_p_values)
    return true_causal_snps, identified_causal_snps, snp_log10p_values, snp_p_values, total_num_snps


def make_manchetan_plot(output_path, true_causal_snps, snp_log10p_values, p_value_threshold=5e-8):
    df = snp_log10p_values.copy()
    df['-log10(P_VALUE)'] = df['LOG10P']
    df['color'] = df['ID'].apply(lambda x: 'red' if x in true_causal_snps else 'blue')
    plt.figure(figsize=(10, 6))
    plt.scatter(df['ID'], df['-log10(P_VALUE)'], color=df['color'], s=10)
    plt.axhline(y=-np.log10(p_value_threshold), color='red', linestyle='--')  # Add a horizontal line for the significance threshold
    plt.xlabel('Chromosome 15')
    plt.ylabel('-log10(P_VALUE)')
    plt.title('Manhattan Plot')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close() 



def compare_causal_snps(identified_causal_snps, true_causal_snps, total_num_snps):
    # Compare identified causal SNPs with true causal SNPs
    true_positives = set(identified_causal_snps) & set(true_causal_snps)  
    false_positives = set(identified_causal_snps) - set(true_causal_snps)
    false_negatives = set(true_causal_snps) - set(identified_causal_snps)
    true_negatives = total_num_snps - len(true_positives) - len(false_positives) - len(false_negatives)  # Total SNPs minus the ones we identified as causal or missed


    return true_positives, false_positives, false_negatives, true_negatives

def extract_prefix_from_filename(filename):
    prefix = filename.split('.')[0] # Remove the file extension
    prefix = prefix.split('/')[-1] # Get the last part of the path, which is the filename without extension
    prefix = prefix.replace('_step2_Phenotype_1', '') # Remove the specific suffix to get the condition name
    return prefix


# Lists of paths to the step2 output files
step2_linear_outputs = ['results/regenie/b_lin_ld_step2_Phenotype_1.regenie', 'results/regenie/b_lin_step2_Phenotype_1.regenie',
                 'results/regenie/q_lin_ld_step2_Phenotype_1.regenie', 'results/regenie/q_lin_step2_Phenotype_1.regenie',
                  "results/regenie/addit2_step2_Phenotype_1.regenie", "results/regenie/addit2_ld_step2_Phenotype_1.regenie" ]  


p_value_threshold = 5e-8  # Threshold for identifying causal SNPs
#list of paths to the files containing true causal SNPs for each phenotype, in the same order as the step2 output files

true_causal_snps_files = ['data/runs/b_lin_ld/phenotype_snps.csv',
                          'data/runs/b_lin/phenotype_snps.csv','data/runs/q_lin_ld/phenotype_snps.csv','data/runs/q_lin/phenotype_snps.csv',
                          'data/runs/addit2/phenotype_snps.csv','data/runs/addit2_ld/phenotype_snps.csv']
def main():
    results_list  = []
    for step2_output, true_causal_snps_file in zip(step2_linear_outputs, true_causal_snps_files):
        print(f"Processing file: {step2_output}")
        try:
            true_causal_snps, identified_causal_snps, snp_log10p_values, snp_p_values, total_num_snps = extract_data_from_files([step2_output, true_causal_snps_file])
            true_positives, false_positives, false_negatives, true_negatives = compare_causal_snps(identified_causal_snps, true_causal_snps, total_num_snps)
            prefix = extract_prefix_from_filename(step2_output)
            make_manchetan_plot(f"results/manchetan_plots/{prefix}_manchetan_plot.png", true_causal_snps, snp_log10p_values)
            results_list.append({ 
                'Conditions': prefix,
                "True Positives" :len(true_positives),
                "False Positives": len(false_positives),
                "True Negatives": true_negatives, 
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