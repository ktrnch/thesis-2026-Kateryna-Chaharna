
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

#print("Extracting data from files...") # to see if  it runs
# DATA EXTRACTION PART
def file_list_from_directory(gen2phen_dir, regenie_dir):
    true_causal_snps_files = []
    step2_outputs = []
    #extract file path to phernotype_snps.csv - only for directories with _seed* suffix
    for file in Path(gen2phen_dir).rglob('*seed*/phenotype_snps.csv'):
        true_causal_snps_files.append(str(file))
    #step 2 output files - only for files with _seed* suffix
    for file in Path(regenie_dir).rglob('*seed*_step2_Phenotype_1.regenie'):
        step2_outputs.append(str(file))
    return step2_outputs, true_causal_snps_files

def extract_data_from_files(file_names: list, p_value_threshold=5e-8):
    step2_outputs = []
    true_causal_snps = []
    for file_name in file_names:
        if 'results/regenie' in file_name:
            df = pd.read_csv(file_name, sep=r'\s+', engine='python')
            step2_outputs.append(df)
        elif 'data/runs' in file_name:
            df = pd.read_csv(file_name)
            if 'Selected_SNPs' in df.columns:
                selected = df['Selected_SNPs'].dropna().astype(str).tolist()
                true_causal_snps.extend([snp.strip() for row in selected for snp in str(row).split(',') if snp.strip()])
            else:
                true_causal_snps.extend(df.iloc[:, 0].dropna().astype(str).tolist())
    snp_log10p_values = step2_outputs[0][['ID', 'LOG10P']].copy()
    snp_log10p_values['P_VALUE'] = 10 ** (-snp_log10p_values['LOG10P'])
    snp_p_values = snp_log10p_values[['ID', 'P_VALUE']]
    identified_causal_snps = snp_p_values[snp_p_values['P_VALUE'] < p_value_threshold]['ID'].tolist()
    total_num_snps = len(snp_p_values)
    return true_causal_snps, identified_causal_snps, snp_log10p_values, snp_p_values, total_num_snps

#CALCULATION PART
def compare_causal_snps(identified_causal_snps, true_causal_snps, total_num_snps):
    # Compare identified causal SNPs with true causal SNPs
    true_positives = set(identified_causal_snps) & set(true_causal_snps)  
    false_positives = set(identified_causal_snps) - set(true_causal_snps)
    false_negatives = set(true_causal_snps) - set(identified_causal_snps)
    true_negatives = total_num_snps - len(true_positives) - len(false_positives) - len(false_negatives)  # Total SNPs minus the ones we identified as causal or missed
    return {
        "True Positives": len(true_positives),
        "False Positives": len(false_positives),
        "False Negatives": len(false_negatives),
        "True Negatives": true_negatives
    }
# PLOTTING PART
def make_manchetan_plot(output_path, true_causal_snps, snp_log10p_values,condition, p_value_threshold=5e-8):
    df = snp_log10p_values.copy()
    df['-log10(P_VALUE)'] = df['LOG10P']
    df['color'] = df['ID'].apply(lambda x: 'red' if x in true_causal_snps else 'blue')
    plt.figure(figsize=(10, 6))
    plt.scatter(df['ID'], df['-log10(P_VALUE)'], color=df['color'], s=10)
    plt.axhline(y=-np.log10(p_value_threshold), color='red', linestyle='--')  # Add a horizontal line for the significance threshold
    plt.xlabel('Chromosome 15')
    plt.ylabel('-log10(P_VALUE)')
    plt.title(f'{condition} - Manhattan Plot')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()




data_list = []

try: 
    #extract prefixes from file names
    def regenie_key(path):
        name = Path(path).stem
        return name.replace('_step2_Phenotype_1', '')
    def phenotype_key(path):
        return Path(path).parent.name

    #files to work with
    step2_outputs, true_causal_snps_files = file_list_from_directory("data/runs", "results/regenie")
    # mapping outputs using dictionaries
    regenie_map = {regenie_key(f): f for f in step2_outputs}
    phenotype_map = { phenotype_key(f): f for f in true_causal_snps_files}

    keys = sorted(set(regenie_map) & set(phenotype_map))
    try:
        for k in keys:
            step2_output = regenie_map[k]
            true_causal_snps_file = phenotype_map[k]
            true_causal_snps, identified_causal_snps, snp_log10p_values, snp_p_values, total_num_snps = extract_data_from_files([step2_output, true_causal_snps_file])
            results = compare_causal_snps(identified_causal_snps, true_causal_snps, total_num_snps)
            data_list.append({
                "Prefix": k,
                "True Positives": results["True Positives"],
                "False Positives": results["False Positives"],
                "False Negatives": results["False Negatives"],
                "True Negatives": results["True Negatives"]
            })
            make_manchetan_plot(f'results/manhattan_{k}.png', true_causal_snps, snp_log10p_values, condition=k)

    except Exception as e:
        print(f"Failed on {k}: {e}")

except FileNotFoundError as e:
    print(f"File not found: {e}")




df = pd.DataFrame(data_list)
df.to_csv('results/extracted_data.csv', index=False)


 