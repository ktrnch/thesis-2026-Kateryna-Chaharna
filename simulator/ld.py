import pandas as pd
from typing import Dict, List

def chromosome_metadata(snp_ids, chromosome_file):
    chromosome_data = pd.read_csv(chromosome_file, sep="\t", header=None)
    chromosome_data.columns = ["SNP", "CHR", "POS"]
    chromosome_data = chromosome_data[chromosome_data["SNP"].isin(snp_ids)]
    return chromosome_data


def find_ld_pairs(positional_data: pd.DataFrame, ld_threshold: int) -> Dict[str, List[str]]:
    """
    Identifies pairs of SNPs that are in Linkage Disequilibrium (LD)
    based on chromosome and physical distance.
    Returns a dictionary where keys are SNPs and values are lists of other
    SNPs they are in LD with.
    """
    ld_map = {}
    valid_snps = positional_data.dropna(subset=['CHR', 'POS'])

    if valid_snps.empty:
        print("No SNPs with valid CHR/POS metadata found for LD analysis.")
        return {}

    for chrom, group in valid_snps.groupby('CHR'):
        # Sort by position to find nearby SNPs faster
        sorted_group = group.sort_values('POS').reset_index(drop=True)
        #print(sorted_group)

        snps_in_chrom = sorted_group['SNP'].tolist()
        positions = sorted_group['POS'].tolist()

        for i in range(len(snps_in_chrom)):
            snp_i = snps_in_chrom[i]
            pos_i = int(positions[i])

            if snp_i not in ld_map:
                ld_map[snp_i] = []

            # Check following SNPs in the same chromosome
            for j in range(i + 1, len(snps_in_chrom)):
                snp_j = snps_in_chrom[j]
                pos_j = int(positions[j])

                if abs(pos_i - pos_j) <= ld_threshold:
                    ld_map[snp_i].append(snp_j)
                    if snp_j not in ld_map:
                        ld_map[snp_j] = []
                    ld_map[snp_j].append(snp_i)
                else:
                    # Since the group is sorted by position, if the distance
                    # to snp_j is already too large, subsequent SNPs will also be too far.
                    break

    # Remove duplicates from the lists and sort for consistency
    for snp, linked_snps in ld_map.items():
        ld_map[snp] = sorted(list(set(linked_snps)))

    return ld_map
