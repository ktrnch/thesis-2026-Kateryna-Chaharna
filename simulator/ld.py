from typing import Dict, List, Optional, Set
import pandas as pd

def chromosome_metadata(snp_ids, chromosome_file):
    if chromosome_file is None:
        raise ValueError("chromosome_file is required when ld_mode=True")

    chromosome_data = pd.read_csv(chromosome_file, sep="\t", header=None)
    chromosome_data.columns = ["SNP", "CHR", "POS"]
    chromosome_data = chromosome_data[chromosome_data["SNP"].isin(snp_ids)].copy()
    chromosome_data["POS"] = pd.to_numeric(chromosome_data["POS"], errors="coerce")

    return chromosome_data


def find_ld_pairs(positional_data: pd.DataFrame, ld_threshold: int) -> Dict[str, List[str]]:
    """
    Identifies SNP pairs that are close on the same chromosome.
    Here LD is approximated by physical distance <= ld_threshold.
    Returns:
        {snp: [nearby_snp_1, nearby_snp_2, ...]}
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

            ld_map.setdefault(snp_i, [])

            # Check following SNPs in the same chromosome
            for j in range(i + 1, len(snps_in_chrom)):
                snp_j = snps_in_chrom[j]
                pos_j = int(positions[j])

                if abs(pos_i - pos_j) <= ld_threshold:
                    ld_map.setdefault(snp_j, [])

                    ld_map[snp_i].append(snp_j)
                    ld_map[snp_j].append(snp_i)
                else:
                    # Since the group is sorted by position, if the distance
                    # to snp_j is already too large, subsequent SNPs will also be too far.
                    break

    # Remove duplicates from the lists and sort for consistency
    for snp, linked_snps in ld_map.items():
        ld_map[snp] = sorted(set(linked_snps))

    return ld_map


def prune_snps_by_ld(
    snps: List[str],
    ld_map: Dict[str, List[str]],
    protected_snps: Optional[Set[str]] = None
) -> List[str]:
    """
    Greedy LD pruning.

    Keeps one SNP from each LD/proximity group.
    SNPs in protected_snps are kept preferentially.
    """
    if protected_snps is None:
        protected_snps = set()

    kept = []
    removed = set()

    # Protected SNPs first, so user-defined SNPs are kept preferentially.
    ordered_snps = (
        [s for s in snps if s in protected_snps]
        + [s for s in snps if s not in protected_snps]
    )

    for snp in ordered_snps:
        if snp in removed:
            continue

        kept.append(snp)

        for linked_snp in ld_map.get(snp, []):
            if linked_snp not in protected_snps:
                removed.add(linked_snp)

    return kept
