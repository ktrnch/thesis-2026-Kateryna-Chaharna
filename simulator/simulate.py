#!/usr/bin/env python3
from typing import Dict, Optional
import numpy as np
import pandas as pd
import random
import json
import shutil
from pathlib import Path

from simulator.sampling import (
    set_seeds, parse_csv_names, ensure_unique,
    sample_without, split_second, split_third,
)
from simulator.io_utils import load_raw_header, selective_read_raw, full_read_raw
from simulator.ld import chromosome_metadata, find_ld_pairs
from simulator.transformer import recode_to_minus101

def simulate(
    raw_path: str, #path to the file with genotype
    output_prefix: str = "sim",
    output_dir: str = "data/runs/sim",
    config_source_path: Optional[str] = None,
    config_archive_dir: str = "data/configs",

    seed: int = 7, #sets seed of random generator fixed seed -> same parameeters give same results
    phenotype_type: str = "quantitative",  # "quantitative" or "binary"
    phenotype_number: int = 1,
    n_indiv_max: Optional[int] = None,     # trim to first N individuals for speed/memory (None = all)

    # Causal pool configuration
    causal_var_pool_num: int = 1000,       # how many SNPs form the causal pool
    causal_var_min: int = 100,              # min #causal SNPs per phenotype
    causal_var_max: int = 100,              # max #causal SNPs per phenotype
    causal_snps: Optional[str] = None,     # optional fixed list of SNP names to include in causal pool

    #LD mode
    ld_mode:  bool = False,                 #if true compute
    ld_threshold: int = 1000,               # i am not sure if it should be optional
    chromosome_file: str = None,             # path to the file with chromosome information

    # Nonlinear effects configuration
    second_inter_num: int = 0,             # number of 2nd-order pairs
    third_inter_num: int = 0,              # number of 3rd-order triples
    dom_num: int = 0,                      # count of SNPs encoded as dominant (presence >= 1)
    rec_num: int = 0,                      # count of SNPs encoded as recessive (homozygous alt == 2)
    int_var_max: int = 0,                  # max # interaction initiators included per phenotype (random mode)
    second_inter_snps: Optional[str] = None,  # explicit names for 2nd-order: A1..Ak,B1..Bk
    third_inter_snps: Optional[str]  = None,  # explicit names for 3rd-order: A1..Ak,B1..Bk,C1..Ck
    dom_snps: Optional[str] = None,           # explicit dominant SNP names
    rec_snps: Optional[str] = None,           # explicit recessive SNP names
    force_defined_interactions: bool = False, # if True: include ALL defined interaction initiators in every phenotype

    # Weights
    weight_renewal: bool = False,               # if True: draw weights per-phenotype; else one shared vector
    separate_interaction_weights: bool = False, # if True: draw separate weights for interaction terms
    interaction_weight_scale: float = 1.0,      # multiplier applied to interaction weights
    interactions_mode: str = "additive",        # "additive" = linear + interaction; "pure" = interaction only for initiators

    # Noise / heritability
    target_h2: Optional[float] = None,  # if set for quantitative phenotypes, we enforce exact h² on observed scale
    noise_factor: float = 0.0,          # std of Gaussian noise (ignored if target_h2 is not None)

    # Output mapping
    quant_range_min: float = 0.0,       # min of final quantitative scale
    quant_range_max: float = 1.0,       # max of final quantitative scale
    percentile_or_threshold: float = 80.0  # binary threshold percentile
) -> Dict[str, pd.DataFrame]:
    """
    Returns a dict of DataFrames and writes files to disk:
      - phenotype_without_noise.txt (genetic component)
      - phenotype_with_noise.txt    (genetic + noise / h²-enforced)
      - final_phenotype.txt         (mapped to range / thresholded)
      - phenotype_snps.csv          (selected SNPs per phenotype)
      - {prefix}_SNP_ASSIGNMENTS.txt (role lists)
      - {prefix}_weights*.tsv        (weights)
      - {prefix}_h2_report.tsv       (realized h²) for quantitative
      - {prefix}_meta.json           (run metadata)
    """

    # ---------------- validations ----------------
    set_seeds(seed)
    if causal_var_min > causal_var_max:
        raise ValueError("causal_var_min <= causal_var_max required")
    if causal_var_max > causal_var_pool_num:
        raise ValueError("causal_var_max <= causal_var_pool_num required")
    if phenotype_type == "quantitative" and quant_range_min >= quant_range_max:
        raise ValueError("quant_range_min < quant_range_max required")
    if target_h2 is not None:
        # Target heritability is only meaningful for quantitative outputs
        if not (0.0 < target_h2 < 1.0):
            raise ValueError("target_h2 in (0,1)")
        if phenotype_type != "quantitative":
            raise ValueError("target_h2 for quantitative only")
        if noise_factor != 0.0:
            raise ValueError("set noise_factor=0 when using target_h2")

    # Parse name lists (if user provided explicit SNP names)
    specific_causal = ensure_unique(parse_csv_names(causal_snps)) if causal_snps else []
    second_list = ensure_unique(parse_csv_names(second_inter_snps)) if second_inter_snps else []
    third_list  = ensure_unique(parse_csv_names(third_inter_snps))  if third_inter_snps  else []
    dom_list    = ensure_unique(parse_csv_names(dom_snps)) if dom_snps else []
    rec_list    = ensure_unique(parse_csv_names(rec_snps)) if rec_snps else []

    # Basic sanity checks for provided interaction lists
    if second_list and len(second_list) < 2*second_inter_num:
        raise ValueError("Too few 2nd-order names")
    if third_list  and len(third_list)  < 3*third_inter_num:
        raise ValueError("Too few 3rd-order names")
    if (second_inter_num + third_inter_num) < int_var_max and not force_defined_interactions:
        raise ValueError("int_var_max > total interactions; set force_defined_interactions or lower int_var_max")

    # ---------------- load genotype table ----------------

    hdr = load_raw_header(raw_path)
    first6 = hdr[:6]       # FID, IID, PAT, MAT, SEX, PHENOTYPE
    all_snp_cols = hdr[6:] # all genotype columns

    proximity_map = {}
    # If LD=True
    if ld_mode:
        snp_metadata = chromosome_metadata(all_snp_cols, chromosome_file)
        proximity_map = find_ld_pairs(snp_metadata, ld_threshold)

    # If user specified causal SNP names, we can selectively read just t    hose (plus any interaction/dom/rec names).
    if specific_causal:
        miss = [s for s in specific_causal if s not in all_snp_cols]
        if miss:
            raise ValueError(f"--causal_snps missing (e.g. {miss[:5]})")
        need = set(specific_causal) | set(second_list) | set(third_list) | set(dom_list) | set(rec_list)
        df_raw = selective_read_raw(raw_path, sorted(need))
    else:
        df_raw = full_read_raw(raw_path)

    # Ensure first 6 columns are aligned and in place
    if list(df_raw.columns[:6]) != first6:
        miss6 = [c for c in first6 if c not in df_raw.columns]
        if miss6:
            raise ValueError(f"Missing first6: {miss6}")
        df_raw = df_raw[first6 + [c for c in df_raw.columns if c not in first6]]

    # Genotypes as numbers; fill NAs to 1 (heterozygote) as a neutral-ish default
    df_012 = df_raw.iloc[:, 6:].apply(pd.to_numeric, errors="coerce").fillna(1)

    # Optional: restrict to first N individuals for speed/memory (useful in Colab)
    if n_indiv_max is not None:
        df_raw = df_raw.iloc[:n_indiv_max].reset_index(drop=True)
        df_012 = df_012.iloc[:n_indiv_max].reset_index(drop=True)

    # Recode to {-1,0,1} and switch orientation: rows=SNPs, cols=individuals
    geno = recode_to_minus101(df_012)
    fids = df_raw["FID"].reset_index(drop=True)
    iids = df_raw["IID"].reset_index(drop=True)
    all_snp_names = geno.index.tolist()

    # ---------------- causal SNP pool ----------------
    # Build the set of SNPs from which each phenotype will draw its causal subset.
    if specific_causal:
        # Keep provided names that exist; optionally top up to causal_var_pool_num with random extras.
        causal_pool = [s for s in specific_causal if s in all_snp_names]
        if len(causal_pool) < len(specific_causal):
            miss = [s for s in specific_causal if s not in all_snp_names]
            raise ValueError(f"Causal SNPs not present: {miss[:5]}")
        need_extra = max(0, causal_var_pool_num - len(causal_pool))
        if need_extra > 0:
            extra = sample_without(all_snp_names, set(causal_pool), need_extra)
            causal_pool = ensure_unique(causal_pool + extra)
    else:
        # Fully random causal pool from available SNPs
        causal_pool = random.sample(all_snp_names, min(causal_var_pool_num, len(all_snp_names)))
    causal_pool = sorted(causal_pool)

    # Submatrix of the causal pool only
    df_causal = geno.loc[causal_pool, :]
    X = df_causal.to_numpy()       # rows=SNPs, cols=individuals
    n_snps, n_ind = X.shape

    # ---------------- define interactions & special encodings ----------------
    # We select names for 2nd/3rd-order interactions either from user lists or randomly from the pool.
    inter2_total = second_inter_num * 2
    inter3_total = third_inter_num  * 3

    if second_list:
        sec_all = [s for s in second_list if s in causal_pool]
        if len(sec_all) < inter2_total:
            raise ValueError("2nd-order names not in pool / too few")
        sec_all = sec_all[:inter2_total]
    else:
        sec_all = random.sample(causal_pool, inter2_total) if inter2_total > 0 else []

    if third_list:
        thr_all = [s for s in third_list if s in causal_pool]
        if len(thr_all) < inter3_total:
            raise ValueError("3rd-order names not in pool / too few")
        thr_all = thr_all[:inter3_total]
    else:
        # Avoid reusing 2nd-order picks if possible (purely aesthetic; not required)
        remaining_for_third = [s for s in causal_pool if s not in sec_all]
        thr_all = random.sample(remaining_for_third, inter3_total) if inter3_total > 0 else []

    # Split into initiators/counterparts (2nd) and triplets (3rd)
    init2, pair2 = ([], [])
    init3, b3, c3 = ([], [], [])
    if inter2_total > 0:
        init2, pair2 = split_second(sec_all, second_inter_num)
    if inter3_total > 0:
        init3, b3, c3 = split_third(thr_all,  third_inter_num)

    # Keep track of names already reserved for interactions (so dominant/recessive avoid collisions)
    used_for_inter = set(init2 + pair2 + init3 + b3 + c3)

    def choose_special(n, given_list, forbidden):
        """
        Pick n special SNPs (dominant/recessive).
        If names provided, use them (and fill if too few). Otherwise, sample from pool.
        """
        if n <= 0:
            return []
        if given_list:
            src = [s for s in ensure_unique(given_list) if s in causal_pool and s not in forbidden]
            if len(src) < n:
                src += sample_without(causal_pool, forbidden.union(src), n - len(src))
            return src[:n]
        return sample_without(causal_pool, forbidden, n)

    dominant_pos  = choose_special(dom_num, parse_csv_names(dom_snps) if dom_snps else [], used_for_inter)
    used_inter_dom = used_for_inter.union(dominant_pos)
    recessive_pos = choose_special(rec_num, parse_csv_names(rec_snps) if rec_snps else [], used_inter_dom)

    # All SNPs that participate in any interaction (initiator or counterpart)
    interacting_all = sorted(set(init2 + pair2 + init3 + b3 + c3))
    snp_to_row = {s:i for i,s in enumerate(df_causal.index)}  # map SNP name -> row index in X

    # ---------------- draw weights ----------------
    # If weight_renewal=False: one common weight vector across phenotypes
    # If weight_renewal=True: a separate weight vector per phenotype (stored in a dict)
    if weight_renewal:
        weights_linear = {p: np.random.normal(0,1,n_snps) for p in range(phenotype_number)}
        weights_inter  = ({p: np.random.normal(0,1,n_snps) for p in range(phenotype_number)}
                          if separate_interaction_weights else weights_linear)
    else:
        w_lin = np.random.normal(0,1,n_snps)
        w_int = np.random.normal(0,1,n_snps) if separate_interaction_weights else w_lin
        weights_linear, weights_inter = w_lin, w_int

    # ommit snps that are in LD with selected ones from the pool
    if proximity_map:
        linked_snps = set()
        for s in causal_pool:
            if s in proximity_map:
                linked_snps.update(proximity_map[s])
        causal_pool = [s for s in causal_pool if s not in linked_snps]     

    # ---------------- choose causal SNPs per phenotype ----------------
    # If force_defined_interactions=True, we include *all* defined initiators + their counterparts in each phenotype.
    # Otherwise, we randomly choose up to `int_var_max` initiators (and include their pairs), then fill with random SNPs.
    #account for ld
    def pick_pheno_snps(force_defined: bool,):
        out = []
        forced_inits = set(init2 + init3)
        forced_pairs = set(pair2 + b3 + c3)
        for _ in range(phenotype_number):
            k = random.randint(causal_var_min, causal_var_max)  # total causal SNP count for this phenotype
            if force_defined:
                sel_inits = list(forced_inits)
                pair_set = set(forced_pairs)
            else:
                num_inter = min(int_var_max, len(forced_inits))
                choose_k  = random.randint(0, num_inter)
                sel_inits = random.sample(list(forced_inits), choose_k) if choose_k > 0 else []
                # add their counterparts
                pair_set  = set()
                for s in sel_inits:
                    if s in init2:
                        idx = init2.index(s); pair_set.add(pair2[idx])
                    if s in init3:
                        idx = init3.index(s); pair_set.add(b3[idx]); pair_set.add(c3[idx])

            used = set(sel_inits) | pair_set

            # if proximity_map:
            #     to_check = list(used)
            #     for i in range(len(to_check)):
            #         snp = to_check[i]
            #         if snp in proximity_map:
            #             for linked_snp in proximity_map[snp]:
            #                 if linked_snp in causal_pool and linked_snp not in used:
            #                     used.add(linked_snp)
            #                     to_check.append(linked_snp)

            remaining = [s for s in causal_pool if s not in used]
            need = max(0, k - len(used))
            picked_rest = remaining if need >= len(remaining) else random.sample(remaining, need)
            out.append(sorted(list(used) + picked_rest))
        return out

    phenotype_snps = pick_pheno_snps(force_defined=force_defined_interactions,)

    # ---------------- build genetic component G (P x N) ----------------
    # For each phenotype:
    #   - build an "interaction map" for initiator rows (accumulate pairwise/triplewise products on initiator rows)
    #   - sum up linear terms and (optionally) interaction terms using weight vectors
    G = []
    for p_idx in range(phenotype_number):
        phset = set(phenotype_snps[p_idx])

        # inter_map: initiator-row-index -> interaction vector across individuals
        inter_map = {}

        # 2nd-order interactions: X[a] * X[b] added to initiator 'a'
        for i in range(len(init2)):
            a, b = init2[i], pair2[i]
            if a in phset and b in phset:
                ra, rb = snp_to_row[a], snp_to_row[b]
                inter_map[ra] = inter_map.get(ra, np.zeros(n_ind)) + (X[ra,:] * X[rb,:])

        # 3rd-order interactions: X[a] * X[b] * X[c] added to initiator 'a'
        for i in range(len(init3)):
            a, bb, cc = init3[i], b3[i], c3[i]
            if a in phset and bb in phset and cc in phset:
                ra, rb, rc = snp_to_row[a], snp_to_row[bb], snp_to_row[cc]
                inter_map[ra] = inter_map.get(ra, np.zeros(n_ind)) + (X[ra,:] * X[rb,:] * X[rc,:])

        # Sum contributions for each SNP included in this phenotype
        sums = np.zeros(n_ind)
        for s in phenotype_snps[p_idx]:
            r = snp_to_row[s]
            x = X[r,:]

            # pick appropriate weights (shared or per-phenotype)
            if weight_renewal:
                w_lin = weights_linear[p_idx][r]
                w_int = interaction_weight_scale * weights_inter[p_idx][r]
            else:
                w_lin = weights_linear[r]
                w_int = interaction_weight_scale * weights_inter[r]

            # Dominant/Recessive encodings override the usual genotype*x form
            add_linear = True
            if s in dominant_pos:
                sums += (x >= 1).astype(float) * w_lin
                add_linear = False
            elif s in recessive_pos:
                sums += (x == 2).astype(float) * w_lin
                add_linear = False

            # Interaction term (if this SNP is an initiator and we precomputed a vector for it)
            if s in interacting_all and (r in inter_map):
                sums += inter_map[r] * w_int
                # If 'pure', we remove the linear term for initiators (keep only interaction effect)
                if interactions_mode == "pure":
                    add_linear = False

            # Ordinary additive linear effect (unless suppressed)
            if add_linear:
                sums += x * w_lin

        G.append(sums)

    G = np.vstack(G)  # Shape: (phenotypes, individuals)

    # ---------------- add noise / enforce heritability ----------------
    # If target_h2 is given (quantitative), we construct y = g + e with Var(e) chosen so that:
    #      h² = Var(g) / Var(y)  =>  Var(e) = Var(g) * (1 - h²) / h²
    if phenotype_type == "quantitative" and target_h2 is not None:
        Y = np.zeros_like(G)
        realized = []
        for i in range(G.shape[0]):
            g = G[i,:].astype(float)
            vg = np.var(g, ddof=1)
            if vg == 0:
                raise ValueError(f"Phenotype {i+1}: zero genetic variance")
            # unit Gaussian noise, then scale to achieve target Var(e)
            eps = np.random.randn(g.shape[0])
            ve_target = vg * (1.0 - target_h2) / target_h2
            scale = np.sqrt(ve_target / np.var(eps, ddof=1))
            y = g + scale * eps
            Y[i,:] = y
            realized.append(vg / np.var(y, ddof=1))
        G_with_noise = Y
        print(f"[info] target_h2={target_h2} | realized={ [round(x,4) for x in realized] }")
    else:
        # Simpler path: add iid Gaussian noise with std=noise_factor (or none)
        G_with_noise = G + (np.random.randn(*G.shape) * noise_factor if noise_factor > 0 else 0.0)

    # ---------------- map to requested output type ----------------
    if phenotype_type == "quantitative":
        # Scale all phenotypes together, using global min/max over the whole G_with_noise matrix
        mn, mx = float(np.min(G_with_noise)), float(np.max(G_with_noise))
        if mx > mn:
            F = (G_with_noise - mn) / (mx - mn) * (quant_range_max - quant_range_min) + quant_range_min
        else:
            F = np.full_like(G_with_noise, quant_range_min)  # edge case: constant vector
    else:
        # Binary mapping via percentile threshold over the entire matrix
        thr = np.percentile(G_with_noise, percentile_or_threshold)
        F = (G_with_noise > thr).astype(int)

    # ---------------- assemble outputs as DataFrames (N x P) ----------------
    pheno_labels = [f"Phenotype_{i+1}" for i in range(phenotype_number)]
    df_wo = pd.DataFrame(G.T, columns=pheno_labels)
    df_w  = pd.DataFrame(G_with_noise.T, columns=pheno_labels)
    df_f  = pd.DataFrame(F.T, columns=pheno_labels)
    for d in (df_wo, df_w, df_f):
        d.insert(0, "FID", iids)
        d.insert(1, "IID", iids)

    # ---------------- write files ----------------
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg_dir = Path(config_archive_dir)
    cfg_dir.mkdir(parents=True, exist_ok=True)

    # Core phenotype outputs
    df_wo.to_csv(out_dir / "phenotype_without_noise.txt", sep="\t", index=False)
    df_w.to_csv(out_dir / "phenotype_with_noise.txt", sep="\t", index=False)
    df_f.to_csv(out_dir / "final_phenotype.txt", sep="\t", index=False)

    # SNP selection per phenotype (names)
    pd.DataFrame({
        "Phenotype": pheno_labels,
        "Selected_SNPs": [",".join(lst) for lst in phenotype_snps]
    }).to_csv(out_dir / "phenotype_snps.csv", index=False)

    # Role assignments (causal pool, interaction sets, dominant/recessive)
    with open(out_dir / f"{output_prefix}_SNP_ASSIGNMENTS.txt", "w") as f:
        def w(label, lst):
            f.write(f"{label}\t" + ("\t".join(lst) if lst else "") + "\n")
        w("Causal_Positions", causal_pool)
        w("Init_Second_Order", init2);  w("Second_Second_Order", pair2)
        w("Init_Third_Order", init3);   w("Second_Third_Order", b3); w("Third_Third_Order", c3)
        w("Dominant_Positions", dominant_pos); w("Recessive_Positions", recessive_pos)

    # Weights
    if weight_renewal:
        wlin = pd.DataFrame({"SNP": df_causal.index})
        wint = pd.DataFrame({"SNP": df_causal.index})
        for p in range(phenotype_number):
            wlin[f"lin_w_p{p+1}"] = weights_linear[p]
            wint[f"int_w_p{p+1}"] = interaction_weight_scale * weights_inter[p]
        wlin.to_csv(out_dir / f"{output_prefix}_weights_linear.tsv", sep="\t", index=False, float_format="%.6f")
        wint.to_csv(out_dir / f"{output_prefix}_weights_interaction.tsv", sep="\t", index=False, float_format="%.6f")
    else:
        pd.DataFrame({
            "SNP": df_causal.index,
            "linear_weight": weights_linear,
            "interaction_weight": interaction_weight_scale * weights_inter
        }).to_csv(out_dir / f"{output_prefix}_weights.tsv", sep="\t", index=False, float_format="%.6f")

    # Heritability summary
    h2_df = None
    if phenotype_type == "quantitative":
        rows = []
        for i in range(G.shape[0]):
            vg = float(np.var(G[i, :], ddof=1))
            vy = float(np.var(G_with_noise[i, :], ddof=1))
            rows.append({
                "phenotype_id": i + 1,
                "var_genetic": vg,
                "var_total": vy,
                "h2_realized": (vg / vy) if vy > 0 else np.nan
            })
        h2_df = pd.DataFrame(rows)
        h2_df.to_csv(out_dir / f"{output_prefix}_h2_report.tsv", sep="\t", index=False)

    # Metadata
    meta = dict(
        seed=seed, phenotype_type=phenotype_type, phenotype_number=phenotype_number,
        causal_var_pool_num=causal_var_pool_num, causal_var_min=causal_var_min, causal_var_max=causal_var_max,
        second_inter_num=second_inter_num, third_inter_num=third_inter_num, dom_num=dom_num, rec_num=rec_num,
        int_var_max=int_var_max, interactions_mode=interactions_mode,
        separate_interaction_weights=separate_interaction_weights,
        interaction_weight_scale=interaction_weight_scale, target_h2=target_h2,
        noise_factor=noise_factor, quant_range=[quant_range_min, quant_range_max],
        percentile_or_threshold=percentile_or_threshold,
        force_defined_interactions=force_defined_interactions,
        n_individuals=int(len(iids)), n_causal_pool=int(len(causal_pool))
    )
    with open(out_dir / f"{output_prefix}_meta.json", "w") as jf:
        json.dump(meta, jf, indent=2)

    # Save a copy of the config separately
    if config_source_path:
        src = Path(config_source_path)
        if src.exists():
            shutil.copy2(src, cfg_dir / f"{output_prefix}_config.json")

    print(f"[done] outputs written to {out_dir}")
    if phenotype_type == "quantitative":
        print(f"[done] {out_dir / f'{output_prefix}_h2_report.tsv'}")
    print(f"[done] config snapshot: {cfg_dir / f'{output_prefix}_config.json'}")

    # Return key DataFrames 
    return {
        "phenotype_without_noise": df_wo,
        "phenotype_with_noise":    df_w,
        "final_phenotype":         df_f,
        "h2_report":               h2_df
    }

