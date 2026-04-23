
import random
import numpy as np
from typing import List, Optional

def set_seeds(seed:int):
    """Fix Python + NumPy RNGs(random generators) for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)

def parse_csv_names(s:Optional[str]):
    """Turn a comma-separated string into a list of stripped tokens; handles None."""
    return [t.strip() for t in s.split(",") if t.strip()] if s else []

def ensure_unique(xs):
    """Preserve order while removing duplicates."""
    return list(dict.fromkeys(xs))

def sample_without(seq, forbidden:set, k:int):
    """Sample k distinct elements from seq excluding any in `forbidden` set."""
    cand = [x for x in seq if x not in forbidden]
    return cand if k >= len(cand) else random.sample(cand, k)

def split_second(lst:List[str], n_pairs:int):
    """
    Split a list representing 2nd-order interactions into (initiators, counterparts),
    e.g., [A1..Ak, B1..Bk] -> ([A1..Ak], [B1..Bk])
    """
    if len(lst) < 2*n_pairs:
        raise ValueError("Not enough names for 2nd-order pairs")
    return lst[:n_pairs], lst[n_pairs:2*n_pairs]

def split_third(lst:List[str], n_triples:int):
    """
    Split a list representing 3rd-order interactions into (A, B, C) thirds:
    e.g., [A1..Ak, B1..Bk, C1..Ck] -> (A[...], B[...], C[...])
    """
    if len(lst) < 3*n_triples:
        raise ValueError("Not enough names for 3rd-order triples")
    return lst[:n_triples], lst[n_triples:2*n_triples], lst[2*n_triples:3*n_triples]
