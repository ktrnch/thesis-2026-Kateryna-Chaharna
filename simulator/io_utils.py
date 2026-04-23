
import pandas as pd
from typing import List

def load_raw_header(raw_path:str):
    """Read only the header row of a PLINK .raw (whitespace-separated)."""
    hdr = pd.read_csv(raw_path, sep=r"\s+", engine="python", nrows=0)
    return hdr.columns.tolist()

def selective_read_raw(raw_path:str, snp_names:List[str]):
    """
    Memory-friendly loader: keep first 6 meta columns + only requested SNPs.
    - The first 6 columns are the PLINK meta columns: FID, IID, PAT, MAT, SEX, PHENOTYPE
    - SNPs are expected to be in {0,1,2}
    """
    cols = load_raw_header(raw_path)
    first6 = cols[:6]
    missing = [s for s in snp_names if s not in cols]
    if missing:
        raise ValueError(f"SNPs not in RAW header (e.g. {missing[:5]}); total missing={len(missing)}")
    usecols = first6 + snp_names
    # Use small integers to reduce RAM footprint
    dtype_map = {s: "Int8" for s in snp_names}
    return pd.read_csv(raw_path, sep=r"\s+", engine="python", usecols=usecols, dtype=dtype_map)

def full_read_raw(raw_path:str):
    """Load the full .raw (all SNP columns). More RAM, but simple."""
    return pd.read_csv(raw_path, sep=r"\s+", engine="python")