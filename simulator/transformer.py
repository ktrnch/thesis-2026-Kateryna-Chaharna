import pandas as pd

def recode_to_minus101(df012:pd.DataFrame):
    """
    Convert SNP dosages (0,1,2) to centered codes (-1,0,1) and transpose so that:
      - Rows = SNPs
      - Cols = individuals
    """
    return (df012.T - 1).astype("int16")