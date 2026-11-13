"""Expected edge cases of small modalities for run_tensorqtl.py (R33).

Each case is detected before the tensorQTL call that would fail on it, and
leads to valid output plus a one-line note. Anything else still fails.

Imported by run_tensorqtl.py from its own directory; the QTL rules register
this file with script() so that it sits beside the cached copy of the script.
Pure pandas/numpy/scipy, so it is unit-tested without tensorQTL.
"""

import numpy as np
import pandas as pd
from scipy import stats

# Result columns of tensorQTL 1.0.10 (map_cis, then calculate_qvalues; and
# map_independent), after the phenotype_id index/column. Checked against real
# outputs in pheast/tests/test_qtl_edge_cases.py.
CIS_COLUMNS = [
    "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant_id",
    "start_distance", "end_distance", "ma_samples", "ma_count", "af",
    "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta",
]
GROUP_COLUMNS = ["group_id", "group_size"]


def empty_cis_table(grouped: bool) -> pd.DataFrame:
    """A cis table with no phenotypes, written with its phenotype_id index.

    Columns are those of a table with no significant phenotype (q-values but
    no pval_nominal_threshold, which calculate_qvalues adds only if some
    phenotype is significant).
    """
    columns = CIS_COLUMNS + (GROUP_COLUMNS if grouped else []) + ["qval"]
    return pd.DataFrame(columns=columns, index=pd.Index([], name="phenotype_id"))


def empty_independent_table(grouped: bool) -> pd.DataFrame:
    """A cis_independent table with no rows, written without an index."""
    columns = ["phenotype_id"] + CIS_COLUMNS + (GROUP_COLUMNS if grouped else []) + ["rank"]
    return pd.DataFrame(columns=columns)


def qvalue_pval_column(res_df: pd.DataFrame) -> str:
    """The p-value column tensorQTL's calculate_qvalues uses."""
    return "pval_beta" if not res_df["pval_beta"].isnull().all() else "pval_perm"


def pi0_fallback_reason(res_df: pd.DataFrame, qvalue_lambda: float):
    """Why tensorQTL's calculate_qvalues would fail on res_df, or None.

    - Fewer than 2 phenotypes: it correlates beta and permutation p-values,
      and scipy's pearsonr needs at least 2 values.
    - No (non-missing) p-value >= qvalue_lambda: qvalue estimates pi0 with a
      single lambda as mean(p >= lambda) / (1 - lambda), which is then 0, and
      stops with "The estimated pi0 <= 0".
    """
    n = len(res_df)
    if n < 2:
        return f"{n} phenotype tested, too few to estimate pi0"
    pvals = res_df[qvalue_pval_column(res_df)]
    if not (pvals >= qvalue_lambda).any():
        return f"no p-value >= lambda {qvalue_lambda} among {n} phenotypes, so the pi0 estimate is 0"
    return None


def calculate_qvalues_pi0_one(res_df: pd.DataFrame, qvalue_pi0_one, fdr: float) -> None:
    """tensorQTL's calculate_qvalues with pi0 = 1 (Benjamini-Hochberg q-values).

    qvalue_pi0_one(p) returns q-values for p with pi0 fixed at 1. Beyond the
    q-values, this mirrors calculate_qvalues of tensorQTL 1.0.10: the
    nominal p-value threshold column is added only if some phenotype is
    significant. The beta/permutation p-value correlation it logs is skipped.
    """
    pval_col = qvalue_pval_column(res_df)
    res_df["qval"] = np.asarray(qvalue_pi0_one(res_df[pval_col].to_numpy()), dtype=float)
    print(f"  * QTL phenotypes @ FDR {fdr:.2f}: {(res_df['qval'] <= fdr).sum()}")
    if pval_col == "pval_beta":
        lb = res_df.loc[res_df["qval"] <= fdr, "pval_beta"].sort_values()
        ub = res_df.loc[res_df["qval"] > fdr, "pval_beta"].sort_values()
        if len(lb) > 0:
            lb = lb.iloc[-1]
            if len(ub) > 0:
                pthreshold = (lb + ub.iloc[0]) / 2
            else:
                pthreshold = lb
            print(f"  * min p-value threshold @ FDR {fdr}: {pthreshold:.6g}")
            res_df["pval_nominal_threshold"] = stats.beta.ppf(
                pthreshold, res_df["beta_shape1"], res_df["beta_shape2"]
            )


def has_significant(cis_df: pd.DataFrame, fdr: float, fdr_col: str) -> bool:
    """Whether map_independent would find any significant phenotype in cis_df."""
    return bool((cis_df[fdr_col] <= fdr).any())
