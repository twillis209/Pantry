"""Tests for pheast/scripts/qtl_edge_cases.py (R33)."""

import gzip
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SCRIPTS = Path(__file__).parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("qtl_edge_cases", SCRIPTS / "qtl_edge_cases.py")
edge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(edge)

# Header lines of real tensorQTL 1.0.10 outputs (module-v0.1.6 github() host run,
# host-github-v016, geuvadis): expression (ungrouped) and isoforms (grouped).
_BASE = (
    "phenotype_id,num_var,beta_shape1,beta_shape2,true_df,pval_true_df,variant_id,"
    "start_distance,end_distance,ma_samples,ma_count,af,pval_nominal,slope,slope_se,"
    "pval_perm,pval_beta"
)
REAL_CIS = {
    False: _BASE + ",qval,pval_nominal_threshold",
    True: _BASE + ",group_id,group_size,qval,pval_nominal_threshold",
}
REAL_INDEPENDENT = {
    False: _BASE + ",rank",
    True: _BASE + ",group_id,group_size,rank",
}


def header_line(df, path, index):
    df.to_csv(path, sep="\t", index=index)
    with gzip.open(path, "rt") as handle:
        return handle.read()


def bh(p):
    """Benjamini-Hochberg adjusted p-values (= qvalue with pi0 = 1), NaN kept."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    order = np.argsort(pv)[::-1]
    ranked = np.minimum.accumulate(pv[order] * m / np.arange(m, 0, -1))
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def cis_table(pval_beta, pval_perm=None):
    n = len(pval_beta)
    return pd.DataFrame(
        {
            "beta_shape1": [1.03] * n,
            "beta_shape2": [380.0] * n,
            "pval_perm": pval_perm if pval_perm is not None else [min(1, p + 1e-4) for p in pval_beta],
            "pval_beta": pval_beta,
        },
        index=pd.Index([f"P{i}" for i in range(n)], name="phenotype_id"),
    )


class TestHeaders(unittest.TestCase):
    def test_empty_cis_header_matches_real_minus_threshold(self):
        # With no significant phenotype, calculate_qvalues adds no
        # pval_nominal_threshold column.
        with tempfile.TemporaryDirectory() as tmp:
            for grouped in (False, True):
                line = header_line(edge.empty_cis_table(grouped), Path(tmp) / "c.txt.gz", True)
                expected = REAL_CIS[grouped].removesuffix(",pval_nominal_threshold")
                self.assertEqual(line, expected.replace(",", "\t") + "\n")

    def test_empty_independent_header_matches_real(self):
        with tempfile.TemporaryDirectory() as tmp:
            for grouped in (False, True):
                line = header_line(edge.empty_independent_table(grouped), Path(tmp) / "i.txt.gz", False)
                self.assertEqual(line, REAL_INDEPENDENT[grouped].replace(",", "\t") + "\n")

    def test_empty_cis_table_reads_back_with_no_significant(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.txt.gz"
            edge.empty_cis_table(False).to_csv(path, sep="\t")
            cis_df = pd.read_csv(path, sep="\t", index_col=0)
            self.assertEqual(len(cis_df), 0)
            self.assertFalse(edge.has_significant(cis_df, 0.05, "qval"))


class TestPi0FallbackReason(unittest.TestCase):
    def test_one_phenotype(self):
        self.assertIn("1 phenotype tested", edge.pi0_fallback_reason(cis_table([0.9]), 0.85))

    def test_zero_phenotypes(self):
        self.assertIsNotNone(edge.pi0_fallback_reason(cis_table([]), 0.85))

    def test_no_pvalue_at_or_above_lambda(self):
        reason = edge.pi0_fallback_reason(cis_table([0.2, 0.5, 0.6, 0.4]), 0.85)
        self.assertIn("pi0 estimate is 0", reason)

    def test_strong_signal_only(self):
        self.assertIsNotNone(edge.pi0_fallback_reason(cis_table([1e-6, 1e-5, 0.01]), 0.85))

    def test_normal_case_uses_tensorqtl(self):
        self.assertIsNone(edge.pi0_fallback_reason(cis_table([0.2, 0.5, 0.9, 0.4]), 0.85))
        self.assertIsNone(edge.pi0_fallback_reason(cis_table([0.2, 0.85]), 0.85))

    def test_nan_pvalues_do_not_count(self):
        self.assertIsNotNone(edge.pi0_fallback_reason(cis_table([0.2, np.nan, 0.3]), 0.85))

    def test_permutation_pvalues_when_no_beta(self):
        df = cis_table([np.nan, np.nan], pval_perm=[0.1, 0.9])
        self.assertEqual(edge.qvalue_pval_column(df), "pval_perm")
        self.assertIsNone(edge.pi0_fallback_reason(df, 0.85))
        df = cis_table([np.nan, np.nan], pval_perm=[0.1, 0.2])
        self.assertIsNotNone(edge.pi0_fallback_reason(df, 0.85))


class TestQvaluesPi0One(unittest.TestCase):
    def test_bh_qvalues_and_threshold(self):
        p = [1e-6, 1e-5, 0.01, 0.3]
        df = cis_table(p)
        edge.calculate_qvalues_pi0_one(df, bh, 0.05)
        np.testing.assert_allclose(df["qval"], [4e-6, 2e-5, 0.01 * 4 / 3, 0.3])
        # As tensorQTL: midpoint of the largest significant and smallest
        # non-significant pval_beta, mapped through each phenotype's beta.
        pthreshold = (0.01 + 0.3) / 2
        np.testing.assert_allclose(
            df["pval_nominal_threshold"], stats.beta.ppf(pthreshold, 1.03, 380.0)
        )
        self.assertEqual(list(df.columns)[-2:], ["qval", "pval_nominal_threshold"])

    def test_all_significant_threshold_is_largest(self):
        df = cis_table([1e-6, 1e-5])
        edge.calculate_qvalues_pi0_one(df, bh, 0.05)
        np.testing.assert_allclose(df["pval_nominal_threshold"], stats.beta.ppf(1e-5, 1.03, 380.0))

    def test_single_phenotype(self):
        df = cis_table([0.003])
        edge.calculate_qvalues_pi0_one(df, bh, 0.05)
        self.assertEqual(df["qval"].tolist(), [0.003])
        self.assertIn("pval_nominal_threshold", df.columns)

    def test_none_significant_adds_no_threshold(self):
        df = cis_table([0.2, 0.5])
        edge.calculate_qvalues_pi0_one(df, bh, 0.05)
        self.assertNotIn("pval_nominal_threshold", df.columns)
        self.assertFalse(edge.has_significant(df, 0.05, "qval"))

    def test_permutation_pvalues_add_no_threshold(self):
        df = cis_table([np.nan, np.nan], pval_perm=[1e-4, 0.2])
        edge.calculate_qvalues_pi0_one(df, bh, 0.05)
        np.testing.assert_allclose(df["qval"], [2e-4, 0.2])
        self.assertNotIn("pval_nominal_threshold", df.columns)


class TestHasSignificant(unittest.TestCase):
    def test_boundary_and_nan(self):
        df = pd.DataFrame({"qval": [0.05, np.nan]})
        self.assertTrue(edge.has_significant(df, 0.05, "qval"))
        df = pd.DataFrame({"qval": [0.0501, np.nan]})
        self.assertFalse(edge.has_significant(df, 0.05, "qval"))
        self.assertFalse(edge.has_significant(pd.DataFrame({"qval": [math.nan]}), 0.05, "qval"))


if __name__ == "__main__":
    unittest.main()
