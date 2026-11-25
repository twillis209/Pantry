"""Tests for pheast/scripts/susie_tables.py (R36)."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("susie_tables", SCRIPTS / "susie_tables.py")
st = importlib.util.module_from_spec(spec)
spec.loader.exec_module(st)


def alleles_for(variant_ids, flip=None):
    flip = flip or [False] * len(variant_ids)
    return pd.DataFrame(
        {
            "chrom": "1",
            "pos": [10 * (i + 1) for i in range(len(variant_ids))],
            "ref": "A",
            "alt": "G",
            "ot_variant_id": [f"1_{10 * (i + 1)}_A_G" for i in range(len(variant_ids))],
            "flip": flip,
        },
        index=pd.Index(variant_ids, name="variant_id"),
    )


def make_fit(L=2, p=5, cs=None, seed=0):
    rng = np.random.default_rng(seed)
    lbf = rng.normal(size=(L, p))
    alpha = st.alpha_from_lbf(lbf)
    mu = rng.normal(size=(L, p))
    mu2 = mu**2 + 0.01
    cs = cs if cs is not None else {"L2": np.array([3, 4])}
    purity = pd.DataFrame({"min_abs_corr": [0.9] * len(cs)}, index=list(cs))
    return {
        "alpha": alpha,
        "mu": mu,
        "mu2": mu2,
        "lbf_variable": lbf,
        "pip": 1 - np.prod(1 - alpha, axis=0),
        "sets": {"cs": cs, "purity": purity},
    }


def test_alpha_from_lbf_rows_sum_to_one():
    lbf = np.random.default_rng(1).normal(scale=20, size=(3, 50))
    alpha = st.alpha_from_lbf(lbf)
    np.testing.assert_allclose(alpha.sum(1), 1)
    assert (alpha.argmax(1) == lbf.argmax(1)).all()


def test_cs_rows_one_set():
    vids = pd.Index([f"chr1:{10 * (i + 1)}" for i in range(5)])
    fit = make_fit()
    af = np.full(5, 0.3)
    out = st.cs_rows("P1", "G1", vids, fit, af, alleles_for(vids))
    assert list(out.columns) == st.CS_COLUMNS
    assert len(out) == 2
    assert (out.cs_id == 2).all()
    np.testing.assert_allclose(out.alpha, fit["alpha"][1, [3, 4]])
    np.testing.assert_allclose(out.cs_coverage, fit["alpha"][1, [3, 4]].sum())
    assert (out.cs_size == 2).all()
    np.testing.assert_allclose(out.beta, fit["mu"][1, [3, 4]])
    np.testing.assert_allclose(out.beta_sd, 0.1)
    np.testing.assert_allclose(out.lbf, fit["lbf_variable"][1, [3, 4]])
    assert (out.cs_purity == 0.9).all()
    assert list(out.variant_id) == ["chr1:40", "chr1:50"]
    assert list(out.ot_variant_id) == ["1_40_A_G", "1_50_A_G"]


def test_cs_rows_flip_negates_beta_and_af():
    vids = pd.Index(["chr1:10"])
    fit = make_fit(L=1, p=1, cs={"L1": np.array([0])})
    fit["mu"][:] = 0.5
    fit["mu2"][:] = 0.26
    out = st.cs_rows("P1", "P1", vids, fit, np.array([0.2]), alleles_for(vids, flip=[True]))
    assert out.beta.iloc[0] == pytest.approx(-0.5)
    assert out.af.iloc[0] == pytest.approx(0.8)
    assert out.beta_sd.iloc[0] == pytest.approx(0.1)


def test_cs_rows_ungrouped_group_id_is_phenotype():
    vids = pd.Index([f"chr1:{10 * (i + 1)}" for i in range(5)])
    out = st.cs_rows("P1", None, vids, make_fit(), np.full(5, 0.3), alleles_for(vids))
    assert (out.group_id == "P1").all()


def test_cs_rows_two_sets_in_signal_order():
    vids = pd.Index([f"chr1:{10 * (i + 1)}" for i in range(5)])
    fit = make_fit(cs={"L2": np.array([4]), "L1": np.array([0, 1])})
    out = st.cs_rows("P1", "P1", vids, fit, np.full(5, 0.3), alleles_for(vids))
    assert list(out.cs_id) == [1, 1, 2]


def test_cs_rows_missing_allele_keeps_row_with_empty_ot_id():
    vids = pd.Index([f"chr1:{10 * (i + 1)}" for i in range(5)])
    alleles = alleles_for(vids)
    alleles.loc["chr1:40", "ot_variant_id"] = ""
    out = st.cs_rows("P1", "P1", vids, make_fit(), np.full(5, 0.3), alleles)
    assert list(out.ot_variant_id) == ["", "1_50_A_G"]


def test_empty_cs_table_header():
    t = st.empty_cs_table()
    assert list(t.columns) == st.CS_COLUMNS
    assert len(t) == 0


def test_summary_row_ok_and_without_fit():
    fit = make_fit()
    fit.update({"converged": True, "niter": 7, "elbo": np.array([-10.0, -5.0]), "sigma2": 0.8})
    row = st.summary_row("P1", "G1", 5, "ok", fit)
    assert list(row) == st.SUMMARY_COLUMNS
    assert row["n_cs"] == 1 and row["elbo"] == -5.0 and row["converged"] is True
    empty = st.summary_row("P2", "G2", 0, "no_variants")
    assert list(empty) == st.SUMMARY_COLUMNS
    assert empty["n_cs"] == 0 and empty["status"] == "no_variants"
