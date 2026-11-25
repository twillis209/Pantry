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


# Task 2: orient_alleles against a 30-base reference for chromosome 1.
REF_SEQ = "ACGTACGTAATGCATGCAGGTTCCAAGGTT"  # position 1 = 'A'


def ref_base(chrom, pos, length):
    assert chrom == "1"
    return REF_SEQ[pos - 1 : pos - 1 + length]


def bim(rows):
    return pd.DataFrame(rows, columns=["chrom", "snp", "pos", "a1", "a2"])


def test_orient_snp_ok():
    # position 10 = 'A'; counted a1 = G (ALT), a2 = A (REF)
    out = st.orient_alleles(bim([["1", "v1", 10, "G", "A"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert (r.ref, r.alt, bool(r.flip), r.status) == ("A", "G", False, "ok")
    assert r.ot_variant_id == "1_10_A_G"


def test_orient_snp_swapped():
    # counted allele is the reference base: flip
    out = st.orient_alleles(bim([["1", "v1", 10, "A", "G"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert (r.ref, r.alt, bool(r.flip), r.status) == ("A", "G", True, "flipped")
    assert r.ot_variant_id == "1_10_A_G"


def test_orient_indel():
    # positions 10-11 = 'AT': deletion REF 'AT', counted ALT 'A'
    out = st.orient_alleles(bim([["1", "v1", 10, "A", "AT"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert (r.ref, r.alt, bool(r.flip), r.status) == ("AT", "A", False, "ok")
    # swapped: bim says REF 'A' (matches too), keep the bim's orientation
    out = st.orient_alleles(bim([["1", "v1", 10, "AT", "A"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert (r.ref, r.alt, bool(r.flip)) == ("A", "AT", False)
    # bim REF 'AG' does not match 'AT', counted 'A' does: flip
    out = st.orient_alleles(bim([["1", "v1", 10, "A", "AG"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert (r.ref, r.alt, bool(r.flip), r.status) == ("A", "AG", True, "flipped")


def test_orient_mismatch():
    out = st.orient_alleles(bim([["1", "v1", 10, "C", "T"]]), ref_base, "a1")
    r = out.loc["v1"]
    assert r.status == "mismatch" and r.ot_variant_id == "" and not bool(r.flip)


def test_chr_prefix_stripped():
    out = st.orient_alleles(bim([["chr1", "chr1:10", 10, "G", "A"]]), ref_base, "a1")
    assert out.loc["chr1:10"].ot_variant_id.startswith("1_")
    assert out.loc["chr1:10"].chrom == "1"


def test_read_fasta_by_chrom(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">1 dna:chromosome\nACGT\nACGT\n>2\nTTTT\n>X\nGGGG\n")
    got = dict(st.read_fasta(fa, keep={"1", "2"}))
    assert got == {"1": "ACGTACGT", "2": "TTTT"}


def test_variant_alleles_cli(tmp_path):
    import subprocess
    import sys as _sys

    (tmp_path / "ref.fa").write_text(">chr1\n" + REF_SEQ + "\n>2\nAAAA\n")
    (tmp_path / "g.bim").write_text("1\tchr1:10\t0\t10\tG\tA\n1\tchr1:12\t0\t12\tC\tT\n")
    out = tmp_path / "alleles.tsv.gz"
    r = subprocess.run(
        [_sys.executable, str(SCRIPTS / "variant_alleles.py"), str(tmp_path / "g.bim"),
         str(tmp_path / "ref.fa"), str(out), "a1"],
        capture_output=True, text=True, check=True,
    )
    assert "ok 1, flipped 0, mismatch 1" in r.stdout
    t = pd.read_csv(out, sep="\t", index_col=0, keep_default_na=False)
    assert list(t.index) == ["chr1:10", "chr1:12"]
    assert t.loc["chr1:10", "ot_variant_id"] == "1_10_A_G"
    assert t.loc["chr1:12", "status"] == "mismatch"


def test_cs_rows_rescales_effects_to_genotype_units():
    # tensorQTL's SuSiE works on standardised genotypes: mu is per genotype SD.
    vids = pd.Index(["chr1:10", "chr1:20"])
    fit = make_fit(L=1, p=2, cs={"L1": np.array([0, 1])})
    fit["mu"][:] = [[0.2, 0.4]]
    fit["mu2"][:] = [[0.05, 0.17]]
    fit["x_scale"] = np.array([0.5, 2.0])
    out = st.cs_rows("P1", "P1", vids, fit, np.full(2, 0.3), alleles_for(vids))
    np.testing.assert_allclose(out.beta, [0.4, 0.2])
    np.testing.assert_allclose(out.beta_sd, [0.2, 0.05])
