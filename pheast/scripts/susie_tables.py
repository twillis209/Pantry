"""Tables for run_tensorqtl.py --mode susie (R36), without torch.

One row per variant per credible set (signal), in the form coloc_pip takes:
alpha is the signal's posterior over variants, beta its posterior effect of
alt, and ot_variant_id the GRCh38 chrom_pos_ref_alt ID Open Targets uses.
"""

import numpy as np
import pandas as pd

CS_COLUMNS = [
    "phenotype_id", "group_id", "cs_id", "variant_id", "chrom", "pos", "ref", "alt",
    "ot_variant_id", "alpha", "pip", "lbf", "beta", "beta_sd", "af",
    "cs_coverage", "cs_purity", "cs_size",
]
SUMMARY_COLUMNS = [
    "phenotype_id", "group_id", "n_variants", "converged", "niter", "elbo", "n_cs",
    "sigma2", "status",
]


def alpha_from_lbf(lbf):
    """Each signal's posterior over variants under a uniform prior: row softmax."""
    lbf = np.asarray(lbf, dtype=float)
    w = np.exp(lbf - lbf.max(axis=1, keepdims=True))
    return w / w.sum(axis=1, keepdims=True)


def empty_cs_table():
    return pd.DataFrame(columns=CS_COLUMNS)


def cs_rows(phenotype_id, group_id, variant_ids, fit, af, alleles):
    """Credible-set rows of one phenotype's SuSiE fit.

    fit: numpy alpha, mu, mu2, lbf_variable (L x p), pip (p), and tensorQTL's
    sets ({'cs': {'L<k>': indices}, 'purity': frame indexed by 'L<k>'}).
    alleles: indexed by variant_id, with chrom, pos, ref, alt, ot_variant_id
    and flip (the counted allele is GRCh38's ref: negate beta, af -> 1 - af).
    """
    group_id = phenotype_id if group_id is None else group_id
    cs = fit["sets"]["cs"] or {}
    out = []
    for name in sorted(cs, key=lambda c: int(c[1:])):
        l = int(name[1:]) - 1
        ix = np.asarray(cs[name], dtype=int)
        vids = variant_ids[ix]
        a = alleles.loc[vids]
        flip = np.where(a["flip"].to_numpy(dtype=bool), -1.0, 1.0)
        mu = fit["mu"][l, ix]
        alpha = fit["alpha"][l, ix]
        afs = np.asarray(af, dtype=float)[ix]
        out.append(pd.DataFrame({
            "phenotype_id": phenotype_id,
            "group_id": group_id,
            "cs_id": l + 1,
            "variant_id": list(vids),
            "chrom": a["chrom"].to_numpy(),
            "pos": a["pos"].to_numpy(),
            "ref": a["ref"].to_numpy(),
            "alt": a["alt"].to_numpy(),
            "ot_variant_id": a["ot_variant_id"].to_numpy(),
            "alpha": alpha,
            "pip": np.asarray(fit["pip"])[ix],
            "lbf": fit["lbf_variable"][l, ix],
            "beta": flip * mu,
            "beta_sd": np.sqrt(np.maximum(fit["mu2"][l, ix] - mu**2, 0)),
            "af": np.where(flip < 0, 1 - afs, afs),
            "cs_coverage": alpha.sum(),
            "cs_purity": float(fit["sets"]["purity"].loc[name, "min_abs_corr"]),
            "cs_size": len(ix),
        }))
    if not out:
        return empty_cs_table()
    return pd.concat(out, ignore_index=True)[CS_COLUMNS]


def summary_row(phenotype_id, group_id, n_variants, status, fit=None):
    """One phenotype's line of the summary table; fit is absent when SuSiE did not run."""
    group_id = phenotype_id if group_id is None else group_id
    row = dict.fromkeys(SUMMARY_COLUMNS)
    row.update(phenotype_id=phenotype_id, group_id=group_id, n_variants=int(n_variants),
               status=status, n_cs=0)
    if fit is not None:
        elbo = np.asarray(fit.get("elbo", []), dtype=float)
        cs = fit.get("sets", {}).get("cs") if fit.get("converged") else None
        row.update(
            converged=bool(fit.get("converged")),
            niter=int(fit["niter"]) if fit.get("niter") is not None else None,
            elbo=float(elbo[-1]) if elbo.size else None,
            n_cs=len(cs) if cs else 0,
            sigma2=float(fit["sigma2"]) if fit.get("sigma2") is not None else None,
        )
    return row
