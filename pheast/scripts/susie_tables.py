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
    sets ({'cs': {'L<k>': indices}, 'purity': frame indexed by 'L<k>'}), and
    optionally x_scale (p): the SD of each residualised genotype, since
    tensorQTL's SuSiE fits standardised genotypes and mu is per SD; beta and
    beta_sd are per allele.
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
        scale = np.asarray(fit.get("x_scale", np.ones(fit["mu"].shape[1])), dtype=float)[ix]
        mu = fit["mu"][l, ix]
        sd = np.sqrt(np.maximum(fit["mu2"][l, ix] - mu**2, 0))
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
            "beta": flip * mu / scale,
            "beta_sd": sd / scale,
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


def _matches(ref_base, chrom, pos, allele):
    return ref_base(chrom, int(pos), len(allele)).upper() == allele.upper()


def orient_alleles(bim, ref_base, counted):
    """GRCh38 ref/alt and Open Targets IDs for a bim's variants.

    bim: chrom, snp, pos, a1 (column 5), a2 (column 6). counted names the
    column whose allele tensorQTL's dosage counts; the other is taken as REF.
    ref_base(chrom, pos, length) gives the reference sequence (chrom without
    'chr'). The bim's orientation is kept when its REF matches the reference,
    so indels where both alleles match stay as the bim has them; otherwise the
    alleles are swapped if the counted one matches (flip: tensorQTL's dosage
    counts ref), else the variant is a mismatch with no ot_variant_id.
    """
    other = "a2" if counted == "a1" else "a1"
    rows = []
    for chrom, snp, pos, c, o in zip(
        bim["chrom"].astype(str).str.removeprefix("chr"), bim["snp"], bim["pos"],
        bim[counted], bim[other],
    ):
        if _matches(ref_base, chrom, pos, o):
            ref, alt, flip, status = o, c, False, "ok"
        elif _matches(ref_base, chrom, pos, c):
            ref, alt, flip, status = c, o, True, "flipped"
        else:
            ref, alt, flip, status = o, c, False, "mismatch"
        ot = "" if status == "mismatch" else f"{chrom}_{int(pos)}_{ref}_{alt}"
        rows.append((snp, chrom, int(pos), ref, alt, ot, flip, status))
    return pd.DataFrame(
        rows, columns=["variant_id", "chrom", "pos", "ref", "alt", "ot_variant_id", "flip", "status"]
    ).set_index("variant_id")


def read_fasta(path, keep):
    """Yield (name, sequence) for the FASTA records named in keep, one at a time.

    The name is the header's first word without 'chr'. No index is read or
    written, so the FASTA may sit in a read-only directory.
    """
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name in keep:
                    yield name, "".join(parts)
                name = line[1:].split()[0].removeprefix("chr")
                parts = []
            elif name in keep:
                parts.append(line.strip())
    if name in keep:
        yield name, "".join(parts)
