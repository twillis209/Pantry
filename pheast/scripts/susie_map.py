"""SuSiE fine-mapping for run_tensorqtl.py --mode susie (R36).

The loop of tensorqtl.susie.map (tensorQTL 1.0.10, susie.py lines 627-738)
with its defaults, keeping what susie.map drops: each signal's alpha and
posterior moments mu, mu2 (for the effect direction), and a summary row for
every phenotype in the BED, including those not fine-mapped.
"""

import numpy as np
import pandas as pd
import torch
from tensorqtl import genotypeio, susie
from tensorqtl.core import Residualizer, calculate_maf, impute_mean

import susie_tables  # beside this script

# susie.map's defaults (max_iter is 200 there, 100 in susie()).
SUSIE_ARGS = dict(
    L=10, scaled_prior_variance=0.2, estimate_residual_variance=True,
    estimate_prior_variance=True, tol=1e-3, coverage=0.95, min_abs_corr=0.5,
    max_iter=200,
)


def _numpy(x):
    return x.cpu().numpy() if torch.is_tensor(x) else np.asarray(x)


def fine_map(genotype_df, variant_df, pheno, pheno_pos, covar, groups, alleles,
             maf_threshold, window):
    """(credible-set table, summary table) for every phenotype of pheno."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    group_of = (lambda p: p) if groups is None else (lambda p: groups.get(p, p))
    constant = pheno.index[np.all(pheno.values == pheno.values[:, [0]], 1)]
    summary = {p: susie_tables.summary_row(p, group_of(p), 0, "constant") for p in constant}
    cs_parts = []

    if len(constant) < len(pheno):
        residualizer = Residualizer(torch.tensor(covar.values, dtype=torch.float32).to(device))
        genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in pheno.columns])
        genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
        igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, pheno, pheno_pos, window=window)
        for phenotype, genotypes, genotype_range, phenotype_id in igc.generate_data(verbose=True):
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:, genotype_ix_t]
            impute_mean(genotypes_t)
            variant_ids = variant_df.index[genotype_range[0]:genotype_range[-1] + 1].rename("variant_id")
            mask_t = ~(genotypes_t == genotypes_t[:, [0]]).all(1)
            if maf_threshold > 0:
                mask_t &= calculate_maf(genotypes_t) >= maf_threshold
            if not mask_t.any():
                summary[phenotype_id] = susie_tables.summary_row(
                    phenotype_id, group_of(phenotype_id), 0, "no_variants")
                continue
            genotypes_t = genotypes_t[mask_t]
            variant_ids = variant_ids[mask_t.cpu().numpy().astype(bool)]

            phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
            genotypes_res_t = residualizer.transform(genotypes_t)
            phenotype_res_t = residualizer.transform(phenotype_t.reshape(1, -1))
            res = susie.susie(genotypes_res_t.T, phenotype_res_t.T, **SUSIE_ARGS)
            af = _numpy(genotypes_t.sum(1) / (2 * genotypes_t.shape[1]))
            fit = {
                "alpha": _numpy(res["alpha"]), "mu": _numpy(res["mu"]), "mu2": _numpy(res["mu2"]),
                "lbf_variable": _numpy(res["lbf_variable"]), "pip": _numpy(res["pip"]),
                "sets": res["sets"], "converged": bool(res["converged"]), "niter": res["niter"],
                "elbo": _numpy(res["elbo"]), "sigma2": float(_numpy(res["sigma2"])),
            }
            status = "ok" if fit["converged"] else "not_converged"
            summary[phenotype_id] = susie_tables.summary_row(
                phenotype_id, group_of(phenotype_id), len(variant_ids), status, fit)
            if fit["converged"] and fit["sets"]["cs"] is not None:
                cs_parts.append(susie_tables.cs_rows(
                    phenotype_id, group_of(phenotype_id), variant_ids, fit, af, alleles))

    # Phenotypes InputGeneratorCis dropped (no variants in the window or on
    # the chromosome).
    for p in pheno.index:
        summary.setdefault(p, susie_tables.summary_row(p, group_of(p), 0, "no_variants"))
    summary_df = pd.DataFrame([summary[p] for p in pheno.index], columns=susie_tables.SUMMARY_COLUMNS)
    cs_df = pd.concat(cs_parts, ignore_index=True) if cs_parts else susie_tables.empty_cs_table()
    return cs_df, summary_df
