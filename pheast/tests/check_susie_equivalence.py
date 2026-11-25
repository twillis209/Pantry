"""Check susie_map.fine_map against tensorqtl.susie.map (R36).

Run inside the tensorQTL container (it needs torch and tensorqtl), e.g.
  apptainer exec --nv <image> python3 pheast/tests/check_susie_equivalence.py \
      GENO_PREFIX BED COVAR N [MAF]
On the first N phenotypes of BED, both must give the same credible sets
(phenotype, signal, variant), the same PIPs, and alpha must equal the
softmax of susie.map's lbf_variable. Prints EQUIVALENT on success.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, susie

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
import susie_map  # noqa: E402
import susie_tables  # noqa: E402


def main(geno_prefix, bed, covar_path, n, maf="0"):
    n, maf = int(n), float(maf)
    pheno, pheno_pos = tensorqtl.read_phenotype_bed(bed)
    pheno, pheno_pos = pheno.iloc[:n], pheno_pos.iloc[:n]
    covar = pd.read_csv(covar_path, sep="\t", index_col=0).T
    pr = genotypeio.PlinkReader(geno_prefix)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    samples = genotype_df.columns.intersection(pheno.columns).intersection(covar.index)
    genotype_df, pheno, covar = genotype_df[samples], pheno[samples], covar.loc[samples]
    alleles = pd.DataFrame(
        {"chrom": "1", "pos": 0, "ref": "A", "alt": "G", "ot_variant_id": "", "flip": False},
        index=variant_df.index,
    )

    cs_df, summary_df = susie_map.fine_map(
        genotype_df, variant_df, pheno, pheno_pos, covar, None, alleles, maf, 1000000)
    ref_summary, ref_res = susie.map(
        genotype_df, variant_df, pheno, pheno_pos, covar, summary_only=False,
        maf_threshold=maf, window=1000000)

    ours = set(zip(cs_df.phenotype_id, cs_df.cs_id.astype(str), cs_df.variant_id))
    theirs = set(zip(ref_summary.phenotype_id, ref_summary.cs_id.astype(str), ref_summary.variant_id)) \
        if len(ref_summary) else set()
    assert ours == theirs, f"credible sets differ: {len(ours - theirs)} ours only, {len(theirs - ours)} theirs only"
    if len(cs_df):
        m = cs_df.merge(ref_summary, on=["phenotype_id", "variant_id"], suffixes=("", "_ref"))
        np.testing.assert_allclose(m.pip, m.pip_ref, atol=1e-5)
    for pid, r in ref_res.items():
        alpha_ref = susie_tables.alpha_from_lbf(r["lbf_variable"])
        mine = cs_df[cs_df.phenotype_id == pid]
        for cs_id in sorted(set(mine.cs_id)):
            sub = mine[mine.cs_id == cs_id]
            ix = r["pip"].index.get_indexer(sub.variant_id)
            # lbf_variable rows are the credible-set signals in cs_index order
            l = list(np.asarray(r["sets"]["cs_index"])).index(cs_id - 1)
            np.testing.assert_allclose(sub.alpha, alpha_ref[l, ix], atol=1e-5)
    print(f"EQUIVALENT n_phenotypes={len(pheno)} n_cs={len({(p, c) for p, c, _ in ours})} "
          f"statuses={summary_df.status.value_counts().to_dict()}")


if __name__ == "__main__":
    main(*sys.argv[1:])
