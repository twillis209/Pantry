# HPC baseline commands

H=/data/hps/assoc/private/gdml/user/twil15. All via `sjob`; ssh prefix: `source /etc/profile; module load slurm; source ~/.dotfiles/.bash_functions;`

- Test data: tarballs from README.md (dropbox, `dl=1`) extracted locally, then `tar | ssh` to `$H/pantry-testdata/{phenotyping,pheast}/input`.
- `git clone https://github.com/twillis209/Pantry.git $H/Pantry -b module` (profile at `$H/Pantry/tests/profile`);
  `git clone ... $H/pantry-baseline && git checkout 8da7787`.
- Envs: committed files `tests/hpc/envs/baseline-{pheno,pheast}.yml` = upstream `environment.yml` (identical at 8da7787) with the
  `defaults` channel removed (SSL errors on repo.anaconda.com, job 8387946) and, for pheast, `pandas<3` (name changed to `pantry-baseline-*`).
  Create (inside an sjob, e.g. `-c 4 -m 7500M -t 0-02:00:00 -e base`; job 8391825 took 1h07):
  `mamba env create -y -f tests/hpc/envs/baseline-pheno.yml --override-channels -c conda-forge -c bioconda` and likewise `baseline-pheast.yml`.
  Why pandas<3: the unpinned env resolved pandas 3.0.6, and tensorqtl 1.0.10 `read_phenotype_bed` then asserts
  "Positions in BED file must be sorted" on the (sorted) test BEDs. The actual HPC env was fixed post hoc with
  `mamba install -n pantry-baseline-pheast --override-channels -c conda-forge -c bioconda "pandas<3"` (job 8392577 -> pandas 2.3.3);
  the committed yml encodes that, so a fresh create should match.
  Exact snapshots of the real envs: `tests/hpc/envs/baseline-{pheno,pheast}.lock.txt` (`mamba list --explicit`; the pheast file also
  carries the pip packages as trailing comments, since `--explicit` omits them, e.g. tensorqtl 1.0.10, torch 2.3.1).
- Stage dirs: `input -> $H/pantry-testdata/<stage>/input`, `config.yml` copied from `tests/configs/*_baseline.yml`;
  `pheast_run2` = `cp -a` of `pheast`.
- Run (in each stage dir): `sjob -c 2 -m 7500M -t 0-12:00:00 -n <name> -e <env> snakemake --profile $H/Pantry/tests/profile`

## Jobs
| Job | Name | Result | Elapsed | MaxRSS (batch) |
|---|---|---|---|---|
| 8391926 | base_pheno | FAILED in 10s (my config: laddr.regress_modalities still listed intron_retention; fixed in tests/configs) | | |
| 8391935 | base_pheno | COMPLETED | 00:13:51 | 333 MB (ReqMem 15000M) |
| 8391927 / 8391928 | base_pheast / base_pheast_run2 | tensorqtl_cis failed (pandas 3), keep-going | ~11 min | |
| 8392578 | base_pheast_b (run1, resumed) | COMPLETED | 01:58:34 | 171 MB |
| 8392579 | base_pheast_run2_b (run2, resumed) | COMPLETED | 01:54:44 | 172 MB |

Archives: `$H/pantry-baseline/archive/{phenotyping,pheast_run1,pheast_run2,phenotyping_reference_files.txt}`.
Pheast modalities in test data: alt_polyA, alt_TSS, expression, intron_retention, isoforms, RNA_editing, splicing, stability (all 8 kept).

## Resume after the pandas failure
The first pheast attempts (8391927 / 8391928) failed only in `tensorqtl_cis` (8 jobs); with keep-going everything else finished, and the
resumed jobs (8392578 / 8392579) re-ran only `all`, `tensorqtl_cis` x8, `tensorqtl_cis_independent` x8. Reused from the first attempt (identically in
both copies): `qtl_genotype_qc`, `prune_for_covar`, `plink_covariates` x8, `covariates` x8, `calculate_heritability_chr_batch` x80,
`combine_heritability_chr` x8, `twas_geno`, `twas_pos_file` x8, `twas_compute_weights_batch` x37, `twas_assemble_summary` x8,
`twas_compress_output` x8 (all heritability, TWAS and covariate outputs; only QTL outputs are from the post-fix env). None of these depend on pandas-sensitive code paths that failed.

## Run1 vs run2
heritability/*.hsq.tsv and qtl/genotype_qc.tsv are identical. `qtl/*.cis_qtl.txt.gz` and `cis_independent` differ
(permutation-based tensorQTL, no seed: e.g. expression cis_qtl differs in beta_shape/pval_perm/qval and even variant_id;
cis_independent row counts 119 vs 120). Later comparisons of QTL outputs need tolerance / structural checks, not exact diffs.

## Task 2: tests/compare_outputs.py on run1 vs run2 (job 8392771, no --noise)
`cis_qtl` files: only the declared noisy columns differ, plus "tie-break?" variant_id differences (equal pval_nominal), 13-24 rows per file.
`cis_independent` files: exact columns that cannot be classified as tie-break/noise (so the comparison design needs a ruling; see
.superpowers task-2 report): (a) same-variant pval_nominal/slope/slope_se float-level differences (rel <= 2.2e-5) in RNA_editing, alt_TSS,
splicing, stability, also at rank 1; (b) row-set differences (alt_polyA 12 keys, isoforms 8, intron_retention 2, expression 1) with large
differences in pval_nominal/slope/num_var for phenotypes whose signal set changed. Heritability, genotype_qc and TWAS are not affected by this.

### Task 2 update: cis_independent judged by summary statistics (controller ruling)
With `--noise`, each cis_independent file is judged on summary stats (rows, distinct phenotype_id, key overlap, same-variant max/median relative
diffs of pval_nominal/slope/slope_se/num_var) against 2x the noise-pair deviation (count floor 1, rel floor 1e-6); without `--noise` rows must match exactly.
- Job 8392773: run1 vs run2 `--noise run1 run2` -> MATCH (0 mismatches), `archive/pheast_noise_calibrated.md`.
- Job 8392774: run2 vs run1, no `--noise` -> MISMATCH (exit 1, expected; raw differences), `archive/pheast_noise.md` (replaces the earlier raw report).

## Spike result (Task 4): checkpoint inside a prefixed module
Snakemake 9.27.0, local run on macOS (throwaway code in /tmp/ckpt-spike, not committed). Host imports the module with
`use rule * from m as m_*` and module-level `pathvars:`. All checks PASS:
- PASS: checkpoint `summary` in the module, read by input function `checkpoints.summary.get(dataset=w.dataset)` (original name, NOT `m_summary`), works under renaming; `snakemake -c1` completes, `res/x/d1/part{0,1,2}.txt` and `gathered.txt` created.
- PASS: checkpoint output path contains a pathvar (`<pantry_output>/{dataset}/n.txt`).
- PASS: module loaded by relative path (`../mod/Snakefile`) and by absolute path; host pathvars referencing other pathvars (`pantry_output="<results>/x"`, `results="res"`) land outputs under `res/x/`.
- PASS: per-rule `localrule: True` inside a rule survives renaming. Checked by `print(r.name, r.name in workflow._localrules)` over `workflow.rules` in the host: only `m_note` is local. (The "localrule" label in `-c1` run logs is NOT discriminating: every rule is labelled local under the local executor.)
- Gotchas: `pathvars:` entries need trailing commas; a top-level `localrule: True` is a silent no-op (parsed as a Python annotation; re-checked: no rule became local), so it must sit inside each rule; the first rule defined is the default target, so with `rule all` after `use rule *` the first imported rule becomes the default; put `default_target: True` on `rule all` as the robust fix.

```
# host/Snakefile
pathvars:
    results="res",
    pantry_output="<results>/x",
rule all:
    input: "<pantry_output>/d1/gathered.txt"
module m:
    snakefile: "/abs/path/mod/Snakefile"
    pathvars:
        pantry_output="<results>/x",
use rule * from m as m_*
# mod/Snakefile
checkpoint summary:
    output: "<pantry_output>/{dataset}/n.txt"
def parts(w):
    n = int(open(checkpoints.summary.get(dataset=w.dataset).output[0]).read())
    return expand("<pantry_output>/{dataset}/part{i}.txt", dataset=w.dataset, i=range(n))
rule note:
    localrule: True
```

## Task 2b: second phenotyping baseline and numeric BED comparison (R12)
kallisto quant is nondeterministic, so with `--noise A B` a BED-like file (`*.bed`, `*.bed.gz`) that is not byte-identical is compared numerically:
header and row count identical, first 4 columns identical row by row, NaN pattern identical, max |diff| over the rest <= 2 x the noise pair's max |diff|
(0 if the noise pair is byte-identical). Without `--noise`: byte equality. Other text files stay byte-exact.
- Second baseline: `$H/pantry-baseline/phenotyping_run2` (rsync copy of `phenotyping` without output/intermediate/logs, input symlink, same config),
  job 8393120 (`base_pheno_run2`, env pantry-baseline-pheno, 00:14:55), archived to `archive/phenotyping_run2` (21 files, same as `archive/phenotyping`).
- Job 8393446: `archive/phenotyping` vs `archive/phenotyping_run2` (no noise) -> 6 files differ: {alt_polyA,expression,isoforms}.bed.gz and
  unnorm/{alt_polyA,expression,isoforms}.bed (alt_TSS, stability, splicing, RNA_editing, phenotype_groups all identical). Report: `archive/phenotyping_noise.md`.
- Job 8393447: `archive/phenotyping` vs gate1-pheno `output` with `--noise archive/phenotyping archive/phenotyping_run2` -> exit 0, MATCH
  (max |diff| / allowed: alt_polyA 0.171/0.328, expression 0.366/0.732, isoforms 0.0835/0.167, unnorm alt_polyA 8.4e-4/2.8e-3, unnorm expression 338/676, unnorm isoforms 3.0e-3/6.1e-3).
  Report: `$H/pantry-runs/gate1-pheno/compare_noise.md`.
