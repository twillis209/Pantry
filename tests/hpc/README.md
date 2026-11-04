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

## Gates 1-3 (CPU, controller env, standalone, old-format config)
Run records from the Task 5/7/11 reports, so every result cited in `docs/module-status.md` is traceable here. All in `$H/pantry-runs/<dir>`,
`config.yml` = `tests/configs/{phenotyping,pheast}_baseline.yml` (old format), `input` -> the test data, no `--sdm` (every rule in the
controller env, as the baseline), `sjob -c 2 -m 7500M -t 0-12:00:00 -n <name> -e <env> snakemake -s $H/Pantry/<stage>/Snakefile --profile $H/Pantry/tests/profile`.
- Gate 1, `gate1-pheno` (Task 5, `$H/Pantry` at 8aba4e6; pathvars, outputs still at `output/`), env pantry-baseline-pheno: dry run 8392776 (341 jobs),
  run 8392777 (341/341, ~16 min). Compare 8393107 (no `--noise`): exit 1, 15/21 byte-identical, the 6 kallisto-derived files differ (`compare.md`).
  Substitution check (copy `gate1-pheno-subst`, since removed): the baseline's kallisto outputs for the 5 divergent samples swapped in, run 8393112
  reran only the 11 downstream jobs, compare 8393118 -> MATCH (accepted under R11). With `--noise`: job 8393447 (Task 2b above) -> MATCH (`compare_noise.md`).
- Gate 2, `gate2-pheno` (Task 7, `$H/Pantry` at eb9862e; `{dataset}` layout, outputs under `output/default/`), env pantry-baseline-pheno: run 8393449
  (343/343, 12m25s). Compare 8393775 (`compare.sh`: `compare_outputs.py archive/phenotyping output/default --noise archive/phenotyping
  archive/phenotyping_run2`) -> exit 0, MATCH (kallisto files within the noise envelope: expression 0.366/0.732, isoforms 0.0835/0.167,
  unnorm alt_polyA 7e-4/2.8e-3, unnorm expression 338/676, unnorm isoforms 3.0e-3/6.1e-3; all other files identical). Report `compare.md`.
- Gate 3, `gate3-pheast` (Task 11, `$H/Pantry` at 0d21fc6; `{dataset}`, checkpoints, `install_fusion`), env pantry-baseline-pheast: run 8394390
  (195/195, 1h55m59s; cold dry run 69 jobs, 195 = upstream's 185 + check_genotype_samples, 8 phenotype_summary, install_fusion).
  Rerun check 8394625 (`rerun_check.sh`, all rerun triggers): nothing to be done. Compare 8394624 (`compare.sh`, `archive/pheast_run1` vs
  `output/default`, `--noise run1 run2`): exit 1, 10 mismatches, all in `*.cis_qtl`/`*.cis_independent_qtl`; heritability x8, TWAS x8
  (`.wgt.RDat` equal via Rscript), genotype_qc identical (`compare.md`; raw `compare_raw.md`; vs run2 8394627 `compare_vs_run2.md`: 17, again
  tensorQTL only). Attribution (R21): covar inputs byte-identical to `$H/pantry-baseline/pheast/intermediate/covar` (only `geno_pruned.log`
  differs); job 8394629 (`cmds.sh`, `-n -p --forcerun tensorqtl_cis tensorqtl_cis_independent qtl_genotype_qc` in both dirs): the 17
  run_tensorqtl.py/genotype_qc.py command lines identical after normalising paths (`cmds_baseline.norm` vs `cmds_gate3.norm`).

## Task 12: per-rule conda envs, tensorQTL in a GPU container
(The francois4 image below was replaced in Task 16 by `ghcr.io/twillis209/pantry-tensorqtl`, tensorqtl 1.0.10.)
Envs in `phenotyping/envs/` and `pheast/envs/` pin every tool/library to its version in the baseline lock files; the tensorQTL rules
use `container: docker://francois4/tensorqtl@sha256:f6efb9e592eb32c46cb75070be2769b34381d60cbb2709d2885771324abfe32a`
(genex's `:latest`, last pushed 2024-01-06, resolved via the Docker Hub API). The test profile sets `software-deployment-method: [conda, apptainer]`,
`apptainer-args: "--bind=/data/hps/assoc/private/gdml --nv"`, shared `conda-prefix: $H/pantry-runs/conda-envs` and
`apptainer-prefix: $H/pantry-runs/apptainer` (so every run dir reuses the envs and image), and `set-resources` for
tensorqtl_cis/_cis_independent/_trans: `slurm_partition: gpu-core`, `slurm_account: gdml`, `gpu: 1`. genex's `slurm_extra: "'--gpus=1'"` is rejected
by the current slurm executor plugin ("The --GPU-allocation option is not allowed in the 'slurm_extra' parameter").
`sjob` has no GPU option; for a manual GPU job use `SBATCH_GPUS=1 sjob -p gpu-core -a gdml ...` (sbatch reads `SBATCH_GPUS`).

- Image contents (job 8394915, A100 node, driver 610.43.02): Ubuntu 22.04, Python 3.10.12, tensorqtl 1.0.9, torch 2.1.2+cu121, pandas 2.1.4,
  numpy 1.26.3, scipy 1.11.4, pandas_plink 2.2.9, rpy2 3.5.15 (R qvalue loads), dask 2023.12.1; `torch.cuda.is_available()` True.
  Baseline env: tensorqtl 1.0.10, torch 2.3.1, pandas 2.3.3, Python 3.12. Local copy: `$H/pantry-runs/containers/tensorqtl_f6efb9e5.sif`.
- Env builds (`--conda-create-envs-only` via the profile): 8394943 pheno 7m (9 envs), 8394944 pheast 6m (5 envs, MaxRSS 11 GB with `-c 4`).
  Built envs vs baseline lock (`t12_envdiff.sh`, job 8394950): all tools and pinned libs identical; only transitive/system libs differ.
- Gate 2 conda (`$H/pantry-runs/gate2-pheno-conda`, job 8394948, 17m52s, 343/343): compare job 8395584
  `--noise archive/phenotyping archive/phenotyping_run2` -> exit 0, MATCH (same notes as gate 1/2: kallisto files within the noise envelope).
- Gate 3 conda+GPU (`$H/pantry-runs/gate3-pheast-conda`): 8394949 ran everything but tensorQTL (GPU submission rejected, see above), resumed
  as 8395487 (4m25s) after the profile fix. tensorqtl_cis jobs 1m15-1m37 and cis_independent 26-55 s on A100s (whole tensorQTL step ~4 min vs
  ~2 h on CPU in the baseline); a manual run of the image (job 8395587) logs `* using GPU (NVIDIA A100-SXM4-80GB)` and computes q-values.
  Compare (job 8395585): heritability x8, TWAS x8, genotype_qc identical (all blocking files pass); 85 mismatches, all in
  `*.cis_qtl`/`*.cis_independent_qtl` (non-blocking per R22). R22(a): covar inputs (8 covar.tsv, 8 covar.plink.tsv,
  geno_pruned.{bed,bim,fam,prune.in,prune.out}) byte-identical to the baseline intermediates; the 17 run_tensorqtl.py/genotype_qc.py command
  lines identical to the baseline's (`cmds_gate3.norm` vs `gate3-pheast/cmds_baseline.norm`).
- How far GPU tensorQTL is from the CPU baseline (compare.md, `--noise run1 run2`):
  - cis_qtl (8 files): same top variant in most rows, with float-level differences: pval_nominal differs in ~120-130 rows/file (max rel
    ~1e-4), slope/slope_se in the 6th-7th significant digit; 4-18 rows/file pick a different top variant with a different pval_nominal
    (so start/end_distance, and af once, change), besides 4-10 rows/file of equal-pval tie-breaks (also in the noise pair). Permutation
    columns stay within the run1-run2 noise, except intron_retention qval (77 rows above the 0.0134 noise, e.g. 0.286 vs 0.313).
  - cis_independent (8 files): median same-variant rel diff of pval_nominal 3-5e-6 (noise 0); row-set changes of 1 key (RNA_editing,
    alt_TSS, splicing) and 7 rows / 5 phenotypes (intron_retention; noise 2 / 1); expression and intron_retention show the large
    pval/slope/num_var changes of phenotypes whose signal set changed, as in the CPU noise pair (Task 2).
  - Version vs GPU (expression cis only, job 8395593: the image on CPU, `* WARNING: using CPU!`, 1.7 min; `containers/cmp3.py`, job 8395594):
    baseline run1 vs run2: pval_nominal/slope/slope_se/af identical where the variant agrees (19/200 variant_id differences = ties);
    image-on-CPU vs baseline: pval_nominal differs in 104 rows (max rel 4.4e-5), slope 11, slope_se 4; image-on-GPU vs baseline: 116 / 18 / 6
    rows (max rel 5.9e-5); GPU vs image-on-CPU: 102 rows (8.2e-5). So the image's software stack (tensorqtl 1.0.9 / torch 2.1.2 / pandas 2.1.4)
    already moves results at float level on CPU, and the GPU adds a comparable float-level difference; neither changes af or the analysis.
    (Task 16 correction: that CPU check ran on 4 cores; the baseline's tensorQTL jobs had 1. With 1 core the 1.0.10 image matches the
    baseline exactly, and 8 cores give differences like these, so thread count, not the version, likely explains them.)

## Task 13: multi-dataset acceptance (criterion 2)
Configs `tests/configs/{phenotyping,pheast}_multi.yml`; derived inputs from `tests/fixtures/make_multi_inputs.py`, run on the HPC
against the test data into the run dirs (job 8395740, env pantry-baseline-pheno for bgzip/tabix):
`make_multi_inputs.py phenotyping input multi_inputs/half_se` (first 8 samples, R1-only map) and
`make_multi_inputs.py pheast input multi_inputs/subset` (first 300 samples; 8 BEDs cut to them, bgzip + tabix; 6 groups files).
Run dirs `$H/pantry-runs/multi-pheno` and `multi-pheast` (`input` link, config.yml = the multi config), same profile, prefixes and
controller envs as Task 12: `sjob -c 2 -m 7500M -t 0-12:00:00 -e <env> bash t13_run.sh <stage> <rundir>`.
`tests/fixtures/check_targets.py <stage> config.yml <dataset> output/<dataset>` checks every `rule all` target of a dataset
(exists, non-empty, >= 1 data row / archive member).
- Phenotyping 8395741: 19m15s, 468/468, EXIT=0. Job stats: `star_index 1`, `kallisto_index 1` (alt_TSS_polyA_kallisto_index 4 = 2 groups x
  2 positions, reference-level), `write_samples 2`; `intermediate/reference/grch38_chr1/` has one `star_index_75/` and one `kallisto.idx`.
  Compare 8396520: `full` vs `archive/phenotyping` `--noise archive/phenotyping archive/phenotyping_run2` -> exit 0, MATCH (kallisto files within
  the noise envelope: expression 0.366/0.732, isoforms 0.149/0.167, unnorm expression 440/676). check_targets: full 20/20, half_se 20/20
  (half_se rows: alt_TSS 166, alt_polyA 137, expression 106, isoforms 129, RNA_editing 11, splicing 72, stability 15).
- Pheast 8395742: 11m02s, 388/388, EXIT=0 (tensorQTL on GPU). Compare 8396510: `full` vs `archive/pheast_run1` `--noise run1 run2` ->
  76 mismatches, all in `*.cis_qtl`/`*.cis_independent_qtl` (non-blocking, R22, as in gate 3 conda); heritability x8, TWAS x8, genotype_qc pass.
  R22(a): `intermediate/full/covar` (covar.tsv x8, covar.plink.tsv x8, geno_pruned.{bed,bim,fam,prune.in,prune.out}) byte-identical to the
  baseline's; the 17 `full` run_tensorqtl.py/genotype_qc.py command lines identical to `gate3-pheast/cmds_baseline.norm`.
  check_targets: full 33/33, subset 33/33 (subset covariates have 300 samples; genotype_qc threshold 20).

## Task 14: host workflow test (criterion 3)
`tests/host/Snakefile` imports both stages by local path (spec §1 pathvars under `results/`, `resources/`); `tests/host/config.yaml`
has phenotyping `full` (+ cross_modality over the 7 modalities) and pheast `geuvadis` (test data) and `chained`
(`phenotype_dir: results/pantry/phenotypes/full`, 16 phenotyping samples, 7 modalities, all analyses; dry run only). `chained` has no
cross_modality: pheast allows that modality only with `qtl` alone, and the checkpoint is used only by heritability/TWAS.
`tests/host/profile` = `tests/profile` with set-resources keyed by the host rule names (`pantry_pheast_tensorqtl_*`).
Run dir `$H/pantry-runs/host` (`input/{phenotyping,pheast}` links, config.yaml copied), controller env `pantry-baseline-pheno`
(Snakemake 9.27.1; the phenotyping module needs gtfparse at parse time): `sjob -c 2 -m 7500M -t 0-12:00:00 -e pantry-baseline-pheno
$H/pantry-runs/t14_run.sh <rundir> <host Snakefile> [targets/args]`.
- Dry run of `all` (8396530): 472 jobs, EXIT=0. `pantry_pheast_covariates` for `chained` takes `results/pantry/phenotypes/full/expression.bed.gz`
  produced by `pantry_pheno_normalize_expression` in the same DAG; 7 `chained` `phenotype_summary` checkpoints pending; heritability/TWAS
  batch jobs deferred on them.
- Run of `tested` (= phenotyping `full` + pheast `geuvadis`) 8396531: 17m44s, 539/539, EXIT=0; 16 tensorQTL jobs on `gpu-core` with
  `gres/gpu:a100=1` (the prefixed set-resources keys apply).
- Compare 8398057 (`$H/pantry-runs/host/compare.sh`): phenotyping `full` (symlink view without `cross_modality.*`, absent from the baseline)
  vs `archive/phenotyping` `--noise archive/phenotyping archive/phenotyping_run2` -> exit 0, MATCH; check_targets 23/23.
  `cross_modality.bed.gz` (decompressed) and `.phenotype_groups.txt` byte-identical to upstream `scripts/combine_modalities.sh` (8da7787)
  run on the same 7 BEDs (771 phenotypes). Pheast `geuvadis` vs `archive/pheast_run1` `--noise run1 run2`: 79 mismatches, all in
  `*.cis_qtl`/`*.cis_independent_qtl` (non-blocking, R22); heritability x8, TWAS x8, genotype_qc pass; check_targets 33/33. R22(a): covar
  inputs byte-identical to the baseline's (only `geno_pruned.log` differs); the 17 run_tensorqtl.py/genotype_qc.py command lines identical
  to the baseline's (`-p` prints each command twice in this log, so `sort -u` first).
- github(), first attempt (tag `module-rc1` then 56c8c22): dry run 8398062 failed at parse with the Snakemake bug below (pheast's
  install_fusion called `source_path` on `gcta_nr_robust` at parse time). Phenotyping alone via github() worked (scratch
  `host-github/Snakefile.pheno_only`, job 8398063: dry run 345 jobs; `preprocess_gtf_for_txrevise`, SLURM job 8398064, EXIT=0, output
  byte-identical to the local-path host's).
- Fix (ruling R24, a140598): gcta is resolved only in install_fusion's run block: configured `gcta` -> nothing installed for it; module
  loaded from a local path (its source dir exists on disk) -> bundled file copied; else downloaded from `gcta_url` (default: upstream
  8da7787's file on raw.githubusercontent.com) and checked against `gcta_sha256` (7b51e51d...93e9); then chmod +x.
  Tag `module-rc1` moved to a140598 (force-pushed). Stale fork entries removed from `~/.cache/snakemake/snakemake/source-cache/`
  ({https/raw.githubusercontent.com,hosted-git/githubfile/github.com,snakemake-git-cache/github.com}/twillis209): the cache is keyed by tag.
- Re-verification, job 8398071 (`$H/pantry-runs/t14_r24.sh`): local-path host dry run of `all` EXIT=0, 60 jobs, all `chained` (nothing of
  `full`/`geuvadis` reruns); `pantry_pheast_install_fusion` in a fresh dir (`host-fusion-local`) EXIT=0, 4 files byte-identical to the
  host run's, gcta executable. github() host (`Snakefile.github`): dry run of `all` EXIT=0 (470 jobs); run of `pantry_pheast_install_fusion`
  + `resources/pantry/reference/grch38_chr1/exonic.gtf` (`exonic_intronic_from_GTF`, script via source_path) EXIT=0: gcta_nr_robust
  downloaded (absent from the source cache), sha256 7b51e51d...93e9, executable, identical to the bundled one; exonic.gtf identical to the
  local host's. The 3 FUSION R scripts differ from the local copies only by the missing final newline (second bug below).

### Snakemake source-cache bugs (to report upstream; not filed)
Snakemake 9.27.1 (HPC env pantry-baseline-pheno; same code in 9.27.0 locally), `snakemake/sourcecache.py`, `HostingProviderFile.open`
(line 595; `_open` for `LocalGitFile` at line 883 has the same pattern): `io.BytesIO(self.hosted_repo.repo.git.show(f"{ref}:{path}").encode())`.
GitPython returns the blob as str (decoded with surrogateescape, trailing newline stripped), so
1. binary files fail: `UnicodeEncodeError: 'utf-8' codec can't encode character '\udc87' in position 42: surrogates not allowed`
   (any `workflow.source_path` on a binary in a github()/gitlab() module);
2. text files lose their final newline (FUSION R scripts installed from github() are 1 byte shorter).
Likely fix: `git.show(..., stdout_as_string=False, strip_newline_in_stdout=False)` and no `.encode()`.

## Task 16: tensorQTL 1.0.10 image on GHCR
`containers/tensorqtl/` (Dockerfile on `mambaorg/micromamba` 2.10.0 by digest; `environment.yml` conda pins and `constraints.txt` = `pip list
--format=freeze` of `pantry-baseline-pheast`; `check_image.py`) is built by `.github/workflows/tensorqtl-image.yml` (linux/amd64, on push to
`containers/tensorqtl/**`) and pushed to `ghcr.io/twillis209/pantry-tensorqtl` (`:1.0.10`, `:sha-<short>`; public).
Digest `sha256:65a42f9fd74a4bc5f0f6a1a2a88a19e1140b62bbef5b3f2e6fc9de4330bde9b8` (run 41716028853, from 206c892), pinned in `pheast/steps/qtl.smk`.
As in the baseline env, pip replaces conda's numpy 1.26.4 with the numpy 2.5.4 wheel and adds scipy 1.18.1. Only pyarrow differs
(25.0.1 vs 26.0.0: the lock's pyarrow 26 needs numpy>=2, dask 2024.1.1 numpy<2; used only by cis_nominal).
- Pull: this account's stored ghcr.io credentials (`~/.apptainer/remote.yaml`, `~/.docker/config.json`) make apptainer get 403 on the
  bearer token even though the image is public (jobs 8399091, 8399106). `$H/pantry-runs/t16/pull.sh` (job 8399094) pulls anonymously with a
  clean HOME/DOCKER_CONFIG into the profile's apptainer-prefix as `7b93bcb6ac4e77c33b84b5d040ff650a.simg` (= md5 of the URI, Snakemake's
  cache name), so runs never pull. Fix for fresh pulls: `apptainer registry logout docker://ghcr.io` and drop ghcr.io from `~/.docker/config.json`.
- `check_image.py` on an A100 (`SBATCH_GPUS=1 sjob -p gpu-core -a gdml`, `--nv`; job 8399095): Python 3.12.16, tensorqtl 1.0.10,
  torch 2.3.1+cu121, pandas 2.3.3, numpy 2.5.4, scipy 1.18.1; rpy2 + qvalue load; `torch.cuda.is_available()` True.
- Gate 3 GPU (`$H/pantry-runs/gate3-pheast-gpu1010`, job 8399104, 17m15s, 195/195; tensorQTL jobs 28-55 s on A100s). Compare 8399295:
  81 mismatches, all in `*.cis_qtl`/`*.cis_independent_qtl` (non-blocking, R22); heritability x8, TWAS x8, genotype_qc pass; check_targets
  33/33. R22(a): covar inputs byte-identical to the baseline's (only `geno_pruned.log` differs); 17 command lines identical. Same pattern
  as with the old image: pval_nominal differs in 103-117 rows/cis file (max rel ~1e-4).
- Image on CPU vs baseline CPU (`$H/pantry-runs/t16/cpucis*.sh`, `cmp16.py`; same command lines and covar inputs; nominal stats on
  same-variant rows):
  - 1 core, as the baseline's tensorQTL jobs (AllocCPUS 1; jobs 8399314, 8399317-8399323; compare 8399443): pval_nominal, slope, slope_se,
    af identical in all 8 modalities (0 rows differ). Top-variant ties and permutation columns vary as between baseline run1 and run2.
  - 8 cores (job 8399105): pval_nominal differs in 35-90 rows/file (max rel 1e-5..1e-4), about as much as on the GPU.
  - So the remaining GPU-vs-baseline differences come from parallel reduction order (GPU or multithreaded CPU), not the software stack.

## Task 17: MAJIQ (intron_retention) at run time
- (Task 20 moved `tests/hpc/envs/majiq3{.yml,-install.sh,.lock.txt}` to `containers/majiq/`.)
- Env `majiq3` (named, as `majiq_env` defaults): `tests/hpc/envs/majiq3.yml` (conda part: python 3.12, htslib, zstd, xz/liblzma-devel,
  bzip2, zlib, cmake, compilers) + `tests/hpc/envs/majiq3-install.sh` (clones the public https://bitbucket.org/biociphers/majiq_academic at
  tag v3.0.22, commit d841baad, sets `HTSLIB_*` and `CMAKE_PREFIX_PATH` to the env, `pip install .` as in its README). No credentials needed.
  Run in `$H/pantry-runs/t17`: `sjob -c 4 -m 7500M -t 0-02:00:00 -e base "bash majiq3-install.sh $PWD/build"` (job 8399463, 5 min; 8399456
  failed: without the conda xz headers CMake picked the system liblzma, which has no `lzma.h`). Lock: `tests/hpc/envs/majiq3.lock.txt`
  (pip packages as trailing comments: rna_majiq/rna_voila/rna_moccasin 3.0.22, numpy 2.0.1, zarr 2.18.3, xarray 2024.6.0, dask 2024.7.1).
- Licence: the public build sets `LICENSE_CHECKING_TYPE = "None"` (`voila/rna_voila/api_static/licensing.py`), so `--license` is accepted
  but not read. The bundled `majiq_license_academic_official.lic` decodes to an academic licence without expiry.
- Config `tests/configs/phenotyping_baseline_ir.yml` = upstream `phenotyping/config.yml` (intron_retention kept, LaDDR off). Targets only
  `intron_retention.bed.gz{,.tbi}` and `intron_retention.phenotype_groups.txt` (refactor: under `output/default/`).
- Baselines `$H/pantry-baseline/phenotyping_ir{,_run2}` (`git archive 8da7787 phenotyping`, `input` link), env pantry-baseline-pheno,
  `--profile tests/profile` unchanged: upstream has `conda:` only on the MAJIQ rules, so the named env is activated for them alone and every
  other rule runs in the controller env, as in the first baseline. `$H/pantry-runs/t17/run.sh <rundir> <Snakefile> [targets]`.
  Dry runs 8399457 (59 jobs) / 8399458 (gate-ir: same + validate_annotation, write_samples; MAJIQ command lines equal up to path prefixes).
  STAR/BAM steps 8399464-6, IR steps 8399661 / 8399662 (baselines), 8399663 (gate-ir, module 9464ad2): all EXIT=0, MAJIQ at 8 threads.
  Archived to `archive/phenotyping_ir{,_run2}` (bed.gz, .tbi, phenotype_groups.txt, unnorm/intron_retention.bed).
- MAJIQ noise: `majiq.psi.tsv` differs run to run in the bootstrap columns only (`bootstrap_psi_std`, `psi_quantile_*`; up to ~0.1),
  in run1 vs run2 and in gate-ir vs run1 alike. `extract_ir_psi.py` uses `raw_psi_mean`, which is identical, so every final file is
  byte-identical across the two baselines and gate-ir.
- Compare 8399692 (`gate-ir/compare.sh`): run1 vs run2 MATCH; gate-ir vs run1 byte-exact (no `--noise`) MATCH; with `--noise` MATCH.
- Symlinked BAMs (fix 21ec32d): `$H/pantry-runs/gate-ir-symlink` = `cp -a` of gate-ir with `bam/{s}.bam{,.bai}` -> `../bam_real/{s}.renamed.bam{,.bai}`
  and the IR outputs removed. Dry run 8399693: 9 IR jobs, `.sj` paths `build/{s}.renamed.sj`. Run 8399696 (9/9, EXIT=0): MAJIQ wrote
  `build/{s}.renamed.sj`, psi-coverage read them, and the 3 final files are byte-identical to gate-ir's.

## Task 18: github() full host run at module-v0.1.1; dry-run scaling
- Tag `module-v0.1.1` (annotated) = befa72c. Fork entries cleared from `~/.cache/snakemake/snakemake/source-cache/` first
  (`hosted-git/githubfile/github.com/twillis209`, `snakemake-git-cache/github.com/twillis209`).
- Run dir `$H/pantry-runs/host-github-v011`: `config.yaml` = `tests/host/config.yaml` plus an `intron_retention` group for `full`
  (`psicov_batch_size: 25`, `min_experiments: 0.01`, as `tests/configs/phenotyping_baseline_ir.yml`; `$H/pantry-runs/t18/mk_gh_config.py`), env
  pantry-baseline-pheno. Job 8399793 (`$H/pantry-runs/t18/gh_run.sh`, `sjob -c 2 -m 7500M -t 0-12:00:00`): fetched befa72c; dry run of `all`
  482 jobs (Task 14's 470 + intron_retention), EXIT=0; run of `tested` EXIT=0, 549/549 steps, 1h06m (incl. 13 conda envs rebuilt: env files
  fetched by github() lack the final newline, so their hashes differ from the local-path ones; source-cache bug 2 above), 16 tensorQTL jobs
  on `gpu-core` (`gres/gpu:a100=1`). gcta downloaded from the full-SHA `gcta_url`, sha256 7b51e51d...93e9.
  - Job counts: `tested` started at 423 jobs = 482 minus the 59 `chained` jobs (7 modalities x 8 per-modality pheast rules + check_genotype_samples,
    prune_for_covar, qtl_genotype_qc); the checkpoints added 126 (geuvadis' 8 modalities: heritability batches 80, TWAS weight batches 37,
    twas_geno 1, plink_covariates 8, as gate 3's 69 -> 195), giving 549. The in-script dry run of `tested` failed (exit 2) because
    `--quiet rules tested` takes `tested` as a `--quiet` choice; after the run, `t18/dry_tested.sh` (job 8400337, target before `--quiet`): EXIT=0,
    nothing to be done.
- Compare 8400338 (`host-github-v011/compare.sh`, env pantry-baseline-pheast; Task 14's script plus intron_retention):
  - phenotyping `full` (view without `cross_modality.*`, `intron_retention.*`) vs `archive/phenotyping` `--noise`: exit 0, MATCH (alt_polyA
    0.213/0.328, expression 0.366/0.732, isoforms 0.0835/0.167, unnorm alt_polyA 4.8e-4/2.8e-3, unnorm expression 338/676, unnorm isoforms 3.0e-3/6.1e-3);
    check_targets 26/26; `cross_modality` bed and groups byte-identical to upstream `combine_modalities.sh` (771 rows).
  - intron_retention vs `archive/phenotyping_ir` without `--noise`: exit 0, MATCH; all 4 files byte-identical (`cmp`).
  - pheast `geuvadis` vs `archive/pheast_run1` `--noise run1 run2`: 80 mismatches, all in the 16 `*.cis_qtl`/`*.cis_independent_qtl` files
    (non-blocking, R22); heritability x8, TWAS x8, genotype_qc pass; check_targets 33/33. Covar inputs: 21 byte-identical, only `geno_pruned.log` differs.
  - Command lines (job 8400342, `cmds.sh` in env pantry-baseline-pheno; the same step in compare.sh failed for lack of gtfparse in the pheast env):
    17/17 run_tensorqtl.py/genotype_qc.py lines identical to the baseline's.
  - vs the Task 14 local-path host run: every phenotyping file byte-identical except `unnorm/alt_polyA.bed` (kallisto noise), every pheast
    non-tensorQTL file byte-identical.
- Scaling (dry runs, `$H/pantry-runs/t18/scale.sh <name> <datasets> <samples> <mode>`, `tests/fixtures/make_scaling_config.py`, env
  pantry-baseline-pheno, Snakemake 9.27.1, `-n --quiet rules` (no target; a target must come before `--quiet rules`, which otherwise takes it as a choice), `/usr/bin/time -v`, `sjob -c 2`). Synthetic datasets = copies of the host's `full`
  (7 modalities + cross_modality) with inline sample lists and FASTQ maps of empty files; pheast datasets chained to them (cold, so the
  checkpoint-dependent jobs are not counted). Phenotyping-only 1-dataset runs use the old-format baseline config with file maps.

  | Case | Job | Jobs | Wall | Max RSS |
  |---|---|---|---|---|
  | upstream 8da7787 standalone phenotyping, 1 x 1000 | 8399787 | 12149 | 32.8 s | 393 MB |
  | module standalone phenotyping, 1 x 1000 | 8399786 | 12151 | 38.0 s | 401 MB |
  | host, 17 x 100 (phenotyping + pheast) | 8399785 | 22053 | 63.6 s | 628 MB |
  | host, 17 x 300, under cProfile | 8399788 | 62853 | 6m43 | 1.75 GB |
  | host, 17 x 1000 (`-m 16000M`) | 8399789 | 205653 | 12m20 | 4.38 GB |

  About 3 ms per job throughout (2.7 upstream, 3.1 ours at 1 x 1000, 3.6 at 17 x 1000), so time and memory grow roughly linearly with
  jobs. Profile (`t18/p17x300/prof.txt`, `t18/ours.py`): Snakemake internals dominate: dry-run job printing evaluates resources
  (`expand_resources` 170 s of 393 s, much of it `inspect.signature`), DAG building 125 s, input sizes 63 s, conda env specs 38 s. All of
  the workflow's own functions together take under 6 s (largest: `parent_dir` 1.6 s over 56k calls); no helper of ours dominates.

## Task 19: FUSION image replaces install_fusion (R26, module-v0.1.2)
- Image `ghcr.io/twillis209/pantry-fusion@sha256:87335916d0395fe07cb6abe01e6cfa52838a85cd8bf538f32bfa11c17fdb4009` (tag `sha-11d208d`),
  built by `.github/workflows/fusion-image.yml` run 41755433329 from `containers/fusion/` (context = repo root; `Dockerfile.dockerignore`
  sends only `pheast/envs/fusion.yaml`, `pheast/scripts/fusion_twas/`, `check_image.py`). First run 41755267647 failed in `check_image.py`
  (micromamba 2's `list --json` format); fixed by reading `conda-meta`. Build gate: python 3.12.16, r-base 4.5.3, r-glmnet 5.0, r-optparse 1.8.2,
  r-matrix 1.7_6, plink 1.90, gemma 0.98.5; gcta sha256 OK, `ldd`: not a dynamic executable (static), runs; R packages load; plink, gemma run.
- Anonymous pull from the HPC (token + manifest by digest): HTTP 200, so public without a visibility change. `$H/pantry-runs/t19/pull.sh`
  (job 8400427, as `t16/pull.sh`) pulled it into the profile's apptainer-prefix as `46cbe05b95ea22cb3a829882a86c1c08.simg` and ran
  `check_image.py` in it: same versions, CHECK_EXIT=0.
- Cold dry runs of `all` (`t19/dry.sh`, env pantry-baseline-pheno; standalone pheast on `tests/configs/pheast_baseline.yml`, local-path host
  `tests/host/Snakefile` on `tests/host/config.yaml`): before (`$H/Pantry` at 474d3ef, code = 27860e1; job 8400376) 69 and 472 jobs, each with
  one `install_fusion`; after (ffe95c5; job 8400428) 68 and 471. Exactly the one `install_fusion` job fewer.
- Gate-3 config run `$H/pantry-runs/gate3-pheast-fusion` (`t12_run.sh pheast`, `sjob -c 2 -m 7500M`, env pantry-baseline-pheno, `tests/profile`:
  `--sdm conda apptainer`, tensorQTL on GPU). The first attempt (8400482, at ffe95c5; logs in `t19/failed1`) failed in all 8 `covariates` jobs:
  `covariates.R` sources `fusion_twas/utils/plink_utils.R` relative to itself, and only `install_fusion`'s params had put that file into the
  source cache. Fixed in 6c5aaf6 (`covariates` params `source_path` it). Run 8400497 at 6c5aaf6: 194/194 (195 minus `install_fusion`), EXIT=0;
  TWAS jobs ran in the image (`Activating singularity image …/46cbe05b….simg`), the cached scripts under `$HOME` being bound by apptainer by default.
- Compare 8400498 (`gate3-pheast-fusion/compare.sh`, env pantry-baseline-pheast; `compare_outputs.py archive/pheast_run1 output/default --noise
  run1 run2`): exit 1, 81 mismatches, all in the 16 `*.cis_qtl`/`*.cis_independent_qtl` files (non-blocking, R22); heritability x8,
  TWAS x8 (every tar member byte-identical except `.profile.err`, `.wgt.RDat` all.equal at 1e-8) and genotype_qc pass; check_targets 33/33.
  TWAS intermediates `twas/{modality}.{profile,pos,list}` vs gate 3: 24/24 byte-identical vs `gate3-pheast-gpu1010` and 24/24 vs `gate3-pheast` (CPU).
  Covar inputs: 21 byte-identical to the baseline's, only `geno_pruned.log` differs. Job 8400678 (`cmds.sh`): 17/17 run_tensorqtl.py/genotype_qc.py
  command lines identical to the baseline's; rerun check: nothing to be done.
- Tag `module-v0.1.2` (annotated) = 0d30215 (code and docs 8c43615 + `Snakefile.github` re-pointed). `github()` host dry run of `all` at the tag
  (`t19/gh_dry.sh`, job 8400681; fork entries cleared from the source cache first; run dir `t19/dry-host-github`): fetched 0d30215, EXIT=0,
  471 jobs (Task 18's 482 at v0.1.1 included 12 intron_retention jobs from its own config; the local-path host on the same config: 471).
  No `UnicodeEncodeError`; the 62 fetched files (`hosted-git/…/module-v0.1.2/`) are all text except 3 `__pycache__/*.pyc` written locally
  by Python imports; from `fusion_twas/` only `utils/plink_utils.R` (for `covariates`) is fetched, no gcta.
- Full `github()` host run at `module-v0.1.2` (= 0d30215; local clone f9c2345 for the host Snakefile, profile and tools). Run dir
  `$H/pantry-runs/host-github-v012`, config as `host-github-v011` (phenotyping `full` + intron_retention, pheast `geuvadis`). Job 8402594
  (`t19/gh_run.sh`, `sjob -c 2 -m 7500M -t 0-12:00:00`, env pantry-baseline-pheno; targets now before `--quiet rules`): DRY_ALL_EXIT=0,
  DRY_TESTED_EXIT=0, RUN_EXIT=0, 548/548 steps (Task 18's 549 minus `install_fusion`).
- Compare 8403790 (`host-github-v012/compare.sh`, as Task 18's, env pantry-baseline-pheast): PHENO_COMPARE_EXIT=0, PHENO_CHECK_EXIT=0;
  intron_retention IR_COMPARE_EXIT=0 and IR_BYTES_SAME x4 vs upstream; XM_BED_SAME, XM_GROUPS_SAME (771 rows); PHEAST_COMPARE_EXIT=1 with
  84 mismatches, all in the 16 `*.cis_qtl`/`*.cis_independent_qtl` files (non-blocking, R22; cis_independent row counts within 2x noise);
  heritability x8, TWAS x8 (tar members, RDat 1e-8) and genotype_qc pass; PHEAST_CHECK_EXIT=0; covar inputs 21 SAME, `geno_pruned.log` DIFF
  (as before). The in-script command check fails in this env (no gtfparse, as Task 18); job 8403792 (`cmds.sh`, env pantry-baseline-pheno):
  CMDS_EXIT=0, 17/17, CMDS_IDENTICAL. HOSTDIFF vs the Task 14 local-path host run: `alt_polyA.bed.gz(.tbi)`, `unnorm/alt_polyA.bed`,
  `cross_modality.bed.gz(.tbi)` (kallisto noise; cross_modality includes alt_polyA; both within the noise envelope in PHENO_COMPARE).
