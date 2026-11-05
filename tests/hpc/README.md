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
  TWAS jobs ran in the image (`Activating singularity image …/46cbe05b….simg`), the cached scripts under `$HOME` being bound by Snakemake by default.
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

## Task 20: user-built MAJIQ image, `majiq_container` (R27, module-v0.1.3)
- Scratch `$H/pantry-runs/t20`; `t20/repo` = the working tree by `tar` (code files checked identical to the commits, by sha1).
- Spike 1 (job 8403888, cpu-20, apptainer 1.5.2-1.el9; `t20/spike1/mini.def` from the pinned micromamba base): `apptainer build plain.sif
  mini.def` PLAIN_EXIT=0, `apptainer build --fakeroot` FAKEROOT_EXIT=0. Both log "User not listed in /etc/subuid, trying root-mapped
  namespace" and run `%post` under the fakeroot command, so the plain form is enough.
- Build (job 8403896, `t20/build/build.sh`, `sjob -c 4 -m 7500M -t 0-02:00:00`): `cd containers/majiq && apptainer build majiq.sif majiq.def`
  -> `t20/majiq.sif`, BUILD_EXIT=0, 4:59 wall, 3.6 GB max RSS; `%test` passed; in-build `majiq --version: majiq 3.0.22`, `import rna_majiq,
  rna_voila: ok`. `apptainer exec`: `/opt/conda/bin/majiq`; `/opt/majiq` holds NOTICE (Penn copyright, creators) and the clone's 10 licence
  files (boost, bxzstr, moccasin, eigen); labels `majiq.copyright`, `majiq.creators`, `majiq.licence`, `majiq.commit`. The def's header
  comments were reworded after this build (no other change). Image env vs `containers/majiq/majiq3.lock.txt`: conda packages identical
  except the build-only kernel-headers 6.12.0 (lock 5.14.0) and sysroot 2.39 (lock 2.34); pip numpy 2.0.1, zarr 2.18.3, xarray 2024.6.0,
  dask 2024.7.1, pandas 2.2.2, scipy 1.14.0, numcodecs 0.13.1, h5py 3.11.0 as the lock (`t20/build/image_*.txt`, `lock_vs_image.diff`).
- Docker route: not built. Local Docker (colima) cannot pull from Docker Hub here (TLS interception, x509 unknown authority); the HPC has no
  Docker or podman.
- Spike 2 (job 8403890 dry runs `t20/spike2/spike2.sh`; job 8403893 `rules.py`, Snakemake API): `conda: None if MAJIQ_CONTAINER else MAJIQ_ENV`
  and `container: MAJIQ_CONTAINER` parse; dry runs EXIT=0 standalone (gate-ir config) and local-path host (host config + intron_retention),
  each with `majiq_container` unset / absolute `.sif` / relative `majiq.sif`. Parsed MAJIQ rules: unset `conda='majiq3' container=None`; set
  `conda=None container=<path as given>`, also as `pantry_pheno_majiq_*` in the host. Local image paths are used as given (relative = against
  the job's working directory).
- Unset route unchanged (job 8403902, `t20/before/dry.sh`, `-n -p`): before = `$H/Pantry` f9c2345 (code = ccc56d6), after = `t20/repo`.
  Standalone 61 / 61 jobs, host 481 / 481; job sets identical (rule, wildcards, inputs, outputs, resources, shell commands, after
  normalising repo paths; `t20/before/jobsets.py`, run from that dir).
- Licence file in the container: job 8403903 (`t20/lic/lic.sh`, Snakemake's form `apptainer exec --home <run dir>` + the profile's
  `--bind=/data/hps/assoc/private/gdml --nv`): the source cache `$HOME/.cache/snakemake/...` is **not** visible (`ls`: No such file;
  `mount home = yes`, but `--home` replaces it). Job 8403905 (`lic2.sh`, `majiq.sif`): with `--bind=/data/hps/home/twil15/.cache/snakemake`
  the bundled `.lic` is visible, md5 99b38b20... as on the host; PATH in the container `/opt/conda/bin:/usr/local/sbin:...` (no host PATH).
  So `tests/profile` and `tests/host/profile` now also bind the source cache (d95e2d4).
  **Correction (Task 20 review, I1):** this test left out Snakemake's own `--bind <source cache>` (Snakemake >= 9.14.0, `runtime_paths`), so the cache is
  visible in real jobs and the bind is redundant on 9.14 and later (all runs here used 9.27.1). It is needed only on 9.12-9.13, for every containerised
  rule (tensorQTL, TWAS, MAJIQ), since `script()` paths point into the source cache for local-path and `github()` modules alike.
- Acceptance (job 8403914, `t20/accept.sh`, env pantry-baseline-pheno, `tests/profile` with `--sdm conda apptainer`,
  `--slurm-keep-successful-logs`; the first attempt 8403906 was cancelled because the slurm plugin deletes successful job logs).
  Run dirs `$H/pantry-runs/gate-ir-sif` and `gate-ir-symlink-sif` = `cp -a` of Task 17's `gate-ir` / `gate-ir-symlink` with the
  intron_retention intermediates and outputs removed and `majiq_container: $H/pantry-runs/t20/majiq.sif` appended to `config.yml`.
  No `majiq` on the controller's PATH.
  - gate-ir-sif: DRY 9 jobs; RUN_EXIT=0, 9 of 9 steps; MAJIQ jobs 8403915 (build), 8403917 (psi-coverage), 8403918 (psi): each log has
    `Activating singularity image /data/hps/assoc/private/gdml/user/twil15/pantry-runs/t20/majiq.sif`, 0 lines matching conda/majiq3.
    `cmp` vs `archive/phenotyping_ir`: SAME `intron_retention.bed.gz`, `.bed.gz.tbi`, `.phenotype_groups.txt`, `unnorm/intron_retention.bed`.
    `.sj` = `build/HG00108.sj` etc. Rerun dry run: nothing to be done.
  - gate-ir-symlink-sif: DRY 9 jobs; RUN_EXIT=0, 9 of 9; MAJIQ jobs 8403925 / 8403927 / 8403928, same image line, 0 conda; `.sj` =
    `build/HG00108.renamed.sj` etc. (symlink targets); `cmp` SAME for all 4 files. Rerun: nothing to be done.
- `github()` host dry run at the tag (job 8403936, `t20/gh_dry.sh`, fork entries cleared from the source cache): tag `module-v0.1.3` =
  da54ff9 locally and on the remote; host config + intron_retention, `majiq_container` unset: GH_DRY_unset_EXIT=0, 481 jobs; set: GH_DRY_set_EXIT=0,
  481 jobs (5 MAJIQ-related jobs each); no UnicodeEncodeError.
- `pytest tests/unit -q`: 151 passed (local).

## Task 21: pins, Snakemake floor 9.14, image-pin check, legacy env headers (R28, module-v0.1.4)
- Scratch `$H/pantry-runs/t21`; `t21/repo` = the working tree by `tar`. Explicit-list check for each env: every package basename in the
  built env's `conda list --explicit` must appear verbatim in the lock (`comm -23` on sorted basenames).
- txrevise (`t21/tx/tx.sh`: gate-2 `config.yml`, `--profile tests/profile`, `--conda-create-envs-only`, `conda list --explicit`, then the
  target `intermediate/reference/default/txrevise/completed.txt`). With only the yaml's 15 direct packages pinned, the existing env
  (`conda-envs/7710e88d…`) had 47 of 261 packages not in `baseline-pheno.lock.txt` (python 3.15, r-base build hc28aa72_21, …). Every package
  now pinned `name=version=build` from the lock: first build (job 8404058, cancelled) 271 of 272 in the lock (`packaging` 26.3); with it
  pinned, job 8404170 (env `8abe62ff…`): 272 of 272 in the lock, RUN_EXIT=0. `t21/tx/cmp.sh` (8404628), `cmp2.sh` (8404697), `rds2.R`
  (8404733) vs `$H/pantry-baseline/phenotyping/intermediate/reference/txrevise`: `gene_annotations.gtf`, `transcript_tags.txt` and the
  4 merged `txrevise.grp_{1,2}.{upstream,downstream}.gff3` byte-identical; the 600 `batch/*.gff3`: 225 byte-identical, the other 375
  identical except the `##date 2026-11-01`/`2026-11-04` header line; `txrevise_annotations.rds` differs in bytes (as gate 2's does, same env
  as the baseline): `transcript_metadata` identical, `exons`/`cdss` identical once `metadata()$genomeInfo` `Data source` (the GTF path,
  `reference/txrevise/` vs `reference/default/txrevise/`) and `Creation time` are dropped.
- LaDDR (`t21/laddr/laddr.sh`, `mamba env create -p`). First build (8404059), LaDDR v1.0.0 with the yaml's version pins: CREATE_EXIT=0 but
  pip had replaced conda's scipy 1.17.1 with 1.16.3 (LaDDR requires scipy<1.17) and `laddr --help` failed (scikit-fda 0.10.1 with
  multimethod 2.1: `TypeError: metaclass conflict`). Final yaml (conda stack pinned `name=version=build` from the lock, scipy
  1.16.3=py312h54fa4ab_2, `multimethod==1.12`), job 8404179: 105 conda packages, all in the lock except scipy; pip scipy 1.16.3 (no
  uninstall), laddr 1.0.0, multimethod 1.12; `laddr --help` and `import laddr, skfda, pyBigWig, bx` OK. LaDDR v1.0.0 = 2addf6fe; the only
  later commit, fae1b0d (2026-04-29, what an unpinned install at 8da7787 got), adds a stranded-bigWig option Pantry does not use.
- MAJIQ image rebuild (job 8404060, `t21/build/build.sh` as Task 20's; `majiq3.yml` pinned to `majiq3.lock.txt`, pip by `constraints.txt`):
  BUILD_EXIT=0, 4:59, 3.6 GB; `majiq --version: majiq 3.0.22`. Image vs lock: conda 64/64 identical (incl. kernel-headers 5.14.0 and
  sysroot 2.34, which drifted in Task 20), pip list identical to the lock's pip lines. IR check with this image and the profile without the
  source-cache bind (job 8404171, `t21/accept.sh` = Task 20's on `$H/pantry-runs/gate-ir{,-symlink}-sif21`): 9/9 each, all 4 IR files `cmp`
  SAME as upstream, every MAJIQ job log `Activating singularity image …/t21/majiq.sif`, 0 conda, `.sj` after the symlink targets.
- Floor: `gh api repos/snakemake/snakemake/compare/f1e8b629...<tag>`: v9.13.0-v9.13.7 behind (v9.13.7 by 10), v9.14.0 ahead 5, behind 0;
  v9.14.0 published 2025-11-27 (the commit is dated 2025-11-27 15:16). v9.12.0/v9.13.7 `singularity.py` already bind
  `<user_cache_dir>/snakemake/source-cache`, which is `workflow.source_cache_path` under the default shared-FS settings.
- Dry runs before (`$H/Pantry` at da54ff9, code = 73eb0e6) and after (`t21/repo`), job 8404177 (`t21/dry/dry.sh`) + 8404187
  (`jobsets.py`, as Task 20's): gate-2 phenotyping 345/345, `dryrun_majiq_laddr.yml` 166/166, `pheast_baseline.yml` 68/68, local-path
  host 481/481; identical job sets in all four.
- Image pins: `containers/pins/{tensorqtl,fusion}`, seeded from the Actions logs (`gh run view <id> --log`, "exporting manifest"):
  tensorqtl run 41716028853 at 206c8920 pushed sha256:65a42f9f… (= `qtl.smk`); run 41737670211 at d05c168f rebuilt the same tree (34f2a4f3)
  and pushed sha256:111e79ff…; fusion run 41755433329 at 11d208d4 pushed sha256:87335916… (= `twas.smk`). Trees at those commits equal
  HEAD's. Workflow `image-pins` run 41908908129 (push of 01b8623): PINS_OK, both manifests HTTP 200. Local failure demos: a pins copy with
  a wrong tree hash and a wrong digest (`check_pins.sh --offline /tmp/pins_bad`), and a scratch worktree commit editing
  `containers/tensorqtl/Dockerfile`: PINS_FAIL, exit 1.
- Full `github()` host run, `$H/pantry-runs/host-github-v014` (config as `host-github-v012` plus `majiq_container:
  $H/pantry-runs/t21/majiq.sif` under `pantry.phenotyping`; `tests/host/profile` without the source-cache bind; `Snakefile.github` at
  1e5c858 fetching commit 01b8623; fork entries cleared from the source cache first). Job 8404292 (`t21/gh_run.sh`, env
  pantry-baseline-pheno, Snakemake 9.27.1, `--slurm-keep-successful-logs`): DRY_ALL_EXIT=0 (481), DRY_TESTED_EXIT=0 (422), RUN_EXIT=0,
  548/548. MAJIQ jobs 8404709 (build), 8404722 (psi-coverage), 8404730 (psi): `Activating singularity image
  /data/hps/assoc/private/gdml/user/twil15/pantry-runs/t21/majiq.sif`, 0 conda lines; all 16 tensorQTL and 45 TWAS-image job logs show
  their image activation. The run built its txrevise env from the fetched 01b8623 yaml (`conda-envs/c60512a8…`, a different hash from the
  local-path `8abe62ff…`: source-cache bug 2); its explicit list is the same 272 packages, all in the lock.
- Compare 8404824 (`host-github-v014/compare.sh`, env pantry-baseline-pheast): PHENO_COMPARE_EXIT=0 (alt_polyA 0.222/0.328, expression
  0.366/0.732, isoforms 0.0835/0.167, unnorm 4.6e-4/2.8e-3, 338/676, 3.0e-3/6.1e-3), PHENO_CHECK_EXIT=0; IR_COMPARE_EXIT=0, IR_BYTES_SAME x4;
  XM_BED_SAME, XM_GROUPS_SAME (771 rows); PHEAST_COMPARE_EXIT=1 with 78 mismatches, all in the 16 `qtl/*.cis_qtl`/`*.cis_independent_qtl`
  files (non-blocking, R22); heritability x8, TWAS x8, genotype_qc no mismatch; PHEAST_CHECK_EXIT=0; covar inputs 21 SAME, `geno_pruned.log`
  DIFF. Job 8404825 (`cmds.sh`, env pantry-baseline-pheno): CMDS_EXIT=0, 17/17, CMDS_IDENTICAL. HOSTDIFF vs the Task 14 local-path run:
  `alt_polyA.bed.gz(.tbi)`, `unnorm/alt_polyA.bed`, `cross_modality.bed.gz(.tbi)` (kallisto noise, as in Task 19).
- Tag `module-v0.1.4` (annotated) = f72dc9d (stage code = 01b8623, which the run above fetched; then docs and `Snakefile.github` re-pointed;
  the tag was held until the run passed). `github()` dry run of `all` at the tag (job 8404833, `t21/gh_dry.sh`, fork entries cleared from the
  source cache, `host-github-v014` config): fetched f72dc9d (= remote), GH_DRY_TAG_EXIT=0, 481 jobs, no UnicodeEncodeError.
- `pytest tests/unit -q`: 151 passed (local).

## Task 22: LaDDR at run time against upstream 8da7787 (R29, module-v0.1.5)
- Scratch `$H/pantry-runs/t22` (scripts `clone.sh`, `menvs.sh`, `sklearn.sh`, `scipy.sh`, `imports.sh`, `run.sh`, `subst.sh`, `compare.sh`,
  `compare_seed.sh`, `cmds.sh` + `cmds_norm.py`, `cmp_laddr_interm.py`, `latdiff.py`, `normcheck.sh`, `seed/sitecustomize.py`, `gh_run.sh`,
  `compare_m2.sh`, `subst_m2.sh`). Config `tests/configs/phenotyping_baseline_laddr.yml` = `phenotyping_baseline_ir.yml` with
  `laddr.enabled: true` (upstream's default LaDDR settings, `dataset_name: pantry`, all 8 `regress_modalities`); nothing adjusted for the
  small data. The module reads it through the old-format shim (`project_dir` ignored with the note).
- Upstream env: upstream's `environment.yml` has no LaDDR, so `pantry-baseline-pheno-laddr` = `conda create --clone pantry-baseline-pheno
  --offline` (the original is untouched) + `pip install` of LaDDR v1.0.0 (`2addf6fe`) and `multimethod==1.12`, with the module LaDDR env's
  pip versions as constraints for packages the clone lacked (job 8405129). Its `conda list --explicit` equals the lock; pip added only the
  LaDDR stack (bx-python 0.16.0, scikit-fda 0.10.1, statsmodels 0.15.0, tqdm 4.70.1, numba 0.68.1, …) and replaced scipy 1.17.1 by a 1.16.3 wheel.
- Module envs (`--conda-create-envs-only` on the m1 config): first build (job 8405131) `deeptools.yaml` (only `deeptools=3.5.5` pinned)
  resolved python 3.13, numpy 2.5.4, pyBigWig 0.3.25, pysam 0.24.0, scipy 1.18.1; `laddr.yaml` left pyBigWig 0.3.29, matplotlib 3.11.2,
  pillow 13.0.0, wrapt 2.5.2 to pip although the lock has 0.3.22, 3.10.1, 10.3.0, 1.17.3. Fixed in 35ebcd7 (deeptools: its dependency
  closure in `pantry-baseline-pheno`'s `conda-meta`, `t22/closure.py`, plus pip/setuptools/wheel, all `name=version=build` from the lock;
  laddr: those four and their closure from the lock). Job 8405161: deeptools 123/123 conda packages in the lock; laddr all in the lock but scipy.
- Import-level check (`imports.sh`, jobs 8405430, 8406305): the top-level modules actually loaded by `laddr.{cli,binning,coverage,models,setup,init}`
  and by `deeptools.bamCoverage`, with versions, in the clone and in the module envs: bamCoverage identical (13 distributions: deeptools 3.5.5,
  pyBigWig 0.3.22, pysam 0.22.1, numpy 1.26.4, matplotlib 3.10.1, …); LaDDR identical (39 distributions in the module env: laddr 1.0.0, numpy 1.26.4, scipy 1.16.3,
  pandas 2.3.3, scikit-learn, scikit-fda 0.10.1, statsmodels 0.15.0, pyBigWig 0.3.22, bx-python 0.16.0, polars 0.20.31, …) except `psutil`
  7.2.2, which the clone has (from the lock) and imports opportunistically, and the module env lacks.
- First runs (fresh dirs `t22/{u1,u2}` = `git archive 8da7787 phenotyping` + config + `input` link, env `pantry-baseline-pheno-laddr`, upstream's
  own Snakefile; `t22/m1` = config + `input` link, `$H/Pantry` at 35ebcd7, env pantry-baseline-pheno; all `run.sh` = `snakemake --profile
  tests/profile --slurm-keep-successful-logs`, `sjob -c 2 -m 7500M -t 0-12:00:00`): U1 8405163, U2 8405164, M1 8405700 all failed in
  `laddr_coverage` only (U 370/378, M1 372/380 steps; keep-going finished everything else): `narwhals.exceptions.DuplicateError: Expected
  unique column names` from `sklearn.decomposition.PCA.fit_transform` in LaDDR's `remove_pcs_from_phenos` (the phenotype frame is a
  `pd.concat` of per-file `reset_index()` frames, so its transpose has duplicate integer column names, which scikit-learn 1.9.1's feature-name
  check rejects). **Upstream 8da7787 + LaDDR v1.0.0 cannot run in upstream's own env.** LaDDR's own `requirements.txt` pins scikit-learn
  1.5.2, so both sides now use conda-forge `scikit-learn=1.5.2=py312h7a48858_1` (54423dc for `laddr.yaml`; job 8406302 for the clone, `--no-deps`)
  and, to make the binaries equal too, conda-forge `scipy=1.16.3=py312h54fa4ab_2` in the clone instead of the pip wheel (job 8406307; the scipy
  trees and `_flapack` are byte-identical to the module env's). The clone's explicit list now differs from the lock in exactly scikit-learn and
  scipy, as the module's `laddr.yaml` env (job 8406304: 135 conda packages, the rest in the lock; `laddr --help` OK).
- Resumed U1 8406309, U2 8406310, M1 8406311 (M1 now at 54423dc, env `conda-envs/640155e0…`) with `--forcerun laddr_setup` (10/10 each,
  about 8 min): every LaDDR step from `laddr_setup` on ran under the final envs; bigWigs and all non-LaDDR outputs come from the first attempt.
- Compare 8406340 (`compare.sh m1 … u2`), 8406341 (`cmds.sh`):
  - bigWigs (bamCoverage, 16), `info/` (5) and `gene_bins/` (3) byte-identical in U1, U2 and M1; `coverage_manifest.tsv` identical.
  - `covg_norm/pantry/batch_*.npy` (3 of 7 files; the `bins.tsv.gz` and `samples.txt` are identical), `models/` (3) and the latent table differ U1~U2
    and U1~M1. The final BEDs differ too: U1 624 latent phenotypes, U2 625 (extra `ENSG00000197530__PC8`), M1 624 (same IDs and groups file as U1);
    max |diff| on common rows 3.13 (normalised) and 41.9 (unnorm) U1~U2, 3.13 and 39.4 U1~M1. `compare_outputs.py U1 M1 --noise U1 U2`: exit 1,
    "noise pair not numerically comparable: row count differs (624 vs 625)" for both BEDs. The U1/U2 envelope is not usable (a 16-sample
    rank-normal spans about ±2 so 3.13 is the whole range) and the tool is right to refuse it; it was not extended.
  - Command lines (`-n -p --forcerun bam_to_bigwig … pheno_groups_latent`, paths normalised): 28 upstream vs 29 module lines; the differences are
    the module's `mv` of the latent table to `latent_phenos.tsv.gz`, the `assemble_latent_bed` input path that follows from it, and
    `normalize_latent --samples intermediate/samples.txt` (`cmp`-identical to `input/samples.txt`; as every module normalize rule). With those
    three design mappings (`cmds_norm.py --design`): CMDS_IDENTICAL_AFTER_DESIGN (26 commands: 16 bamCoverage, `write_laddr_config.py`, `laddr
    setup/binning/coverage/fit/transform`, assemble, normalize, bgzip, groups; plus two `output:`/`reason:` lines the filter also keeps).
  - LaDDR config YAML written by `write_laddr_config.py`: identical after replacing the run dirs and `output/default/` → `output/` (CONFIG_SAME_NORM).
- Determinism with identical inputs: U1b = `cp -a` of U1, LaDDR rerun from `laddr_setup` (job 8406344): the same three `covg_norm` batches, models,
  latent table and final BEDs differ from U1, again 624 vs 625 rows. Cause: `remove_pcs_from_phenos` calls `PCA(n_components=4)` without
  `random_state` on a 16 x thousands matrix, for which scikit-learn's `auto` solver is `randomized` (sklearn `_pca.py`: `1 <= n_components <
  0.8 * min(X.shape)`), seeded from numpy's global RNG, which LaDDR never seeds. LaDDR v1.0.0 is nondeterministic run to run (upstream
  nondeterminism, not the refactor); the substitution run M1s (`cp -a` of M1 with U1's 8 regressed BEDs + `.tbi`, job 8406345) differs from U1 the same way.
- Seeded runs (R29 addendum): `seed/sitecustomize.py` seeds numpy's global RNG (`np.random.seed(0)`) in every Python process when `T22_SEED` is set
  (and logs it in the `laddr` processes). The seed applies to every Python process of the rerun jobs, on both sides alike; `PYTHONPATH=t22/seed T22_SEED=0`, all rerun with `--forcerun laddr_config`: U1b 8406372, U1c (another
  `cp -a` of U1) 8406373, M1s 8406432. Snakemake drops `PYTHONPATH` for conda jobs, so M1s needed `--conda-not-block-search-path-envvars` (the
  first seeded M1s, 8406374, had no seed line in its `laddr_coverage` log; discarded). Compare 8406545 (`compare_seed.sh m1s`):
  - U1b vs U1c: every LaDDR file byte-identical (bigWigs 16, info 5, gene_bins 3, covg_norm 7, models 3, latent table, config, manifest) and the 4 finals
    (SEED_u1b_u1c_EXIT=0 without `--noise`). Seeded LaDDR is deterministic.
  - **M1s vs U1b: byte-identical**: the 4 finals (`latent_residual.bed.gz`, `.tbi`, `.phenotype_groups.txt`, `unnorm/latent_residual.bed`, 624
    phenotypes), `covg_norm` 7/7, `models` 3/3, the latent table (module `latent_phenos.tsv.gz` vs upstream `phenotypes/latent_phenos.pantry.tsv.gz`),
    `info`, `gene_bins`, bigWigs; `compare_outputs.py u1b m1s --noise u1b u1c`: exit 0, MATCH; config CONFIG_SAME_NORM. Seed lines present in all
    three `laddr_coverage` logs.
- Non-LaDDR sanity (`compare.sh`, U1, U2, M1 vs `archive/phenotyping` with `--noise archive/phenotyping archive/phenotyping_run2`; IR vs
  `archive/phenotyping_ir`): exit 1 for all three, only in kallisto-derived files. U1: `expression.bed.gz` 0.744 > 0.732 and `isoforms.bed.gz`
  0.213 > 0.167 (their `unnorm/` files within the envelope: 338/676, 3.0e-3/6.1e-3); U2 and M1: `alt_TSS.bed.gz` and `unnorm/alt_TSS.bed`
  differ although the archive pair had them identical (noise 0). RNA_editing, splicing, stability and the 4 intron_retention files are
  byte-identical to the archive in all three. The single-pair envelope broke here (the concern of Task 13); the three runs were concurrent.
  Env attribution (job 8406357, `normcheck.sh`): `normalize_phenotypes.py` on U1's 9 `unnorm/*.bed` (8 modalities + latent_residual) in the clone
  (scipy 1.16.3, scikit-learn 1.5.2), in pantry-baseline-pheno (1.17.1, 1.9.1) and in the module's `python.yaml` env: all 27 outputs
  byte-identical, and equal to U1's `output/*.bed.gz` (KNN imputation ran on 5 of them). So the clone changes no normalisation, and
  `normalize_latent` (clone upstream, `python.yaml` in the module) is env-independent; `unnorm/latent_residual.bed` has no missing values.
- M2, full `github()` host run, `$H/pantry-runs/host-github-v015` (`t22/mk_m2.sh`: `host-github-v014`'s config, i.e. phenotyping `full` +
  intron_retention from `t21/majiq.sif` and pheast `geuvadis`, plus a `laddr:` block for `full` with upstream's defaults written out,
  `dataset_name: pantry` and all 8 `regress_modalities`; committed as `tests/host/config_laddr.yaml`, which also drives the dry run at the tag below;
  `majiq_container` there is the HPC path of Task 21's `.sif`); `Snakefile.github` at 7baa456 fetching commit 54423dc (the code M1 ran); fork entries
  cleared from the source cache first. Job 8406551 (`t22/gh_run.sh`, env pantry-baseline-pheno, `sjob -c 2 -m 7500M -t 0-12:00:00`):
  DRY_ALL_EXIT=0 (508 = Task 21's 481 + 27 LaDDR jobs), DRY_TESTED_EXIT=0 (449), RUN_EXIT=0, 575/575.
- Compare 8407162 (`t22/compare_m2.sh`, env pantry-baseline-pheno-laddr; Task 21's compare plus LaDDR):
  - phenotyping `full` (view without cross_modality, intron_retention, latent_residual) vs `archive/phenotyping` `--noise`: PHENO_COMPARE_EXIT=1, one
    mismatch, `unnorm/alt_TSS.bed` 1.2e-3 > 0 (the archive pair has it identical; `alt_TSS.bed.gz` itself identical; kallisto, as U2/M1 above);
    others within the envelope (alt_polyA 0.222/0.328, expression 0.366/0.732, isoforms 0.108/0.167, unnorm 5.8e-4/2.8e-3, 338/676, 3.0e-3/6.1e-3).
    PHENO_CHECK_EXIT=0 (29/29 targets incl. `latent_residual.bed.gz` 624 rows, `.tbi`, groups). IR_COMPARE_EXIT=0, IR_BYTES_SAME x4. XM_BED_SAME,
    XM_GROUPS_SAME (771 rows).
  - pheast `geuvadis`: PHEAST_COMPARE_EXIT=1, 89 mismatches, all in `*.cis_qtl`/`*.cis_independent_qtl` (non-blocking, R22); 0 others;
    PHEAST_CHECK_EXIT=0 (33/33); covar inputs 21 SAME, `geno_pruned.log` DIFF (as before).
  - LaDDR envs built from the fetched yamls (`conda-envs/ca351d83…` laddr, `85691e07…` deeptools; new hashes, source-cache bug 2): `conda list
    --explicit` and `pip list` equal to the local-path envs (job 8406304) for both.
  - LaDDR raw vs U1 (informational, kallisto noise in the inputs): 624 rows, same IDs, max |diff| 3.13 / 37.9; bigWigs, info, gene_bins, manifest
    identical; config CONFIG_SAME_NORM. Command lines (`-n -p --forcerun pantry_pheno_<rule>`, host paths mapped to upstream's): the same three
    design differences as M1; CMDS_IDENTICAL_AFTER_DESIGN.
- M2 substitution (`t22/subst_m2.sh`, job 8407165): M2's own LaDDR outputs, LaDDR dir and the 8 regressed BEDs kept in
  `host-github-v015/snap_unseeded/`; U1's 8 BEDs (+ `.tbi`) copied into `results/pantry/phenotypes/full/`, then `tested` with
  `PYTHONPATH=t22/seed T22_SEED=0 --conda-not-block-search-path-envvars --forcerun pantry_pheno_laddr_config`: RUN_EXIT=0, 15/15 (the LaDDR
  chain, `combine_modalities` and the groups files of the substituted BEDs; nothing in pheast). Compare 8407187 (`compare_seed.sh m2s`): seed
  line in the `laddr_coverage` log; **M2s vs U1b byte-identical**: the 4 finals (`compare_outputs.py --noise u1b u1c` exit 0, MATCH), `covg_norm`
  7/7, models 3/3, latent table, info, gene_bins, bigWigs, manifest; config CONFIG_SAME_NORM.
- Tag `module-v0.1.5` (annotated) = a8720eb (stage code = 54423dc, which M1 and the M2 run fetched: `git diff 54423dc a8720eb -- phenotyping
  pheast` is empty; then docs, `tests/configs/dryrun_majiq_laddr.yml` header and `Snakefile.github` re-pointed; the tag was held until M2
  passed). `github()` dry run of `all` at the tag (job 8407193, `t22/gh_dry.sh`, fork entries cleared from the source cache, `host-github-v015`
  config = `tests/host/config_laddr.yaml`): fetched a8720eb (= remote), GH_DRY_TAG_EXIT=0, 508 jobs, no UnicodeEncodeError.
- `pytest tests/unit -q`: 151 passed (local).
- Review fixes (after the tag; docs, tests and the compare tool only, tag unmoved):
  - pheast command lines of M2 (`t22/cmds_pheast_m2.sh` = `host-github-v014/cmds.sh` for v015, job 8407214, env pantry-baseline-pheno; `Snakefile.github`
    now fetches the tag, whose pheast equals 54423dc's): CMDS_EXIT=0, 17 vs 17 lines, CMDS_IDENTICAL. With PHENO_COMPARE_EXIT=1 (`unnorm/alt_TSS.bed`,
    1.2e-3 > 0) as the one exception, M2 met the whole Task 21 bar.
  - R12 amendment (`docs/module-status.md`): `compare_outputs.py --noise` takes two or more runs; envelope = 2× max |diff| over all pairs. Job 8407236
    (`t22/r12a.sh`, `t22/pairtable.py`; reports and `pairs_{rest,ir}.tsv` in `t22/cmp/r12a/`): noise `{archive/phenotyping, archive/phenotyping_run2,
    views/u1_rest, views/u2_rest}` and, for intron_retention, `{archive/phenotyping_ir, archive/phenotyping_ir_run2, views/u1_ir, views/u2_ir}`;
    candidates U1, U2, M1, M2 (view `views/m2_rest`/`m2_ir` = `host-github-v015/results/.../full` with the 8 normalised BEDs + `.tbi` from
    `snap_unseeded/orig_beds`, since job 8407165 had put U1's there) and `host-github-v014/cmp_pheno_full`/`cmp_ir` (Task 21). The old (a65ac3e) and new
    tool give identical two-run reports on M2 (TWO_RUN_REPORT_IDENTICAL), whose lines equal job 8407162's. `{u1,u2,m1,m2,v014}_REST_R12A_EXIT=0`,
    `_IR_R12A_EXIT=0` for all five.
