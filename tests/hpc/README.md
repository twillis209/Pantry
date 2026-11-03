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

## Task 12: per-rule conda envs, tensorQTL in a GPU container
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
