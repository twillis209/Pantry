# Pantry multi-dataset module refactor — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Also invoke `writing-snakemake-pipelines` before touching any `.smk`/Snakefile.

**Goal:** Make `phenotyping/` and `pheast/` importable once each as Snakemake modules that run many datasets (`{dataset}` wildcard) with host-set output roots and shared references, without changing any analysis result.

**Architecture:** Keep both stages and every rule body; replace parse-time globals with per-dataset lookups backed by a small pure-Python config library per stage (`<stage>/lib/config.py`, loaded via `workflow.source_path`, unit-tested); move paths onto `pantry_*` pathvars; turn parse-time file reads into rules/checkpoints; add cross-modality and FUSION-install rules; add per-rule conda envs. Verify against an upstream `8da7787` baseline on the HPC.

**Tech Stack:** Snakemake ≥ 9.12 (pathvars), Python 3 (pandas, pytest), JSON-schema config validation, SLURM via `sjob`, conda.

**Spec:** `docs/superpowers/specs/2026-10-30-module-refactor-design.md` (read it first; this plan does not repeat its rationale). Requirements source: `~/projects/scPantry/docs/handover/2026-10-30-pantry-fork-handover.md`.

## Global Constraints

- No change to analysis logic, default parameters or output formats. Rule bodies change only where a path/parameter must come from somewhere new.
- Never push to `upstream` (PejLab/Pantry). Push `module` and tags to `origin` (`twillis209/Pantry`) only.
- `min_version("9.12")` in both Snakefiles.
- Pathvars: phenotyping `pantry_output` (default `output`), `pantry_intermediate` (`intermediate`), `pantry_reference` (`intermediate/reference`); pheast `pantry_output` (`output`), `pantry_intermediate` (`intermediate`), `pantry_resources` (`intermediate/resources`).
- Layout: `<pantry_output>/{dataset}/…`, `<pantry_intermediate>/{dataset}/…`, `<pantry_reference>/{reference}/…` (STAR: `star_index_{read_length}`), `<pantry_resources>/fusion/…`.
- `wildcard_constraints`: `dataset`, `reference`, `sample_id`, `modality`, `pheno`, `type`, `group`, `position`, `chrom`, `batch`, `feature_type` = `[^/]+`; `read_length` = `\d+`.
- No `localrules:` blocks; use `localrule: True` per rule.
- Scripts referenced only through `workflow.source_path(...)`; nothing depends on the CWD except standalone `configfile: 'config.yml'` (which Snakemake skips inside modules).
- MAJIQ (`intron_retention`) and LaDDR: restructured, verified by dry run only; excluded from run-time test configs.
- Every Snakemake invocation on the HPC (including `-n`) runs via `sjob` with `--profile tests/profile`; non-interactive SSH starts with `source /etc/profile; module load slurm; source ~/.dotfiles/.bash_functions`. `sjob -m` is per CPU. Bulk transfer via `tar` over one `ssh`.
- HPC locations: `H=/data/hps/assoc/private/gdml/user/twil15`; `$H/pantry-baseline` (upstream `8da7787` checkout), `$H/pantry-testdata/{phenotyping,pheast}/input`, `$H/Pantry` (fork clone on `module`), runs under `$H/pantry-runs/<name>`.
- `snakefmt` on every new/changed Snakemake file; small single-purpose commits ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **Two datasets sharing one reference at different read lengths** — both STAR indexes build, no `AmbiguousRuleException`. Test: Task 7 Step "dry-run two-reference fixture" (`tests/configs/dryrun_multi.yml` has datasets with read_length 75 and 100 on one reference).
2. **Sample IDs/dataset names containing `.` or `_`** (`HG00108`, `onek1k_CD4_TCM`) — must not be split across wildcards (e.g. `{sample_id}.{feature_type}.counts.txt`). Test: Task 7 dry-run fixture uses dataset `ct_CD4.TCM` and checks the job list resolves.
3. **Old-format `config.yml` (no `datasets:` key) with `intermediate_dir` set** — runs, writes intermediates there. Test: Task 6 unit test `test_shim_maps_intermediate_dir` + Task 7 standalone dry-run.
4. **Host DAG where pheast's BEDs don't exist yet** — checkpoint-derived chroms/TWAS batches resolve after phenotyping finishes. Test: Task 14 chained host dry-run + Task 4 spike.
5. **Genotype `.fam` missing samples** — fails at job time with the upstream message `Samples missing from genotypes: [...]`, not at parse time. Test: Task 11 step with a fixture `.fam` lacking one sample.

---

## Phase A — baseline, tooling, spike

### Task 1: Test data, profile and upstream baseline on HPC

**Files:**
- Create: `tests/profile/config.yaml`, `tests/configs/phenotyping_baseline.yml`, `tests/configs/pheast_baseline.yml`, `tests/hpc/README.md` (commands used), `.gitignore` additions (`phenotyping/input`, `pheast/input` already ignored; add `tests/**/input`, `tests/**/results`, `tests/**/resources`).

**Interfaces:**
- Produces: `$H/pantry-baseline/archive/phenotyping/` (contents of upstream `phenotyping/output/`), `$H/pantry-baseline/archive/pheast_run1/`, `…/pheast_run2/` (two independent pheast runs), test data at `$H/pantry-testdata/…/input`, local copies in `phenotyping/input` and `pheast/input` (for local dry runs). Profile `tests/profile`.

- [ ] **Step 1:** Write `tests/profile/config.yaml` modelled on `~/projects/genex/smk/profile/config.yaml`: `executor: slurm`, `jobs: 100`, `default-resources: {slurm_account: gdml, slurm_partition: cpu-core, runtime: 60, mem_mb: threads * 7500}`, `rerun-triggers: mtime`, `keep-going: True`, `rerun-incomplete: True`, `nolock: True`, `max-status-checks-per-second: 1`. No `software-deployment-method` (baseline uses monolithic envs; Task 12 adds `--sdm conda` on the command line).
- [ ] **Step 2:** Write `tests/configs/phenotyping_baseline.yml` = upstream `phenotyping/config.yml` with the `intron_retention` group deleted and `laddr.enabled: false`. Write `tests/configs/pheast_baseline.yml` = upstream `pheast/config.yml`, with `modalities` reduced to those present in the downloaded `input/phenotypes/` (check after Step 3; record which).
- [ ] **Step 3:** Download both test tarballs (URLs in `README.md` lines 73, 87) on HPC into `$H/pantry-testdata/{phenotyping,pheast}/` and extract; also locally into `phenotyping/input`, `pheast/input`. Expected: `phenotyping/input/{fastq_map.txt,samples.txt,ref/}`, `pheast/input/{phenotypes/,GEUVADIS.445_samples.GRCh38.chr1_120Mb.*,samples.txt,ld_snps_chr1.txt}`.
- [ ] **Step 4:** On HPC: `git clone` the fork to `$H/pantry-baseline`, `git checkout 8da7787`; create conda envs `pantry-baseline-pheno` and `pantry-baseline-pheast` from the two upstream `environment.yml` files (via `sjob`, `-c 4 -t 0-02:00:00`). In each stage dir symlink `input -> $H/pantry-testdata/<stage>/input` and copy in the baseline config as `config.yml`.
- [ ] **Step 5:** Run phenotyping: `sjob -c 2 -m 7500M -t 0-12:00:00 -n base_pheno -e pantry-baseline-pheno snakemake --profile <fork>/tests/profile` in `$H/pantry-baseline/phenotyping`. Expected: exit 0; `output/` contains alt_TSS, alt_polyA, expression, isoforms, RNA_editing, splicing, stability BED/tbi/groups files. Copy `output/` → `archive/phenotyping/` (and `intermediate/reference/` file listing → `archive/phenotyping_reference_files.txt`).
- [ ] **Step 6:** Run pheast twice in two separate copies (`pheast` and `pheast_run2`, each a fresh copy of the stage dir), same command with `-e pantry-baseline-pheast`. Expected: `output/{qtl,heritability,twas}/…` per modality, `output/qtl/genotype_qc.tsv`. Archive to `archive/pheast_run1/`, `archive/pheast_run2/`. Record controller job IDs and `sacct -o JobID,ReqMem,MaxRSS,Elapsed` in `tests/hpc/README.md`.
- [ ] **Step 7:** Commit `tests/` and `.gitignore` (`git commit -m "Add HPC test profile and baseline configs"`).

### Task 2: Output comparison tool

**Files:**
- Create: `tests/compare_outputs.py`, `tests/unit/test_compare_outputs.py`

**Interfaces:**
- Produces: CLI `python tests/compare_outputs.py <expected_dir> <actual_dir> --report <md> [--noise <dirA> <dirB>]`, exit 0 iff all files match; `compare_tables(a: Path, b: Path, exact_cols, noisy_cols, noise: dict[str,float]) -> list[str]` (list of mismatch descriptions).

Rules the tool applies (relative paths matched between dirs; a file present on one side only is a mismatch, except `*.tbi`, which is compared via the decompressed BED it indexes):
- `*.bed.gz`, `*.txt`, `*.tsv`, `*.pos`, `*.list`: byte-identical after gzip decompression.
- tensorQTL `*.cis_qtl.txt.gz`, `*.cis_independent_qtl.txt.gz`: noisy columns `beta_shape1, beta_shape2, true_df, pval_true_df, pval_perm, pval_beta, qval, pval_nominal_threshold`; all other columns exact (floats compared with `rtol=1e-9`). Rows keyed by `phenotype_id` (+`rank` for independent). Allowed noisy |diff| per column = 2 × max |diff| seen between the `--noise` pair (0 if no noise pair given). Differing `variant_id` rows are reported as "tie-break?" if `pval_nominal` is equal on both sides, otherwise mismatch.
- `*.tar.bz2` (TWAS): extract both; `.list`, `.pos`, `.profile` byte-identical; `.profile.err` ignored (timings); each `*.wgt.RDat` compared with `Rscript -e 'a<-new.env();b<-new.env();load(A,a);load(B,b);stopifnot(isTRUE(all.equal(as.list(a),as.list(b),tolerance=1e-8)))'`.
- `heritability/*.hsq.tsv`: sorted rows, exact.

- [ ] **Step 1:** Write tests with tiny fixture files built in `tmp_path`: `test_identical_bed_passes`, `test_bed_value_change_fails`, `test_noisy_column_within_noise_passes` (noise pair diff 0.01, actual diff 0.015 → pass), `test_noisy_column_outside_noise_fails` (0.03 → fail), `test_exact_column_change_fails` (`slope`), `test_tie_break_reported_not_failed`, `test_missing_file_fails`.
- [ ] **Step 2:** `pytest tests/unit/test_compare_outputs.py -v` → FAIL (module missing).
- [ ] **Step 3:** Implement `tests/compare_outputs.py` (pandas; stdlib `tarfile`, `gzip`).
- [ ] **Step 4:** `pytest tests/unit/test_compare_outputs.py -v` → PASS.
- [ ] **Step 5:** On HPC: `python tests/compare_outputs.py archive/pheast_run1 archive/pheast_run2 --report archive/pheast_noise.md`. Expected: only noisy columns/tie-breaks differ (if any exact column differs between two upstream runs, stop and report — the comparison design must change). Commit tool, tests and `tests/hpc/README.md` note of the result.

### Task 3: Pure formatting commit

- [ ] **Step 1:** `snakefmt phenotyping/Snakefile phenotyping/steps pheast/Snakefile pheast/steps`.
- [ ] **Step 2:** Locally: `cd phenotyping && cp ../tests/configs/phenotyping_baseline.yml config.yml && snakemake -n --quiet rules > /tmp/after.txt`; same on the parent commit → `/tmp/before.txt`; `diff` → empty. Same for pheast. (Local dry runs are allowed: they're not cluster work.) Restore `config.yml` from git.
- [ ] **Step 3:** Commit `"Apply snakefmt (formatting only)"`.

### Task 4: Spike — checkpoint inside a prefixed module (throwaway, not committed)

- [ ] **Step 1:** In `/tmp/ckpt-spike/`: `mod/Snakefile` with `checkpoint summary` (`output: "<pantry_output>/{dataset}/n.txt"`, writes `3`) and `rule gather` whose input function calls `checkpoints.summary.get(dataset=w.dataset).output[0]` and expands `<pantry_output>/{dataset}/part{i}.txt` for i < n; `rule part` makes each part. `host/Snakefile`: `module m: snakefile: "../mod/Snakefile"; pathvars: pantry_output="res/x"`, `use rule * from m as m_*`, `rule all` requesting `res/x/d1/gathered.txt`.
- [ ] **Step 2:** `cd host && snakemake -c1`. Expected: completes; `res/x/d1/part{0,1,2}.txt` exist. Also confirm a module-level `localrule: True` rule survives renaming (`snakemake -n --verbose` shows it as local).
- [ ] **Step 3:** Record outcome in `tests/hpc/README.md` ("spike result"). If it fails, stop and switch Task 11 to the spec's fallback (pheast reads BEDs at DAG time; host runs pheast as a second invocation) and tell the user.

---

## Phase B — phenotyping

### Task 5: Phenotyping relocatable (still single-dataset)

**Files:** Modify `phenotyping/Snakefile`, `phenotyping/steps/*.smk`.

**Interfaces:**
- Produces: `script(name: str) -> str` in `phenotyping/steps/setup.smk` returning `workflow.source_path(f"../scripts/{name}")`; pathvars `pantry_output`, `pantry_intermediate`, `pantry_reference` with standalone defaults; helper `parent_dir(w, output) -> str` (directory of `output[0]`) for `mkdir -p` params.

- [ ] **Step 1:** `min_version("9.12")`; top-level `pathvars:` with the defaults in Global Constraints; global `wildcard_constraints:`; keep `configfile: 'config.yml'`.
- [ ] **Step 2:** Every `localrules:` block → `localrule: True` on those rules.
- [ ] **Step 3:** Every `python3 scripts/…`/`Rscript scripts/…` → `{params.script}` with `params: script=script("…")`.
- [ ] **Step 4:** Replace `output_dir`/`interm_dir`/`ref_dir` in input/output paths with `"<pantry_output>/…"`, `"<pantry_intermediate>/…"`, `"<pantry_reference>/…"` strings (no `{dataset}` yet). Params/shell uses of these dirs → derived from `input`/`output` (`parent_dir`) since pathvars don't reach params. `intermediate_dir`/`intermediate_ref_dir` config keys: if set, emit a conditional top-level `pathvars:` (same level as defaults, later wins).
- [ ] **Step 5:** Generic rules scoped: `index_bam` input `<pantry_intermediate>/bam/{sample_id}.bam`; `index_bed` `<pantry_output>/{pheno}.bed.gz`; `index_fasta` under `<pantry_reference>`.
- [ ] **Step 6:** Local check: standalone `snakemake -n` with baseline config → same rule counts as Task 3's `/tmp/before.txt`; run once from a different CWD with `snakemake -s <abs>/phenotyping/Snakefile --directory <tmp> --configfile <abs config>` → resolves.
- [ ] **Step 7 (gate 1):** HPC: pull branch in `$H/Pantry`, run phenotyping standalone in `$H/pantry-runs/gate1-pheno` (monolithic env `pantry-baseline-pheno`); `compare_outputs.py $H/pantry-baseline/archive/phenotyping output` → exit 0.
- [ ] **Step 8:** Commit in ≥3 commits: localrules; script paths; pathvars + constraints + scoped generic rules.

### Task 6: Phenotyping config library, shim and schema

**Files:**
- Create: `phenotyping/lib/config.py`, `phenotyping/schemas/config.schema.yaml`, `tests/unit/test_phenotyping_config.py`
- Modify: `phenotyping/steps/setup.smk` (load lib via `importlib.util.spec_from_file_location` on `workflow.source_path("../lib/config.py")`; delete the old load-time reads)

**Interfaces (all pure; file reads only where named):**
- `normalise(config: dict) -> dict` — if no `datasets` key: shim to `{"references": {"default": {genome, annotation, edit_sites_bed?, annotation_gff3?}}, "datasets": {"default": rest}}` (moves `ref_genome`→`genome`, `ref_anno`→`annotation`, `modality_groups.RNA_editing.edit_sites_bed`→reference, `modality_groups.intron_retention.annotation_gff3`→reference; renames `samples_file`→`samples`). Deep-merges `defaults` under each dataset. Returns new dict plus `legacy_pathvars: dict` (from `intermediate_dir`, `intermediate_ref_dir`). Raises `ValueError` on: unknown `reference`, `ref_anno`/`annotation` not ending `.gtf`, single-end without fragment lengths, invalid `laddr.regress_modalities`.
- `load_samples(spec: list | str) -> list[str]` — list returned as is; str read as one-column file (same pandas call as upstream).
- `load_fastq_map(spec: dict | str, fastq_dir: str | None, paired_end: bool | None) -> dict[str, tuple]` — same return shape and errors as upstream `load_fastq_map`; dict form `{sample: [[R1...],[R2...]]}` or `{sample: [[R...]]}`, relative paths joined to `fastq_dir`.
- `edit_sites_min_samples(configured: int, fraction: float, n: int) -> int` — upstream formula `min(max(configured, ceil(fraction*n)), n)`.
- `available_modalities(modality_groups: dict) -> list[str]` — unchanged.
- `dataset_targets(ds: dict) -> list[str]` — file names relative to `<pantry_output>/{dataset}/` (all `files` of all groups, LaDDR files if enabled, cross-modality files if `cross_modality` set).
- In `setup.smk`: `CONFIG = normalise(config)`; validate with `snakemake.utils.validate(CONFIG, "../schemas/config.schema.yaml")` (schema carries the upstream defaults for LaDDR keys, `edit_sites_min_samples_fraction: 0.1`, `edit_sites_cluster: true`, `edit_sites_correlation_threshold: 0.9`, `psicov_batch_size: 25`, `min_experiments: 0.01`, MAJIQ `threads: 8`; top-level `majiq_env: majiq3`); helpers `ds(w) -> dict`, `ref(w) -> dict` (by `ds(w)["reference"]` or `w.reference`), `samples(dataset: str) -> list[str]` and `fastq(dataset: str) -> dict` (both `functools.cache`d, calling the loaders lazily).

- [ ] **Step 1:** Tests: `test_shim_wraps_legacy_config` (upstream `phenotyping/config.yml` → one reference `default`, one dataset `default`, `edit_sites_bed` on the reference), `test_shim_maps_intermediate_dir`, `test_defaults_deep_merged`, `test_unknown_reference_raises`, `test_gtf_suffix_enforced`, `test_single_end_requires_fragment_length`, `test_inline_fastq_map_paired_and_single`, `test_file_fastq_map_matches_upstream` (parse `phenotyping/input/fastq_map.txt` if present, else a fixture; compare to upstream function copied into the test), `test_mixed_endedness_raises`, `test_edit_sites_min_samples` (`(10, 0.1, 16) → 10`, `(10, 0.1, 8) → 8`, `(10, 0.5, 30) → 15`), `test_dataset_targets_includes_laddr_and_cross_modality`.
- [ ] **Step 2:** `pytest tests/unit/test_phenotyping_config.py -v` → FAIL.
- [ ] **Step 3:** Implement `lib/config.py`, schema, and `setup.smk` wiring. The transcript-row GTF check is **not** here (Task 7 rule).
- [ ] **Step 4:** Tests PASS; local standalone `snakemake -n` with baseline config still resolves (rules unchanged yet; setup exposes the old global names temporarily from `CONFIG["datasets"]["default"]` so Task 7 can convert rules file by file).
- [ ] **Step 5:** Commit `"Add phenotyping config library, legacy shim and schema"`.

### Task 7: Phenotyping `{dataset}` / `{reference}` wildcards (core modalities)

**Files:** Modify `phenotyping/Snakefile`, `phenotyping/steps/{setup,align,expression,splicing,stability,alt_TSS_polyA,RNA_editing}.smk`. Create `tests/configs/dryrun_multi.yml`, `phenotyping/config.yml` (converted to new format; old format stays accepted).

**Interfaces:**
- Consumes: Task 6 helpers.
- Produces: rule paths per Global Constraints; new rules `write_samples` (`<pantry_intermediate>/{dataset}/samples.txt`, from `samples(dataset)`), `validate_annotation` (`<pantry_reference>/{reference}/annotation.ok`, upstream transcript-row check, same error text); `rule all` (`default_target: True`) over `dataset_targets` for every dataset.

- [ ] **Step 1:** Write `tests/configs/dryrun_multi.yml`: one reference (local test files), datasets `full` (16 paired, read_length 75), `ct_CD4.TCM` (inline 4-sample list, inline single-end fastq map, read_length 100, fragment 200/50), both with all core groups and RNA_editing.
- [ ] **Step 2:** Local `snakemake -n -s phenotyping/Snakefile --configfile tests/configs/dryrun_multi.yml --directory /tmp/dr` → currently FAILS (no dataset support). Keep the failing output.
- [ ] **Step 3:** Convert rules: references under `<pantry_reference>/{reference}/` (cDNA, kallisto.idx, `star_index_{read_length}`, txrevise tree, exonic/intronic GTF, `edit_sites_to_genes.tsv`), reference inputs (`genome`, `annotation`, `edit_sites_bed`) via `lookup`/input functions on `w.reference`; per-dataset rules under `{dataset}`, consuming references through input functions (`ref_path(w, "kallisto.idx")` → `<pantry_reference>/{ds reference}/kallisto.idx`); `expand(..., sample_id=samples)` → input functions over `samples(w.dataset)`; `fastq_inputs`, `fastq_star_param`, `kallisto_single_end_flag`, featureCounts `paired_flag` keyed by `(w.dataset, w.sample_id)`; RNA-editing `min_sam` param via `edit_sites_min_samples(...)` per dataset; scripts' `--samples` take `write_samples` output; reference-building rules take `annotation.ok` as input. `genomeSAindexNbases` computed from `input.ref_genome` at job time (params lambda with `input`).
- [ ] **Step 4:** Dry run from Step 2 → resolves; job list contains `star_index` twice (75, 100), one `kallisto_index`, one `txrevise_prepare_annotations`; no ambiguity warnings. Standalone dry run (old config) also resolves.
- [ ] **Step 5:** Commit per area (references/align; expression; splicing+stability; alt_TSS_polyA; RNA_editing; rule all + config.yml).
- [ ] **Step 6 (gate 2):** HPC standalone run with the *old-format* baseline config in `$H/pantry-runs/gate2-pheno`; compare `output/default/` against `archive/phenotyping` → exit 0.

### Task 8: MAJIQ and LaDDR static rules (dry-run only)

**Files:** Modify `phenotyping/steps/intron_retention.smk`, `phenotyping/steps/laddr.smk`, `phenotyping/Snakefile` (always include both).

- [ ] **Step 1:** MAJIQ: GFF3 input function picks reference `annotation_gff3` or `<pantry_reference>/{reference}/majiq/annotation.gff3`; `majiq_build` output `directory(<pantry_intermediate>/{dataset}/intron_retention/build)` plus `splicegraph.zarr` inside it; `majiq_psi_coverage` params compute `.sj` paths from `Path(bam).resolve().stem` at job time (params function); batches from `samples(w.dataset)` and `psicov_batch_size`; MAJIQ rules use `conda: MAJIQ_ENV` where `MAJIQ_ENV = CONFIG.get("majiq_env", "majiq3")` is a top-level setting (one MAJIQ env per module instance; documented).
- [ ] **Step 2:** LaDDR: project dir `<pantry_intermediate>/{dataset}/laddr`; `laddr_dataset_name(w)` defaults to `w.dataset`; existing vs built bins via distinct outputs `bins_ready.existing`/`bins_ready.built` and an input function in `laddr_coverage`; `laddr_transform` moves LaDDR's `latent_phenos.<name>.tsv.gz` to stable `<…>/laddr/latent_phenos.tsv.gz`.
- [ ] **Step 3:** Dry runs: `tests/configs/dryrun_multi.yml` + a variant `tests/configs/dryrun_majiq_laddr.yml` enabling intron_retention (with and without `annotation_gff3`) and LaDDR (with and without `use_existing_bins`) → all resolve.
- [ ] **Step 4:** Commit (two commits: MAJIQ, LaDDR).

### Task 9: Cross-modality rule

**Files:** Create `phenotyping/steps/cross_modality.smk`; modify `phenotyping/Snakefile`, `phenotyping/scripts/combine_modalities.sh` (header comment pointing at the rule; logic unchanged).

**Interfaces:** Produces `<pantry_output>/{dataset}/cross_modality.bed.gz`, `.bed.gz.tbi`, `.phenotype_groups.txt` for datasets with `cross_modality.modalities`.

- [ ] **Step 1:** Implement rule `combine_modalities` reproducing the script's commands exactly (header `diff`, `modality:` prefix awk, `sort -k1,1 -k2,2n`, bgzip, `tabix -p bed`, groups awk), with temp file under the dataset's intermediate dir.
- [ ] **Step 2:** Equivalence check on HPC: in a scratch copy of `archive/phenotyping` run `bash combine_modalities.sh <7 modalities>`; run the rule on gate-2 outputs (add `cross_modality` to that config); `compare_outputs.py` on the three files → exit 0.
- [ ] **Step 3:** Commit `"Add per-dataset cross-modality assembly rule"`.

---

## Phase C — pheast

### Task 10: Pheast relocatable, config library, shim, schema

**Files:**
- Create: `pheast/lib/config.py`, `pheast/schemas/config.schema.yaml`, `tests/unit/test_pheast_config.py`
- Modify: `pheast/Snakefile`, `pheast/steps/*.smk`

**Interfaces:**
- `normalise(config: dict) -> dict` — shim: no `datasets` → `{"datasets": {"default": rest}}`, `samples_file`→`samples`, `intermediate_dir` → legacy pathvar; deep-merge `defaults`; raises on `cross_modality` with non-`qtl` analyses or trans files (upstream message).
- `load_samples(spec)` as Task 6. `batch_ranges` loaded from `scripts/twas_batches.py` via `source_path` (no `from scripts…` import).
- Schema defaults: `qtl_min_allele_count: 10`, `twas_models_per_job: 200`, `twas_threads: 4`; `qtl_cis_window` (integer, **no schema default** — `run_tensorqtl.py`'s own default 1000000 applies), `fusion_dir`, `gcta` optional.
- `qtl_filter_args(ds: dict) -> str` — upstream string plus ` --window {qtl_cis_window}` **only if** the key was set explicitly (so default commands are byte-identical).

- [ ] **Step 1:** Tests: `test_shim_wraps_legacy_config`, `test_cross_modality_with_twas_raises`, `test_qtl_filter_args_default_matches_upstream` (`'--min_allele_count 10'`), `test_qtl_filter_args_with_maf_and_window`, `test_defaults_deep_merged`.
- [ ] **Step 2:** FAIL → implement → PASS.
- [ ] **Step 3:** Snakefile/steps: `min_version`, pathvars, constraints, `localrule: True`, `source_path` scripts (incl. `twas_assemble_summary`'s `realpath scripts`), `<pantry_*>` paths (no `{dataset}` yet), wired to `normalise`.
- [ ] **Step 4:** Local standalone dry run (baseline pheast config) → same rules as `/tmp/before` for pheast.
- [ ] **Step 5:** Commits: relocatable; config lib + schema.

### Task 11: Pheast `{dataset}`, job-time checks, checkpoint, FUSION install

**Files:** Modify `pheast/steps/{setup,covariates,qtl,heritability,twas}.smk`, `pheast/Snakefile`, `pheast/config.yml` (new format). Create `pheast/steps/fusion.smk`.

**Interfaces:**
- Rules: `check_genotype_samples` → `<pantry_intermediate>/{dataset}/genotype_samples.ok` (reads `.fam`, raises `Exception(f'Samples missing from genotypes: {missing}')`), input to `prune_for_covar` and `qtl_genotype_qc`; `checkpoint phenotype_summary` → `<pantry_intermediate>/{dataset}/summary/{modality}.tsv` (columns `chroms` comma-joined, `n_models`); `install_fusion` → `<pantry_resources>/fusion/{FUSION.compute_weights.R,utils/plink_utils.R,utils/FUSION.profile_wgt.R,gcta_nr_robust}` (`chmod +x` gcta) unless `fusion_dir`/`gcta` configured.
- Phenotype inputs: `pheno_bed(w)` → `Path(ds(w)["phenotype_dir"]) / f"{w.modality}.bed.gz"` (and `.tbi`, groups).

- [ ] **Step 1:** Fixture test (local, dry run + one local job): copy `pheast/input/GEUVADIS…fam` minus one sample to `/tmp/badfam/`, config pointing at it; `snakemake … <intermediate>/default/genotype_samples.ok -c1` → fails with `Samples missing from genotypes: ['<id>']`; parse succeeds.
- [ ] **Step 2:** Convert all pheast paths to `{dataset}`; `combine_heritability_chr` and `twas_batch_status_input` use `checkpoints.phenotype_summary.get(dataset=…, modality=…)`; `twas_pos_file` `n_samples` from `samples(w.dataset)`; FUSION paths from `install_fusion` outputs or config.
- [ ] **Step 3:** Step 1 fixture now passes as specified; standalone dry run resolves (heritability/TWAS jobs appear after checkpoint as expected).
- [ ] **Step 4 (gate 3):** HPC pheast standalone, old-format baseline config, `$H/pantry-runs/gate3-pheast`; `compare_outputs.py archive/pheast_run1 output/default --noise archive/pheast_run1 archive/pheast_run2` → exit 0.
- [ ] **Step 5:** Commits: dataset wildcard; genotype check; checkpoint; FUSION install.

---

## Phase D — environments, multi-dataset, host, finish

### Task 12: Per-rule conda envs

**Files:** Create `phenotyping/envs/*.yaml`, `pheast/envs/*.yaml` (split as in spec §7; carry over upstream constraints, e.g. `kallisto=0.48.0`, `regtools>=1.0.0`, `gemma=0.98.5`, `plink2=2`, `dask<2024.2`, `torch==2.3.1` via pip); add `conda:` to every `shell:`/`script:` rule (MAJIQ keeps named env). Keep both `environment.yml`.

- [ ] **Step 1:** HPC: `snakemake --sdm conda --conda-create-envs-only` for both stages via `sjob` (`-c 4 -t 0-04:00:00`) → all envs build.
- [ ] **Step 2:** Re-run gates 2 and 3 with `--sdm conda` in fresh run dirs → comparisons exit 0. If a tool version drift changes outputs, pin that tool to the version in the baseline env (`conda list -n pantry-baseline-*`) and rerun.
- [ ] **Step 3:** Commit per stage.

### Task 13: Multi-dataset acceptance (criterion 2)

**Files:** Create `tests/configs/phenotyping_multi.yml`, `tests/configs/pheast_multi.yml`, `tests/fixtures/make_multi_inputs.py`.

- [ ] **Step 1:** `make_multi_inputs.py` writes: phenotyping `half_se` sample list (first 8 of `samples.txt`) and single-end fastq map (R1 only); pheast `subset` (first 300 samples): BEDs subset to those columns (bgzip + tabix), groups files copied, samples list.
- [ ] **Step 2:** Phenotyping config: one reference; `full` (baseline settings) + `half_se` (single-end, fragment 200/50, same read_length). Pheast config: `full` (baseline) + `subset` (`qtl_min_allele_count: 20`).
- [ ] **Step 3:** HPC runs (`--sdm conda`) → exit 0; one STAR index and one kallisto index built; `full` matches baseline (phenotyping exact; pheast with noise pair); `half_se` and `subset` outputs exist and are non-empty for every target.
- [ ] **Step 4:** Commit.

### Task 14: Host workflow test (criterion 3)

**Files:** Create `tests/host/Snakefile`, `tests/host/config.yaml`.

- [ ] **Step 1:** Host imports both stages by relative path with the spec §1 snippet (pathvars under `results/`/`resources/`), datasets: phenotyping `full` with `cross_modality`; pheast `geuvadis` (baseline test data) and `chained` whose `phenotype_dir` is `results/pantry/phenotypes/full`.
- [ ] **Step 2:** HPC `snakemake -n` of the whole host (incl. `chained`) → resolves (chained pheast jobs pending on phenotyping + checkpoint).
- [ ] **Step 3:** HPC full run of targets for phenotyping `full` and pheast `geuvadis` → exit 0; outputs under `results/pantry/…` match baseline.
- [ ] **Step 4:** Push `module` to `origin`; tag `module-rc1`; host variant using `github("twillis209/Pantry", path=…, tag="module-rc1")` → dry run + `install_fusion` job succeed (gcta executable).
- [ ] **Step 5:** Commit host test.

### Task 15: Lint, docs, status note, release tag

**Files:** Create `CHANGES-module.md`, `docs/module-status.md`. Modify `README.md` (Snakemake ≥ 9.12, new config format, module usage).

- [ ] **Step 1:** `snakemake --lint` for both stages → fix or list each remaining warning with justification (expected: `run:` blocks) in `CHANGES-module.md`. `snakefmt --check phenotyping pheast tests/host` → clean.
- [ ] **Step 2:** `CHANGES-module.md`: every change, config migration (old → new key table), path layout change (`output/default/…`), deviations, MAJIQ/LaDDR dry-run-only status.
- [ ] **Step 3:** `docs/module-status.md` (handover "Report back"): branch/tag, exact `module` snippet + config shape, acceptance results with run dirs and comparison reports, deviations, follow-ups (MAJIQ/LaDDR runtime verification, env pinning).
- [ ] **Step 4:** Commit, push, tag `module-v0.1.0`.
