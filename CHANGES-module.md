# Pantry as a reusable, multi-dataset Snakemake module: changes

Branch `module` of `twillis209/Pantry`, based on upstream `PejLab/Pantry` `8da7787`. Release tag `module-v0.1.0`.
Design: `docs/superpowers/specs/2026-10-30-module-refactor-design.md`. Status and acceptance results: `docs/module-status.md`.
Test commands and HPC run records: `tests/hpc/README.md`.

**What it does.** `phenotyping/` and `pheast/` can now each be imported once as a Snakemake `module` and run any number of datasets in one DAG. The host sets the output roots with pathvars. Datasets that name the same reference share its reference files. Standalone use inside `phenotyping/` or `pheast/` still works, and the old single-dataset config format is still accepted.

**What it does not change.** Analysis logic, default parameters and output file formats are unchanged. Rule bodies changed only where a path or parameter now comes from somewhere else. `git diff 8da7787 -- pheast/scripts` is empty. In `phenotyping/scripts`, only a usage comment was added to `combine_modalities.sh`. The one exception to "same software" is tensorQTL's runtime: under apptainer it runs in a pinned container (see [Environments](#environments-containers-gpu)).

Requires **Snakemake ≥ 9.12** (pathvars). Both Snakefiles set `min_version("9.12")`. Tested with 9.26–9.27.1.

## Contents
- [Layout and pathvars](#layout-and-pathvars)
- [Config format and migration table](#config-format)
- [Validation](#validation)
- [References](#references)
- [Parse-time reads moved to rules and checkpoints](#parse-time-reads-moved-to-rules-and-checkpoints)
- [New and removed rules](#new-and-removed-rules)
- [Environments, containers, GPU](#environments-containers-gpu)
- [FUSION](#fusion)
- [Tests](#tests)
- [Lint and formatting](#lint-and-formatting)
- [Migrating a single-dataset config](#migrating-a-single-dataset-config)
- [Deviations from the handover contract and the spec](#deviations-from-the-handover-contract-and-the-spec)
- [Known limitations and follow-ups](#known-limitations-and-follow-ups)

## Layout and pathvars

- Every rule path is `<pathvar>/{dataset}/…` or `<pantry_reference>/{reference}/…`. The pathvar defaults reproduce the standalone layout, with one added `{dataset}` level:

  | Pathvar | Stage | Default | Contents |
  |---|---|---|---|
  | `pantry_output` | both | `output` | phenotyping `{dataset}/{modality}.bed.gz`, `{dataset}/unnorm/…`; pheast `{dataset}/{qtl,heritability,twas}/…` |
  | `pantry_intermediate` | both | `intermediate` | `{dataset}/…` (BAMs, kallisto, covariates, TWAS batches, LaDDR, `samples.txt`, `summary/`) |
  | `pantry_reference` | phenotyping | `intermediate/reference` | `{reference}/…` shared reference files |
  | `pantry_resources` | pheast | `intermediate/resources` | `fusion/` (installed FUSION scripts and gcta) |

- **Standalone output paths gain a dataset level.** `output/expression.bed.gz` becomes `output/default/expression.bed.gz`, `output/unnorm/…` becomes `output/default/unnorm/…`, and `pheast/output/qtl/…` becomes `pheast/output/default/qtl/…`. The same applies to `intermediate/`: `intermediate/reference/…` becomes `intermediate/reference/default/…`.
- Pathvars are prefixed (`pantry_*`) so that a host's own `results`/`resources` pathvars cannot capture them. A host sets them in its `module … pathvars:` block, and they may refer to the host's pathvars (`pantry_output="<results>/pantry/phenotypes"`).
- Pathvars don't resolve in `params` or shell strings. Directory params are therefore derived from the resolved outputs or inputs (named helpers `parent_dir`, `bgzip_source`, `laddr_project_dir`, …). Bare `{output_dir}`-style globals in shell bodies were replaced the same way.
- Global `wildcard_constraints`:
  - Phenotyping: `dataset`, `reference`, `sample_id`, `pheno`, `laddr_name` = `[^/]+`; `read_length`, `batch`, `n_batches` = `\d+`; `type` = `TSS|polyA`; `group` = `grp_[12]`; `position` = `upstream|contained|downstream`; `feature_type` = `exonic|intronic`.
  - Pheast: `dataset`, `modality`, `chrom`, `batch` = `[^/]+`; `batch_start`, `batch_end` = `\d+`.
  - The enumerations are narrower than the spec's `[^/]+` (ruling R4). They keep dotted sample IDs like `NA1.exonic` unambiguous.
- Generic rules are scoped: `index_bam` matches only `<pantry_intermediate>/{dataset}/bam/{sample_id}.bam`, and `index_bed` matches only `<pantry_output>/{dataset}/{pheno}.bed.gz`. `index_bed` excludes `cross_modality`, whose `.tbi` `combine_modalities` writes.
- All `localrules:` blocks were replaced by per-rule `localrule: True`. The module's `localrules:` names are unprefixed, so they no longer match after `use rule * … as prefix_*`. pheast's empty `localrules:` block was removed (R10).
- `rule all` has `default_target: True` in both stages.
- Scripts and bundled files are referenced only through `workflow.source_path(...)`. This works for local-path and `github()` modules. `pheast/Snakefile`'s `from scripts.twas_batches import …` became a load by path. `twas_assemble_summary`'s `realpath scripts/fusion_twas` now uses the installed FUSION path. Script paths sit in callable params, so the per-process source-cache path doesn't trigger params reruns.
- Nothing depends on the CWD except the standalone `configfile: "config.yml"`, which Snakemake skips inside a module.

## Config format

Both stages take a multi-dataset config: `datasets: {name: {...}}`, plus an optional `defaults:` block that is deep-merged under every dataset (dataset values win; lists replace). The old single-dataset format is still accepted. A config without `datasets` is wrapped as dataset `default` (and, for phenotyping, reference `default`) by `normalise()` in `{phenotyping,pheast}/lib/config.py`, then validated. The repo's `phenotyping/config.yml` and `pheast/config.yml` are in the new format. Their `normalise()` output equals that of the old files (unit tests). The examples are `phenotyping/config.yml`, `pheast/config.yml`, `tests/host/config.yaml` (host), and `tests/configs/*_multi.yml`.

### Phenotyping: old key → new key

| Old (top level) | New |
|---|---|
| `ref_genome` | `references.<ref>.genome` |
| `ref_anno` | `references.<ref>.annotation` (uncompressed GTF with transcript rows, as before) |
| `modality_groups.RNA_editing.edit_sites_bed` | `references.<ref>.edit_sites_bed` (required if any dataset using the reference has RNA_editing) |
| `modality_groups.intron_retention.annotation_gff3` | `references.<ref>.annotation_gff3` (optional; otherwise converted from the GTF) |
| — | `datasets.<ds>.reference: <ref>` (required) |
| `samples_file` | `datasets.<ds>.samples`: a path to a one-column file (as before) **or** an inline list |
| `fastq_map`, `fastq_dir` | `datasets.<ds>.fastq_map`: a path in the existing map format (with `fastq_dir`) **or** an inline dict `{sample: [[R1…], [R2…]]}` (paired) / `{sample: [[R…]]}` (single); relative paths are under `fastq_dir` |
| `paired_end`, `read_length`, `fragment_length_mean`, `fragment_length_sd` | same keys under `datasets.<ds>` |
| `modality_groups` (incl. RNA-editing knobs `edit_sites_min_coverage`, `edit_sites_min_samples`, `edit_sites_min_samples_fraction`, `edit_sites_correlation_threshold`, `edit_sites_cluster`; MAJIQ `license`, `threads`, `psicov_batch_size`, `min_experiments`, `group`) | same, under `datasets.<ds>.modality_groups` |
| `laddr` | `datasets.<ds>.laddr`. `project_dir` is ignored, with a note: the project dir is always `<pantry_intermediate>/{dataset}/laddr`. `dataset_name` defaults to the dataset key (was `pantry` in the template) |
| `intermediate_dir`, `intermediate_ref_dir` | Still honoured in the old format: they set `pantry_intermediate` / `pantry_reference`. In the new format, use Snakemake's config `pathvars:` (`pathvars: {pantry_intermediate: …, pantry_reference: …}`) or a host's module pathvars |
| (hard-coded conda env `majiq3`) | top-level `majiq_env` (default `majiq3`) |
| — | new `datasets.<ds>.cross_modality.modalities: [...]` → `cross_modality.bed.gz`, `.tbi`, `.phenotype_groups.txt` |

### Pheast: old key → new key

| Old (top level) | New |
|---|---|
| `phenotype_dir`, `geno_prefix`, `modalities`, `analyses`, `qtl_min_allele_count`, `qtl_maf_threshold`, `twas_snps`, `twas_models_per_job`, `twas_threads` | same keys under `datasets.<ds>` (or `defaults`) |
| `samples_file` | `datasets.<ds>.samples`: path or inline list. Meaning unchanged: used for the genotype check and N in the TWAS `.pos` file; it does not subset the analysis |
| `intermediate_dir` | Old format: sets `pantry_intermediate`, and `pantry_resources` = `<intermediate_dir>/resources`. New format: config `pathvars:` or module pathvars |
| — | new `datasets.<ds>.qtl_cis_window` → `run_tensorqtl.py --window` (unset: the script's own default, 1000000, applies as before) |
| (FUSION assumed at `scripts/fusion_twas/`) | top-level `fusion_dir`, `gcta` (overrides), `gcta_url`, `gcta_sha256` (see [FUSION](#fusion)) |
| `modalities.<m>.chroms` (set internally at parse time) | gone; read from the `phenotype_summary` checkpoint |

`phenotype_dir` is read through input functions. A host that points it at the phenotyping module's output therefore links the two stages in one DAG.

## Validation

- JSON schemas `phenotyping/schemas/config.schema.yaml` and `pheast/schemas/config.schema.yaml` are applied with `snakemake.utils.validate()` after `normalise()`. Defaults that upstream set with `setdefault` (LaDDR knobs, RNA-editing knobs, MAJIQ knobs, `qtl_min_allele_count` 10, `twas_models_per_job` 200, `twas_threads` 4, `majiq_env`) moved into the schemas with the same values.
- Cross-field checks a schema can't express stay in Python and run on the config only:
  - fragment lengths are required for single-end data;
  - RNA_editing needs `edit_sites_bed` on the dataset's reference;
  - the reference must exist;
  - `cross_modality.modalities` must be a list of modalities with BED outputs, without duplicates;
  - `laddr.regress_modalities` must be a subset of the available modalities (its default is all of them);
  - existing LaDDR bins need both dirs;
  - pheast `cross_modality` is allowed only with `qtl` and without trans-QTL (upstream's message);
  - dataset and reference names must be single path components.
- **Stricter than upstream (rulings R14, R18):**
  - Integer knobs reject quoted numbers (`"10"`) and non-integral floats. Upstream `int()`-coerced some of them, e.g. `qtl_min_allele_count`. Integral floats such as `40.0` are accepted (JSON Schema treats them as integers).
  - TWAS knobs (`twas_models_per_job`, `twas_threads`) are checked by the schema, so jsonschema's message replaces upstream's "must be at least 1".
- **New:** a config with `datasets:` and also legacy top-level keys (e.g. `read_length`, `samples_file`, `geno_prefix`) is rejected instead of silently ignoring them. Snakemake merges `--configfile` with the workflow's `config.yml`, so an old-format `--configfile` used inside `phenotyping/` or `pheast/` would otherwise have run the example `default` dataset of `config.yml`.
- The phenotyping annotation check (transcript rows present) is now a job, `validate_annotation`. The pheast genotype-vs-samples check is now a job, `check_genotype_samples`, with the same message as upstream. See below.

## References

- References are keyed by **named reference sets**: `<pantry_reference>/{reference}/{cDNA.fa.gz, kallisto.idx, star_index_{read_length}/, txrevise/…, exonic.gtf, intronic.gtf, majiq/annotation.gff3, edit_sites_to_genes.tsv, annotation.ok}`. The STAR index is also keyed by read length. Datasets that name the same reference share all of these. Each file has exactly one producing rule (Task 13 ran 2 datasets on one reference: 1 STAR index job, 1 kallisto index job).
- `edit_sites_bed` and `annotation_gff3` are reference data and moved to the reference set. The shim moves them automatically.

## Parse-time reads moved to rules and checkpoints

- Phenotyping: sample lists and FASTQ maps given as paths are read lazily on first use and cached per dataset, at DAG-build time, not at load. The inline forms read no file. The mixed/inconsistent endedness checks moved into that loader. The GTF transcript-row check is the job `validate_annotation`, an input of the reference-building rules.
- Pheast:
  - the `.fam` sample check is the job `check_genotype_samples`, an input of `prune_for_covar` and `qtl_genotype_qc`;
  - per-modality chromosomes (heritability) and model counts (TWAS batches) come from the **checkpoint** `phenotype_summary` per `{dataset},{modality}`, read by `combine_heritability_chr` and `twas_assemble_summary`. A cold dry run therefore shows fewer jobs (69 vs upstream's 185 for the test config); the rest appear after the checkpoints run (195 = 185 + 10 new jobs);
  - `write_samples` materialises `<pantry_intermediate>/{dataset}/samples.txt` for scripts that take `--samples FILE`. Its contents are identical to the old samples file.

## New and removed rules

- Phenotyping, new:
  - `write_samples`, `validate_annotation` (above);
  - `combine_modalities` and `pheno_groups_cross_modality`: the cross-modality BED/`.tbi`/groups per dataset. The body reproduces `scripts/combine_modalities.sh`, and the outputs are byte-identical to the script's on the same BEDs (Task 14). The script stays for manual use;
  - static versions of rules that upstream defined conditionally at parse time: `majiq_gff3_from_gtf` (used only when no `annotation_gff3` is configured), `laddr_setup`/`laddr_binning` vs `laddr_prepare_existing_bins` (distinct `bins_ready.built` / `bins_ready.existing` outputs; the consumer picks one).
- Phenotyping, changed:
  - `majiq_build` outputs `directory(build)` and the nested `splicegraph.zarr`. `majiq_psi_coverage` derives the per-sample `.sj` paths in params, resolving BAM symlinks at job time (keeps fix `21ec32d`);
  - LaDDR's output, whose name embeds `dataset_name`, is moved to a stable path `<pantry_intermediate>/{dataset}/laddr/latent_phenos.tsv.gz`.
- Phenotyping, **removed: `index_fasta`** (R15). No rule used it, and under the new layout nothing could request it. This deviates from spec Corrections 5, which said it would be kept.
- Pheast, new: `check_genotype_samples`, checkpoint `phenotype_summary`, `install_fusion`.
- All step files are always included. Which outputs a dataset gets depends on its `modality_groups` / `analyses` in `rule all` (spec §6).
- Thread and resource knobs that already existed (`twas_threads`, MAJIQ `threads`, LaDDR `bam_coverage_threads`) are kept for compatibility.

## Environments, containers, GPU

- Per-rule conda envs, split by tool and pinned to the versions in the baseline env lock files (`tests/hpc/envs/baseline-*.lock.txt`):
  - `phenotyping/envs/`: `star`, `samtools`, `kallisto` (0.48.0), `gffread`, `regtools`, `subread`, `python`, `txrevise`, `rna_editing`, `deeptools`, `laddr`;
  - `pheast/envs/`: `plink`, `covariates`, `python`, `gcta`, `fusion`.

  Built envs match the baseline env for every tool and pinned library (Task 12). MAJIQ keeps a named env (`majiq_env`, default `majiq3`). The monolithic `environment.yml` files stay for users not using `--sdm conda`.
- **tensorQTL rules** (`tensorqtl_cis`, `tensorqtl_cis_independent`, `tensorqtl_trans`) have **no `conda:`**, only `container: docker://francois4/tensorqtl@sha256:f6efb9e5…` (R20). This is genex's image pinned by digest; Snakemake runs conda inside a container when a rule has both.
  - **With `--sdm apptainer`**, they run tensorqtl **1.0.9** / torch 2.1.2 / pandas 2.1.4 / Python 3.10. Upstream's pip env had **1.0.10** / 2.3.1 / 2.3.3 / 3.12 (R23). On CPU, the image alone moves `pval_nominal` by up to ~4e-5 relative; the GPU adds a similar float-level difference.
  - **Without apptainer**, these rules run the host `python3`, which must have tensorqtl (e.g. the `pheast/environment.yml` env), as upstream did.
- GPU:
  - The profile must give the tensorQTL rules a GPU partition and `gpu: 1` (the slurm executor plugin rejects `--gpus`/`--gres` in `slurm_extra`), plus `apptainer-args: "--nv --bind=<dirs>"`. `--nv` exposes the GPU. The bind must cover the data, including symlink targets, and for a local-path module the module checkout.
  - In a host, `set-resources` keys are the **prefixed** rule names (e.g. `pantry_pheast_tensorqtl_cis`).
  - Examples: `tests/profile/config.yaml` (standalone) and `tests/host/profile/config.yaml` (host).
- Not pinned yet: the `txrevise` env's python (unused by the scripts), and LaDDR, which is pip-installed from GitHub HEAD.

## FUSION

- FUSION is vendored upstream (`pheast/scripts/fusion_twas/`: `FUSION.compute_weights.R`, `utils/plink_utils.R`, `utils/FUSION.profile_wgt.R`, Linux x86-64 `gcta_nr_robust`; gemma comes from the `fusion` env). `install_fusion` copies the scripts into `<pantry_resources>/fusion/`, keeping the `utils/` layout that `compute_weights` sources relative to itself, and makes gcta executable. The source cache doesn't preserve the exec bit.
- Overrides: `fusion_dir` (a directory with the three R files) and/or `gcta` (a binary path). Each one removes only its own files from `install_fusion`'s outputs.
- **R24 (Task 14):** Snakemake's git source cache can't serve binary files. gcta is therefore resolved only in `install_fusion`'s job:
  - configured `gcta` → nothing installed for it;
  - module loaded from a local path → bundled file copied;
  - otherwise (e.g. `github()`) → downloaded from `gcta_url` (default: upstream `8da7787`'s file on raw.githubusercontent.com), with a 300 s timeout, and checked against `gcta_sha256` (default `7b51e51d…93e9`).

  **Remote modules need network access in that job. Offline sites set `gcta`.** The spec's §7 FUSION bullet carries a dated amendment for this.
- `install_fusion` and the FUSION files are inputs of the TWAS rules, so FUSION can now run from any working directory.

## Tests

- `tests/unit/` (pytest, 146 tests): config shim/normalise/validation, the FASTQ-map loader against an oracle copy of upstream's (R5), samples, targets, phenotype summaries against upstream's reads, FUSION paths, gcta install, and the comparison tool.
- `tests/compare_outputs.py`: compares an output tree against a baseline.
  - Exact for text files and decompressed `.gz`.
  - With `--noise A B`: BED-like files are compared numerically within 2× the noise pair's max |diff| (R12); `cis_qtl` permutation columns and tie-breaks are judged against the noise pair; `cis_independent` is judged on summary statistics (R8, R9).
  - Writes a markdown report.
- `tests/configs/`: baseline (old-format) configs for both stages, multi-dataset configs, and dry-run-only configs (`dryrun_multi.yml` with a dotted dataset name; `dryrun_majiq_laddr.yml` for MAJIQ and LaDDR).
- `tests/fixtures/`: `make_multi_inputs.py` (derives the single-end and 300-sample test datasets), `check_targets.py`.
- `tests/host/`: a host workflow that imports both stages by local path (`Snakefile`) or with `github()` at `module-rc1` (`Snakefile.github`), plus its profile.
- `tests/profile/` (SLURM test profile), `tests/hpc/` (run records, baseline envs).

## Lint and formatting

`snakemake --lint` was run on both stages (phenotyping with `tests/configs/dryrun_majiq_laddr.yml`, which enables every rule incl. MAJIQ and LaDDR; pheast with `tests/configs/pheast_baseline.yml`) and on `tests/host`. Snakemake 9.27.1.

Compared with upstream `8da7787` (default configs), these are gone:
- every "Specify a conda environment or container" (41 phenotyping + 13 pheast);
- "Shell command directly uses variable output_dir" (6);
- "Param … is a prefix of input or output file but hardcoded" (8 + 7);
- "Path composition with '+'" (1).

Remaining warnings, each justified:

| Warning | Where | Justification |
|---|---|---|
| No log directive defined | every rule (66 phenotyping, 17 pheast) | Upstream has no `log:` anywhere. Adding logs means adding `2> {log}` to every rule body, which the no-rule-body-rewrite constraint excludes. Follow-up. |
| Mixed rules and functions in same snakefile | `phenotyping/Snakefile`, `align`, `alt_TSS_polyA`, `intron_retention`, `laddr`, `stability`; pheast `qtl`, `heritability`, `twas` | Per-dataset input/params functions live beside the rules they serve (spec §6). Moving them all to `setup.smk` would separate each function from its only user. |
| Migrate long run directives into scripts | `write_samples`, `make_majiq_experiments_tsv`, `laddr_manifest`; `plink_covariates`, `install_fusion` | `run:` exception (spec §6). The full set of `run:` rules is `write_samples`, `validate_annotation`, `make_majiq_experiments_tsv`, `laddr_manifest`, `plink_covariates`, `check_genotype_samples`, `phenotype_summary`, `install_fusion`; lint flags only the five long ones. They are stdlib/pandas only (Snakemake provides both), are local, and are small. `plink_covariates`, `make_majiq_experiments_tsv` and `laddr_manifest` are upstream's own `run:` blocks. |
| Do not access input and output files individually by index | `combine_heritability_chr` (`head -n 1 {input[0]}`) | Upstream's rule body, unchanged. |
| Absolute path "/transcript_id "";/d" | `extract_cdna` in `phenotyping/Snakefile` | False positive: a `sed` delete expression, not a path. Upstream body. |

`snakefmt` (2.0.3): `snakefmt --check --exclude 'run_twas\.smk' phenotyping pheast tests/host` passes on all 20 files, and `tests/host/Snakefile.github` is clean too. The one exclusion, `pheast/scripts/run_twas.smk`, is upstream's standalone TWAS-association example. The module doesn't load it, and it is left byte-identical to upstream so that `pheast/scripts` has no diff (the tensorQTL attribution evidence in R21/R22 relies on that).

## Migrating a single-dataset config

1. **Nothing is required to keep running.** An old-format `config.yml` in `phenotyping/` or `pheast/` still works and becomes dataset `default` (and reference `default`). Replace the repo's `config.yml` with yours, or run from another directory with `-s …/Snakefile --directory <dir>` where `<dir>/config.yml` is yours. Don't pass an old-format `--configfile` while the repo's new-format `config.yml` is present: Snakemake merges the two, and this is now an error rather than a silent use of the example dataset.
2. **Outputs move**: `output/X` → `output/default/X`, `intermediate/X` → `intermediate/default/X`, references under `intermediate/reference/default/`. Pheast's `phenotype_dir: ../phenotyping/output` therefore becomes `../phenotyping/output/default`. An existing run's intermediates are not picked up at the new paths, so expect a full rerun (or move them yourself).
3. To convert to the new format (recommended; see `phenotyping/config.yml` and `pheast/config.yml`):
   - phenotyping: move `ref_genome`/`ref_anno` to `references: {<ref>: {genome, annotation}}`, and `modality_groups.RNA_editing.edit_sites_bed` (and `intron_retention.annotation_gff3`) to that reference. Put everything else under `datasets: {<name>: {reference: <ref>, …}}`. Rename `samples_file` → `samples`. Drop `laddr.project_dir`. Set `laddr.dataset_name: pantry` if you relied on upstream's template value.
   - pheast: put everything under `datasets: {<name>: {…}}` and rename `samples_file` → `samples`.
   - Both: replace `intermediate_dir` / `intermediate_ref_dir` with top-level `pathvars: {pantry_intermediate: …, pantry_reference: …}` (pheast: `pantry_intermediate`, `pantry_resources`).
   - Integer knobs must be unquoted integers.
4. **tensorQTL**: run with `--sdm conda apptainer` to use the pinned container (tensorqtl 1.0.9). Otherwise these three rules use your `python3`'s tensorqtl, as upstream did. GPU profiles: see [Environments](#environments-containers-gpu).
5. Snakemake ≥ 9.12.

## Deviations from the handover contract and the spec

- **Contract 1 (one `config["pantry"]`)**: the stages are two modules, each with its own config section (`config["pantry"]["phenotyping"]`, `config["pantry"]["pheast"]`) and distinct prefixes (`pantry_pheno_*`, `pantry_pheast_*`). Both define `all`, so one shared prefix would clash.
- **Contract 8 (fetch a pinned FUSION commit and GCTA/GEMMA)**: FUSION is vendored with Pantry's modifications and installed from the module. gcta is copied, or downloaded from a pinned URL with a sha256 check for remote modules (R24). gemma comes from the `fusion` conda env.
- **Contract 11 / criterion 11 (standalone)**: works, but standalone outputs gain the `default/` level.
- **Acceptance criterion 1 for tensorQTL**: tensorQTL is unseeded and nondeterministic upstream, and the container's tensorqtl is 1.0.9, so `*.cis_qtl`/`*.cis_independent_qtl` were judged by attribution rather than equality (R21, R22): identical inputs and command lines, comparison reported but not blocking. kallisto outputs were judged within the run-to-run noise envelope (R11, R12). See `docs/module-status.md`.
- **MAJIQ (`intron_retention`) and LaDDR** were restructured and checked by dry run only (design decision); runtime verification is a follow-up.
- Spec Corrections 5: `index_fasta` removed, not kept (R15).
- Spec §2 constraints: some wildcards are narrower enumerations (R4).
- Spec §6: LaDDR `dataset_name` defaults to the dataset key, so legacy configs that omit it get `default`, not upstream's `pantry` (R16). This affects only LaDDR's internal manifest names.
- Spec §7: gcta download at job time for remote modules (R24; spec amended).
- Validation is stricter on quoted or non-integral numbers (R14, R18).

## Known limitations and follow-ups

- **Runtime verification of MAJIQ and LaDDR** (dry run only). Not covered by any run: the nested directory outputs of `majiq_build`, the `.sj` derivation, and the LaDDR `mv` to a stable name.
- **Env pinning**: the tensorQTL image with tensorqtl 1.0.10 to match upstream (R23 follow-up); the txrevise python; LaDDR from GitHub HEAD; MAJIQ's named env.
- **Snakemake source-cache bugs** (9.27.x, `sourcecache.py` `HostingProviderFile.open` / `LocalGitFile._open`; to report upstream, not filed). `git show` output is decoded as text, so (1) `source_path` on a binary file in a `github()` module fails with `UnicodeEncodeError` (worked around by R24), and (2) text files lose their final newline: FUSION R scripts installed from `github()` are 1 byte shorter. This is harmless. See `tests/hpc/README.md`, Task 14.
- **Nondeterminism**: kallisto 0.48.0 quant varies run to run (node/load dependent; not thread count). tensorQTL permutations and tie-breaks are unseeded. The noise envelopes come from a single pair of baseline runs, and the isoforms kallisto margin in Task 13 was close (max |diff| 0.149 vs allowed 0.167).
- **Legacy pathvars vs host pathvars**: a legacy `intermediate_dir`/`intermediate_ref_dir` in a module's config overrides the host's module pathvars, because it registers later. Hosts should pass new-format configs only.
- `--configfile` with a new-format config inside `phenotyping/` or `pheast/` is merged with the repo's `config.yml`, which adds its `default` dataset. Use `--directory` to point elsewhere, or edit `config.yml`. This merging is standard Snakemake behaviour.
- **Rerun triggers**: script paths are derived params, kept out of the params trigger because the source-cache path changes per process. As upstream (scripts called from shell), editing a bundled script doesn't trigger a rerun; use `-R <rule>`. The configured MAJIQ `license` is a derived param too, so changing it doesn't trigger a rerun. Wildcard-only lambda params remain in the params trigger (stable values). `install_fusion` has no inputs, so editing the vendored FUSION does not reinstall it (delete `<pantry_resources>/fusion/`).
- `combine_modalities` keeps the script's behaviour: `set +o pipefail` (the script runs without pipefail) lets a truncated `.bed.gz` through (a `gzip -t` check would catch it), and the `sort` locale is unpinned. Left as upstream to keep outputs identical to the script's.
- LaDDR: if the final `mv` in `laddr_transform` fails, a stale latent file from an earlier run can remain in the LaDDR project dir.
- The R scripts in `install_fusion` are copied without a temp file. This is harmless: the job output is removed on failure.
- Test tooling: the comparison tool's floors (1 row; 1e-6 relative) tolerate a single variant-changed key or NaN mismatch in `cis_independent`, and RDat pairs after an OSError are not compared. The fixture/checker scripts are not unit-tested; `.tbi` files are checked by size only. The setup.smk injection of `annotation_gff3` is covered by dry run only. Intermediate commits `ad85347..c5c3ab2` don't dry-run on their own (bisect across them with care).
- pheast `samples` still doesn't subset the analysis (upstream semantics): to analyse fewer samples, supply BEDs with those samples.
- Ideas not taken: the standard `workflow/` layout; wrapping `prepare_phenotypes.py` in a rule.
