# Pantry as a reusable, multi-dataset Snakemake module: changes

Branch `module` of `twillis209/Pantry`, based on upstream `PejLab/Pantry` `8da7787`. Release tag `module-v0.1.6` (earlier: `module-v0.1.5`, `module-v0.1.4`, `module-v0.1.3`, `module-v0.1.2`, `module-v0.1.1`, `module-v0.1.0`, `module-rc1`).
Design: `docs/superpowers/specs/2026-10-30-module-refactor-design.md`. Status and acceptance results: `docs/module-status.md`.
Test commands and HPC run records: `tests/hpc/README.md`.

**What it does.** `phenotyping/` and `pheast/` can now each be imported once as a Snakemake `module` and run any number of datasets in one DAG. The host sets the output roots with pathvars. Datasets that name the same reference share its reference files. Standalone use inside `phenotyping/` or `pheast/` still works, and the old single-dataset config format is still accepted.

**What it does not change.** Analysis logic, default parameters and output file formats are unchanged, with one deliberate exception: tensorQTL's cis and cis_independent mapping are seeded by default (R30, `module-v0.1.6`; `qtl_seed: null` runs as upstream). Rule bodies changed only where a path or parameter now comes from somewhere else, plus that seed. `git diff 8da7787 -- pheast/scripts` shows only `run_tensorqtl.py --seed`. In `phenotyping/scripts`, only a usage comment was added to `combine_modalities.sh`. The exceptions to "same software" are the runtimes of tensorQTL and FUSION: under apptainer they run in pinned containers (see [Environments](#environments-containers-gpu) and [FUSION](#fusion)).

Requires **Snakemake ≥ 9.14** (R28, `module-v0.1.4`; pathvars need 9.12). Both Snakefiles set `min_version("9.14")`: 9.14.0 is the first release with Snakemake commit f1e8b629 (PR #3840), which binds the source cache path actually in use and the local storage prefix into container jobs. Tested with 9.27.1 (HPC runs and local dry runs) and 9.27.0 (local spike).

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
  | `pantry_resources` | pheast | `intermediate/resources` | nothing since `module-v0.1.2` (was `fusion/`, the installed FUSION scripts and gcta); kept so hosts that set it still load |

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
- Scripts and bundled files are referenced only through `workflow.source_path(...)`. This works for local-path and `github()` modules. `pheast/Snakefile`'s `from scripts.twas_batches import …` became a load by path. `twas_assemble_summary`'s `realpath scripts/fusion_twas` now uses the FUSION path in the container (`/opt/fusion`). Script paths sit in callable params, so the per-process source-cache path doesn't trigger params reruns.
- Nothing depends on the CWD except the standalone `configfile: "config.yml"`. Snakemake skips it inside a module only when the host's `module` directive passes `config:` (`snakemake/modules.py`: `skip_configfile = config is not None`). **Hosts must pass `config:`**; without it, the module's `configfile:` would load a `config.yml` from the host's working directory, if there is one.

## Config format

Both stages take a multi-dataset config: `datasets: {name: {...}}`, plus an optional `defaults:` block that is deep-merged under every dataset (dataset values win; lists replace). A dataset can add keys to `defaults` or override them, but cannot remove one: a modality group (phenotyping `modality_groups`) or modality/analysis (pheast `modalities`, `analyses`) in `defaults` is run for every dataset. Keep modality sets that vary between datasets out of `defaults` and list them per dataset. (Removal support is deferred.) The old single-dataset format is still accepted. A config without `datasets` is wrapped as dataset `default` (and, for phenotyping, reference `default`) by `normalise()` in `{phenotyping,pheast}/lib/config.py`, then validated. The repo's `phenotyping/config.yml` and `pheast/config.yml` are in the new format. Their `normalise()` output equals that of the old files (unit tests). The examples are `phenotyping/config.yml`, `pheast/config.yml`, `tests/host/config.yaml` (host), and `tests/configs/*_multi.yml`.

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
| (hard-coded conda env `majiq3`) | top-level `majiq_env` (default `majiq3`), or top-level `majiq_container` (a user-built image, which takes precedence; R27) |
| — | new `datasets.<ds>.cross_modality.modalities: [...]` → `cross_modality.bed.gz`, `.tbi`, `.phenotype_groups.txt` |

### Pheast: old key → new key

| Old (top level) | New |
|---|---|
| `phenotype_dir`, `geno_prefix`, `modalities`, `analyses`, `qtl_min_allele_count`, `qtl_maf_threshold`, `twas_snps`, `twas_models_per_job`, `twas_threads` | same keys under `datasets.<ds>` (or `defaults`) |
| `samples_file` | `datasets.<ds>.samples`: path or inline list. Meaning unchanged: used for the genotype check and N in the TWAS `.pos` file; it does not subset the analysis |
| `intermediate_dir` | Old format: sets `pantry_intermediate`, and `pantry_resources` = `<intermediate_dir>/resources`. New format: config `pathvars:` or module pathvars |
| — | new `datasets.<ds>.qtl_cis_window` → `run_tensorqtl.py --window` (unset: the script's own default, 1000000, applies as before) |
| — | new `datasets.<ds>.qtl_seed` (integer 0 to 2³²−1, or `null`; **default 0**) → `run_tensorqtl.py --seed` for `tensorqtl_cis` and `tensorqtl_cis_independent` (`module-v0.1.6`, R30). A top-level `qtl_seed` in the old format goes to dataset `default` |
| (FUSION assumed at `scripts/fusion_twas/`) | no key: FUSION is in the TWAS rules' container (see [FUSION](#fusion)). `fusion_dir`, `gcta`, `gcta_url`, `gcta_sha256` (added in `module-v0.1.0`/`v0.1.1`) were removed in `module-v0.1.2` and are rejected with an error |
| `modalities.<m>.chroms` (set internally at parse time) | gone; read from the `phenotype_summary` checkpoint |

`phenotype_dir` is read through input functions. A host that points it at the phenotyping module's output therefore links the two stages in one DAG. The link is by path text: `phenotype_dir` must equal the phenotyping output path `<pantry_output>/{dataset}` as the host resolves it (e.g. `results/pantry/phenotypes/full`). A relative and an absolute spelling of the same directory do not match, and pheast then fails with missing input files instead of waiting on the phenotyping rules.

## Validation

- JSON schemas `phenotyping/schemas/config.schema.yaml` and `pheast/schemas/config.schema.yaml` are applied with `snakemake.utils.validate()` after `normalise()`. Defaults that upstream set with `setdefault` (LaDDR knobs, RNA-editing knobs, MAJIQ knobs, `qtl_min_allele_count` 10, `twas_models_per_job` 200, `twas_threads` 4, `majiq_env`) moved into the schemas with the same values.
- Cross-field checks a schema can't express stay in Python and run on the config only:
  - fragment lengths are required for single-end data: with `paired_end: false` at config time, and with `paired_end` omitted when the loaded FASTQ map has a single-end sample (an error naming the dataset, raised when the map is read; before, a `KeyError` in a params function);
  - RNA_editing needs `edit_sites_bed` on the dataset's reference;
  - the reference must exist;
  - `cross_modality.modalities` must be a list of modalities with BED outputs from `modality_groups`, without duplicates. LaDDR's `latent_residual` is not accepted (upstream's manual script would take any BED in `output/`);
  - `laddr.regress_modalities` must be a subset of the available modalities (its default is all of them);
  - existing LaDDR bins need both dirs;
  - pheast `cross_modality` is allowed only with `qtl` and without trans-QTL (upstream's message);
  - dataset and reference names must be single path components.
- **Stricter than upstream (rulings R14, R18):**
  - Integer knobs reject quoted numbers (`"10"`) and non-integral floats. Upstream `int()`-coerced some of them, e.g. `qtl_min_allele_count`. Integral floats such as `40.0` are accepted (JSON Schema treats them as integers).
  - TWAS knobs (`twas_models_per_job`, `twas_threads`) are checked by the schema, so jsonschema's message replaces upstream's "must be at least 1".
- **New:** a config with `datasets:` and also legacy top-level keys (e.g. `read_length`, `samples_file`, `geno_prefix`) is rejected instead of silently ignoring them. Snakemake merges `--configfile` with the workflow's `config.yml`, so an old-format `--configfile` used inside `phenotyping/` or `pheast/` would otherwise have run the example `default` dataset of `config.yml`.
  The same applies to `--config key=value` with a legacy key (e.g. `--config read_length=100`) inside a stage directory: it is now an error. The new-format equivalent is a `--configfile` holding just the override, which Snakemake deep-merges into `config.yml`:

  ```yaml
  # override.yml, used as: snakemake --configfile override.yml
  datasets:
    default:
      read_length: 100
  ```

  Or edit `config.yml`, or run an old-format config from another directory (`-s <path>/Snakefile --directory <dir>`, with `<dir>/config.yml` yours). `--config 'datasets={default: {read_length: 100}}'` works only for string values: Snakemake parses nested `--config` values as strings, so `100` arrives as `'100'` and fails the integer check.
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
  - `combine_modalities` and `pheno_groups_cross_modality`: the cross-modality BED/`.tbi`/groups per dataset. The body reproduces `scripts/combine_modalities.sh`, and the outputs are byte-identical to the script's on the same BEDs (Task 14). One difference: if the BED headers (sample columns) differ, the job fails, as its shell runs with `set -e`; the script printed the `diff` and went on. The script stays for manual use;
  - static versions of rules that upstream defined conditionally at parse time: `majiq_gff3_from_gtf` (used only when no `annotation_gff3` is configured), `laddr_setup`/`laddr_binning` vs `laddr_prepare_existing_bins` (distinct `bins_ready.built` / `bins_ready.existing` outputs; the consumer picks one).
- Phenotyping, changed:
  - `majiq_build` outputs `directory(build)` and the nested `splicegraph.zarr`. `majiq_psi_coverage` derives the per-sample `.sj` paths in params, resolving BAM symlinks at job time (keeps fix `21ec32d`);
  - LaDDR's output, whose name embeds `dataset_name`, is moved to a stable path `<pantry_intermediate>/{dataset}/laddr/latent_phenos.tsv.gz`.
- Phenotyping, **removed: `index_fasta`** (R15). No rule used it, and under the new layout nothing could request it. This deviates from spec Corrections 5, which said it would be kept.
- Pheast, new: `check_genotype_samples`, checkpoint `phenotype_summary`. (`install_fusion`, added in `module-v0.1.0`, was removed in `module-v0.1.2`; R26.)
- All step files are always included. Which outputs a dataset gets depends on its `modality_groups` / `analyses` in `rule all` (spec §6).
- Thread and resource knobs that already existed (`twas_threads`, MAJIQ `threads`, LaDDR `bam_coverage_threads`) are kept for compatibility.

## Environments, containers, GPU

- Per-rule conda envs, split by tool and pinned to the versions in the baseline env lock files (`tests/hpc/envs/baseline-*.lock.txt`):
  - `phenotyping/envs/`: `star`, `samtools`, `kallisto` (0.48.0), `gffread`, `regtools`, `subread`, `python`, `txrevise`, `rna_editing`, `deeptools`, `laddr`;
  - `pheast/envs/`: `plink`, `covariates`, `python`, `gcta`, `fusion` (`fusion` is also the env of the FUSION image; as a conda env it is used only by `twas_compress_output`).

  Built envs match the baseline env for every tool and pinned library (Task 12). MAJIQ uses an image the user builds (`majiq_container`, R27) or else a named env (`majiq_env`, default `majiq3`). The monolithic `environment.yml` files (upstream's, unchanged) stay for users not using `--sdm conda`, but they cover only the conda rules: the tensorQTL and TWAS rules, and MAJIQ with `majiq_container`, need apptainer. Since `module-v0.1.4` each file has a header saying so and that the per-rule `envs/` are the maintained route.
- **tensorQTL rules** (`tensorqtl_cis`, `tensorqtl_cis_independent`, `tensorqtl_trans`) have **no `conda:`**, only `container: docker://ghcr.io/twillis209/pantry-tensorqtl@sha256:65a42f9f…` (R20), pinned by digest. Snakemake runs conda inside a container when a rule has both.
  - The image is built from `containers/tensorqtl/` by `.github/workflows/tensorqtl-image.yml` (public). Pull it by digest only (`@sha256:65a42f9f…`, the tested build): image builds are not byte-reproducible, and every rebuild moves the `1.0.10` tag, which now points at an untested rebuild (`111e79ff…`, run 41737670211). It has upstream's pip env versions: tensorqtl **1.0.10**, torch 2.3.1+cu121, pandas 2.3.3, Python 3.12, numpy 2.5.4, scipy 1.18.1, rpy2 + R qvalue (Task 16; supersedes R23, genex's 1.0.9 image). Only pyarrow differs (25.0.1 vs 26.0.0), used only by `cis_nominal`.
  - The image is public; apptainer pulls it anonymously. Stale stored ghcr.io credentials (`apptainer registry login`, `~/.docker/config.json`) make the pull fail with 403; log out of ghcr.io or refresh them.
  - **With `--sdm apptainer`** on 1 CPU core, as upstream's jobs ran, cis nominal statistics (`pval_nominal`, `slope`, `slope_se`, `af`) are identical to upstream's on rows with the same top variant, at written precision (Task 16). On a GPU, or on several CPU cores, `pval_nominal` differs by up to ~1e-4 relative (reduction order).
  - **Without apptainer**, these rules run the host `python3`, which must have tensorqtl (e.g. the `pheast/environment.yml` env), as upstream did.
- GPU:
  - The profile must give the tensorQTL rules a GPU partition and `gpu: 1` (the slurm executor plugin rejects `--gpus`/`--gres` in `slurm_extra`), plus `apptainer-args: "--nv --bind=<dirs>"`. `--nv` exposes the GPU. The bind must cover the data, including symlink targets, and for a local-path module the module checkout.
  - In a host, `set-resources` keys are the **prefixed** rule names (e.g. `pantry_pheast_tensorqtl_cis`).
  - Examples: `tests/profile/config.yaml` (standalone) and `tests/host/profile/config.yaml` (host).
- **Full pins (Task 21, `module-v0.1.4`).** `phenotyping/envs/txrevise.yaml` pins every package of the env, transitive ones and python included, as `name=version=build` from `tests/hpc/envs/baseline-pheno.lock.txt`; the built env's `conda list --explicit` is a subset of that lock (272 packages, none differing). On the gate-2 config its txrevise reference files are equivalent to the baseline's: the GTF, tags and the four merged `txrevise.{grp_1,grp_2}.{upstream,downstream}.gff3` (all the baseline has) byte for byte; the 600 batch GFF3s except their `##date` line; `txrevise_annotations.rds` except its provenance metadata (`Data source` path and `Creation time`), which differ between any two runs (Task 21). `phenotyping/envs/laddr.yaml` pins LaDDR to release v1.0.0 (commit `2addf6fe`), and every conda package from the same lock except scipy: LaDDR requires scipy < 1.17, so the env has conda-forge scipy 1.16.3 (the lock's 1.17.1 was silently replaced by pip before). Its other pip dependencies (bx-python, scikit-fda, statsmodels, tqdm and theirs) are in no lock and resolve at build time, except `multimethod==1.12`: with multimethod 2.x, scikit-fda 0.10.1 fails to import, so `laddr` failed to start (Task 22 moved pyBigWig, matplotlib, pillow and wrapt from pip to lock pins, below). The image-pins tree lines do not cover the image workflow files (`.github/workflows/{tensorqtl,fusion}-image.yml`), which set labels and build context; both changed after the pinned builds, so recording them would need a rebuild. The 9.14 floor is checked against Snakemake's source; every run used 9.27.1. `containers/majiq/majiq3.yml` pins every conda package (compilers, sysroot and kernel headers included) to `majiq3.lock.txt`, and `containers/majiq/constraints.txt` (used by `majiq3-install.sh`) the lock's pip packages.
- **LaDDR verified at run time (Task 22, R29, `module-v0.1.5`).**
  - Env fixes:
    - `deeptools.yaml` pinned only `deeptools=3.5.5` and resolved python 3.13, numpy 2.5.4, pyBigWig 0.3.25, pysam 0.24.0 and scipy 1.18.1. It now pins deeptools' dependency closure in the baseline env (plus pip, setuptools, wheel), all `name=version=build` from the lock; the built env has 123 conda packages, all in the lock.
    - `laddr.yaml` now pins `pybigwig`, `matplotlib-base`, `pillow`, `wrapt`, `cffi` and their closure from the lock, instead of letting pip install pyBigWig 0.3.29, matplotlib 3.11.2, pillow 13.0.0 and wrapt 2.5.2.
    - `laddr.yaml` pins **scikit-learn 1.5.2** (`py312h7a48858_1`) instead of the lock's 1.9.1, because 1.9.1 was seen to fail: `laddr coverage` passes PCA a frame whose transpose has duplicate integer column names (a `pd.concat` of per-file `reset_index()` frames), and scikit-learn 1.9.1 rejects duplicate column names (narwhals `DuplicateError`). 1.5.2 is a known-good version, chosen because LaDDR's `requirements.txt` names it; that file's other pins (pandas, scikit-fda, statsmodels, bx-python, pyBigWig) are not followed, and LaDDR's `pyproject.toml` leaves scikit-learn unbounded. **Upstream 8da7787 cannot run LaDDR v1.0.0 in its own env** (`pip install` of LaDDR into it, as upstream's README says, keeps scikit-learn 1.9.1).
  - Method:
    - Upstream ran twice in a clone of the baseline env plus LaDDR (with the same scikit-learn and scipy builds).
    - Every module that LaDDR's CLI and bamCoverage import has the same version on both sides.
    - LaDDR v1.0.0 is nondeterministic: its PCA in `remove_pcs_from_phenos` uses scikit-learn's randomized solver unseeded. Two runs on identical inputs differ in coverage residuals, models and latent phenotypes, even in the phenotype count. Kallisto noise in the regressed modalities adds to this.
    - Acceptance was therefore on identical inputs (U1's regressed BEDs substituted), with numpy's global RNG seeded by a test-only `sitecustomize.py`.
  - Results:
    - Seeded upstream is byte-identical to itself.
    - The module standalone and the `github()` host match it in every LaDDR file: the finals byte-identical (`cmp` on the raw files), the intermediates byte-identical, except the `.gz` ones (`gene_bins/*.bed.gz`, `covg_norm/*.bins.tsv.gz`, `info/exons.tsv.gz`, the latent table), which are content-identical (compared decompressed).
    - The LaDDR config and command lines are identical after path normalisation and three design mappings: the latent table is moved to a stable name, so the assemble input path changes, and `--samples` reads `intermediate/samples.txt`.
    - The LaDDR `mv` to a stable name is now covered by these runs. See `tests/hpc/README.md`, Task 22, and `docs/module-status.md`.
- **Seeded tensorQTL (Task 23, ruling R30, `module-v0.1.6`).**
  - **Ruling R30** (asked for by the user, 2026-11-05): `datasets.<ds>.qtl_seed`, default 0, seeds `tensorqtl_cis` and `tensorqtl_cis_independent` through `run_tensorqtl.py --seed`. This is the one deliberate change of a default and of a script against upstream, which is unseeded; `qtl_seed: null` gives upstream's command lines exactly. The comparison with upstream stays attribution (R22), since upstream's outputs vary run to run and seeding the module does not tighten that bar. A seed gives one fixed draw of the permutations. All runs with the same seed and sample count use the same permutation indices, including across modalities, datasets and cis/cis_independent, so their Monte Carlo errors are not independent. This does not affect per-modality FDR; if an analysis needs independent draws between datasets, set different seeds. Choose the seed before looking at results; changing it shows the Monte Carlo spread.
  - `--seed` does two things, because tensorQTL has two random sources: `seed=` to `map_cis`/`map_independent` (`np.random.seed`, the permutations) and `torch.manual_seed` (the random tie-break, `torch.randint` at `cis.py:66`, which draws from torch's CPU generator and which tensorQTL's own `seed` does not touch). `trans` has no randomness and gets no seed.
  - Verified on the gate-3 config: two runs give the same decompressed files on GPU (A100, 16 of 16 files) and on 1 CPU core (16 of 16); `.gz` files differ only in the gzip header's write time (13 of 16 on GPU, 16 of 16 on CPU). Seed 1 changes all 16 files. GPU and CPU results differ from each other, as expected. The `github()` host at bb0b690 gave the same decompressed files as the standalone GPU run, 16 of 16.
  - Reproducibility is shown for one GPU model and for 1 CPU core: other GPU models and more CPU threads are untested. See `tests/hpc/README.md`, Task 23, and `docs/module-status.md`.
- **Image pin check** (`.github/workflows/image-pins.yml`, `containers/check_pins.sh`, Task 21). Image builds are not byte-reproducible (the same `containers/tensorqtl/` tree built twice gave two digests), so instead of rebuilding, the check runs on every push and pull request (`contents: read`, no docker) and fails if (1) a recipe input's git tree hash differs from the one recorded in `containers/pins/<image>` for the commit that built the pinned image, (2) the step file (`pheast/steps/qtl.smk`, `twas.smk`) references another digest, or (3) the pinned manifest is not anonymously pullable (HTTP 200). Recipe inputs: tensorQTL `containers/tensorqtl/`; FUSION `containers/fusion/`, `pheast/scripts/fusion_twas/`, `pheast/envs/fusion.yaml`. The pins files sit outside the recipe directories, so editing them changes no recorded tree and starts no image build. The MAJIQ image is never published and is not checked. **Updating an image:** change the recipe and push to `module`; the image workflow builds and pushes a new digest (image-pins fails on that push, as the tree no longer matches; to avoid that, put the new `tree` lines, known before the push, in the recipe commit, so only the digest check fails until the update). Then put the new digest in the step file and update `containers/pins/<image>`: `digest` (from the build log), `run`, `commit` (the run's head commit) and each `tree` line (`git rev-parse <commit>:<path>`), and push.

## FUSION

- **R26 (Task 19, supersedes R24): FUSION runs in a digest-pinned container.** `twas_compute_weights_batch` and `twas_assemble_summary` have no `conda:`, only `container: docker://ghcr.io/twillis209/pantry-fusion@sha256:87335916…4009`.
  - The image is built from `containers/fusion/` by `.github/workflows/fusion-image.yml` (public; tag `sha-<commit>`), from the repository root. It holds the vendored FUSION of `pheast/scripts/fusion_twas/` (`FUSION.compute_weights.R`, `utils/plink_utils.R`, `utils/FUSION.profile_wgt.R`, Linux x86-64 `gcta_nr_robust`, sha256-checked at build time and made executable) under `/opt/fusion`, keeping the `utils/` layout that `compute_weights` sources relative to itself, and the env of `pheast/envs/fusion.yaml`: R 4.5.3, r-glmnet 5.0, r-optparse 1.8.2, r-matrix 1.7_6, plink 1.90, gemma 0.98.5, Python 3.12.16, the baseline env's versions. `check_image.py` gates the build on these versions, on gcta running (statically linked), on the R packages loading, and on plink and gemma running.
  - The rules pass `/opt/fusion/...` paths as params, not inputs (they exist only in the container). `pheast/scripts/fusion_twas/` stays the image's source of truth; changing it (or the env file) rebuilds the image, and the new digest must be pinned in `pheast/steps/twas.smk`.
  - `twas_geno`, `twas_pos_file` and `twas_compress_output` keep their conda envs.
- **Behaviour differences from `module-v0.1.1`:**
  - **TWAS needs apptainer** (or another container runtime Snakemake supports): run with `--sdm conda apptainer`. A conda-only run of the TWAS rules is no longer possible (without a container they would run the host's `python`/`Rscript` with `/opt/fusion` paths that don't exist).
  - No network access is needed at run time beyond the image pull (no gcta download).
  - Removed config keys `fusion_dir`, `gcta`, `gcta_url` and `gcta_sha256`. `normalise()` rejects them with an error at the top level, under `defaults` or under a dataset (the schema allows unknown top-level keys, so the check is in Python). Upstream `8da7787` had none of them.
  - The `install_fusion` rule and `<pantry_resources>/fusion/` are gone; the `pantry_resources` pathvar stays (unused) so hosts that set it still load.
  - The Snakemake source-cache bug with binary files (below) no longer affects Pantry: the module reads no binary file.
- The image is public; apptainer pulls it anonymously. As for the tensorQTL image, stale stored ghcr.io credentials make the pull fail with 403.
- History: `module-v0.1.0` installed FUSION into `<pantry_resources>/fusion/` with an `install_fusion` job (overrides `fusion_dir`, `gcta`); R24 (`module-rc1` to `module-v0.1.1`) downloaded gcta from a pinned URL for `github()` modules, which needed network access in that job.

## MAJIQ

- **R27 (Task 20, `module-v0.1.3`): MAJIQ in an image each user builds; not published.** MAJIQ's [Penn academic licence](https://majiq.biociphers.org/Academic_License.txt) forbids redistribution (5(d)), so unlike the tensorQTL and FUSION images no MAJIQ image is built or published by any workflow (no `.github/workflows` `paths:` filter matches `containers/majiq/**`; `*.sif` is gitignored). The repository ships only recipes, `containers/majiq/majiq.def` (apptainer) and `Dockerfile`, which fetch MAJIQ at build time from https://bitbucket.org/biociphers/majiq_academic, tag `v3.0.22`, and stop unless it is commit `d841baad`. MAJIQ is not patched.
  - Both recipes start from the same digest-pinned micromamba base as the other images, create the env from `containers/majiq/majiq3.yml` and run `containers/majiq/majiq3-install.sh`, the one set of install steps also used by the conda route (`bash majiq3-install.sh <workdir>` creates the named env `majiq3`). The script checks `majiq --version` (3.0.22) and imports `rna_majiq`, `rna_voila`. The image puts `/opt/conda/bin` on PATH and keeps MAJIQ's notices (licence Section 7): the licence files from the clone and a `NOTICE` with "Copyright © 2016 University of Pennsylvania; All Rights Reserved" and the creators Yoseph Barash, Jorge V. Garcia, Alejandro Barrera and Matthew Gazzara, in `/opt/majiq` and in labels. The recipe headers say that users must accept the licence first, that use is non-commercial or educational only, and that built images must not be redistributed.
  - `apptainer build majiq.sif majiq.def` works without root or a subuid entry (apptainer 1.5.2 falls back to a root-mapped namespace; Task 20). The Docker route (`docker build`, then `apptainer build majiq.sif docker-archive://…`) is documented but untested: it was not built in Task 20.
  - Config: optional top-level `majiq_container` (a local `.sif` path, best absolute, or a `docker://`/`oras://` URI), in the schema and lifted to the top level by the legacy shim like `majiq_env`. When set, `majiq_build`, `majiq_psi_coverage` and `majiq_psi` have `container: <majiq_container>` and no `conda:` (a rule with both would run conda inside the container, as for R20). Unset, they keep `conda: <majiq_env>` and no container, so the `majiq_env` route's jobs and commands are unchanged.
  - Run with `--sdm conda apptainer`. Snakemake runs `apptainer exec --home <workdir>`, so the real home directory is not mounted: `apptainer-args` must bind the data, the targets of symlinked BAMs (MAJIQ names `.sj` files after the resolved BAM) as the test profiles do. Snakemake binds its own source cache (`~/.cache/snakemake`, where `script()` paths and the bundled licence file are) into container jobs, so no profile bind is needed; the test profiles' bind was removed in `module-v0.1.4` (R28). (Correction to the Task 20 docs, which said 9.12-9.13 needed it: those releases already bind the default cache path `<user cache>/snakemake/source-cache`, which is the cache in use under default settings; from 9.14.0 Snakemake binds whichever cache path is in use.) The Docker route and the `docker://`/`oras://` URI route are untested.
  - Verified (Task 20): the Task 17 configs rerun with the image give intron_retention outputs byte-identical to upstream 8da7787's, with plain and with symlinked BAMs. See `tests/hpc/README.md`, Task 20.
  - Moved: `tests/hpc/envs/majiq3.yml`, `majiq3-install.sh`, `majiq3.lock.txt` → `containers/majiq/`. `majiq3-install.sh` gained `--prefix <env>` (install into an existing env, for the images) and `--notices <dir>`, and now also copies the notices and runs the version and import check.

## Tests

- `tests/unit/` (pytest, 174 tests): config shim/normalise/validation, the FASTQ-map loader against an oracle copy of upstream's (R5), samples, targets, phenotype summaries against upstream's reads, rejection of the removed FUSION keys, `majiq_container`, and the comparison tool.
- `tests/compare_outputs.py`: compares an output tree against a baseline.
  - Exact for text files and decompressed `.gz`.
  - With `--noise A B [C ...]`: BED-like files are compared numerically within 2× the max |diff| over all pairs of noise runs (R12; two runs = the single pair; more than two: R12 amendment, Task 22 review, ratified by the user on 2026-11-05); `cis_qtl` permutation columns and tie-breaks are judged against the noise pair; `cis_independent` is judged on summary statistics (R8, R9).
  - Writes a markdown report.
- `tests/configs/`: baseline (old-format) configs for both stages, multi-dataset configs, and dry-run-only configs (`dryrun_multi.yml` with a dotted dataset name; `dryrun_majiq_laddr.yml` for MAJIQ and LaDDR), `phenotyping_baseline_ir.yml` (upstream's config with intron_retention, for the MAJIQ run) and `phenotyping_baseline_laddr.yml` (the same with LaDDR on, for the LaDDR run).
- `tests/fixtures/`: `make_multi_inputs.py` (derives the single-end and 300-sample test datasets), `check_targets.py`, `make_scaling_config.py` (synthetic many-dataset configs for dry-run timing, Task 18).
- `tests/host/`: a host workflow that imports both stages by local path (`Snakefile`) or with `github()` at `module-v0.1.6` (`Snakefile.github`), plus its profile.
- `tests/profile/` (SLURM test profile), `tests/hpc/` (run records, baseline envs); `containers/majiq/` (MAJIQ recipes and the `majiq3` env, R27).

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
| No log directive defined | every rule (66 phenotyping, 17 pheast; 16 since `install_fusion` was removed) | Upstream has no `log:` anywhere. Adding logs means adding `2> {log}` to every rule body, which the no-rule-body-rewrite constraint excludes. Follow-up. |
| Mixed rules and functions in same snakefile | `phenotyping/Snakefile`, `align`, `alt_TSS_polyA`, `intron_retention`, `laddr`, `stability`; pheast `qtl`, `heritability`, `twas` | Per-dataset input/params functions live beside the rules they serve (spec §6). Moving them all to `setup.smk` would separate each function from its only user. |
| Migrate long run directives into scripts | `write_samples`, `make_majiq_experiments_tsv`, `laddr_manifest`; `plink_covariates` | `run:` exception (spec §6). The full set of `run:` rules is `write_samples`, `validate_annotation`, `make_majiq_experiments_tsv`, `laddr_manifest`, `plink_covariates`, `check_genotype_samples`, `phenotype_summary`; lint flags only the four long ones (`install_fusion`, the fifth, was removed in `module-v0.1.2`). They are stdlib/pandas only (Snakemake provides both), are local, and are small. `plink_covariates`, `make_majiq_experiments_tsv` and `laddr_manifest` are upstream's own `run:` blocks. |
| Do not access input and output files individually by index | `combine_heritability_chr` (`head -n 1 {input[0]}`) | Upstream's rule body, unchanged. |
| Absolute path "/transcript_id "";/d" | `extract_cdna` in `phenotyping/Snakefile` | False positive: a `sed` delete expression, not a path. Upstream body. |

`snakefmt` (2.0.3): `snakefmt --check --exclude 'run_twas\.smk' phenotyping pheast tests/host` passes on all 20 files, and `tests/host/Snakefile.github` is clean too. The one exclusion, `pheast/scripts/run_twas.smk`, is upstream's standalone TWAS-association example. The module doesn't load it, and it is left byte-identical to upstream so that `pheast/scripts` differs only by `run_tensorqtl.py --seed` (R30) (the tensorQTL attribution evidence in R21/R22 relies on that).

## Migrating a single-dataset config

1. **Nothing is required to keep running.** An old-format `config.yml` in `phenotyping/` or `pheast/` still works and becomes dataset `default` (and reference `default`). Replace the repo's `config.yml` with yours, or run from another directory with `-s …/Snakefile --directory <dir>` where `<dir>/config.yml` is yours. Don't pass an old-format `--configfile` while the repo's new-format `config.yml` is present: Snakemake merges the two, and this is now an error rather than a silent use of the example dataset.
2. **Outputs move**: `output/X` → `output/default/X`, `intermediate/X` → `intermediate/default/X`, references under `intermediate/reference/default/`. Pheast's `phenotype_dir: ../phenotyping/output` therefore becomes `../phenotyping/output/default`. An existing run's intermediates are not picked up at the new paths, so expect a full rerun (or move them yourself).
3. To convert to the new format (recommended; see `phenotyping/config.yml` and `pheast/config.yml`):
   - phenotyping: move `ref_genome`/`ref_anno` to `references: {<ref>: {genome, annotation}}`, and `modality_groups.RNA_editing.edit_sites_bed` (and `intron_retention.annotation_gff3`) to that reference. Put everything else under `datasets: {<name>: {reference: <ref>, …}}`. Rename `samples_file` → `samples`. Drop `laddr.project_dir`. Set `laddr.dataset_name: pantry` if you relied on upstream's template value.
   - pheast: put everything under `datasets: {<name>: {…}}` and rename `samples_file` → `samples`.
   - Both: replace `intermediate_dir` / `intermediate_ref_dir` with top-level `pathvars: {pantry_intermediate: …, pantry_reference: …}` (pheast: `pantry_intermediate`, `pantry_resources`).
   - Integer knobs must be unquoted integers.
4. **tensorQTL and TWAS**: run with `--sdm conda apptainer` to use the pinned containers (tensorqtl 1.0.10; FUSION). Without apptainer the three tensorQTL rules use your `python3`'s tensorqtl, as upstream did, but the TWAS rules need the FUSION container. Drop any `fusion_dir`, `gcta`, `gcta_url`, `gcta_sha256` keys (added by `module-v0.1.0`/`v0.1.1`). GPU profiles: see [Environments](#environments-containers-gpu).
5. Snakemake ≥ 9.14.

## Deviations from the handover contract and the spec

- **Contract 1 (one `config["pantry"]`)**: the stages are two modules, each with its own config section (`config["pantry"]["phenotyping"]`, `config["pantry"]["pheast"]`) and distinct prefixes (`pantry_pheno_*`, `pantry_pheast_*`). Both define `all`, so one shared prefix would clash.
- **Contract 8 (fetch a pinned FUSION commit and GCTA/GEMMA)**: FUSION is vendored with Pantry's modifications and shipped, with gcta and gemma, in a digest-pinned container built from the module's own files (R26; supersedes R24's install job and gcta download).
- **Contract 11 / criterion 11 (standalone)**: works, but standalone outputs gain the `default/` level.
- **Default behaviour differs from upstream for tensorQTL (R30, `module-v0.1.6`)**: cis and cis_independent mapping are seeded by default (`qtl_seed: 0`); upstream is unseeded. `qtl_seed: null` restores upstream's command lines. Seeded outputs are still judged by R22 (identical inputs and command lines apart from `--seed`), since upstream's outputs vary run to run.
- **Acceptance criterion 1 for tensorQTL**: tensorQTL is unseeded and nondeterministic upstream (and the module's default seed changes neither that nor the baseline), so `*.cis_qtl`/`*.cis_independent_qtl` were judged by attribution rather than equality (R21, R22): identical inputs and command lines, comparison reported but not blocking. kallisto outputs were judged within the run-to-run noise envelope (R11, R12). See `docs/module-status.md`.
- **LaDDR** and MAJIQ (`intron_retention`) were first checked by dry runs (design decision); both are now verified at run time (Tasks 17 and 22). The LaDDR env deviates from the baseline lock in scikit-learn (1.5.2) and scipy (1.16.3), since LaDDR v1.0.0 cannot run with the lock's versions (R29).
- Spec Corrections 5: `index_fasta` removed, not kept (R15).
- Spec §2 constraints: some wildcards are narrower enumerations (R4).
- Spec §6: LaDDR `dataset_name` defaults to the dataset key, so legacy configs that omit it get `default`, not upstream's `pantry` (R16). This affects only LaDDR's internal manifest names.
- Spec §7: FUSION and gcta come from a digest-pinned container (R26), superseding the job-time gcta download for remote modules (R24; spec amended).
- Validation is stricter on quoted or non-integral numbers (R14, R18).

## Known limitations and follow-ups

- **LaDDR upstream issues** (to report; Task 22): v1.0.0 fails in `laddr coverage` under scikit-learn 1.9.1 (duplicate column names in `remove_pcs_from_phenos`) and is nondeterministic (unseeded randomized PCA there). The module's outputs equal upstream's only with the RNG seeded; unseeded, two runs of either differ.
- **MAJIQ verified at run time** (Task 17): MAJIQ 3.0.22 (public academic build, which does not check the licence file) in a named env `majiq3` (then `tests/hpc/envs/majiq3.yml`, `majiq3-install.sh`, lock; now in `containers/majiq/`), and again in a user-built image (Task 20, R27). The intron_retention outputs of the module are byte-identical to upstream 8da7787's, including the nested `majiq_build` outputs and the `.sj` derivation; with symlinked, renamed BAMs the `.sj` files follow the targets and outputs are unchanged. `majiq.psi.tsv` bootstrap columns vary run to run; the extracted `raw_psi_mean` does not.
- **Env pinning**: LaDDR's pip-only dependencies (bx-python, scikit-fda, statsmodels, tqdm, numba and theirs) other than `multimethod` resolve at build time (no lock covers them). MAJIQ's named env and image are reproducible from the pinned `containers/majiq/` recipes, but the rules name the env or take the image path from config rather than pin them, and MAJIQ's build-isolation packages (pip's build env) are not pinned.
- **Snakemake source-cache bugs** (9.27.x, `sourcecache.py` `HostingProviderFile.open` / `LocalGitFile._open`; to report upstream, not filed). `git show` output is decoded as text, so (1) `source_path` on a binary file in a `github()` module fails with `UnicodeEncodeError` (worked around by R24 up to `module-v0.1.1`; since R26 the module reads no binary file), and (2) text files lose their final newline: FUSION R scripts installed from `github()` were 1 byte shorter (up to `module-v0.1.1`), and env files hash differently from the local-path ones, so a `github()` module builds its own conda envs (13 envs, about 15 min in Task 18). Outputs are unaffected. See `tests/hpc/README.md`, Task 14.
- **Nondeterminism**: kallisto 0.48.0 quant varies run to run (node/load dependent; not thread count). tensorQTL permutations and tie-breaks are unseeded upstream; the module seeds them by default (R30), so its reruns are reproducible on one device type, but the comparison with upstream's unseeded runs stays attribution (R22). The noise envelopes came from a single pair of baseline runs, and the isoforms kallisto margin in Task 13 was close (max |diff| 0.149 vs allowed 0.167). **R12 amendment** (Task 22 review; methodology ruling, ratified by the user on 2026-11-05): phenotyping envelopes are now 2× the max over all pairs of the four upstream runs {archive/phenotyping, archive/phenotyping_run2, U1, U2} (job 8407236); every phenotyping candidate (U1, U2, M1, M2, Task 21's v0.1.4 run; later Task 23's run, job 8433373) passes, where U1, U2, M1 and M2 failed the single-pair envelope on kallisto files. Normalised noise is quantised in rank steps (0.148 to 0.378 at 16 samples); `alt_TSS.bed.gz`'s amended envelope (0.216) is still below one extreme-rank step. tensorQTL envelopes are still from one pair. See `docs/module-status.md`, "R12 amendment".
- **Legacy pathvars vs host pathvars**: a legacy `intermediate_dir`/`intermediate_ref_dir` in a module's config overrides the host's module pathvars, because it registers later. Hosts should pass new-format configs only.
- `--configfile` with a new-format config inside `phenotyping/` or `pheast/` is merged with the repo's `config.yml`, which adds its `default` dataset. Use `--directory` to point elsewhere, or edit `config.yml`. This merging is standard Snakemake behaviour.
- **Rerun triggers**: script paths are derived params, kept out of the params trigger because the source-cache path changes per process. As upstream (scripts called from shell), editing a bundled script doesn't trigger a rerun; use `-R <rule>`. The configured MAJIQ `license` is a derived param too, so changing it doesn't trigger a rerun. Wildcard-only lambda params remain in the params trigger (stable values). Under `rerun-triggers: mtime` (as in `tests/profile`), a changed `qtl_seed` reruns nothing, and a run directory from before `module-v0.1.6` keeps its unseeded tensorQTL results; force them with `-R tensorqtl_cis tensorqtl_cis_independent` (in a host, the prefixed rule names). Changing the FUSION image digest doesn't trigger a TWAS rerun either (Snakemake has no container rerun trigger), nor does switching MAJIQ between `majiq_env` and `majiq_container` (the profiles use `rerun-triggers: mtime`). The FUSION paths (`/opt/fusion/...`) are plain params, so TWAS outputs from v0.1.1 or earlier (which took them as inputs from `<pantry_resources>/fusion`) rerun once under v0.1.2.
- `pheast/scripts/covariates.R` sources `fusion_twas/utils/plink_utils.R` relative to itself (upstream). The `covariates` rule's otherwise unused `plink_utils` param keeps that file beside it in Snakemake's source cache (before Task 19, `install_fusion`'s params did); keep both if either moves.
- `combine_modalities` keeps the script's behaviour: `set +o pipefail` (the script runs without pipefail) lets a truncated `.bed.gz` through (a `gzip -t` check would catch it), and the `sort` locale is unpinned. Left as upstream to keep outputs identical to the script's.
- LaDDR's `use_existing_bins` route (`laddr_prepare_existing_bins`) was not run in Task 22; `tests/configs/dryrun_majiq_laddr.yml` covers it in a dry run.
- LaDDR: if the final `mv` in `laddr_transform` fails, a stale latent file from an earlier run can remain in the LaDDR project dir.
- Test tooling: the comparison tool's floors (1 row; 1e-6 relative) tolerate a single variant-changed key or NaN mismatch in `cis_independent`, and RDat pairs after an OSError are not compared. The fixture/checker scripts are not unit-tested; `.tbi` files are checked by size only. The setup.smk injection of `annotation_gff3` is covered by dry run only. Intermediate commits `ad85347..c5c3ab2` don't dry-run on their own (bisect across them with care).
- **`defaults` cannot remove keys**: a dataset can add to or override `defaults` but not drop a modality group, modality or analysis from it. Keep varying modality sets out of `defaults` (see [Config format](#config-format)). Removal support is deferred.
- **Stage linking is by path text**: pheast `phenotype_dir` must spell the phenotyping `<pantry_output>/{dataset}` path exactly as the host resolves it (relative vs absolute matters); otherwise pheast reports missing input files.
- **Hosts must pass `config:`** in each `module` directive; otherwise the module's standalone `configfile: config.yml` can load a `config.yml` from the host's working directory.
- pheast `samples` still doesn't subset the analysis (upstream semantics): to analyse fewer samples, supply BEDs with those samples.
- **Scaling** (dry runs, Task 18): about 3 ms and ~20 KB per job, roughly linear; 17 datasets x 1000 samples (205,653 jobs) dry-run in 12m20 with 4.4 GB, against 33 s for upstream's 1 x 1000 (38 s here). Snakemake's own job handling dominates, not Pantry's helpers. See `tests/hpc/README.md`, Task 18.
- Ideas not taken: the standard `workflow/` layout; wrapping `prepare_phenotypes.py` in a rule.
