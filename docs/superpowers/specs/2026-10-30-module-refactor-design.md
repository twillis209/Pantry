# Pantry as a reusable, multi-dataset Snakemake module — design

Date: 2026-10-30
Branch: `module` on `twillis209/Pantry` (base: upstream `8da7787`)
Source of requirements: `~/projects/scPantry/docs/handover/2026-10-30-pantry-fork-handover.md` (the "handover"). This document does not repeat the handover's contract; it says how the fork meets it.

## Intent

Pantry's `phenotyping/` and `pheast/` are standalone projects run once per dataset from inside their directory. scPantry needs to import each stage **once** as a Snakemake `module` and run ~17 datasets through one DAG, with host-controlled output locations and shared references. The fork must stay upstreamable: same analysis, same defaults, same output formats, standalone use still works.

## Goals, non-goals, anti-goals

Goals: handover contract items 1–11 and acceptance criteria 1–6.

Non-goals (from the handover; the one later exception is tensorQTL's default seed, `qtl_seed: 0`, ruling R30, asked for by the user on 2026-11-05): changing analysis logic, defaults or output formats; TWAS association/coloc; anything scPantry-specific. Also out of scope here: wrapping `pheast/scripts/prepare_phenotypes.py` (sample-ID remapping) in a rule — a host that needs it runs it itself; moving to the standard `workflow/` layout (would rewrite every path and make the upstream diff unreadable).

Anti-goals (a technically working result that would still be a failure):
- Any output differing from upstream `8da7787` on the test data beyond documented nondeterminism.
- Rule bodies rewritten "while we're there". Shell/script logic changes only where a path or parameter must come from somewhere new.
- Swallowing errors (cf. the previous author's `run_tensorqtl.py`).
- A refactor that only works via `github()` or only locally; both must work.

## Corrections to the handover (verified at `8da7787`)

1. **FUSION is vendored**: `pheast/scripts/fusion_twas/` contains `FUSION.compute_weights.R`, `utils/plink_utils.R`, `utils/FUSION.profile_wgt.R` and a statically linked **Linux x86-64** `gcta_nr_robust`. `gemma` comes from the conda env. So "fetch a pinned FUSION commit" becomes "install the vendored FUSION into a resources directory" (below). Consequence: pheast TWAS cannot run on macOS; tests run on Linux/HPC.
2. **`localrules:` breaks under rule renaming.** The module's `localrules: index_bed, …` registers unprefixed names, which no longer exist after `use rule * … as pantry_*`. All `localrules:` blocks become per-rule `localrule: True`.
3. **Cross-modality**: upstream has only the manual `phenotyping/scripts/combine_modalities.sh`; no rule. Confirmed.
4. **tensorQTL is nondeterministic**: `run_tensorqtl.py` uses permutations and `random_tiebreak=True` with no seed. Permutation-derived columns cannot be byte-compared (see Verification).
5. `index_fasta` is not used by any rule in `phenotyping/steps/`; kept (scoped) for compatibility.

## Snakemake facts the design relies on (checked in Snakemake source, v9.27)

- Pathvars (`<name>`) resolve in input/output/log/benchmark paths, **including paths returned by input functions** (`IOFile` applies `rule.pathvars`), but **not in `params`** or shell strings. Pathvar values may reference other pathvars. Precedence: rule > module config `pathvars:` > `module … pathvars:` > global config > Snakefile `pathvars:`.
- Checkpoints inside a module get their own namespace (`checkpoints.spawn_new_namespace()`), registered under the original rule name, so `checkpoints.<orig_name>.get(...)` works inside a prefixed module. **Still to be proved by the host test** (highest-risk assumption; see Risks).
- Pathvars need Snakemake ≥ 9.12: both Snakefiles get `min_version("9.12")` (upstream README says v8; README updated).

## Design

### 1. Module shape and host usage

Two independent modules, unchanged entry points (`phenotyping/Snakefile`, `pheast/Snakefile`). Host:

```python
module pantry_pheno:
    snakefile: github("twillis209/Pantry", path="phenotyping/Snakefile", tag="<tag>")
    config: config["pantry"]["phenotyping"]
    pathvars:
        pantry_output="<results>/pantry/phenotypes",
        pantry_intermediate="<results>/pantry/intermediate/phenotyping",
        pantry_reference="<resources>/pantry/reference",
use rule * from pantry_pheno as pantry_pheno_*

module pantry_pheast:
    snakefile: github("twillis209/Pantry", path="pheast/Snakefile", tag="<tag>")
    config: config["pantry"]["pheast"]
    pathvars:
        pantry_output="<results>/pantry/pheast",
        pantry_intermediate="<results>/pantry/intermediate/pheast",
        pantry_resources="<resources>/pantry",
use rule * from pantry_pheast as pantry_pheast_*
```

Distinct prefixes are required (both stages define `all`). Neither Snakefile has `configfile:` any more in module mode; standalone mode is handled in §8.

### 2. Path layout and pathvars

Prefixed pathvar names (`pantry_*`) rather than stock `results`/`resources`, so module defaults can reproduce the standalone layout and cannot be silently captured by a host's global `results`. Snakefile-level defaults = standalone layout:

| Pathvar | Stage | Standalone default | Contents |
|---|---|---|---|
| `pantry_output` | both | `output` | phenotyping: `{dataset}/{modality}.bed.gz` etc.; pheast: `{dataset}/{analysis}/…` |
| `pantry_intermediate` | both | `intermediate` | `{dataset}/…` per-dataset intermediates (STAR, kallisto, covar, TWAS, LaDDR…) |
| `pantry_reference` | phenotyping | `intermediate/reference` | `{reference}/…` shared references |
| `pantry_resources` | pheast | `intermediate/resources` | `fusion/` installed FUSION |

Every rule path is `<pathvar>/{dataset}/…` or `<pantry_reference>/{reference}/…`. **Deviation (criterion 11):** standalone outputs gain a dataset level (`output/default/expression.bed.gz`).

Because pathvars don't reach `params`, directory params are derived from resolved outputs/inputs, e.g. `params: outdir=lambda w, output: Path(output[0]).parent` (named helpers, not inline lambdas, where repeated). Bare `{output_dir}`/`{majiq_psicov_dir}` in shell bodies are replaced the same way.

`wildcard_constraints` (global in each Snakefile): `dataset`, `reference`, `sample_id`, `modality`, `pheno`, `type`, `group`, `position`, `chrom`, `batch` = `[^/]+`; `read_length` = `\d+`. Generic rules are scoped: `index_bam` only matches `<pantry_intermediate>/{dataset}/bam/{sample_id}.bam`; `index_bed` only `<pantry_output>/{dataset}/{pheno}.bed.gz`; `index_fasta` only under `<pantry_reference>`.

### 3. Config shape (phenotyping)

```yaml
references:
  grch38_e106:
    genome: /path/genome.fa            # uncompressed
    annotation: /path/annotation.gtf   # uncompressed, must contain transcript rows
    edit_sites_bed: /path/sites.bed    # optional; needed if any dataset uses RNA_editing
    annotation_gff3: /path/x.gff3      # optional, MAJIQ; else converted from GTF
defaults:            # optional; deep-merged under every dataset
  modality_groups: {...}
datasets:
  onek1k_all:
    reference: grch38_e106
    samples: [s1, s2]                  # list, or path to one-column file
    fastq_map: {s1: [[r1.fq.gz], [r2.fq.gz]]}  # dict, or path to existing map file format
    fastq_dir: /path                   # used with the file form
    paired_end: true                   # optional, as today
    read_length: 75
    fragment_length_mean: 200          # single-end only
    fragment_length_sd: 50
    modality_groups:                   # as today, per dataset; RNA_editing knobs unchanged
      expression: {files: [...]}
      RNA_editing: {edit_sites_min_coverage: 20, ...,  files: [...]}
    cross_modality:                    # optional, §6
      modalities: [alt_polyA, alt_TSS, expression, isoforms, RNA_editing, splicing, stability]
    laddr: {enabled: false, ...}       # as today, per dataset
```

`edit_sites_bed` moves from `modality_groups.RNA_editing` to the reference set (it is reference data, and `edit_sites_to_genes.tsv` is keyed by reference). The shim (§8) moves it automatically.

Dict form of `fastq_map`: paired = `[R1 list, R2 list]`, single = `[list]`; paths absolute or relative to `fastq_dir`.

### 4. Config shape (pheast)

```yaml
defaults: {...}                     # optional, deep-merged
datasets:
  onek1k_all:
    phenotype_dir: results/pantry/phenotypes/onek1k_all   # host points this at phenotyping output
    geno_prefix: /path/plink_prefix
    samples: [...] or path
    modalities: {expression: {grouped: false}, cross_modality: {grouped: true}, ...}
    analyses: {qtl: {files: [...]}, heritability: {...}, twas: {...}}
    qtl_min_allele_count: 10
    qtl_maf_threshold: null
    qtl_cis_window: 1000000         # new knob → run_tensorqtl.py --window (its existing default)
    twas_snps: /path or absent
    twas_models_per_job: 200
    twas_threads: 4
fusion_dir: null                    # optional override, §7
gcta: null                          # optional override of gcta binary
```

pheast `samples` keeps its upstream meaning: it does **not** subset the analysis (scripts use the samples shared by BED and genotypes); it is used only for the genotype check and `N` in the TWAS `.pos` file. A dataset wanting fewer samples supplies BEDs with those samples. (My call: changing this would be an analysis change.)

`phenotype_dir` is consumed through input functions, so if it equals phenotyping's output path the host DAG links the two stages automatically.

### 5. Validation and parse-time reads

- JSON schemas `phenotyping/schemas/config.schema.yaml`, `pheast/schemas/config.schema.yaml`, applied with `snakemake.utils.validate()` **after** the shim. Defaults currently set via `setdefault` (LaDDR, edit-site knobs, qtl/twas knobs) move into the schemas, same values. Cross-field rules that schemas can't express (cross_modality only with `qtl`; `laddr.regress_modalities` ⊆ available modalities; fragment lengths for single-end) stay as Python checks on config only.
- **No file reads at load.** Sample lists and FASTQ maps given as paths are read lazily on first use by an input function (cached per dataset), i.e. at DAG-build time; the inline forms avoid files. Mixed/inconsistent endedness checks move into that loader.
- Reference checks become a job: `validate_annotation` per `{reference}` writes `<pantry_reference>/{reference}/annotation.ok` after the existing transcript-row check; reference-building rules take it as input. Existence of genome/edit-sites files is enforced by them being rule inputs.
- Genotype-vs-sample check becomes `check_genotype_samples` per `{dataset}` (reads `.fam`, fails with the same message), an input of `prune_for_covar` and `qtl_genotype_qc`.
- Things only knowable from phenotype BEDs at job time (pheast `chroms` per modality for heritability; model count for TWAS batches) come from a checkpoint `phenotype_summary` per `{dataset},{modality}` writing chromosomes and row count; `combine_heritability_chr` and `twas_assemble_summary` inputs use `checkpoints.phenotype_summary.get(...)`.
- Each dataset materialises `<pantry_intermediate>/{dataset}/samples.txt` (rule from config) because scripts take `--samples FILE`. Contents identical to the current samples file.

### 6. Rules: static rule set, per-dataset behaviour

- All `steps/*.smk` are always included; which outputs a dataset gets is decided by its `modality_groups`/`analyses` when building `rule all` targets. Globals become per-dataset lookup helpers in `setup.smk` (`ds(w)`, `samples(w)`, `fastq(w)`, `ref(w)`), used via `lookup()` or named input functions; `expand(..., sample_id=samples)` becomes input functions over the dataset's samples.
- Parse-time `if` rule pairs become static rules with input-function branching: MAJIQ GFF3 (configured file vs converted), LaDDR existing-bins vs built bins (distinct outputs `bins_ready.existing` / `bins_ready.built`, consumer picks one).
- `majiq_build`: its output list is currently computed at load from resolved BAM paths. It becomes `directory(<build dir>)`; `majiq_psi_coverage` derives the per-sample `.sj` paths in params (resolving BAM symlinks at job time, preserving fix `21ec32d`). Same commands.
- LaDDR: project dir `<pantry_intermediate>/{dataset}/laddr`; `dataset_name` defaults to the dataset key. Outputs whose names embed `dataset_name` get a stable path (the rule moves LaDDR's file to it) — my call; LaDDR internals unchanged.
- **Cross-modality rule** `combine_modalities` per `{dataset}`: inputs = the listed modalities' normalised BEDs; outputs `<pantry_output>/{dataset}/cross_modality.bed.gz`, `.tbi`, `.phenotype_groups.txt`. Body reproduces `combine_modalities.sh` exactly (header check, `modality:` prefix, `sort -k1,1 -k2,2n`, bgzip, `tabix -p bed`, groups file); the script itself stays, rewritten to call the same logic or left as-is for manual use (builder's call).
- `localrules:` blocks → `localrule: True` per rule.
- `run:` blocks (`make_majiq_experiments_tsv`, `laddr_manifest`, `plink_covariates`) stay as `run:`; they use only stdlib/pandas, which Snakemake itself provides. Justified lint exception.
- Threads/resources knobs that exist today (`twas_threads`, MAJIQ `threads`, LaDDR `bam_coverage_threads`) are kept for compatibility even though they'd normally be profile settings.

### 7. References, scripts, FUSION, environments

- References keyed by **named reference set**: `<pantry_reference>/{reference}/{cDNA.fa.gz, kallisto.idx, star_index_{read_length}/, txrevise/…, exonic.gtf, intronic.gtf, majiq/annotation.gff3, edit_sites_to_genes.tsv}`. STAR also keyed by `read_length`. Datasets sharing a reference name share all of these; each file has exactly one producing rule. Choosing names (not content hashes) is my call: it is explicit and the host already names its references.
- Scripts are referenced via `workflow.source_path("../scripts/…")` in `params`, so they work from local paths and `github()`. `pheast/Snakefile`'s `from scripts.twas_batches import …` becomes loading that file via `source_path` (no CWD import). `twas_assemble_summary`'s `realpath scripts` becomes the installed FUSION path.
- **FUSION**: rule `install_fusion` copies `FUSION.compute_weights.R`, `utils/plink_utils.R`, `utils/FUSION.profile_wgt.R` and `gcta_nr_robust` (via `source_path`) into `<pantry_resources>/fusion/` preserving the `utils/` layout (compute_weights sources it relative to itself) and `chmod +x` the binary (the source cache doesn't preserve the exec bit). Config `fusion_dir` / `gcta` override it. No network fetch of upstream FUSION: the vendored copy carries Pantry's modifications.
  - *Amended 2026-11-02 (ruling R24, commit `a140598`):* Snakemake's git source cache cannot serve binary files (`github()` modules fail at parse time on `source_path` of `gcta_nr_robust`), so gcta is resolved only inside `install_fusion` at job time: a configured `gcta` → nothing installed for it; module loaded from a local path → bundled file copied; otherwise (remote module) downloaded from config `gcta_url` (default: upstream `8da7787`'s file on `raw.githubusercontent.com`) and checked against `gcta_sha256` (default `7b51e51d…93e9`), then `chmod +x`. Remote modules therefore need network access in that job; offline sites set `gcta`. The FUSION R scripts are still copied from the module source (no network).
  - *Amended 2026-11-03 (ruling R26, `module-v0.1.2`; supersedes R24):* `install_fusion` and the `fusion_dir`/`gcta`/`gcta_url`/`gcta_sha256` keys are removed. The FUSION files of `pheast/scripts/fusion_twas/`, gcta and the `fusion` env are shipped in a digest-pinned image (`containers/fusion/`, `ghcr.io/twillis209/pantry-fusion`) that `twas_compute_weights_batch` and `twas_assemble_summary` use via `container:`, with `/opt/fusion` paths as params; TWAS therefore needs apptainer.
- **Envs**: `phenotyping/envs/` and `pheast/envs/`, one per tool group, e.g. phenotyping `star`, `samtools` (samtools+htslib), `kallisto` (0.48.0 as today), `gffread`, `regtools`, `subread`, `python` (pandas, gtfparse, scikit-learn, pysam as needed by scripts), `txrevise` (R txrevise/tidyverse/optparse + python), `rna_editing`, `deeptools`, `laddr`; pheast `plink`, `covariates` (R impute), `tensorqtl` (pip tensorqtl, torch 2.3.1, dask<2024.2, rpy2, qvalue), `gcta`, `fusion` (R glmnet/optparse, gemma 0.98.5, plink). MAJIQ keeps a configurable named env (`majiq_env`, default `majiq3`). Constraints in `environment.yml` carried over; full pinning later. The monolithic `environment.yml` files stay for users not using `--sdm conda`.

### 8. Standalone mode and config shim

`phenotyping/Snakefile` and `pheast/Snakefile` load `config.yml` only when run standalone and no config was supplied (`if not config and Path(workflow.basedir, 'config.yml').exists()` is not enough under `github()`; standalone is detected by the absence of module context — builder to pick the mechanism, e.g. `configfile` guarded by `workflow.modifier` state or a documented `--configfile config.yml`). If config lacks `datasets`, the shim maps it to `references: {default: {genome, annotation, edit_sites_bed}}`, `datasets: {default: {...rest}}`; old `intermediate_dir`/`intermediate_ref_dir`/`phenotype_dir` keys are honoured (the first two by setting `config['pathvars']`, which is how Snakemake reads config pathvars). `config.yml` files in the repo are converted to the new format, with the old format still accepted.

## Verification plan (acceptance criteria)

Run on Linux/HPC (required by `gcta_nr_robust`; tools via conda). Decided with the user:
- The builder submits jobs on the cluster reachable as SSH alias `hpc` (GDML HPC, x86-64, miniforge, Snakemake 9.26). Verification can gate each behaviour-touching commit. Cluster conventions follow `~/projects/genex/CLAUDE.md` "HPC Setup" (detail in `~/projects/genex/docs/context/operations.md`):
  - Non-interactive SSH: `source /etc/profile; module load slurm; source ~/.dotfiles/.bash_functions` first.
  - **Every** Snakemake invocation, controllers and dry runs alike, runs as a SLURM job via `sjob` (e.g. `sjob -c 2 -m 7500M -t 0-12:00:00 -n <name> -e snakemake snakemake --profile <profile> …`; dry runs as a script under `logs/` run by `sjob`). The login node is for `squeue`/`sacct`, `mv` and light inspection only.
  - `sjob -m` is memory **per CPU**; size as 3 × peak ÷ cores and check `sacct -o ReqMem,MaxRSS` after new runners.
  - A SLURM profile for the tests (kept in `tests/profile/`, modelled on `genex/smk/profile/config.yaml`: `executor: slurm`, account `gdml`, partition `cpu-core`, `software-deployment-method: conda`, `rerun-triggers: mtime`) and always passed with `--profile`.
  - Bulk file transfer with `tar` over one `ssh` connection, not per-file `scp`.
  - Working copies on HPC live under `/data/hps/assoc/private/gdml/user/twil15/` (e.g. `…/pantry-baseline` at `8da7787` and `…/Pantry` on `module`).
- **MAJIQ and LaDDR are deferred.** Both are refactored to the new structure but checked by `snakemake -n` only; `intron_retention` is excluded from the baseline and comparison test configs. Their run-time verification is follow-up work, listed in the status note.
- The builder may push `module` and tags to `origin` (`twillis209/Pantry`, never `upstream`), so the host test also runs once via `github(...)`.

1. **Baseline first**: run upstream `8da7787` (unmodified checkout) on both test datasets, archive outputs, before any behaviour-touching commit. Run pheast's QTL twice to measure run-to-run variation.
2. **Comparison**: phenotyping BEDs, phenotype-groups files, RNA-editing site map: identical after decompression. Heritability and TWAS (`set.seed(1)`, `PLINK_LASSO_SEED`): identical, or equal to stated tolerance if thread scheduling changes ordering. tensorQTL: deterministic columns (phenotype, variant, nominal p, slope, se, af…) identical except documented random tie-breaks; permutation columns (`pval_perm`, `pval_beta`, …) within the run-to-run range measured in step 1. A `tests/compare_outputs.py` does this and writes a report.
3. **Multi-dataset**: phenotyping datasets `full` (16 paired samples; must match baseline) and `half_se` (8 samples, R1 only as single-end, same `read_length`) sharing one reference; pheast datasets `full` (445; must match baseline) and `subset` (300 samples, different `qtl_min_allele_count`), whose BEDs and sample list are subset by a test-fixture script in `tests/` (not a module rule).
4. **Host test**: `tests/host/` with a Snakefile importing both stages by local path, custom pathvars, `cross_modality` enabled for phenotyping, and one pheast dataset whose `phenotype_dir` points at phenotyping output (`snakemake -n` only — 16 samples is too few for pheast's 20 phenotype PCs) plus the pheast test-data dataset (full run, proves checkpoints inside a prefixed module).
5. `snakemake --lint` (justified exceptions listed), `snakefmt` on changed files.
6. Standalone: unmodified-format `config.yml` in each stage runs and matches baseline.

## Sequencing (for the plan)

1. Produce baseline (§Verification 1).
2. Pure `snakefmt` commit, no behaviour change.
3. Small single-purpose commits in roughly this order: localrules→localrule; script paths via `source_path`; pathvars for current single-dataset paths; setup → lookup helpers; `{dataset}` wildcard + config shape + shim; references keyed by `{reference}`; parse-time reads → rules/checkpoints; MAJIQ/LaDDR static rules; cross-modality rule; FUSION install; per-rule envs; schemas; tests; `CHANGES-module.md` and status note.
4. Spike early (right after step 2): a toy host importing a prefixed module with a checkpoint, to retire the main risk before building on it.

## Risks

- **Checkpoints in prefixed modules** (needed for pheast chroms/TWAS batches when BEDs are produced in the same DAG). Source suggests it works; spike proves it. Fallback: host builds phenotypes and pheast in separate invocations, with phenotype-derived values read at DAG time.
- `github()` source caching of a 9.5 MB binary and of `directory`-style script trees; covered by `install_fusion` and the host test using `github()` once the branch is pushed (local path before).
- Upstream may not want the `{dataset}` level in standalone output paths.

## Left to the builder

Helper names and file organisation within `setup.smk`; exact env split boundaries; the standalone-detection mechanism; whether `combine_modalities.sh` stays or wraps the rule; tolerance values (documented in the comparison report).

## My calls (not discussed with the user)

Prefixed `pantry_*` pathvars; named reference sets; `edit_sites_bed` moved to references; `defaults:` deep-merge block; materialised per-dataset `samples.txt`; `run:` blocks kept; vendored FUSION installed rather than fetched; new `qtl_cis_window` knob; standard `workflow/` layout not adopted; LaDDR stable output names.
