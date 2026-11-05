# Pantry module refactor: status note

Date: 2026-11-04. This is the report back requested by the scPantry handover (`~/projects/scPantry/docs/handover/2026-10-30-pantry-fork-handover.md`, "Report back").
Changes in full: `CHANGES-module.md`. Design: `docs/superpowers/specs/2026-10-30-module-refactor-design.md`. Run records: `tests/hpc/README.md`.

## Branch and tag

- Fork `twillis209/Pantry`, branch **`module`**, release tag **`module-v0.1.5`** (earlier `module-v0.1.4`, `module-v0.1.3`, `module-v0.1.2`, `module-v0.1.1`, `module-v0.1.0`, `module-rc1`). Base: upstream `PejLab/Pantry` `8da7787`. Nothing was pushed to upstream.
- What ran at which commit:
  - gates 1–3 and the multi-dataset runs ran at earlier commits of the branch (see `tests/hpc/README.md`);
  - the full local-path host run (job 8396531) ran at `56c8c22`;
  - the R24 re-verification (job 8398071: local-path host dry run, `install_fusion`, and the `github()` host dry run plus two jobs) ran at `a140598` = `module-rc1`;
  - R24 changed only how `install_fusion` resolves gcta;
  - the tensorqtl 1.0.10 image (Task 16) came after `module-v0.1.0`: gate 3 GPU (job 8399104) and the CPU checks ran at `f709952`.
- `module-v0.1.0` adds no-rerun changes on top of `module-rc1`, checked locally by unit tests and dry runs. None of them can change any output:
  - a 300 s download timeout for gcta;
  - two validations: legacy top-level keys beside `datasets:` are rejected, and duplicate `cross_modality` modalities are rejected;
  - a comment, the wording of the LaDDR `project_dir` note, and docs.
- `module-v0.1.1` adds, on top of `module-v0.1.0` and the Task 16/17 commits (tensorQTL image, MAJIQ env), changes that cannot alter outputs:
  - a clear error naming the dataset when single-end samples lack `fragment_length_mean`/`_sd` and `paired_end` is omitted;
  - the RNA editing threshold note is printed once per dataset;
  - the default `gcta_url` names upstream's full commit SHA (same file, same sha256);
  - the tensorQTL image workflow builds only on pushes to `module`; `majiq3-install.sh` checks the tag's commit;
  - docs.
- The full `github()` host run (Task 18, job 8399793) fetched `module-v0.1.1` = `befa72c`.
- `module-v0.1.2` (Task 19, R26, superseding R24) moves FUSION (R scripts, gcta_nr_robust, R, plink, gemma) into the digest-pinned image `ghcr.io/twillis209/pantry-fusion`, removes `install_fusion` and the `fusion_dir`/`gcta`/`gcta_url`/`gcta_sha256` keys, and was re-run on the gate-3 config (below). TWAS now needs apptainer.
- `module-v0.1.3` (Task 20, R27) adds the optional `majiq_container`: MAJIQ in an image each user builds from `containers/majiq/` (not published: Penn's academic licence forbids redistribution). With it unset, the MAJIQ rules and every other job dry-run unchanged. It was re-run on the Task 17 MAJIQ configs (below).
- `module-v0.1.4` (Task 21, R28) raises the floor to Snakemake 9.14, pins the txrevise and LaDDR envs (LaDDR v1.0.0) and the MAJIQ recipes to the lock files, removes the test profiles' source-cache bind, adds the `image-pins` CI check (`containers/check_pins.sh`) and headers on the legacy `environment.yml` files. Dry-run job sets are unchanged. A full `github()` host run with `majiq_container` set and no cache bind passed the Task 18/19 bar (below).
- `module-v0.1.5` (Task 22, R29) verifies LaDDR at run time against upstream (below). It fully pins `deeptools.yaml` and more of `laddr.yaml` to the lock, and sets scikit-learn 1.5.2 in `laddr.yaml`, because LaDDR v1.0.0 fails in `laddr coverage` with the lock's 1.9.1. Rule code is unchanged. Later commits on `module` only add run records to the docs.

## How a host imports it

Snakemake **≥ 9.14** (R28; pathvars need 9.12). Two modules with distinct prefixes are needed, because both stages define `all`:

```python
configfile: "config.yaml"

pathvars:
    results="results",
    resources="resources",

module pantry_pheno:
    pathvars:
        pantry_output="<results>/pantry/phenotypes",
        pantry_intermediate="<results>/pantry/intermediate/phenotyping",
        pantry_reference="<resources>/pantry/reference",
    snakefile:
        github("twillis209/Pantry", path="phenotyping/Snakefile", tag="module-v0.1.5")
    config:
        config["pantry"]["phenotyping"]

use rule * from pantry_pheno as pantry_pheno_*

module pantry_pheast:
    pathvars:
        pantry_output="<results>/pantry/pheast",
        pantry_intermediate="<results>/pantry/intermediate/pheast",
        pantry_resources="<resources>/pantry",
    snakefile:
        github("twillis209/Pantry", path="pheast/Snakefile", tag="module-v0.1.5")
    config:
        config["pantry"]["pheast"]

use rule * from pantry_pheast as pantry_pheast_*

rule all:
    default_target: True
    input:
        rules.pantry_pheno_all.input,
        rules.pantry_pheast_all.input,
```

- **Pathvars.** Phenotyping takes `pantry_output`, `pantry_intermediate` and `pantry_reference`; pheast takes `pantry_output`, `pantry_intermediate` and `pantry_resources`. Each has a standalone default (`output`, `intermediate`, `intermediate/reference` or `intermediate/resources`), so set all of them to control the layout. Outputs land in `<pantry_output>/{dataset}/…` and references in `<pantry_reference>/{reference}/…`.
- **Config shape.** The working example is `tests/host/config.yaml`; the annotated templates are `phenotyping/config.yml` and `pheast/config.yml`. In outline:

  ```yaml
  pantry:
    phenotyping:
      references: {<ref>: {genome, annotation, edit_sites_bed?, annotation_gff3?}}
      defaults: {...}        # optional, deep-merged under every dataset
      datasets:
        <name>:
          reference: <ref>
          samples: [..] | path
          fastq_map: {sample: [[R1..], [R2..]]} | path   # + fastq_dir
          paired_end: true
          read_length: 75
          modality_groups: {...}                          # as upstream
          cross_modality: {modalities: [...]}             # optional
          laddr: {enabled: false, ...}                    # optional
    pheast:
      defaults: {...}
      datasets:
        <name>:
          phenotype_dir: results/pantry/phenotypes/<name>  # = phenotyping output links the stages
          geno_prefix: /path/plink_prefix
          samples: [..] | path
          modalities: {expression: {grouped: false}, ...}
          analyses: {qtl: {files: [...]}, heritability: {...}, twas: {...}}
          qtl_min_allele_count: 10    # qtl_maf_threshold, qtl_cis_window optional
  ```

  Pass new-format configs only. A legacy `intermediate_dir` in a module config would override the host's pathvars.
- **Pass `config:`** in each `module` directive, as above. Snakemake skips a module's standalone `configfile: "config.yml"` only when `config:` is given; otherwise it would load a `config.yml` from the host's working directory.
- **Linking the stages is by path text.** A pheast `phenotype_dir` must equal the phenotyping `<pantry_output>/{dataset}` path as the host resolves it (`results/pantry/phenotypes/<name>` above). A relative vs absolute spelling of the same directory does not link, and pheast then fails with missing input files.
- **`defaults:`** can be added to or overridden per dataset, but a dataset cannot remove a key from it. Keep modality groups, modalities and analyses that vary between datasets out of `defaults`. (Removal support is deferred.)
- **Stage directories and `--config`.** Inside `phenotyping/` or `pheast/`, `--config key=value` with a legacy key (e.g. `read_length=100`) is an error next to the new-format `config.yml`. Use `--configfile override.yml` with `datasets: {default: {read_length: 100}}` (deep-merged), edit `config.yml`, or run an old-format config from another directory (`--directory`). Nested `--config` values arrive as strings, so `--config 'datasets={...}'` works only for string-valued keys.
- **TWAS needs apptainer.** The FUSION rules (`twas_compute_weights_batch`, `twas_assemble_summary`) run only in the pinned image `ghcr.io/twillis209/pantry-fusion` (FUSION scripts and gcta under `/opt/fusion`), so the profile needs `software-deployment-method: [conda, apptainer]` (below); no network is needed at run time beyond the image pulls. A `github()` module fetches no binary file. A `github()` module builds its own conda envs: env files fetched from GitHub lose their final newline (Snakemake source-cache bug), so their hashes differ from the local-path ones (13 envs, about 15 min in Task 18). Outputs are unaffected.
- **Profile for the containers and GPU tensorQTL.**

  ```yaml
  software-deployment-method: [conda, apptainer]
  apptainer-args: "--nv --bind=<data root>[,<module checkout if loaded by local path>]"
  set-resources:
    pantry_pheast_tensorqtl_cis:             {slurm_partition: <gpu partition>, gpu: 1}
    pantry_pheast_tensorqtl_cis_independent: {slurm_partition: <gpu partition>, gpu: 1}
    pantry_pheast_tensorqtl_trans:           {slurm_partition: <gpu partition>, gpu: 1}
  ```

  - `set-resources` keys are the **prefixed** rule names.
  - Use `gpu: 1`, not `slurm_extra: --gpus`/`--gres`; the current slurm executor plugin rejects those.
  - `--nv` exposes the GPU to the container. The bind must cover every data path, including symlink targets, and the module checkout when the module is loaded by local path.
  - Without apptainer, the tensorQTL rules run the host `python3`, which must provide tensorqtl; the TWAS rules cannot run.
  - Full example: `tests/host/profile/config.yaml`.

## Acceptance results

Everything ran on the GDML HPC (SLURM, `sjob`, profile `tests/profile` or `tests/host/profile`), with Snakemake 9.27.1. `H=/data/hps/assoc/private/gdml/user/twil15`.

**Baseline.** Upstream `8da7787` was run on both test datasets. Phenotyping ran twice and pheast twice, to measure run-to-run noise. Outputs are in `$H/pantry-baseline/archive/{phenotyping,phenotyping_run2,pheast_run1,pheast_run2}`, and the noise reports are `archive/phenotyping_noise.md`, `archive/pheast_noise.md` and `archive/pheast_noise_calibrated.md`. Two sources of nondeterminism showed up:
- **kallisto 0.48.0 quant** varies run to run. It depends on node/load, not on thread count. This affects 6 of 21 phenotyping files: `{alt_polyA,expression,isoforms}.bed.gz` and their `unnorm/` sources.
- **tensorQTL** permutations and tie-breaks are unseeded. Even `cis_independent` row sets differ between the two baseline runs.

**Comparison method.** `tests/compare_outputs.py <baseline> <new> --noise <run1> <run2>`:
- All other files must be byte-identical (`.gz` compared decompressed).
- Kallisto-derived BEDs must fall within 2× the baseline pair's max |diff| (R12). Since the R12 amendment (Task 22 review, below), `--noise` takes two or more runs and the envelope is 2× the max over all pairs.
- For tensorQTL files, the judgement is attribution (R21, R22): the tensorQTL inputs (covariates, pruned genotypes, BEDs, groups) must be byte-identical to the baseline's, and the `run_tensorqtl.py`/`genotype_qc.py` command lines identical. The comparison of the `*.cis_qtl`/`*.cis_independent_qtl` files is reported, but it does not block.

| Gate (criterion) | What | Run dir (`$H/pantry-runs/`) | Result | Report |
|---|---|---|---|---|
| Gate 1 (1, 5) | phenotyping standalone, pathvars, old-format config | `gate1-pheno` | 15/21 files byte-identical; 6 kallisto files differ. Accepted on substitution evidence (R11): with the baseline's kallisto outputs swapped in → MATCH. Later re-compared with `--noise` → MATCH (run 8392777, compares 8393107, 8393118, 8393447) | `gate1-pheno/compare.md`, `gate1-pheno/compare_noise.md`, `tests/hpc/README.md` (Gates 1-3) |
| Gate 2 (1, 5) | phenotyping `{dataset}` layout, old-format config, standalone (job 8393449) | `gate2-pheno` | MATCH (kallisto files within noise; compare job 8393775) | `gate2-pheno/compare.md`, `tests/hpc/README.md` (Gates 1-3) |
| Gate 3 (1, 5) | pheast `{dataset}`, checkpoints, FUSION install, old-format config, CPU (job 8394390, 195/195) | `gate3-pheast` | heritability ×8, TWAS ×8, genotype_qc identical. 10 mismatches, all tensorQTL; inputs and commands identical → accepted (R21); compare job 8394624, commands job 8394629 | `gate3-pheast/compare.md`, `cmds_{baseline,gate3}.norm`, `tests/hpc/README.md` (Gates 1-3) |
| Gate 2 conda (1) | per-rule conda envs (job 8394948) | `gate2-pheno-conda` | MATCH | `gate2-pheno-conda/compare.md` |
| Gate 3 conda+GPU (1) | conda envs + tensorQTL container on A100 (jobs 8394949, 8395487) | `gate3-pheast-conda` | all non-tensorQTL files identical. 85 tensorQTL mismatches (non-blocking, R22); inputs and commands identical | `gate3-pheast-conda/compare.md` |
| Gate 3 GPU, tensorqtl 1.0.10 image (1) | `ghcr.io/twillis209/pantry-tensorqtl@sha256:65a42f9f…` on A100 (job 8399104, 195/195); image on 1 CPU for all 8 cis runs (jobs 8399314, 8399317–8399323) | `gate3-pheast-gpu1010`, `t16` | non-tensorQTL files pass; 81 tensorQTL mismatches (non-blocking, R22); inputs and commands identical. On 1 CPU core, cis nominal stats identical to the baseline on same-top-variant rows, at written precision | `gate3-pheast-gpu1010/compare.md`, `tests/hpc/README.md` (Task 16) |
| Gate 3, FUSION image (1; Task 19, R26) | `ghcr.io/twillis209/pantry-fusion@sha256:87335916…` for the TWAS rules, `install_fusion` removed, conda + apptainer, tensorQTL on GPU (job 8400497, 194/194 = 195 − `install_fusion`; cold dry runs 69 → 68 standalone, 472 → 471 host) | `gate3-pheast-fusion`, `t19` | heritability ×8, TWAS ×8 (tar members byte-identical except `.profile.err`, RDat all.equal 1e-8), genotype_qc pass; `twas/*.{profile,pos,list}` 24/24 byte-identical to gate 3's (GPU 1.0.10 and CPU); 81 tensorQTL mismatches (non-blocking, R22); inputs and 17 command lines identical | `gate3-pheast-fusion/compare.md`, `tests/hpc/README.md` (Task 19) |
| Multi-dataset (2) | phenotyping `full` + `half_se` (8 samples, R1 as single-end) on one reference (job 8395741, 468/468) | `multi-pheno` | `full` MATCH vs baseline; 1 STAR index, 1 kallisto index shared; `half_se` 20/20 targets | `multi-pheno/compare.md` |
| Multi-dataset (2) | pheast `full` + `subset` (300 samples, MAC 20) (job 8395742, 388/388, GPU) | `multi-pheast` | `full`: heritability/TWAS/genotype_qc pass, 76 tensorQTL mismatches (non-blocking); `subset` 33/33 targets | `multi-pheast/compare.md` |
| Host, local path (3) | `tests/host/Snakefile`: dry run of `all` (472 jobs, incl. pheast `chained` on phenotyping output in the same DAG); full run of phenotyping `full` + pheast `geuvadis` (job 8396531, 539/539, tensorQTL on GPU with prefixed set-resources) | `host` | phenotyping MATCH; `cross_modality` byte-identical to upstream `combine_modalities.sh`; pheast: blocking files pass, 79 tensorQTL mismatches (non-blocking) | `host/compare_pheno.md`, `host/compare_pheast.md` |
| Host, `github()` (3) | At `module-v0.1.2` (= `0d30215`): dry run of `all` (job 8400681, 471 jobs, no binary fetched) and full run as Task 18's (job 8402594, 548/548; `host-github-v012`: phenotyping MATCH, intron_retention byte-identical, `cross_modality` identical, pheast blocking files incl. TWAS pass, 84 tensorQTL mismatches (R22), 17 command lines identical). Earlier, `tests/host/Snakefile.github` at `module-v0.1.1` (= `befa72c`), full run of phenotyping `full` (+ intron_retention) and pheast `geuvadis` (job 8399793, 549/549, tensorQTL on GPU); earlier, at `module-rc1`, dry run and two jobs only (job 8398071) | `host-github-v011` | phenotyping MATCH; intron_retention byte-identical to upstream; `cross_modality` identical to the script; pheast: blocking files pass, 80 tensorQTL mismatches (non-blocking), inputs and 17 command lines identical; gcta downloaded, sha256 OK | `host-github-v011/compare_{pheno,ir,pheast}.md`, `tests/hpc/README.md` (Task 18) |
| LaDDR (1, 3; Task 22, R29) | upstream 8da7787 ×2 + seeded reruns in a LaDDR clone of the baseline env; module standalone (jobs 8405700, 8406311) and `github()` host at 54423dc (job 8406551, 575/575), each rerun seeded on U1's regressed BEDs (jobs 8406432, 8407165) | `t22`, `host-github-v015` | seeded, identical inputs: module and host byte-identical to upstream in all LaDDR files (finals and intermediates); config YAML and 26 command lines identical after path and design mappings; unseeded LaDDR is nondeterministic upstream (R29) | `tests/hpc/README.md` (Task 22) |
| Lint, format (4) | `snakemake --lint`, `snakefmt` | local | remaining warnings justified in `CHANGES-module.md` (no log directives, as upstream; justified `run:` blocks; mixed rules/functions; two upstream-body items) | `CHANGES-module.md` |

How far the tensorQTL results are from upstream's:
- On the cis_qtl files, the GPU + container results match the baseline at float level. `pval_nominal` max relative difference is ~1e-4; 4–18 rows per file pick a different top variant, besides equal-p tie-breaks.
- Permutation columns are within the run-to-run noise, except `intron_retention` qval.
- `cis_independent` row sets change by about as much as between the two baseline runs.
- The tensorQTL image now has upstream's versions (tensorqtl 1.0.10, torch 2.3.1, pandas 2.3.3, Python 3.12; Task 16, superseding R23). On 1 CPU core, as the baseline ran, its nominal statistics are identical to the baseline's on rows with the same top variant, at written precision, in all 8 modalities. The GPU, or several CPU cores, shift `pval_nominal` by up to ~1e-4 relative.

**MAJIQ (`intron_retention`)** is now verified at run time (Task 17): with MAJIQ 3.0.22 in the named env `majiq3`, the refactored rules reproduce upstream 8da7787 byte for byte on the test data (`intron_retention.bed.gz`, `.phenotype_groups.txt`, `unnorm/intron_retention.bed`), as do two upstream runs of each other. A run with renamed, symlinked BAMs names the `.sj` files after the targets and gives the same outputs. See `tests/hpc/README.md`, Task 17.

**MAJIQ image (R27, Task 20).** MAJIQ's [academic licence](https://majiq.biociphers.org/Academic_License.txt) forbids redistribution, so there is no published MAJIQ image: users accept the licence and build one from `containers/majiq/majiq.def` (`apptainer build majiq.sif majiq.def`, unprivileged, about 5 minutes) or the `Dockerfile`, which fetch MAJIQ v3.0.22 (commit `d841baad`) at build time, then set `majiq_container` (README, "MAJIQ (intron_retention)"). The MAJIQ rules then run in the image with no conda env. Rerunning the Task 17 configs with the image under `--sdm conda apptainer` gave intron_retention outputs byte-identical to upstream's, with plain and with symlinked BAMs. Hosts must bind their data and symlinked-BAM targets in `apptainer-args`. Snakemake binds its own source cache (`~/.cache/snakemake`, which holds `script()` paths and the bundled licence file), so no profile bind is needed (Task 20's claim that 9.12-9.13 need one was wrong: they bind the default cache path; R28). The Docker and `docker://`/`oras://` routes are untested. See `tests/hpc/README.md`, Task 20.

**Pins and the full `github()` run (Task 21, `module-v0.1.4`).** Each pinned env was built on the HPC and its `conda list --explicit` checked against the lock: txrevise 272/272 packages in the lock; LaDDR 105 conda packages, all in the lock but scipy (1.16.3, as LaDDR requires < 1.17); the rebuilt MAJIQ image's conda (64) and pip lists equal `majiq3.lock.txt`. The txrevise reference files from the pinned env match the baseline's (merged GFF3s, GTF, tags byte-identical; batch GFF3s except a `##date` line; the `.rds` except its provenance metadata). The `github()` host run (job 8404292; phenotyping `full` + intron_retention from the rebuilt image, pheast `geuvadis`; profile without the source-cache bind) ran 548/548; phenotyping within the noise envelope, intron_retention byte-identical (4 files), cross_modality identical to upstream's script, pheast mismatches only in tensorQTL files (78, R22), heritability, TWAS and genotype_qc passing, tensorQTL command lines 17/17 identical, covariate inputs identical. LaDDR was not part of that run (Task 22, below). See `tests/hpc/README.md`, Task 21.

**HPC note: pulling the images.** The images `ghcr.io/twillis209/pantry-tensorqtl` and `ghcr.io/twillis209/pantry-fusion` are public, but stale stored ghcr.io credentials make apptainer fail with 403 on the token request (seen on this account, Task 16). Run `apptainer registry logout docker://ghcr.io` (and drop `ghcr.io` from `~/.docker/config.json`) or refresh the login. The runs here used an anonymous pre-pull into the profile's `apptainer-prefix` (`tests/hpc/README.md`, Tasks 16 and 19).

**LaDDR at run time (Task 22, R29, `module-v0.1.5`).** Upstream `8da7787` was run twice (U1, U2) on `tests/configs/phenotyping_baseline_laddr.yml` (upstream's config with intron_retention and LaDDR on, upstream's LaDDR defaults, `dataset_name: pantry`, all 8 `regress_modalities`), then the module standalone (M1) and in a `github()` host (M2). Envs:
- Upstream's env has no LaDDR, and upstream 8da7787 with LaDDR v1.0.0 cannot run in it: `laddr coverage` fails under the lock's scikit-learn 1.9.1 (PCA on a frame with duplicate integer column names). The upstream side therefore ran in a clone of the baseline env (the original untouched) with LaDDR v1.0.0 (`2addf6fe`), and scikit-learn 1.5.2 (needed because 1.9.1 fails, narwhals `DuplicateError`; 1.5.2 is the known-good version LaDDR's `requirements.txt` names) and scipy 1.16.3 (LaDDR requires < 1.17) from conda-forge. Its explicit list differs from the lock in those two packages only; pip added the LaDDR stack.
- The module's `laddr.yaml` has the same two versions and builds. `deeptools.yaml`, which had pinned only deeptools and resolved python 3.13, numpy 2.5 and pyBigWig 0.3.25, now pins its dependency closure from the lock. `laddr.yaml` now also pins pyBigWig, matplotlib, pillow and wrapt from the lock rather than leaving them to pip.
- Every module that LaDDR's CLI and bamCoverage actually import has the same version on both sides (scipy as the same conda build). The one exception is `psutil`, which only the clone has and which is imported opportunistically.

Results:
- **LaDDR v1.0.0 is not deterministic.** `remove_pcs_from_phenos` runs an unseeded randomized PCA (scikit-learn's `auto` solver for 16 × thousands). A rerun on identical inputs gives different coverage residuals, models and latent phenotypes, and even a different phenotype count (624 vs 625). Kallisto noise in the regressed modalities adds to this. The U1/U2 noise envelope therefore spans the whole rank-normal range (max |diff| 3.13) and differs in row set, so it is no usable bar; `compare_outputs.py` refuses it.
- **Acceptance (R29) was on identical inputs with numpy's global RNG seeded.** This used a test-only `sitecustomize.py` with `T22_SEED=0`, and for conda jobs `--conda-not-block-search-path-envvars`. Two seeded upstream reruns were byte-identical. With U1's eight regressed BEDs substituted, the module standalone (M1s) and the `github()` host (M2s) were **byte-identical to seeded upstream**: the final `latent_residual.bed.gz`, `.tbi`, `.phenotype_groups.txt` and `unnorm/latent_residual.bed` (624 phenotypes; `cmp` on the raw files). Every LaDDR intermediate (bigWigs, `info/`, `gene_bins/`, `covg_norm/`, `models/`, latent table) is byte-identical too, except the gzip files among them (`gene_bins/*.bed.gz`, `covg_norm/*.bins.tsv.gz`, `info/exons.tsv.gz`, the latent table), which are content-identical (compared decompressed). bigWigs, `info/` and `gene_bins/` are identical in this sense across all runs, seeded or not.
- The LaDDR config written by `write_laddr_config.py` is identical after path normalisation. The 26 LaDDR command lines (16 bamCoverage, config, `laddr setup/binning/coverage/fit/transform`, assemble, normalize, bgzip, groups) are identical once three documented design differences are mapped: the module's `mv` of the latent table to a stable name, the assemble input path that follows from it, and `--samples intermediate/samples.txt` (same content as `input/samples.txt`).
- `normalize_phenotypes.py` gives byte-identical output in the clone, the baseline env and the module's `python.yaml` env on all 9 unnormalised BEDs, so the clone changed no normalisation.
- U1's other outputs vs the archive: RNA_editing, splicing, stability and intron_retention byte-identical. Kallisto files exceeded the single-pair envelope (expression 0.744 vs 0.732, isoforms 0.213 vs 0.167). In U2, M1 and M2, `alt_TSS` varied for the first time (the archive pair had it identical), with three concurrent runs. U1's unnormalised kallisto files stayed within it; `unnorm/alt_TSS.bed` in U2, M1 and M2 did not (9.4e-4, 1.0e-3, 1.2e-3 against an envelope of 0). Under the R12 amendment (below) all of these are within the envelope.
- M2 (`host-github-v015`, job 8406551, 575/575 = Task 21's 548 + 27 LaDDR jobs; host config `tests/host/config_laddr.yaml`) met the Task 21 bar with one exception: PHENO_COMPARE_EXIT=1, on `unnorm/alt_TSS.bed` alone (1.2e-3 against an envelope of 0; within the amended envelope, 1.9e-3). The rest held: intron_retention byte-identical, cross_modality identical to upstream's script, pheast mismatches only in tensorQTL files (89, R22), the pheast tensorQTL/genotype_qc command lines identical to the baseline's (CMDS_IDENTICAL 17/17, job 8407214), and the fetched LaDDR and deeptools envs equal to the local-path ones.

See `tests/hpc/README.md`, Task 22.

**R12 amendment (Task 22 review; methodology ruling, awaiting the user's ratification).** The single archive pair made a fragile envelope for kallisto files: `alt_TSS` had an envelope of 0 because the pair happened to agree, and the rank-normalised files have quantised noise (below), so "2× one pair" had no margin. `compare_outputs.py --noise` now takes two or more runs, and the envelope per file is 2× the max |diff| over **all pairs** of them (with two runs, nothing changes). The noise runs for phenotyping are the four upstream-8da7787 runs `{archive/phenotyping, archive/phenotyping_run2, U1, U2}` (non-LaDDR outputs; intron_retention: `{archive/phenotyping_ir, archive/phenotyping_ir_run2, U1, U2}`). U1 and U2 qualify: the clone differs from the baseline env only in scipy and scikit-learn, and the only non-LaDDR upstream scripts that import either are `normalize_phenotypes.py`, whose output is byte-identical in all three envs (job 8406357), and the RNA_editing preparation, whose output is byte-identical to the archive's. No third baseline run was made. Re-judged (job 8407236; `t22/r12a.sh`, `t22/pairtable.py`):

| file | archive pair | max over 6 pairs (pair) | old allowed | new allowed | U1 | U2 | M1 | M2 | Task 21 (v0.1.4) |
|---|---|---|---|---|---|---|---|---|---|
| `alt_TSS.bed.gz` | 0 | 0.108 (arch~U2) | 0 | 0.216 | 0 | 0.108 | 0.151 | 0 | 0 |
| `alt_polyA.bed.gz` | 0.164 | 0.222 (arch~U2) | 0.328 | 0.444 | 0.108 | 0.222 | 0 | 0.222 | 0.222 |
| `expression.bed.gz` | 0.366 | 0.744 (arch~U1) | 0.732 | 1.488 | **0.744** | 0.366 | 0.366 | 0.366 | 0.366 |
| `isoforms.bed.gz` | 0.0835 | 0.213 (arch~U1) | 0.167 | 0.426 | **0.213** | 0.0835 | 0.0835 | 0.108 | 0.0835 |
| `unnorm/alt_TSS.bed` | 0 | 9.4e-4 (arch~U2) | 0 | 1.9e-3 | 0 | **9.4e-4** | **1.0e-3** | **1.2e-3** | 0 |
| `unnorm/alt_polyA.bed` | 1.4e-3 | 1.4e-3 (arch~run2) | 2.8e-3 | 2.8e-3 | 7.8e-4 | 3.2e-4 | 5.9e-4 | 5.8e-4 | 4.6e-4 |
| `unnorm/expression.bed` | 338 | 338 (arch~run2) | 676 | 676 | 338 | 338 | 338 | 338 | 338 |
| `unnorm/isoforms.bed` | 3.0e-3 | 3.0e-3 (arch~run2) | 6.1e-3 | 6.1e-3 | 3.0e-3 | 3.0e-3 | 3.0e-3 | 3.0e-3 | 3.0e-3 |

Bold: failed the old envelope. Every candidate now passes (REST exit 0 for U1, U2, M1, M2 and the Task 21 run; IR exit 0, all byte-identical). RNA_editing, splicing, stability and intron_retention (normalised and `unnorm/`) are byte-identical in all runs and keep an envelope of 0, as intended for deterministic files; no kallisto file has an envelope of 0 any more. Caveats:
- U1 and U2 pass by construction: archive~U1 and archive~U2 are among the pairs, so their deviation is at most half the envelope. The substantive re-judgements are M1, M2 and the Task 21 run, and no file that passed before can fail now (the envelope only grows).
- **Rank quantisation.** `normalize_phenotypes.py` maps ranks to `norm.ppf(rank / 17)` (16 samples, ties averaged), so every normalised deviation is a difference of those quantiles: all eight values observed here match one to within 4e-5 (0.366 = rank 2 → 3.5, 0.744 = rank 1 → 3.5). The adjacent-rank step is 0.378 at the extremes and 0.148 in the middle. The amended `alt_TSS.bed.gz` envelope (0.216) is still below one extreme-rank step, so a single swap of the lowest or highest sample in `alt_TSS` would fail it; `alt_polyA` (0.444), `isoforms` (0.426) and `expression` (1.488) admit one extreme step or more. The bar is kept at 2× the max over all upstream pairs and not widened further.
- M2's normalised BEDs were taken from `host-github-v015/snap_unseeded/orig_beds`, since the seeded rerun (job 8407165) replaced them in `results/` with U1's. The old two-run compare on that view reproduces job 8407162's lines exactly, and the old (a65ac3e) and new tool give identical two-run reports.

**Scaling (dry runs, Task 18).** Synthetic host configs with N datasets (copies of the test `full` dataset, 7 modalities + cross_modality, chained pheast datasets) on the HPC, Snakemake 9.27.1:

| Case | Jobs | Wall | Max RSS |
|---|---|---|---|
| upstream 8da7787 standalone phenotyping, 1 dataset x 1000 samples | 12149 | 32.8 s | 393 MB |
| this module, standalone phenotyping, 1 x 1000 | 12151 | 38.0 s | 401 MB |
| host, 17 x 100 | 22053 | 64 s | 628 MB |
| host, 17 x 1000 | 205653 | 12m20 | 4.4 GB |

About 3 ms per job in all cases, so cost grows roughly linearly with the job count. A profile of 17 x 300 shows Snakemake's own work (resource evaluation while printing jobs, DAG building, input sizes, conda env specs) dominating; Pantry's helper functions take under 2% of the time. Details: `tests/hpc/README.md`, Task 18.

**Margin to watch.** The isoforms kallisto margin in the multi-dataset run was close: max |diff| 0.149 vs allowed 0.167. The envelope comes from a single baseline pair.

## Deviations from the handover contract

- **Item 1**: each stage is its own module with its own config section and prefix (`config["pantry"]["phenotyping"]`, `config["pantry"]["pheast"]`), not one `config["pantry"]` with prefix `pantry_*`. Both define `all`.
- **Item 8**: FUSION is vendored upstream with modifications; rather than fetching a FUSION commit, the module's copy, gcta and gemma ship in a digest-pinned image built from the repository (R26).
- **Item 11**: standalone works, with the old config format too, but outputs gain a `default/` level (`output/default/…`).
- **Criterion 1**: kallisto- and tensorQTL-derived files are judged within measured noise or by attribution (R11, R12, R21, R22).
- MAJIQ: no published image (licence); each user builds one (R27).
- LaDDR env: scikit-learn 1.5.2 and scipy 1.16.3, not the lock's 1.9.1 and 1.17.1, because LaDDR v1.0.0 cannot run with them (R29).
- Also: `index_fasta` removed (R15). LaDDR `dataset_name` defaults to the dataset key (R16). Validation is stricter on quoted numbers (R14, R18).

## Follow-ups

1. Report upstream: LaDDR v1.0.0 fails under scikit-learn ≥ 1.9 (`laddr coverage`), its PCA in `remove_pcs_from_phenos` is unseeded (run-to-run differences, Task 22), and upstream's env cannot run LaDDR as documented. (Runtime verification of MAJIQ and LaDDR: done, Tasks 17 and 22.)
2. Env pinning (mostly done in Tasks 21 and 22): LaDDR's pip-only dependencies (bx-python, scikit-fda, statsmodels, tqdm, numba and theirs; no lock covers them) are unpinned except `multimethod`; the Task 22 versions are recorded in `tests/hpc/README.md`. MAJIQ's named env and image are built from the pinned `containers/majiq/` recipes, but the rules name the env or take the image path from config rather than pin them.
3. Report the Snakemake source-cache bugs upstream: binary files fail (no longer affects Pantry since R26), and text files lose their final newline under `github()` (which also makes a `github()` module rebuild its conda envs).
4. Optional: `log:` directives for all rules. A third baseline run is no longer needed for phenotyping: the R12 amendment uses U1 and U2 as further noise runs (its ratification by the user is pending). The `alt_TSS` envelope stays below one extreme-rank step; tensorQTL envelopes are still from one pair.
5. Full list: `CHANGES-module.md`, "Known limitations and follow-ups".
