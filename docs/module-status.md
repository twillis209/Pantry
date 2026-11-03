# Pantry module refactor: status note

Date: 2026-11-02. This is the report back requested by the scPantry handover (`~/projects/scPantry/docs/handover/2026-10-30-pantry-fork-handover.md`, "Report back").
Changes in full: `CHANGES-module.md`. Design: `docs/superpowers/specs/2026-10-30-module-refactor-design.md`. Run records: `tests/hpc/README.md`.

## Branch and tag

- Fork `twillis209/Pantry`, branch **`module`**, release tag **`module-v0.1.0`**. Base: upstream `PejLab/Pantry` `8da7787`. Nothing was pushed to upstream.
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

## How a host imports it

Snakemake **≥ 9.12** (pathvars). Two modules with distinct prefixes are needed, because both stages define `all`:

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
        github("twillis209/Pantry", path="phenotyping/Snakefile", tag="module-v0.1.0")
    config:
        config["pantry"]["phenotyping"]

use rule * from pantry_pheno as pantry_pheno_*

module pantry_pheast:
    pathvars:
        pantry_output="<results>/pantry/pheast",
        pantry_intermediate="<results>/pantry/intermediate/pheast",
        pantry_resources="<resources>/pantry",
    snakefile:
        github("twillis209/Pantry", path="pheast/Snakefile", tag="module-v0.1.0")
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
      # fusion_dir / gcta: optional overrides; gcta_url / gcta_sha256: download source
  ```

  Pass new-format configs only. A legacy `intermediate_dir` in a module config would override the host's pathvars.
- **`github()` and gcta.** A `github()` module downloads `gcta_nr_robust` in the `pantry_pheast_install_fusion` job (pinned URL plus sha256; Snakemake's source cache can't fetch binaries). That job needs network access. **Offline sites set `gcta`** (or `fusion_dir` for the R scripts as well). Local-path modules copy the bundled binary.
- **Profile for GPU tensorQTL.**

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
  - Without apptainer, the tensorQTL rules run the host `python3`, which must provide tensorqtl.
  - Full example: `tests/host/profile/config.yaml`.

## Acceptance results

Everything ran on the GDML HPC (SLURM, `sjob`, profile `tests/profile` or `tests/host/profile`), with Snakemake 9.27.1. `H=/data/hps/assoc/private/gdml/user/twil15`.

**Baseline.** Upstream `8da7787` was run on both test datasets. Phenotyping ran twice and pheast twice, to measure run-to-run noise. Outputs are in `$H/pantry-baseline/archive/{phenotyping,phenotyping_run2,pheast_run1,pheast_run2}`, and the noise reports are `archive/phenotyping_noise.md`, `archive/pheast_noise.md` and `archive/pheast_noise_calibrated.md`. Two sources of nondeterminism showed up:
- **kallisto 0.48.0 quant** varies run to run. It depends on node/load, not on thread count. This affects 6 of 21 phenotyping files: `{alt_polyA,expression,isoforms}.bed.gz` and their `unnorm/` sources.
- **tensorQTL** permutations and tie-breaks are unseeded. Even `cis_independent` row sets differ between the two baseline runs.

**Comparison method.** `tests/compare_outputs.py <baseline> <new> --noise <run1> <run2>`:
- All other files must be byte-identical (`.gz` compared decompressed).
- Kallisto-derived BEDs must fall within 2× the baseline pair's max |diff| (R12).
- For tensorQTL files, the judgement is attribution (R21, R22): the tensorQTL inputs (covariates, pruned genotypes, BEDs, groups) must be byte-identical to the baseline's, and the `run_tensorqtl.py`/`genotype_qc.py` command lines identical. The comparison of the `*.cis_qtl`/`*.cis_independent_qtl` files is reported, but it does not block.

| Gate (criterion) | What | Run dir (`$H/pantry-runs/`) | Result | Report |
|---|---|---|---|---|
| Gate 1 (1, 5) | phenotyping standalone, pathvars, old-format config | `gate1-pheno` | 15/21 files byte-identical; 6 kallisto files differ. Accepted on substitution evidence (R11): with the baseline's kallisto outputs swapped in → MATCH. Later re-compared with `--noise` → MATCH | `gate1-pheno/compare.md`, `gate1-pheno/compare_noise.md` |
| Gate 2 (1, 5) | phenotyping `{dataset}` layout, old-format config, standalone (job 8393449) | `gate2-pheno` | MATCH (kallisto files within noise) | `gate2-pheno/compare.md` |
| Gate 3 (1, 5) | pheast `{dataset}`, checkpoints, FUSION install, old-format config, CPU (job 8394390, 195/195) | `gate3-pheast` | heritability ×8, TWAS ×8, genotype_qc identical. 10 mismatches, all tensorQTL; inputs and commands identical → accepted (R21) | `gate3-pheast/compare.md`, `cmds_{baseline,gate3}.norm` |
| Gate 2 conda (1) | per-rule conda envs (job 8394948) | `gate2-pheno-conda` | MATCH | `gate2-pheno-conda/compare.md` |
| Gate 3 conda+GPU (1) | conda envs + tensorQTL container on A100 (jobs 8394949, 8395487) | `gate3-pheast-conda` | all non-tensorQTL files identical. 85 tensorQTL mismatches (non-blocking, R22); inputs and commands identical | `gate3-pheast-conda/compare.md` |
| Gate 3 GPU, tensorqtl 1.0.10 image (1) | `ghcr.io/twillis209/pantry-tensorqtl@sha256:65a42f9f…` on A100 (job 8399104, 195/195); image on 1 CPU for all 8 cis runs (jobs 8399314, 8399317–8399323) | `gate3-pheast-gpu1010`, `t16` | non-tensorQTL files pass; 81 tensorQTL mismatches (non-blocking, R22); inputs and commands identical. On 1 CPU, cis nominal stats identical to the baseline | `gate3-pheast-gpu1010/compare.md`, `tests/hpc/README.md` (Task 16) |
| Multi-dataset (2) | phenotyping `full` + `half_se` (8 samples, R1 as single-end) on one reference (job 8395741, 468/468) | `multi-pheno` | `full` MATCH vs baseline; 1 STAR index, 1 kallisto index shared; `half_se` 20/20 targets | `multi-pheno/compare.md` |
| Multi-dataset (2) | pheast `full` + `subset` (300 samples, MAC 20) (job 8395742, 388/388, GPU) | `multi-pheast` | `full`: heritability/TWAS/genotype_qc pass, 76 tensorQTL mismatches (non-blocking); `subset` 33/33 targets | `multi-pheast/compare.md` |
| Host, local path (3) | `tests/host/Snakefile`: dry run of `all` (472 jobs, incl. pheast `chained` on phenotyping output in the same DAG); full run of phenotyping `full` + pheast `geuvadis` (job 8396531, 539/539, tensorQTL on GPU with prefixed set-resources) | `host` | phenotyping MATCH; `cross_modality` byte-identical to upstream `combine_modalities.sh`; pheast: blocking files pass, 79 tensorQTL mismatches (non-blocking) | `host/compare_pheno.md`, `host/compare_pheast.md` |
| Host, `github()` (3) | `tests/host/Snakefile.github` at `module-rc1` (job 8398071): dry run of `all` (470 jobs); ran `install_fusion` (gcta downloaded, sha256 OK, identical to bundled) and one script rule (`exonic_intronic_from_GTF`, output identical to the local host's) | `host-github` | pass. No full pipeline run via `github()` | `tests/hpc/README.md` (Task 14) |
| Lint, format (4) | `snakemake --lint`, `snakefmt` | local | remaining warnings justified in `CHANGES-module.md` (no log directives, as upstream; justified `run:` blocks; mixed rules/functions; two upstream-body items) | `CHANGES-module.md` |

How far the tensorQTL results are from upstream's:
- On the cis_qtl files, the GPU + container results match the baseline at float level. `pval_nominal` max relative difference is ~1e-4; 4–18 rows per file pick a different top variant, besides equal-p tie-breaks.
- Permutation columns are within the run-to-run noise, except `intron_retention` qval.
- `cis_independent` row sets change by about as much as between the two baseline runs.
- The tensorQTL image now has upstream's versions (tensorqtl 1.0.10, torch 2.3.1, pandas 2.3.3, Python 3.12; Task 16, superseding R23). On 1 CPU, as the baseline ran, it reproduces the baseline's nominal statistics exactly in all 8 modalities. The GPU, or several CPU threads, shift `pval_nominal` by up to ~1e-4 relative.

**MAJIQ (`intron_retention`)** is now verified at run time (Task 17): with MAJIQ 3.0.22 in the named env `majiq3`, the refactored rules reproduce upstream 8da7787 byte for byte on the test data (`intron_retention.bed.gz`, `.phenotype_groups.txt`, `unnorm/intron_retention.bed`), as do two upstream runs of each other. A run with renamed, symlinked BAMs names the `.sj` files after the targets and gives the same outputs. See `tests/hpc/README.md`, Task 17.

**LaDDR** was restructured but checked by dry run only (`tests/configs/dryrun_majiq_laddr.yml`), as agreed. It is excluded from every run above.

**Margin to watch.** The isoforms kallisto margin in the multi-dataset run was close: max |diff| 0.149 vs allowed 0.167. The envelope comes from a single baseline pair.

## Deviations from the handover contract

- **Item 1**: each stage is its own module with its own config section and prefix (`config["pantry"]["phenotyping"]`, `config["pantry"]["pheast"]`), not one `config["pantry"]` with prefix `pantry_*`. Both define `all`.
- **Item 8**: FUSION is vendored upstream with modifications, and `install_fusion` installs it rather than fetching a FUSION commit. gcta is downloaded from a pinned URL with a sha256 check only for remote modules. gemma comes from conda.
- **Item 11**: standalone works, with the old config format too, but outputs gain a `default/` level (`output/default/…`).
- **Criterion 1**: kallisto- and tensorQTL-derived files are judged within measured noise or by attribution (R11, R12, R21, R22).
- Also: `index_fasta` removed (R15). LaDDR `dataset_name` defaults to the dataset key (R16). Validation is stricter on quoted numbers (R14, R18).

## Follow-ups

1. Runtime verification of LaDDR (MAJIQ done, Task 17).
2. Env pinning: the txrevise python; LaDDR from a release instead of GitHub HEAD. MAJIQ's named env is built from `tests/hpc/envs/majiq3.yml` + `majiq3-install.sh` (v3.0.22), but the rules still name it rather than pin it.
3. Report the Snakemake source-cache bugs upstream: binary files fail, and text files lose their final newline under `github()`.
4. Optional: `log:` directives for all rules; a third baseline run to widen the kallisto and tensorQTL noise envelopes.
5. Full list: `CHANGES-module.md`, "Known limitations and follow-ups".
