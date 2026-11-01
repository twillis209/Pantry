import copy
import importlib.util
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
PHENO = REPO / "phenotyping"
SCHEMA = PHENO / "schemas" / "config.schema.yaml"

_spec = importlib.util.spec_from_file_location(
    "pantry_phenotyping_config", PHENO / "lib" / "config.py"
)
cfg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cfg)


# Oracle: upstream load_fastq_map, copied verbatim from
# git show 8da7787:phenotyping/steps/setup.smk
def upstream_load_fastq_map(map_file: Path, fastq_dir: Path, paired_end=None) -> dict:
    """Load the FASTQ file paths for each sample.

    Returns a dictionary where each value is a tuple (paths, is_paired).
    paths is either a list (single-end) or tuple of two lists (paired-end).
    is_paired is a boolean indicating if the sample is paired-end.
    """
    paths = {}
    with open(map_file, 'r') as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            sample_id = fields[-1]
            is_paired = len(fields) == 3

            if sample_id not in paths:
                paths[sample_id] = (([], []), is_paired) if is_paired else ([], is_paired)
            elif paths[sample_id][1] != is_paired:
                raise ValueError(f"Sample {sample_id} has mixed paired-end and single-end reads")

            if is_paired:
                paths[sample_id][0][0].append(str(fastq_dir / fields[0]))
                paths[sample_id][0][1].append(str(fastq_dir / fields[1]))
            else:
                paths[sample_id][0].append(str(fastq_dir / fields[0]))

    # Validate against config if specified
    if paired_end is not None:
        for sample_id, (_, is_paired) in paths.items():
            if is_paired != paired_end:
                raise ValueError(f"Sample {sample_id} has {'paired' if is_paired else 'single'}-end reads, but config specifies {'paired' if paired_end else 'single'}-end")

    return paths


def legacy_config() -> dict:
    with open(PHENO / "config.yml") as f:
        return yaml.safe_load(f)


def new_config(**dataset_over) -> dict:
    dataset = {
        "reference": "grch38",
        "samples": ["s1", "s2"],
        "fastq_map": {"s1": [["s1_1.fq.gz"], ["s1_2.fq.gz"]]},
        "fastq_dir": "/data",
        "read_length": 75,
        "modality_groups": {
            "expression": {"files": ["expression.bed.gz", "expression.bed.gz.tbi"]}
        },
    }
    dataset.update(dataset_over)
    return {
        "references": {
            "grch38": {"genome": "/ref/genome.fa", "annotation": "/ref/anno.gtf"}
        },
        "datasets": {"d1": dataset},
    }


def test_shim_wraps_legacy_config():
    raw = legacy_config()
    before = copy.deepcopy(raw)
    out = cfg.normalise(raw)
    assert raw == before  # input not mutated
    assert out["references"] == {
        "default": {
            "genome": raw["ref_genome"],
            "annotation": raw["ref_anno"],
            "edit_sites_bed": raw["modality_groups"]["RNA_editing"]["edit_sites_bed"],
        }
    }
    assert list(out["datasets"]) == ["default"]
    ds = out["datasets"]["default"]
    assert ds["reference"] == "default"
    assert ds["samples"] == raw["samples_file"]
    for gone in ("samples_file", "ref_genome", "ref_anno"):
        assert gone not in ds
    assert "edit_sites_bed" not in ds["modality_groups"]["RNA_editing"]
    assert ds["modality_groups"]["RNA_editing"]["edit_sites_min_samples"] == 10
    assert ds["fastq_map"] == raw["fastq_map"]
    assert ds["laddr"]["dataset_name"] == "pantry"
    assert out["legacy_pathvars"] == {}


def test_shim_moves_annotation_gff3():
    raw = legacy_config()
    raw["modality_groups"]["intron_retention"]["annotation_gff3"] = "x.gff3"
    out = cfg.normalise(raw)
    assert out["references"]["default"]["annotation_gff3"] == "x.gff3"
    groups = out["datasets"]["default"]["modality_groups"]
    assert "annotation_gff3" not in groups["intron_retention"]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda c: c.update(datasets=["d1"]),
        lambda c: c["datasets"].update(d1="oops"),
        lambda c: c["references"].update(grch38="oops"),
        lambda c: c["datasets"]["d1"].update(reference=["grch38"]),
        lambda c: c["datasets"]["d1"].update(modality_groups=["expression"]),
        lambda c: c["datasets"]["d1"]["modality_groups"].update(expression="x"),
        lambda c: c["datasets"]["d1"]["modality_groups"].update(expression={}),
    ],
)
def test_malformed_config_raises_value_error(mutate):
    c = new_config()
    mutate(c)
    with pytest.raises(ValueError):
        cfg.normalise(c)


def test_shim_maps_intermediate_dir():
    raw = legacy_config()
    raw["intermediate_dir"] = "foo"
    out = cfg.normalise(raw)
    assert out["legacy_pathvars"] == {
        "pantry_intermediate": "foo",
        "pantry_reference": "foo/reference",
    }
    assert "intermediate_dir" not in out["datasets"]["default"]

    raw["intermediate_ref_dir"] = "~/bar"
    out = cfg.normalise(raw)
    assert out["legacy_pathvars"] == {
        "pantry_intermediate": "foo",
        "pantry_reference": str(Path("~/bar").expanduser()),
    }
    assert "intermediate_ref_dir" not in out["datasets"]["default"]

    del raw["intermediate_dir"]
    out = cfg.normalise(raw)
    assert out["legacy_pathvars"]["pantry_intermediate"] == "intermediate"


def test_new_format_has_no_legacy_pathvars():
    assert cfg.normalise(new_config())["legacy_pathvars"] == {}


def test_defaults_deep_merged():
    c = new_config(read_length=100)
    c["datasets"]["d1"]["modality_groups"] = {
        "RNA_editing": {"edit_sites_min_coverage": 30}
    }
    c["datasets"]["d2"] = {"reference": "grch38", "samples": ["a"], "fastq_map": {}}
    c["references"]["grch38"]["edit_sites_bed"] = "/ref/sites.bed"
    c["defaults"] = {
        "read_length": 75,
        "modality_groups": {
            "RNA_editing": {
                "edit_sites_min_coverage": 20,
                "edit_sites_min_samples": 10,
                "files": ["RNA_editing.bed.gz"],
            }
        },
    }
    out = cfg.normalise(c)
    d1, d2 = out["datasets"]["d1"], out["datasets"]["d2"]
    assert d1["read_length"] == 100
    assert d2["read_length"] == 75
    assert d1["modality_groups"]["RNA_editing"] == {
        "edit_sites_min_coverage": 30,
        "edit_sites_min_samples": 10,
        "files": ["RNA_editing.bed.gz"],
    }
    assert d2["modality_groups"]["RNA_editing"]["edit_sites_min_coverage"] == 20
    # Datasets get independent copies of the defaults.
    d2["modality_groups"]["RNA_editing"]["files"].append("x")
    assert d1["modality_groups"]["RNA_editing"]["files"] == ["RNA_editing.bed.gz"]
    assert c["defaults"]["modality_groups"]["RNA_editing"]["files"] == [
        "RNA_editing.bed.gz"
    ]


def test_dataset_lists_replace_default_lists():
    c = new_config()
    c["defaults"] = {"modality_groups": {"expression": {"files": ["a", "b"]}}}
    out = cfg.normalise(c)
    assert out["datasets"]["d1"]["modality_groups"]["expression"]["files"] == [
        "expression.bed.gz",
        "expression.bed.gz.tbi",
    ]


def test_unknown_reference_raises():
    with pytest.raises(ValueError, match="reference"):
        cfg.normalise(new_config(reference="nope"))


def test_gtf_suffix_enforced():
    raw = legacy_config()
    raw["ref_anno"] = "anno.gtf.gz"
    with pytest.raises(ValueError, match="GTF"):
        cfg.normalise(raw)
    c = new_config()
    c["references"]["grch38"]["annotation"] = "anno.gff3"
    with pytest.raises(ValueError, match="GTF"):
        cfg.normalise(c)
    c["references"]["grch38"]["annotation"] = "ANNO.GTF"
    cfg.normalise(c)


def test_single_end_requires_fragment_length():
    with pytest.raises(ValueError, match="fragment_length"):
        cfg.normalise(new_config(paired_end=False))
    with pytest.raises(ValueError, match="fragment_length_sd"):
        cfg.normalise(new_config(paired_end=False, fragment_length_mean=200))
    cfg.normalise(
        new_config(paired_end=False, fragment_length_mean=200, fragment_length_sd=50)
    )
    raw = legacy_config()
    raw["paired_end"] = False
    with pytest.raises(ValueError, match="fragment_length"):
        cfg.normalise(raw)


def test_rna_editing_requires_edit_sites_bed():
    c = new_config()
    c["datasets"]["d1"]["modality_groups"]["RNA_editing"] = {"files": []}
    with pytest.raises(ValueError, match="edit_sites_bed"):
        cfg.normalise(c)


def test_regress_modalities_validated_and_defaulted():
    c = new_config(laddr={"enabled": True})
    out = cfg.normalise(c)
    assert out["datasets"]["d1"]["laddr"]["regress_modalities"] == ["expression"]
    with pytest.raises(ValueError, match="regress_modalities"):
        cfg.normalise(new_config(laddr={"regress_modalities": ["splicing"]}))
    with pytest.raises(ValueError, match="regress_modalities"):
        cfg.normalise(new_config(laddr={"regress_modalities": "expression"}))


def test_existing_bins_require_dirs():
    with pytest.raises(ValueError, match="existing_info_dir"):
        cfg.normalise(new_config(laddr={"use_existing_bins": True}))


def test_load_samples(tmp_path):
    assert cfg.load_samples(["b", "a"]) == ["b", "a"]
    f = tmp_path / "samples.txt"
    f.write_text("s1\n0012\n")
    assert cfg.load_samples(str(f)) == ["s1", "0012"]


def test_inline_fastq_map_paired_and_single():
    spec = {
        "p": [["p_a_1.fq.gz", "p_b_1.fq.gz"], ["p_a_2.fq.gz", "p_b_2.fq.gz"]],
        "s": [["/abs/s.fq.gz"]],
    }
    assert cfg.load_fastq_map(spec, "/data", None) == {
        "p": (
            (["/data/p_a_1.fq.gz", "/data/p_b_1.fq.gz"], ["/data/p_a_2.fq.gz", "/data/p_b_2.fq.gz"]),
            True,
        ),
        "s": (["/abs/s.fq.gz"], False),
    }
    assert cfg.load_fastq_map({"s": [["rel.fq.gz"]]}, None, False) == {
        "s": (["rel.fq.gz"], False)
    }
    with pytest.raises(ValueError):
        cfg.load_fastq_map({"p": [["a"], ["b", "c"]]}, None, None)


def test_file_fastq_map_matches_upstream(tmp_path):
    map_file = PHENO / "input" / "fastq_map.txt"
    if not map_file.exists():
        map_file = tmp_path / "fastq_map.txt"
        map_file.write_text(
            "a_1.fq.gz\ta_2.fq.gz\tS1\nb_1.fq.gz\tb_2.fq.gz\tS1\nc_1.fq.gz\tc_2.fq.gz\tS2\n"
        )
    expected = upstream_load_fastq_map(map_file, Path("input/fastq"), True)
    got = cfg.load_fastq_map(str(map_file), "input/fastq", True)
    assert got == expected
    assert [type(v[0]) for v in got.values()] == [type(v[0]) for v in expected.values()]

    single = tmp_path / "single.txt"
    single.write_text("a.fq.gz\tS1\nb.fq.gz\tS1\nc.fq.gz\tS2\n")
    assert cfg.load_fastq_map(str(single), "fq", None) == upstream_load_fastq_map(
        single, Path("fq"), None
    )


def test_mixed_endedness_raises(tmp_path):
    mixed = tmp_path / "mixed.txt"
    mixed.write_text("a_1.fq.gz\ta_2.fq.gz\tS1\nb.fq.gz\tS1\n")
    with pytest.raises(ValueError, match="mixed"):
        cfg.load_fastq_map(str(mixed), "fq", None)
    single = tmp_path / "single.txt"
    single.write_text("a.fq.gz\tS1\n")
    with pytest.raises(ValueError, match="config specifies paired-end"):
        cfg.load_fastq_map(str(single), "fq", True)
    with pytest.raises(ValueError, match="config specifies single-end"):
        cfg.load_fastq_map({"S1": [["a_1"], ["a_2"]]}, None, False)


def test_edit_sites_min_samples():
    assert cfg.edit_sites_min_samples(10, 0.1, 16) == 10
    assert cfg.edit_sites_min_samples(10, 0.1, 8) == 8
    assert cfg.edit_sites_min_samples(10, 0.5, 30) == 15


def test_available_modalities():
    groups = legacy_config()["modality_groups"]
    assert cfg.available_modalities(groups) == [
        "RNA_editing",
        "alt_TSS",
        "alt_polyA",
        "expression",
        "intron_retention",
        "isoforms",
        "splicing",
        "stability",
    ]


def test_dataset_targets_includes_laddr_and_cross_modality():
    ds = cfg.normalise(new_config())["datasets"]["d1"]
    assert cfg.dataset_targets(ds) == ["expression.bed.gz", "expression.bed.gz.tbi"]
    ds["laddr"] = {"enabled": True}
    ds["cross_modality"] = {"modalities": ["expression"]}
    assert cfg.dataset_targets(ds) == [
        "expression.bed.gz",
        "expression.bed.gz.tbi",
        "latent_residual.bed.gz",
        "latent_residual.bed.gz.tbi",
        "latent_residual.phenotype_groups.txt",
        "cross_modality.bed.gz",
        "cross_modality.bed.gz.tbi",
        "cross_modality.phenotype_groups.txt",
    ]
    ds["laddr"] = {"enabled": False}
    assert "latent_residual.bed.gz" not in cfg.dataset_targets(ds)


@pytest.fixture
def validate():
    utils = pytest.importorskip("snakemake.utils")
    return lambda c: utils.validate(c, str(SCHEMA))


def test_schema_accepts_normalised_legacy_config(validate):
    out = cfg.normalise(legacy_config())
    validate(out)
    assert out["majiq_env"] == "majiq3"


def test_schema_defaults_land_per_dataset(validate):
    c = new_config()
    c["references"]["grch38"]["edit_sites_bed"] = "/ref/sites.bed"
    groups = {
        "RNA_editing": {
            "edit_sites_min_coverage": 20,
            "edit_sites_min_samples": 10,
            "files": ["RNA_editing.bed.gz"],
        },
        "intron_retention": {"files": ["intron_retention.bed.gz"]},
    }
    c["datasets"]["d1"]["modality_groups"] = copy.deepcopy(groups)
    c["datasets"]["d1"]["laddr"] = {"enabled": True}
    c["datasets"]["d2"] = copy.deepcopy(c["datasets"]["d1"])
    del c["datasets"]["d2"]["laddr"]
    out = cfg.normalise(c)
    validate(out)
    for name in ("d1", "d2"):
        groups = out["datasets"][name]["modality_groups"]
        assert groups["RNA_editing"]["edit_sites_min_samples_fraction"] == 0.1
        assert groups["RNA_editing"]["edit_sites_cluster"] is True
        assert groups["RNA_editing"]["edit_sites_correlation_threshold"] == 0.9
        assert groups["intron_retention"]["threads"] == 8
        assert groups["intron_retention"]["psicov_batch_size"] == 25
        assert groups["intron_retention"]["min_experiments"] == 0.01
    laddr = out["datasets"]["d1"]["laddr"]
    assert laddr == {
        "enabled": True,
        "regress_modalities": ["RNA_editing", "intron_retention"],
        "dataset_name": "pantry",
        "project_dir": "<pantry_intermediate>/laddr",
        "use_existing_bins": False,
        "min_samples_expressed": 0.5,
        "binning_method": "adaptive_diffvar",
        "batch_size": 40,
        "max_bin_width": 1024,
        "adaptive_max_samples": 256,
        "adaptive_bins_per_gene": 256,
        "adaptive_min_mean_total_covg": 128,
        "adaptive_max_corr": 0.8,
        "model_var_explained_max": 0.8,
        "model_n_pcs_max": 16,
        "bam_coverage_threads": 4,
        "bam_coverage_bin_size": 1,
    }
    assert "laddr" not in out["datasets"]["d2"]


def test_schema_rejects_bad_types(validate):
    out = cfg.normalise(new_config(read_length="75"))
    with pytest.raises(Exception):
        validate(out)
    out = cfg.normalise(new_config(fastq_map="map.txt"))
    del out["datasets"]["d1"]["fastq_dir"]
    with pytest.raises(Exception, match="fastq_dir"):
        validate(out)


# check_annotation: upstream validate_reference's transcript-row check.
def _gtf(tmp_path, rows):
    path = tmp_path / "anno.gtf"
    path.write_text("".join("\t".join(["1", "src", f, "1", "2", ".", "+", ".", "x"]) + "\n" for f in rows))
    return path


def test_check_annotation_accepts_transcript_rows(tmp_path):
    path = _gtf(tmp_path, ["gene", "transcript", "exon"])
    path.write_text("#!comment\n" + path.read_text())
    cfg.check_annotation(path)


def test_check_annotation_rejects_missing_transcripts(tmp_path):
    path = _gtf(tmp_path, ["gene", "exon"])
    with pytest.raises(ValueError) as e:
        cfg.check_annotation(path)
    assert str(e.value) == (
        "Reference annotation file does not contain transcript entries "
        f"(at least not in the first 1000 lines): {path}"
    )


def test_check_annotation_looks_at_first_1000_rows_only(tmp_path):
    with pytest.raises(ValueError):
        cfg.check_annotation(_gtf(tmp_path, ["exon"] * 1000 + ["transcript"]))
    cfg.check_annotation(_gtf(tmp_path, ["exon"] * 999 + ["transcript"]))
