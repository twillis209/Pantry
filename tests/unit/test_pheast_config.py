import copy
import importlib.util
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
PHEAST = REPO / "pheast"
SCHEMA = PHEAST / "schemas" / "config.schema.yaml"

_spec = importlib.util.spec_from_file_location(
    "pantry_pheast_config", PHEAST / "lib" / "config.py"
)
cfg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cfg)

UPSTREAM_CROSS_MODALITY_MESSAGE = (
    'The "cross_modality" modality should only be used with analysis "qtl" '
    "for cross-modality cis-QTL mapping."
)


# Oracle: upstream qtl_filter_args with its process_config conversions, from
# git show 8da7787:pheast/steps/{qtl,setup}.smk.
def upstream_qtl_filter_args(config: dict) -> str:
    if "qtl_min_allele_count" not in config:
        config["qtl_min_allele_count"] = 10
    qtl_min_allele_count = int(config["qtl_min_allele_count"])
    qtl_maf_threshold = config.get("qtl_maf_threshold")
    if qtl_maf_threshold is not None:
        qtl_maf_threshold = float(qtl_maf_threshold)
    args = f"--min_allele_count {qtl_min_allele_count}"
    if qtl_maf_threshold is not None:
        args += f" --maf_threshold {qtl_maf_threshold}"
    return args


def legacy_config() -> dict:
    with open(PHEAST / "config.yml") as f:
        return yaml.safe_load(f)


def new_config(**dataset_over) -> dict:
    dataset = {
        "phenotype_dir": "pheno",
        "geno_prefix": "geno/prefix",
        "samples": ["s1", "s2"],
        "modalities": {"expression": {"grouped": False}},
        "analyses": {"qtl": {"files": ["{modality}.cis_qtl.txt.gz"]}},
    }
    dataset.update(dataset_over)
    return {"datasets": {"d1": dataset}}


def test_shim_wraps_legacy_config():
    legacy = legacy_config()
    legacy["intermediate_dir"] = "~/interm"
    legacy["fusion_dir"] = "/opt/fusion"
    legacy["gcta"] = "/opt/gcta64"
    before = copy.deepcopy(legacy)
    out = cfg.normalise(legacy)
    assert legacy == before  # input not mutated
    assert list(out["datasets"]) == ["default"]
    ds = out["datasets"]["default"]
    assert ds["samples"] == "input/samples.txt"
    assert "samples_file" not in ds
    assert ds["phenotype_dir"] == "input/phenotypes"
    assert ds["modalities"] == legacy["modalities"]
    assert ds["analyses"] == legacy["analyses"]
    for key in ("intermediate_dir", "fusion_dir", "gcta"):
        assert key not in ds
    assert out["fusion_dir"] == "/opt/fusion"
    assert out["gcta"] == "/opt/gcta64"
    interm = str(Path("~/interm").expanduser())
    assert out["legacy_pathvars"] == {
        "pantry_intermediate": interm,
        "pantry_resources": f"{interm}/resources",
    }


def test_shim_without_intermediate_dir_has_no_legacy_pathvars():
    assert cfg.normalise(legacy_config())["legacy_pathvars"] == {}
    assert cfg.normalise(new_config())["legacy_pathvars"] == {}


def test_cross_modality_with_twas_raises():
    c = new_config(
        modalities={"cross_modality": {"grouped": True}},
        analyses={
            "qtl": {"files": ["{modality}.cis_qtl.txt.gz"]},
            "twas": {"files": ["{modality}.tar.bz2"]},
        },
    )
    with pytest.raises(ValueError) as e:
        cfg.normalise(c)
    assert UPSTREAM_CROSS_MODALITY_MESSAGE in str(e.value)


def test_cross_modality_with_trans_qtl_raises():
    c = new_config(
        modalities={"cross_modality": {"grouped": True}},
        analyses={"qtl": {"files": ["{modality}.trans_qtl.txt.gz"]}},
    )
    with pytest.raises(ValueError, match="cross_modality"):
        cfg.normalise(c)


def test_cross_modality_with_cis_qtl_only_ok():
    c = new_config(modalities={"cross_modality": {"grouped": True}})
    cfg.normalise(c)


def test_qtl_filter_args_default_matches_upstream():
    assert cfg.qtl_filter_args({"qtl_min_allele_count": 10}) == "--min_allele_count 10"
    assert cfg.qtl_filter_args(
        {"qtl_min_allele_count": 10}
    ) == upstream_qtl_filter_args({})
    ds = {"qtl_min_allele_count": 10, "qtl_maf_threshold": None}
    assert cfg.qtl_filter_args(ds) == "--min_allele_count 10"


@pytest.mark.parametrize(
    "ds",
    [
        {"qtl_min_allele_count": 5, "qtl_maf_threshold": 0.05},
        {"qtl_min_allele_count": 0, "qtl_maf_threshold": 1},
        {"qtl_min_allele_count": 20},
    ],
)
def test_qtl_filter_args_matches_upstream_without_window(ds):
    assert cfg.qtl_filter_args(ds) == upstream_qtl_filter_args(copy.deepcopy(ds))


def test_qtl_filter_args_with_maf_and_window():
    ds = {
        "qtl_min_allele_count": 10,
        "qtl_maf_threshold": 0.05,
        "qtl_cis_window": 500000,
    }
    assert (
        cfg.qtl_filter_args(ds)
        == "--min_allele_count 10 --maf_threshold 0.05 --window 500000"
    )
    ds = {"qtl_min_allele_count": 10, "qtl_cis_window": None}
    assert cfg.qtl_filter_args(ds) == "--min_allele_count 10"


def test_defaults_deep_merged():
    c = new_config(analyses={"qtl": {"files": ["{modality}.cis_qtl.txt.gz"]}})
    c["datasets"]["d2"] = new_config(qtl_min_allele_count=5)["datasets"]["d1"]
    del c["datasets"]["d2"]["modalities"]
    c["defaults"] = {
        "qtl_min_allele_count": 20,
        "twas_threads": 2,
        "modalities": {"stability": {"grouped": False}},
        "analyses": {"heritability": {"files": ["{modality}.hsq.tsv"]}},
    }
    before = copy.deepcopy(c)
    out = cfg.normalise(c)
    assert c == before
    d1, d2 = out["datasets"]["d1"], out["datasets"]["d2"]
    assert d1["qtl_min_allele_count"] == 20
    assert d2["qtl_min_allele_count"] == 5  # dataset wins
    assert d1["twas_threads"] == d2["twas_threads"] == 2
    assert d1["modalities"] == {
        "stability": {"grouped": False},
        "expression": {"grouped": False},
    }
    assert d2["modalities"] == {"stability": {"grouped": False}}
    assert d1["analyses"] == {
        "heritability": {"files": ["{modality}.hsq.tsv"]},
        "qtl": {"files": ["{modality}.cis_qtl.txt.gz"]},
    }
    # Merged values are copies, not shared with defaults or other datasets.
    d1["modalities"]["stability"]["grouped"] = True
    assert d2["modalities"]["stability"]["grouped"] is False
    assert c["defaults"]["modalities"]["stability"]["grouped"] is False


def test_dataset_lists_replace_default_lists():
    c = new_config()
    c["defaults"] = {"analyses": {"qtl": {"files": ["{modality}.trans_qtl.txt.gz"]}}}
    out = cfg.normalise(c)
    assert out["datasets"]["d1"]["analyses"]["qtl"]["files"] == [
        "{modality}.cis_qtl.txt.gz"
    ]


@pytest.mark.parametrize("bad", ["a/b", "", ".", "..", 3])
def test_dataset_name_must_be_path_component(bad):
    c = new_config()
    c["datasets"] = {bad: c["datasets"]["d1"]}
    with pytest.raises(ValueError, match="Dataset name"):
        cfg.normalise(c)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda c: c.update(datasets=[]),
        lambda c: c.update(defaults=[]),
        lambda c: c["datasets"].update(d1="x"),
        lambda c: c["datasets"]["d1"].update(modalities=["expression"]),
        lambda c: c["datasets"]["d1"].update(analyses=["qtl"]),
        lambda c: c["datasets"]["d1"].update(analyses={"qtl": "x"}),
        lambda c: c["datasets"]["d1"].update(analyses={"qtl": {"files": "x"}}),
    ],
)
def test_malformed_config_raises_value_error(mutate):
    c = new_config(modalities={"cross_modality": {"grouped": True}})
    mutate(c)
    with pytest.raises(ValueError):
        cfg.normalise(c)


def test_empty_config_raises():
    with pytest.raises(ValueError, match="No config"):
        cfg.normalise({})


def test_load_samples(tmp_path):
    f = tmp_path / "samples.txt"
    f.write_text("s1\ns2\n")
    assert cfg.load_samples(str(f)) == ["s1", "s2"]
    assert cfg.load_samples(["a", "b"]) == ["a", "b"]


@pytest.fixture
def validate():
    utils = pytest.importorskip("snakemake.utils")
    return lambda c: utils.validate(c, str(SCHEMA))


def test_schema_accepts_normalised_legacy_config(validate):
    out = cfg.normalise(legacy_config())
    validate(out)
    ds = out["datasets"]["default"]
    assert ds["qtl_min_allele_count"] == 10
    assert ds["twas_models_per_job"] == 200
    assert ds["twas_threads"] == 4


def test_schema_defaults_land_per_dataset(validate):
    c = new_config()
    c["datasets"]["d2"] = copy.deepcopy(c["datasets"]["d1"])
    c["datasets"]["d2"]["twas_threads"] = 8
    out = cfg.normalise(c)
    validate(out)
    for name in ("d1", "d2"):
        ds = out["datasets"][name]
        assert ds["qtl_min_allele_count"] == 10
        assert ds["twas_models_per_job"] == 200
        assert "qtl_cis_window" not in ds  # run_tensorqtl.py's default applies
        assert "qtl_maf_threshold" not in ds
    assert out["datasets"]["d1"]["twas_threads"] == 4
    assert out["datasets"]["d2"]["twas_threads"] == 8
    assert "fusion_dir" not in out and "gcta" not in out


@pytest.mark.parametrize(
    "over",
    [
        {"twas_models_per_job": 0},
        {"twas_threads": 0},
        {"qtl_min_allele_count": "10"},
        {"qtl_cis_window": 1.5},
        {"modalities": {"expression": {}}},
        {"analyses": {"eqtl": {"files": ["x"]}}},
    ],
)
def test_schema_rejects_bad_values(validate, over):
    out = cfg.normalise(new_config(**over))
    with pytest.raises(Exception):
        validate(out)


def test_schema_requires_dataset_fields(validate):
    for key in ("phenotype_dir", "geno_prefix", "samples", "modalities", "analyses"):
        c = new_config()
        del c["datasets"]["d1"][key]
        with pytest.raises(Exception, match=key):
            validate(cfg.normalise(c))


def test_qtl_filter_args_without_window_for_genotype_qc():
    # genotype_qc.py has no --window option.
    ds = {"qtl_min_allele_count": 10, "qtl_maf_threshold": 0.05, "qtl_cis_window": 5}
    assert (
        cfg.qtl_filter_args(ds, window=False)
        == "--min_allele_count 10 --maf_threshold 0.05"
    )
