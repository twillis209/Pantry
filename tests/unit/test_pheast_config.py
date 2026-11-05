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
    """The old-format pheast config (pheast/config.yml before conversion)."""
    with open(REPO / "tests" / "configs" / "pheast_baseline.yml") as f:
        return yaml.safe_load(f)


def test_converted_config_matches_legacy_config():
    with open(PHEAST / "config.yml") as f:
        converted = yaml.safe_load(f)
    assert "datasets" in converted
    assert cfg.normalise(converted) == cfg.normalise(legacy_config())


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
    assert "intermediate_dir" not in ds
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


def test_legacy_keys_beside_datasets_rejected():
    # e.g. an old-format --configfile merged over the new-format config.yml
    mixed = new_config()
    mixed.update(geno_prefix="other", qtl_min_allele_count=20)
    with pytest.raises(ValueError, match="geno_prefix, qtl_min_allele_count"):
        cfg.normalise(mixed)
    # Other top-level keys stay allowed.
    cfg.normalise({**new_config(), "pathvars": {}})


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
        assert ds["qtl_seed"] == 0
    assert out["datasets"]["d1"]["twas_threads"] == 4
    assert out["datasets"]["d2"]["twas_threads"] == 8


@pytest.mark.parametrize(
    "over",
    [
        {"twas_models_per_job": 0},
        {"twas_threads": 0},
        {"qtl_min_allele_count": "10"},
        {"qtl_cis_window": 1.5},
        {"qtl_seed": 1.5},
        {"qtl_seed": "0"},
        {"qtl_seed": -1},
        {"qtl_seed": 2**32},
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


def test_dataset_targets_match_upstream_outputs():
    ds = cfg.normalise(legacy_config())["datasets"]["default"]
    # Oracle: upstream setup.smk's outputs list (paths relative to output/).
    upstream = []
    for analysis, params in ds["analyses"].items():
        if analysis == "qtl":
            upstream.append(f"{analysis}/genotype_qc.tsv")
        for f in params["files"]:
            upstream.extend(f"{analysis}/{f}".format(modality=m) for m in ds["modalities"])
    targets = cfg.dataset_targets(ds)
    assert targets == upstream
    assert "qtl/genotype_qc.tsv" in targets
    assert "twas/expression.tar.bz2" in targets
    assert len(targets) == 1 + 4 * 8


def test_dataset_targets_without_qtl_has_no_genotype_qc():
    ds = {
        "modalities": {"expression": {"grouped": False}},
        "analyses": {"heritability": {"files": ["{modality}.hsq.tsv"]}},
    }
    assert cfg.dataset_targets(ds) == ["heritability/expression.hsq.tsv"]


def write_fam(path, ids):
    path.write_text("".join(f"0 {i} 0 0 0 -9\n" for i in ids))


def test_check_genotype_samples_passes_when_all_present(tmp_path):
    fam = tmp_path / "g.fam"
    write_fam(fam, ["a", "b", "c"])
    cfg.check_genotype_samples(["a", "c"], fam)


def test_check_genotype_samples_upstream_message(tmp_path):
    fam = tmp_path / "g.fam"
    # Tab-separated, as in the test data .fam.
    fam.write_text("0\ta\t0\t0\t0\t-9\n0\t007\t0\t0\t0\t-9\n")
    with pytest.raises(Exception) as e:
        cfg.check_genotype_samples(["a", "b", "007", "c"], fam)
    assert str(e.value) == "Samples missing from genotypes: ['b', 'c']"


def write_bed_gz(path, rows):
    import gzip

    with gzip.open(path, "wt") as f:
        f.write("#chr\tstart\tend\tphenotype_id\ts1\n")
        for chrom, i in rows:
            f.write(f"{chrom}\t{i}\t{i + 1}\tp{i}\t0.5\n")


def test_phenotype_summary_round_trip(tmp_path):
    bed = tmp_path / "m.bed.gz"
    write_bed_gz(bed, [("chr2", 1), ("chr10", 2), ("chr2", 3), ("1", 4)])
    out = tmp_path / "summary" / "m.tsv"
    cfg.write_phenotype_summary(bed, out)
    assert out.read_text() == "chroms\tn_models\nchr2,chr10,1\t4\n"
    chroms, n_models = cfg.read_phenotype_summary(out)
    assert chroms == ["chr2", "chr10", "1"]  # file order, strings
    assert n_models == 4


@pytest.mark.skipif(
    not list((PHEAST / "input" / "phenotypes").glob("*.bed.gz")),
    reason="pheast test data not downloaded",
)
def test_phenotype_summary_matches_upstream_reads():
    """Oracle: upstream's parse-time chroms read and twas_model_count."""
    import gzip

    import pandas as pd

    for bed in sorted((PHEAST / "input" / "phenotypes").glob("*.bed.gz")):
        chroms, n_models = cfg.phenotype_summary(bed)
        up = pd.read_csv(bed, sep="\t", usecols=[0], dtype=str)
        assert chroms == up.iloc[:, 0].unique().tolist()
        with gzip.open(bed, "rt") as handle:
            assert n_models == sum(1 for _ in handle) - 1


def test_phenotype_summary_empty_bed(tmp_path):
    bed = tmp_path / "m.bed.gz"
    write_bed_gz(bed, [])
    out = tmp_path / "m.tsv"
    cfg.write_phenotype_summary(bed, out)
    assert cfg.read_phenotype_summary(out) == ([], 0)


REMOVED_KEYS = ("fusion_dir", "gcta", "gcta_url", "gcta_sha256")


@pytest.mark.parametrize("key", REMOVED_KEYS)
def test_removed_fusion_keys_rejected(key):
    # Top level of the new format.
    with pytest.raises(ValueError, match=f"{key} are no longer supported"):
        cfg.normalise({**new_config(), key: "/x"})
    # Under defaults and under a dataset.
    c = new_config()
    c["defaults"] = {key: "/x"}
    with pytest.raises(ValueError, match=f"defaults.{key} "):
        cfg.normalise(c)
    c = new_config()
    c["datasets"]["d1"][key] = None
    with pytest.raises(ValueError, match=f"datasets.d1.{key} "):
        cfg.normalise(c)
    # Top level of the legacy format.
    with pytest.raises(ValueError, match="pantry-fusion"):
        cfg.normalise({**legacy_config(), key: "/x"})


def test_covariates_keeps_plink_utils_beside_its_script():
    # covariates.R sources fusion_twas/utils/plink_utils.R relative to itself;
    # the covariates rule's script() param is what puts it in the source cache.
    r_script = (PHEAST / "scripts" / "covariates.R").read_text()
    assert '"fusion_twas", "utils", "plink_utils.R"' in r_script
    rule = (PHEAST / "steps" / "covariates.smk").read_text()
    assert 'script("fusion_twas/utils/plink_utils.R")' in rule


def test_qtl_seed_null_survives_validation(validate):
    # null disables seeding; the schema default must not replace it.
    out = cfg.normalise(new_config(qtl_seed=None))
    validate(out)
    assert out["datasets"]["d1"]["qtl_seed"] is None
    out = cfg.normalise(new_config(qtl_seed=2**32 - 1))
    validate(out)


@pytest.mark.parametrize(
    "ds, arg",
    [
        ({"qtl_seed": 0}, "--seed 0"),
        ({"qtl_seed": 7}, "--seed 7"),
        ({"qtl_seed": None}, ""),
        ({}, ""),
    ],
)
def test_qtl_seed_arg(ds, arg):
    assert cfg.qtl_seed_arg(ds) == arg


def test_qtl_seed_not_in_filter_args():
    # genotype_qc and trans take the filter args but no seed.
    ds = {"qtl_min_allele_count": 10, "qtl_seed": 3}
    assert cfg.qtl_filter_args(ds) == "--min_allele_count 10"
    assert cfg.qtl_filter_args(ds, window=False) == "--min_allele_count 10"


def test_legacy_top_level_qtl_seed(validate):
    c = legacy_config()
    c["qtl_seed"] = 3
    out = cfg.normalise(c)
    validate(out)
    assert out["datasets"]["default"]["qtl_seed"] == 3
    c["qtl_seed"] = None
    out = cfg.normalise(c)
    validate(out)
    assert out["datasets"]["default"]["qtl_seed"] is None
    # Unset in the legacy config: the default applies.
    out = cfg.normalise(legacy_config())
    validate(out)
    assert out["datasets"]["default"]["qtl_seed"] == 0


def test_top_level_qtl_seed_beside_datasets_rejected():
    c = new_config()
    c["qtl_seed"] = 1
    with pytest.raises(ValueError, match="qtl_seed"):
        cfg.normalise(c)
