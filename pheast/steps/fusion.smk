import shutil


rule install_fusion:
    """Install the bundled FUSION scripts and gcta binary for TWAS.

    FUSION.compute_weights.R sources utils/plink_utils.R relative to itself,
    so the utils/ layout is kept. Configured fusion_dir / gcta replace the
    scripts / the binary, which are then not installed. gcta is resolved only
    here, at job time: copied from the bundled file if the module is loaded
    from a local path, else downloaded from gcta_url and checked against
    gcta_sha256 (the source cache cannot fetch binaries from git); either way
    made executable.
    """
    output:
        **pantry_config.install_fusion_outputs(CONFIG),
    localrule: True
    params:
        compute_weights=script("fusion_twas/FUSION.compute_weights.R"),
        plink_utils=script("fusion_twas/utils/plink_utils.R"),
        profile_wgt=script("fusion_twas/utils/FUSION.profile_wgt.R"),
        gcta_url=CONFIG["gcta_url"],
        gcta_sha256=CONFIG["gcta_sha256"],
    run:
        for key in ("compute_weights", "plink_utils", "profile_wgt"):
            if key in output.keys():
                Path(getattr(output, key)).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(getattr(params, key), getattr(output, key))
        if "gcta" in output.keys():
            pantry_config.install_gcta(
                output.gcta,
                src=(
                    workflow.source_path("../scripts/fusion_twas/gcta_nr_robust")
                    if LOCAL_SOURCE
                    else None
                ),
                url=params.gcta_url,
                sha256=params.gcta_sha256,
            )
