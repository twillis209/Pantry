rule install_fusion:
    """Install the bundled FUSION scripts and gcta binary for TWAS.

    FUSION.compute_weights.R sources utils/plink_utils.R relative to itself,
    so the utils/ layout is kept. The source cache (e.g. for a github()
    module) does not keep gcta's exec bit, so it is set here.
    """
    output:
        compute_weights="<pantry_resources>/fusion/FUSION.compute_weights.R",
        plink_utils="<pantry_resources>/fusion/utils/plink_utils.R",
        profile_wgt="<pantry_resources>/fusion/utils/FUSION.profile_wgt.R",
        gcta="<pantry_resources>/fusion/gcta_nr_robust",
    localrule: True
    params:
        compute_weights=script("fusion_twas/FUSION.compute_weights.R"),
        plink_utils=script("fusion_twas/utils/plink_utils.R"),
        profile_wgt=script("fusion_twas/utils/FUSION.profile_wgt.R"),
        gcta=script("fusion_twas/gcta_nr_robust"),
        utils_dir=lambda w, output: str(Path(output.plink_utils).parent),
    shell:
        """
        mkdir -p {params.utils_dir}
        cp {params.compute_weights} {output.compute_weights}
        cp {params.plink_utils} {output.plink_utils}
        cp {params.profile_wgt} {output.profile_wgt}
        cp {params.gcta} {output.gcta}
        chmod +x {output.gcta}
        """
