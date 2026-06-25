from test.common import snakemake_run, check_output


def test_calculate_qc_metrics_for_taxonomic_profiling_simulation(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/qc/qc-simulated_peptides_0_with_1_percent_noise.tsv"]
    snakefile = workflow_path / "workflow/rules/qc_taxonomic_profiling.smk"
    configfile = prepared_workdir.workdir / "config" / "config.yaml"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
            "--configfile",
            configfile.as_posix(),
        ]
    )

    check_output(prepared_workdir, mode="text")
