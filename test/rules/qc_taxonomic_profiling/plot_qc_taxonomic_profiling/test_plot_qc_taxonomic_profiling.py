from test.common import snakemake_run, check_output


def test_plot_qc_taxonomic_profiling(tmpdir, workflow_path, prepared_workdir):
    targets = ["plots/taxonomic_profiling/qc-foo.svg"]
    snakefile = workflow_path / "workflow/rules/qc_taxonomic_profiling.smk"
    configfile = prepared_workdir.workdir / "config" / "config.yaml"
    conda_prefix = workflow_path / ".snakemake" / "conda"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
            "--use-conda",
            "--conda-prefix",
            conda_prefix.as_posix(),
            "--configfile",
            configfile.as_posix(),
        ]
    )

    check_output(prepared_workdir, mode="presence")
