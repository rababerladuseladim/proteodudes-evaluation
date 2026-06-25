rule calculate_qc_metrics_for_taxonomic_profiling_sample:
    input:
        ground_truth=lambda wc: samples.loc[wc.sample_name, "ground_truth"],
        diamond_result="results/diamond/sample_{sample_name}-lineage.tsv",
        unipept_result="results/unipept/sample_{sample_name}.csv",
        megadudes_results=expand(
            "results/megadudes/{method}/sample_{{sample_name}}.out",
            method=config["alignment_methods"],
        ),
    log:
        "logs/qc/calculate_qc_metrics_for_taxonomic_profiling_sample-sample_{sample_name}.txt",
    output:
        report(
            "results/qc/qc-sample_{sample_name}.tsv",
            category="qc",
            subcategory="taxonomic_profiling",
        ),
    script:
        "../scripts/calculate_qc_metrics_for_taxonomic_profiling.py"


rule calculate_qc_metrics_for_taxonomic_profiling_simulation:
    input:
        ground_truth="results/simulation/sample_taxons_lineage_{repeat}.tsv",
        diamond_result="results/diamond/simulated_peptides_{repeat}_with_{percentage}_percent_noise-lineage.tsv",
        unipept_result="results/unipept/simulated_peptides_{repeat}_with_{percentage}_percent_noise.csv",
        megadudes_results=expand(
            "results/megadudes/{method}/simulated_peptides_{{repeat}}_with_{{percentage}}_percent_noise.out",
            method=config["alignment_methods"],
        ),
    log:
        "logs/qc/calculate_qc_metrics_for_taxonomic_profiling_simulation-simulated_peptides_{repeat}_with_{percentage}_percent_noise.txt",
    output:
        report(
            "results/qc/qc-simulated_peptides_{repeat}_with_{percentage}_percent_noise.tsv",
            category="qc",
            subcategory="taxonomic_profiling",
        ),
    script:
        "../scripts/calculate_qc_metrics_for_taxonomic_profiling.py"


rule plot_qc_taxonomic_profiling:
    input:
        qc_data="results/qc/qc-{sample_name}.tsv",
    log:
        "logs/qc/taxonomic_profiling/plot_qc_taxonomic_profiling-{sample_name}.txt",
    output:
        report(
            "plots/taxonomic_profiling/qc-{sample_name}.svg",
            category="qc",
            subcategory="taxonomic_profiling",
        ),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_qc_taxonomic_profiling.py"
