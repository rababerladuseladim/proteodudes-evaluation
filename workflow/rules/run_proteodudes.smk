import os.path


rule build_proteodudes_db:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        uniprot_fastas=config["query_dbs"],
        ncbi_nodes="resources/ncbi/nodes.dmp",
        ncbi_names="resources/ncbi/names.dmp",
    log:
        stdout="logs/proteodudes/build_proteodudes_db-stdout.txt",
        stderr="logs/proteodudes/build_proteodudes_db-stderr.txt",
    benchmark:
        "benchmarks/build_proteodudes_db-benchmark.txt"
    output:
        dudes_db="results/proteodudes/proc/dudes_db.npz",
    conda:
        "../envs/proteodudes.yaml"
    threads: 99
    params:
        db_base_name=lambda wildcards, output: os.path.splitext(output["dudes_db"])[0],
    shell:
        """
        dudesdb -m up \
        -f {input.uniprot_fastas} \
        -g {input.idmap} \
        -n {input.ncbi_nodes} \
        -a {input.ncbi_names} \
        -o {params.db_base_name} \
        -t {threads} > {log.stdout} 2> {log.stderr}
        """


rule run_proteodudes:
    input:
        dudes_db="results/proteodudes/proc/dudes_db.npz",
        custom_blast_format_file="results/{method}/{sample}.tsv",
    log:
        stdout="logs/proteodudes/run_proteodudes-{method}-{sample}-stdout.txt",
        stderr="logs/proteodudes/run_proteodudes-{method}-{sample}-stderr.txt",
    benchmark:
        "benchmarks/run_proteodudes-{method}-{sample}-benchmark.txt"
    output:
        result=report(
            "results/proteodudes/{method}/{sample}.out", category="proteodudes"
        ),
    conda:
        "../envs/proteodudes.yaml"
    threads: 99
    params:
        result_wo_ext=lambda wildcards, output: os.path.splitext(output.result)[0],
    shell:
        """
        dudes \
        -c {input.custom_blast_format_file} \
        -d {input.dudes_db} \
        -t {threads} \
        -o {params.result_wo_ext}\
        --debug > {log.stdout} 2> {log.stderr}
        """
