CHECKM2_SIF = os.path.expanduser("~/lab/software/containers/checkm2.sif")

rule quast:
    input:  os.path.join(RESULTS, "assembly/{sample}/assembly.fasta")
    output: os.path.join(RESULTS, "evaluate/quast/{sample}/report.tsv")
    params: outdir = os.path.join(RESULTS, "evaluate/quast/{sample}")
    log:    os.path.join(LOGS, "quast/{sample}.log")
    threads: THREADS
    conda:  os.path.expanduser("~/lab/software/envs/assembly_core")
    shell:  "quast {input} -o {params.outdir} -t {threads} > {log} 2>&1"

rule busco:
    input:  os.path.join(RESULTS, "assembly/{sample}/assembly.fasta")
    output: os.path.join(RESULTS, "evaluate/busco/{sample}/short_summary.txt")
    params:
        outdir  = os.path.join(RESULTS, "evaluate/busco/{sample}"),
        lineage = config["databases"]["busco_lineage"],
        db      = config["databases"]["busco"],
    log:    os.path.join(LOGS, "busco/{sample}.log")
    threads: THREADS
    conda:  os.path.expanduser("~/lab/software/envs/assembly_core")
    shell:
        "busco -i {input} -o {params.outdir} -l {params.lineage} "
        "--offline --download_path {params.db} -m genome "
        "-c {threads} --force > {log} 2>&1 && "
        "cp {params.outdir}/short_summary*.txt {output}"

rule checkm2:
    input:  os.path.join(RESULTS, "assembly/{sample}/assembly.fasta")
    output: os.path.join(RESULTS, "evaluate/checkm2/{sample}/quality_report.tsv")
    params:
        outdir  = os.path.join(RESULTS, "evaluate/checkm2/{sample}"),
        db      = config["databases"]["checkm2"],
    log:    os.path.join(LOGS, "checkm2/{sample}.log")
    threads: THREADS_H
    container: CHECKM2_SIF
    shell:
        "checkm2 predict --input {input} --output-directory {params.outdir} "
        "--database_path {params.db}/uniref100.KO.1.dmnd "
        "--threads {threads} --force > {log} 2>&1"
