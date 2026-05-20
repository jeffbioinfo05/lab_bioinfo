TETOOLS_SIF = os.path.expanduser("~/lab/software/containers/tetools.sif")

rule repeatmodeler:
    input:  os.path.join(RESULTS, "assembly/{sample}/assembly.fasta")
    output: os.path.join(RESULTS, "repeats/{sample}/repeat_library.fa")
    params:
        outdir = os.path.join(RESULTS, "repeats/{sample}"),
        run_rm = config["repeats"]["run_repeatmodeler"],
        lib    = config["repeats"].get("repeat_library", ""),
    log:    os.path.join(LOGS, "repeatmodeler/{sample}.log")
    threads: THREADS_H
    container: TETOOLS_SIF
    run:
        import subprocess, shutil, os
        os.makedirs(params.outdir, exist_ok=True)
        if params.run_rm:
            cmd = (f"cd {params.outdir} && "
                   f"BuildDatabase -name {wildcards.sample}_db {input} >> {log} 2>&1 && "
                   f"RepeatModeler -database {wildcards.sample}_db -pa {threads} >> {log} 2>&1 && "
                   f"cat {params.outdir}/{wildcards.sample}_db-families.fa > {output}")
            subprocess.run(cmd, shell=True, check=True)
        else:
            shutil.copy(params.lib, output[0])

rule repeatmasker:
    input:
        fasta = os.path.join(RESULTS, "assembly/{sample}/assembly.fasta"),
        lib   = os.path.join(RESULTS, "repeats/{sample}/repeat_library.fa"),
    output: os.path.join(RESULTS, "repeats/{sample}/assembly.masked.fasta")
    params:
        outdir    = os.path.join(RESULTS, "repeats/{sample}"),
        soft_mask = "-xsmall" if config["repeats"]["soft_mask"] else "",
    log:    os.path.join(LOGS, "repeatmasker/{sample}.log")
    threads: THREADS_H
    container: TETOOLS_SIF
    shell:
        "RepeatMasker {params.soft_mask} -lib {input.lib} "
        "-pa {threads} -dir {params.outdir} {input.fasta} > {log} 2>&1 && "
        "mv {params.outdir}/$(basename {input.fasta}).masked {output}"
