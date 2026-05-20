def decontam_input(wc):
    if config["decontam"]["run"]:
        return {"r1": os.path.join(PROC, f"trimmed/{wc.sample}_R1.fastq.gz"),
                "r2": os.path.join(PROC, f"trimmed/{wc.sample}_R2.fastq.gz")}
    return {"r1": os.path.join(PROC, f"trimmed/{wc.sample}_R1.fastq.gz"),
            "r2": os.path.join(PROC, f"trimmed/{wc.sample}_R2.fastq.gz")}

if config["decontam"]["run"]:
    rule decontam:
        input:
            r1  = os.path.join(PROC, "trimmed/{sample}_R1.fastq.gz"),
            r2  = os.path.join(PROC, "trimmed/{sample}_R2.fastq.gz"),
            ref = config["decontam"]["reference"],
        output:
            r1  = os.path.join(PROC, "decontam/{sample}_R1.fastq.gz"),
            r2  = os.path.join(PROC, "decontam/{sample}_R2.fastq.gz"),
            bam = temp(os.path.join(TMP, "decontam/{sample}_contam.bam")),
        log:    os.path.join(LOGS, "decontam/{sample}.log")
        threads: THREADS_H
        conda:  os.path.expanduser("~/lab/software/envs/general")
        shell:
            "mkdir -p $(dirname {output.r1}) $(dirname {output.bam}) && "
            "minimap2 -ax sr -t {threads} {input.ref} {input.r1} {input.r2} 2>>{log} | "
            "samtools view -f 12 -F 256 -b -o {output.bam} 2>>{log} && "
            "samtools sort -n {output.bam} 2>>{log} | "
            "samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null 2>>{log}"

def clean_reads(sample, read):
    if config["decontam"]["run"]:
        return os.path.join(PROC, f"decontam/{sample}_{read}.fastq.gz")
    return os.path.join(PROC, f"trimmed/{sample}_{read}.fastq.gz")
else:
    def clean_reads(sample, read):
        return os.path.join(PROC, f"trimmed/{sample}_{read}.fastq.gz")
