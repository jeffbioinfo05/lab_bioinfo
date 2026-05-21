# decontam.smk — remoção de múltiplos hospedeiros em sequência
# config.yaml:
#   decontam:
#     run: true
#     hosts:
#       - "/mnt/bioinfo/databases/genomes/human/GRCh38.fa"
#       - "/mnt/bioinfo/databases/genomes/cat/Felis_catus.fna"

def get_hosts():
    return config["decontam"].get("hosts", [])

if config["decontam"]["run"] and len(get_hosts()) > 0:

    rule decontam:
        input:
            r1 = os.path.join(PROC, "trimmed/{sample}_R1.fastq.gz"),
            r2 = os.path.join(PROC, "trimmed/{sample}_R2.fastq.gz"),
        output:
            r1     = os.path.join(PROC, "decontam/{sample}_R1.fastq.gz"),
            r2     = os.path.join(PROC, "decontam/{sample}_R2.fastq.gz"),
            report = os.path.join(RESULTS, "qc/decontam/{sample}_decontam_report.tsv"),
        params:
            hosts   = get_hosts(),
            tmp_dir = os.path.join(TMP, "decontam/{sample}"),
        log:     os.path.join(LOGS, "decontam/{sample}.log")
        threads: THREADS_H
        conda:   os.path.expanduser("~/lab/software/envs/general")
        run:
            import subprocess, os, gzip

            os.makedirs(params.tmp_dir, exist_ok=True)
            os.makedirs(os.path.dirname(output.r1), exist_ok=True)
            os.makedirs(os.path.dirname(output.report), exist_ok=True)

            def count_reads(fastq_gz):
                result = subprocess.run(
                    ["bash", "-c", f"zcat {fastq_gz} | wc -l"],
                    capture_output=True, text=True)
                return int(result.stdout.strip()) // 4

            current_r1 = input.r1
            current_r2 = input.r2

            initial = count_reads(current_r1)
            report_rows = [("sample", "step", "reads", "removed", "removed_pct"),
                           (wildcards.sample, "raw", initial, 0, "0.00")]
            prev = initial

            for i, host_ref in enumerate(params.hosts):
                host_name = os.path.basename(host_ref)
                for ext in [".fa", ".fna", ".fasta", ".gz"]:
                    host_name = host_name.replace(ext, "")

                out_r1 = os.path.join(params.tmp_dir, f"step{i+1}_R1.fastq.gz")
                out_r2 = os.path.join(params.tmp_dir, f"step{i+1}_R2.fastq.gz")
                bam    = os.path.join(params.tmp_dir, f"step{i+1}_host.bam")

                with open(log[0], "a") as logf:
                    logf.write(f"\n=== Removendo hospedeiro: {host_name} ===\n")

                cmd = (
                    f"minimap2 -ax sr -t {threads} {host_ref} "
                    f"{current_r1} {current_r2} 2>>{log[0]} | "
                    f"samtools view -f 12 -F 256 -b -o {bam} 2>>{log[0]} && "
                    f"samtools sort -n -@ {threads} {bam} 2>>{log[0]} | "
                    f"samtools fastq "
                    f"-1 {out_r1} -2 {out_r2} "
                    f"-0 /dev/null -s /dev/null 2>>{log[0]}"
                )
                subprocess.run(cmd, shell=True, check=True)

                curr = count_reads(out_r1)
                removed = prev - curr
                pct = f"{removed / prev * 100:.2f}" if prev > 0 else "0.00"
                report_rows.append((wildcards.sample, host_name, curr, removed, pct))

                current_r1 = out_r1
                current_r2 = out_r2
                prev = curr

            subprocess.run(f"cp {current_r1} {output.r1}", shell=True, check=True)
            subprocess.run(f"cp {current_r2} {output.r2}", shell=True, check=True)

            with open(output.report, "w") as f:
                for row in report_rows:
                    f.write("\t".join(str(x) for x in row) + "\n")

    def clean_reads(sample, read):
        return os.path.join(PROC, f"decontam/{sample}_{read}.fastq.gz")

else:
    def clean_reads(sample, read):
        return os.path.join(PROC, f"trimmed/{sample}_{read}.fastq.gz")
