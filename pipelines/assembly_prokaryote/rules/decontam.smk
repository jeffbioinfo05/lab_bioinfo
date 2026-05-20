# decontam.smk — remoção de múltiplos hospedeiros em sequência
# Configuração no config.yaml:
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
            r1 = os.path.join(PROC, "decontam/{sample}_R1.fastq.gz"),
            r2 = os.path.join(PROC, "decontam/{sample}_R2.fastq.gz"),
            report = os.path.join(RESULTS, "qc/decontam/{sample}_decontam_report.tsv"),
        params:
            hosts   = get_hosts(),
            tmp_dir = os.path.join(TMP, "decontam/{sample}"),
        log:     os.path.join(LOGS, "decontam/{sample}.log")
        threads: THREADS_H
        conda:   os.path.expanduser("~/lab/software/envs/general")
        run:
            import subprocess, os

            os.makedirs(params.tmp_dir, exist_ok=True)
            os.makedirs(os.path.dirname(output.r1), exist_ok=True)
            os.makedirs(os.path.dirname(output.report), exist_ok=True)

            # Leitura inicial
            current_r1 = input.r1
            current_r2 = input.r2

            counts = []

            # Conta reads iniciais
            result = subprocess.run(
                f"zcat {current_r1} | wc -l",
                shell=True, capture_output=True, text=True)
            initial = int(result.stdout.strip()) // 4
            counts.append(("raw", initial))

            # Itera sobre cada hospedeiro em sequência
            for i, host_ref in enumerate(params.hosts):
                host_name = os.path.basename(host_ref).replace(".fa","").replace(".fna","")
                out_r1 = os.path.join(params.tmp_dir, f"step{i+1}_R1.fastq.gz")
                out_r2 = os.path.join(params.tmp_dir, f"step{i+1}_R2.fastq.gz")
                bam    = os.path.join(params.tmp_dir, f"step{i+1}_host.bam")

                cmd = (
                    f"minimap2 -ax sr -t {threads} {host_ref} "
                    f"{current_r1} {current_r2} 2>>{log} | "
                    f"samtools view -f 12 -F 256 -b -o {bam} 2>>{log} && "
                    f"samtools sort -n -@ {threads} {bam} 2>>{log} | "
                    f"samtools fastq "
                    f"-1 {out_r1} -2 {out_r2} "
                    f"-0 /dev/null -s /dev/null 2>>{log}"
                )
                subprocess.run(cmd, shell=True, check=True)

                # Conta reads restantes após este hospedeiro
                result = subprocess.run(
                    f"zcat {out_r1} | wc -l",
                    shell=True, capture_output=True, text=True)
                remaining = int(result.stdout.strip()) // 4
                counts.append((host_name, remaining))

                current_r1 = out_r1
                current_r2 = out_r2

            # Copia resultado final
            subprocess.run(f"cp {current_r1} {output.r1}", shell=True, check=True)
            subprocess.run(f"cp {current_r2} {output.r2}", shell=True, check=True)

            # Gera relatório
            with open(output.report, "w") as f:
                f.write("sample\tstep\treads\tremoved\tremoved_pct\n")
                prev = counts[0][1]
                for step, count in counts:
                    removed = prev - count if step != "raw" else 0
                    pct = round(removed / prev * 100, 2) if prev > 0 and step != "raw" else 0
                    f.write(f"{wildcards.sample}\t{step}\t{count}\t{removed}\t{pct}\n")
                    if step != "raw":
                        prev = count

    def clean_reads(sample, read):
        return os.path.join(PROC, f"decontam/{sample}_{read}.fastq.gz")

else:
    def clean_reads(sample, read):
        return os.path.join(PROC, f"trimmed/{sample}_{read}.fastq.gz")
