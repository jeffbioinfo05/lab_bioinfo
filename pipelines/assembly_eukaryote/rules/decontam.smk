GENERAL_BIN = os.path.expanduser("~/lab/software/envs/general/bin")

def get_hosts():
    return config["decontam"].get("hosts", [])

if config["decontam"]["run"] and len(get_hosts()) > 0:

    rule decontam:
        input:
            r1  = os.path.join(PROC, "trimmed/{sample}_R1.fastq.gz"),
            r2  = os.path.join(PROC, "trimmed/{sample}_R2.fastq.gz"),
        output:
            r1     = os.path.join(PROC, "decontam/{sample}_R1.fastq.gz"),
            r2     = os.path.join(PROC, "decontam/{sample}_R2.fastq.gz"),
            report = os.path.join(RESULTS, "qc/decontam/{sample}_decontam_report.tsv"),
        params:
            hosts   = get_hosts(),
            tmp_dir = os.path.join(TMP, "decontam/{sample}"),
            bindir  = GENERAL_BIN,
        log:     os.path.join(LOGS, "decontam/{sample}.log")
        threads: THREADS_H
        shell:
            """
            export PATH={params.bindir}:$PATH
            mkdir -p {params.tmp_dir} $(dirname {output.r1}) $(dirname {output.report})

            current_r1={input.r1}
            current_r2={input.r2}

            prev_count=$(zcat $current_r1 | wc -l | awk "{print $1/4}")
            printf "sample\tstep\treads\tremoved\tremoved_pct\n" > {output.report}
            printf "{wildcards.sample}\traw\t$prev_count\t0\t0.00\n" >> {output.report}

            i=0
            for host_ref in {params.hosts}; do
                i=$((i+1))
                host_name=$(basename $host_ref | sed "s/\.fa.*//;s/\.fna.*//")
                out_r1={params.tmp_dir}/step${i}_R1.fastq.gz
                out_r2={params.tmp_dir}/step${i}_R2.fastq.gz
                bam={params.tmp_dir}/step${i}_host.bam

                minimap2 -ax sr -t {threads} $host_ref $current_r1 $current_r2 2>>{log} | \
                    samtools view -f 12 -F 256 -b -o $bam 2>>{log}
                samtools sort -n -@ {threads} $bam 2>>{log} | \
                    samtools fastq -1 $out_r1 -2 $out_r2 -0 /dev/null -s /dev/null 2>>{log}

                curr_count=$(zcat $out_r1 | wc -l | awk "{print $1/4}")
                removed=$((prev_count - curr_count))
                pct=$(awk "BEGIN {printf \"%.2f\", $removed/$prev_count*100}")
                printf "{wildcards.sample}\t$host_name\t$curr_count\t$removed\t$pct\n" >> {output.report}

                current_r1=$out_r1
                current_r2=$out_r2
                prev_count=$curr_count
            done

            cp $current_r1 {output.r1}
            cp $current_r2 {output.r2}
            """

    def clean_reads(sample, read):
        return os.path.join(PROC, f"decontam/{sample}_{read}.fastq.gz")

else:
    def clean_reads(sample, read):
        return os.path.join(PROC, f"trimmed/{sample}_{read}.fastq.gz")
