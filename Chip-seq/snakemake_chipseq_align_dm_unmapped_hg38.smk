import glob
import os

# === Config ===
HUMAN_REF = "/seqs_for_alignment_pipelines/GRCh38_no_alt_analysis_set_genome"
DM_REF    = "/dm_ref_genome/dm6"
TRIMMED   = "Trimmed"

# Auto-detect samples from Trimmed folder
SAMPLES = [
    os.path.basename(f).replace("_R1.trimmed.fastq.gz", "")
    for f in glob.glob(os.path.join(TRIMMED, "*_R1.trimmed.fastq.gz"))
]

# === Target rule ===
rule all:
    input:
        expand("QC/{sample}.flagstat.txt", sample=SAMPLES),
        expand("Peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),


# === Step 1: Align to dm6, output unmapped reads ===
rule align_dm6:
    input:
        r1 = TRIMMED + "/{sample}_R1.trimmed.fastq.gz",
        r2 = TRIMMED + "/{sample}_R2.trimmed.fastq.gz",
    output:
        bam      = "dm_aligned/{sample}.dm6.bam",
        bai      = "dm_aligned/{sample}.dm6.bam.bai",
        unmap_r1 = "dm_aligned/{sample}.unmapped_to_dm6.1.fq.gz",
        unmap_r2 = "dm_aligned/{sample}.unmapped_to_dm6.2.fq.gz",
    threads: 8
    shell:
        """
        module load bowtie2/2.4.5
        module load samtools/1.20
        bowtie2 -x {DM_REF} \
            -1 {input.r1} -2 {input.r2} \
            --very-sensitive -X 2000 \
            --no-mixed --no-discordant \
            -p {threads} \
            --un-conc-gz dm_aligned/{wildcards.sample}.unmapped_to_dm6.%.fq.gz \
            | samtools view -bS - \
            | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


# === Step 2: Align unmapped reads to hg38 ===
rule align_hg38:
    input:
        r1 = "dm_aligned/{sample}.unmapped_to_dm6.1.fq.gz",
        r2 = "dm_aligned/{sample}.unmapped_to_dm6.2.fq.gz",
    output:
        bam = "hg38_aligned/{sample}.hg38.bam",
        bai = "hg38_aligned/{sample}.hg38.bam.bai",
    threads: 8
    shell:
        """
        module load bowtie2/2.4.5
        module load samtools/1.20
        bowtie2 -x {HUMAN_REF} \
            -1 {input.r1} -2 {input.r2} \
            --very-sensitive -X 2000 \
            --no-mixed --no-discordant \
            -p {threads} \
            | samtools view -bS - \
            | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


# === Step 3: Sort and deduplicate with sambamba ===
rule dedup:
    input:
        bam = "hg38_aligned/{sample}.hg38.bam",
    output:
        sorted_bam = "hg38_aligned/{sample}.hg38.sorted.bam",
        dedup_bam  = "hg38_aligned/{sample}.hg38.sorted.dedup.bam",
        bai        = "hg38_aligned/{sample}.hg38.sorted.dedup.bam.bai",
    threads: 8
    shell:
        """
        module load sambamba/0.7.0
        module load samtools/1.20
        sambamba sort -t {threads} -m 20G \
            -o {output.sorted_bam} {input.bam}
        sambamba markdup -t {threads} --remove-duplicates \
            {output.sorted_bam} {output.dedup_bam}
        samtools index {output.dedup_bam}
        """


# === Step 4: QC - flagstat ===
rule flagstat:
    input:
        bam = "hg38_aligned/{sample}.hg38.sorted.dedup.bam",
    output:
        "QC/{sample}.flagstat.txt",
    shell:
        """
        module load samtools/1.20
        samtools flagstat {input.bam} > {output}
        """


# === Step 5: Peak calling with MACS2 ===
rule macs2:
    input:
        bam = "hg38_aligned/{sample}.hg38.sorted.dedup.bam",
    output:
        "Peaks/{sample}_peaks.narrowPeak",
    shell:
        """
        module load MACS/2.2.7.1
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE -g hs \
            -n {wildcards.sample} \
            -q 0.01 -B \
            --outdir Peaks
        """
