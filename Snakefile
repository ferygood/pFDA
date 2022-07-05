# This pipeline is to conduct pre-alignment QC and read mapping paired end reads to hg19
#
# Created by Eric Juo
# June 2, 2022

# Set global variables
lab = ["LAB" + str(i+1) for i in range(3)]
replicate = ["LIB" + str(i+1) for i in range(4)]
pairA = ["R1", "R3"]
pairB = ["R1", "R2"]

# Output files
rule all:
    input:
        expand("reports/fastqc/pre_alignment/PanelA/PanelA_{lab}_{replicate}_{pair}.fastq.html", lab=lab, replicate=replicate, pair=pairA),
        expand("reports/fastqc/pre_alignment/PanelB/PanelB_{lab}_{replicate}_{pair}.fastq.html", lab=lab, replicate=replicate, pair=pairB),
        "reports/multiqc/pre_alignment/PanelA",
        "reports/multiqc/pre_alignment/PanelB",

# Generate pre-alignment QC reports for each paired reads
rule pre_alignment_qc:
    input: 
        expand("../data/{{panel}}/FASTQ/{{panel}}_{lab}_{replicate}_{{pair}}.fastq.gz", lab=lab, replicate=replicate) 
    output:
        expand("reports/fastqc/pre_alignment/{{panel}}/{{panel}}_{lab}_{replicate}_{{pair}}.fastq.html", lab=lab, replicate=replicate)
    params: 
        outdir = "reports/fastqc/pre_alignment/{panel}/"
    threads: 16
    conda: 
        "envs/bioinfo.yml"
    singularity:
        "docker://staphb/fastqc"
    shell:
        "mkdir -p reports/fastqc/pre_alignment/{wildcards.panel} && fastqc -t {threads} -o {params.outdir} {input}"

# Generate multiqc report
rule pre_alignment_multiqc:
    input:
        "reports/fastqc/pre_alignment/{panel}"
    output:
        "reports/multiqc/pre_alignment/{panel}"
    params:
    conda:
        "envs/bioinfo.yml"
    singularity:
        "docker://staphb/multiqc"
    shell:
        "mkdir -p reports/multiqc/pre_alignment/ && multiqc -f -o {output} {input}"


    