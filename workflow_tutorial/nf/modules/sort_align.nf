#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process sort_align {
    
    conda 'envs/samtools.yml'
    publishDir "results/mapped", mode: 'copy'

    input:
        tuple val(sample_id), file(bam)

    output:
        tuple val(sample_id), file('*.sorted.bam')

    script:
        """
        samtools sort -o ${sample_id}.sorted.bam $bam
        """
}

workflow {
  bam_data = channel.fromPath( 'results/mapped/*.bam' ).map { file -> tuple(file.baseName, file) }
  sort_align( bam_data )
}