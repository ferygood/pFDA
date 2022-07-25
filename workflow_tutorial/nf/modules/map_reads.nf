#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process mapping {
    
    conda 'envs/bwa.yml'
    publishDir 'results/mapped', mode: 'copy'
  
    input:
        tuple val(sample_id), file(fastq)
        //file index
        path index

    output:
        tuple val(sample_id), file('*.bam')

    script:
        """
        bwa $index $fastq | samtools view -b - > ${sample_id}.bam
        """

}

workflow {
    fastq_data = channel.fromFilePairs( 'data/samples/*.fastq', size:1 ) //.map { file -> tuple(file.baseName, file) }
    index = channel.fromPath( 'data/genome.fa' )
    mapping( fastq_data, index.toList() )
}