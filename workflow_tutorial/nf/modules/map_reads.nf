#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process mapping {
    
    conda 'envs/bwa.yml'
    publishDir 'results/mapped', mode: 'copy'
  
    input:
        tuple val(sample_id), file(fastq)
        //file index

    output:
        tuple val(sample_id), file('*.bam')

    script:
        """
        bwa mem data/genome.fa $fastq | samtools view -b - > ${sample_id}.bam
        """

}

workflow {
  fastq_data = channel.fromPath( 'data/samples/*.fastq' ).map { file -> tuple(file.baseName, file) }
  //index = channel.from( 'data/genome.fa' )
  mapping( fastq_data )
}