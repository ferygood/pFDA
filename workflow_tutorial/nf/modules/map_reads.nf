#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process mapping {
    
    conda 'envs/bwa.yml'
  
    input:
        tuple val(sample_id), file(fastq)
        file index

    output:
        file 'results/${sample_id}.bam'

    script:
        """
        bwa mem $index $fastq | samtools view -b - > results/${sample_id}.bam
        """

}

workflow {
  fastq_data = channel.fromPath( 'data/samples/*.fastq' ).map { file -> tuple(file.baseName, file) }
  index = channel.fromPath( 'data/genome.fa' )
  mapping( fastq_data, index )
}