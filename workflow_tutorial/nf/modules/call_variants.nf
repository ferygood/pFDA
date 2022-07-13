#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process call_variants {
    
    conda 'envs/bcftools.yml'
    publishDir "results/calls", mode: 'copy'

    input:
        file all_sorted_bam
        file index

    output:
        file "all.vcf"

    script:
        """
        bcftools mpileup -f $index $all_sorted_bam | bcftools call -mv -> all.vcf
        """
}

workflow {
  all_sorted_bam = channel.fromPath( 'results/mapped/*.sorted.bam' ).collect()
  index = channel.fromPath( 'data/genome.fa' )
  call_variants( all_sorted_bam, index )
}