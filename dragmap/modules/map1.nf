#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DRAGMAP {

    tag "DRAGMAP mapping on $sample_id"

    input:
      tuple val(sample_id), file(reads)
      path hash_reference
      path outdir

    output:
      tuple val(sample_id), file('*.sam')

    script:
      """
      dragen-os -r $hash_reference \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                --output-directory $outdir \
                --output-file-prefix $sample_id
      """
   
}
