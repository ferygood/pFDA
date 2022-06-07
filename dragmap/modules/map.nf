#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DRAGMAP {

    tag "DRAGMAP mapping on $sample_id"
    publishDir "$params.outdir", mode: 'copy'

    input:
      tuple val(sample_id), file(reads)
      path hash_reference

    output:
      tuple val(sample_id), file('*.sam')

    script:
      """
      dragen-os -r $hash_reference \
                -1 $reads \
                --output-directory $outdir \
                --output-file-prefix $sample_id
      """

    workflow {
        params.reads = '/home/azureuser/data/PanelA/FASTQ/*.fastq.gz'
        params.hash_reference = '/home/azureuser/ycchen/hg19_hash_fa'
        params.outdir = '/home/azureuser/ycchen/dragmap/results'
        read_pairs_ch = channel.fromFile( params.reads, checkIfExistts: true )
        hash_reference_ch = channel.fromPath( params.hash_reference )
        DRAGMAP( read_pairs_ch, hash_reference_ch.toList() )
    }
    
}