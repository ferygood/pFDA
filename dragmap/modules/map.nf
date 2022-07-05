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

    workflow {
        params.reads = '/home/azureuser/data/PanelA/FASTQ/*_R{1,3}.fastq.gz'
        params.hash_reference = '/home/azureuser/ycchen/hg19_hash_fa'
        params.outdir = '/home/azureuser/ycchen/dragmap/results'
        read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
        hash_reference_ch = channel.fromPath( params.hash_reference )
        outdir_ch = channel.fromPath( params.outdir )
        DRAGMAP( read_pairs_ch, hash_reference_ch.toList(), outdir_ch )
    }
    
}
