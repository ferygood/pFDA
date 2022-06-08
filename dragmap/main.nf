#!/usr/bin/env nexflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

params.reads = '/home/azureuser/data/PanelA/FASTQ/*_R{1,3}.fastq.gz'
params.hash_reference = '/home/azureuser/ycchen/hg19_hash_fa'
params.outdir = '/home/azureuser/ycchen/dragmap/results'

log.info """\
DRAGMAP - Mapping
====================================
reads               : ${params.reads}
hash_reference      : ${params.hash_reference}
outdir              : ${params.outdir}
"""


//import modules
include { DRAGMAP } from './modules/map1.nf' 

workflow {
        read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
        hash_reference_ch = channel.fromPath( params.hash_reference )
        outdir_ch = channel.fromPath( params.outdir )
        DRAGMAP( read_pairs_ch, hash_reference_ch.toList(), outdir_ch )
}
