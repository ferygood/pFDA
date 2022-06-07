#!/usr/bin/env nexflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

params.reads = '/home/azureuser/data/PanelA/FASTQ/*.fastq.gz'
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
include { DRAGMAP } from './modules/map.nf' 

workflow {
        read_ch = channel.fromFile( params.reads, checkIfExistts: true )
        hash_reference_ch = channel.fromPath( params.hash_reference )
        DRAGMAP( read_ch, hash_reference_ch.toList() )
}
