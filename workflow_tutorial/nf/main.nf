#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */
params.fastq = "data/samples/*.fastq"
params.index = "data/genome.fa"


//import modules
include { mapping } from './modules/map_reads.nf'
include { sort_align } from './modules/sort_align.nf'
include { call_variants } from './modules/call_variants.nf'

//workflow 
workflow {
  fastq_data = channel.fromPath( params.fastq ).map { file -> tuple(file.baseName, file) }
  index = channel.fromPath( params.index )
  mapping( fastq_data, index.toList() )
  sort_align( mapping.out )
  call_variants( sort_align.out, index )
}