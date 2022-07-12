#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

//import modules
include { map_reads; sort_alignments; call_variants } from './modules/'



workflow {

}