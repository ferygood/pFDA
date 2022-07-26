version 1.0

import "subworkflows/bwa.wdl" as mapper
import "subworkflows/mutect2.wdl" as caller

workflow MainWorkflow {
    input {
        String sampleName
        Array[File] inFileFqs
        File refFa
        File refFai
        File refDict
        File refSa
        File refAnn
        File refBwt
        File refPac
        File refAmb
        String bwaCommandline = "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
        Int memory = 16
    }

    call mapper.BwaWorkflow {
        input:
            sampleName = sampleName,
            inFileFqs = inFileFqs,
            refFa = refFa,
            refSa = refSa,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refAmb = refAmb,
            bwaCommandline = bwaCommandline,
            dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest",
    }

    call caller.Mutect2Workflow {
        input:
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            inFileBam = BwaWorkflow.outFileBam,
            sampleName = sampleName,
            memory = memory,
            dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    output {
        File outFileBam = BwaWorkflow.outFileBam
        File outFileSortBam = Mutect2Workflow.outFileSortBam
        File outFileSortBai = Mutect2Workflow.outFileSortBai
        File outFileVcf = Mutect2Workflow.outFileVcf
        File outFileVcfIdx = Mutect2Workflow.outFileVcfIdx
    }
}
