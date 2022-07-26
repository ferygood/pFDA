version 1.0

workflow Mutect2Workflow {
    input {
      File refFa
      File refFai
      File refDict
      File inFileBam
      String sampleName
      Int memory = 16
      String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    call SortSam {
        input:
            inFileUnsortedBam = inFileBam,
            sampleName = sampleName,
            memory = memory,
            dockerImage = dockerImage
    }
    
    call Mutect2 {
        input:
            refFasta = refFa,
            refFai = refFai,
            refDict = refDict,
            inFileBam = SortSam.outFileSortedBam,
            inFileBai = SortSam.outFileSortedBai,
            sampleName = sampleName,
            memory = memory,
            dockerImage = dockerImage
    }

    output {
        File outFileSortBam = SortSam.outFileSortedBam
        File outFileSortBai = SortSam.outFileSortedBai
        File outFileVcf = Mutect2.outFileVcf
        File outFileVcfIdx = Mutect2.outFileVcfIdx
    }
}

task SortSam {
    input {
        File inFileUnsortedBam
        String sampleName
        Int memory
        String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    command <<<
        set -e -o pipefail

        /gatk/gatk-4.2.0.0/gatk --java-options "-Xmx~{memory}G -Xms~{memory}G" \
        SortSam \
        --INPUT ~{inFileUnsortedBam} \
        --OUTPUT ~{sampleName}.sorted.bam \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE false
    >>>

    runtime {
        docker: dockerImage
        cpu: "1"
        memory: "16 GB"
    }

    output {
        File outFileSortedBam = "~{sampleName}.sorted.bam"
        File outFileSortedBai = "~{sampleName}.sorted.bai"
    }
}


task Mutect2 {
    input {
      File refFasta
      File refFai
      File refDict
      File inFileBam
      File inFileBai
      String sampleName
      Float initialTLod = 25
      Float tLodToEmit = 25
      Int callableDepht = 10
      Int memory
      String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    command <<<
        set -e -o pipefail

        gatk --java-options "-Xmx~{memory}G -Xms~{memory}G" \
        Mutect2 \
            --reference ~{refFasta} \
            --input ~{inFileBam} \
            --output ~{sampleName}.vcf.gz \
            --initial-tumor-lod ~{initialTLod} \
            --tumor-lod-to-emit ~{tLodToEmit} \
            --callable-depth ~{callableDepht}

    >>>

    runtime {
        docker: dockerImage
        cpu: "1"
        memory: "16 GB"
    }

    output {
        File outFileVcf = "~{sampleName}.vcf.gz"
        File outFileVcfIdx = "~{sampleName}.vcf.gz.tbi"
    }
}
