version 1.0

workflow DragmapWorkflow {
    input {
        Array[File] inFileFqs
        File refPath
        String outPrefix
        Int numThreads = 16
    }

    call DRAGMAP {
        input: 
            inFileFastqR1 = inFileFqs[0],
            inFileFastqR2 = inFileFqs[1],
            refPath = refPath,
            outPrefix = outPrefix,
            numThreads = numThreads
    }

    call ToBam {
        input:
            inFileSam = DRAGMAP.outFileSam,
            outPrefix = outPrefix
    }

    output {
        File outBam = ToBam.outBam
    }
}

task DRAGMAP {
    input {
        File inFileFastqR1
        File inFileFastqR2
        File refPath
        String outPath = "."
        String outPrefix
        Int numThreads
        String dockerImage = "covid-lineage:1.0"
    }

    command <<<
        set -e -o pipefail

        dragen-os \
        --ref-dir ~{refPath} \
        --fastq-file1 ~{inFileFastqR1} \
        --fastq-file2 ~{inFileFastqR2} \
        --output-directory ~{outPath} \
        --output-file-prefix ~{outPrefix} \
        --num-threads ~{numThreads}
    >>>

    runtime {
        docker: dockerImage
        cpu: "16"
        memory: "16 GB"
    }

    output {
        File outFileSam = "~{outPrefix}.sam"
    }
}

task ToBam {
    input {
        File inFileSam
        String outPrefix
        String dockerImage = "covid-lineage:1.0"
    }

    command <<<
        set -e -o pipefail

        samtools view -@ 16 -1 ~{inFileSam} > ~{outPrefix}.bam
    >>>

    runtime {
        docker: dockerImage
        cpu: "16"
        memory: "16 GB"
    }

    output {
        File outBam = "~{outPrefix}.bam"
    }
}