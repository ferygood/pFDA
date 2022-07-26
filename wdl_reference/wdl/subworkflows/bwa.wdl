version 1.0

workflow BwaWorkflow {
    input {
        String sampleName
        Array[File] inFileFqs
        File refFa
        File refSa
        File refAnn
        File refBwt
        File refPac
        File refAmb

        String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
        String bwaCommandline = "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
    }

    call BwaMem {
        input:
            inFileFastqR1 = inFileFqs[0],
            inFileFastqR2 = inFileFqs[1],
            sampleName = sampleName,

            bwaCommandline = bwaCommandline,
            refFa = refFa,
            refSa = refSa,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refAmb = refAmb,
            dockerImage = dockerImage
    }

    output {
        File outFileBam = BwaMem.outFileBam
    }
}

task BwaMem {
    input {
        File inFileFastqR1
        File inFileFastqR2
        String sampleName
    
        String bwaCommandline
        File refFa
        File refSa
        File refAnn
        File refBwt
        File refPac
        File refAmb

        String dockerImage
        String bwaPath = "/usr/bin/"
    }

    command {
        set -e -o pipefail
    
        bash_ref_fasta=~{refFa}
    
        ~{bwaPath}~{bwaCommandline} \
          -R "@RG\tID:~{sampleName}\tLB:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA" \
          ~{inFileFastqR1} ~{inFileFastqR2} \
        | \
        samtools view -@ 16 -1 - > ~{sampleName}.bam
    }
    runtime {
        docker: dockerImage
        cpu: "16"
        memory: "7 GB"
    }
  
    output {
        File outFileBam = "~{sampleName}.bam"
    }
}