# Pipeline summary
> you may add additional sessions. The following sessions are required. If nothing to enter in a session, put "NA" in that session.
## Team TaipeiBio, Pipeline 041
### Whether AI/ML is involved
Data processing: No
Mapping: No
Variant calling: No
Post processing: No
Other steps involve AI/ML: No
### Oncopanels applied
Oncopanel A: Yes
Oncopanel B: Yes
Oncopanel X: No
### Read processing
Raw reads were processed by FASTQC[1] version 0.11.9 with default parameters.
### Map to the reference genome
The processed sequence data were mapped to the human reference genome hg19 using bwa[2] version 0.7.17 with default parameters.
### Variant calling
we use Strelka2[3], Mutect2[4], Pisces[5], Scalpel[6], VarDict[7] with default parameters.
### Post processing
We then use SomaticSeq[8] to ensemble our VCF result with suggest single option.

[1] FASTQC. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[2] BWA. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
[3] Strelka2. https://github.com/Illumina/strelka
[4] Mutect2. https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
[5] Pisces. https://github.com/Illumina/Pisces
[6] Scalpel. http://scalpel.sourceforge.net/manual.html
[7] VarDict. https://github.com/AstraZeneca-NGS/VarDict
[8] Somaticseq. https://github.com/bioinform/somaticseq