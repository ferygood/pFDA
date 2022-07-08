# Pipeline summary
> you may add additional sessions. The following sessions are required. If nothing to enter in a session, put "NA" in that session.
## Team 0123456, Pipeline P1
### Whether AI/ML is involved
Data processing: No
Mapping: No
Variant calling: Yes
Post processing: No
Other steps involve AI/ML: data manipulation
### Oncopanels applied
Oncopanel A: Yes
Oncopanel B: Yes
Oncopanel X: No
### Read processing
Raw reads were processed by ReadProcess[1] version 1.0 with parameters "-q 30 -t 5".
### Map to the reference genome
The processed sequence data were mapped to the human reference genome GRChNN/hgNN using MapTool[2] version 0.9.2 with default parameters.
### Remove deduplicated reads
Deduplicated reads were detected and removed by SomeTools[3] version 3.6 using function RemoveDuplicates with parameters "-a 0 -b 5 -c 100".
> If UMI tag is used, please describe.
### Local alignment optimization
Local alignment optimization was performed using AnotherTool[4] version 2.10 with parameters "-x 100 -p"
### Variant calling
CallVar[5] version 11.0 was used for variant calling with parameters "-m 0.01 -M 5 -o vcf". A blocklist/allowlist was used in the variant calling. Variant calling was also restrict to our defined regions. These restricted information is provided in a text document (in VCF or BED format) along with the vcf result submission.
### Post processing
The reported VCF files were further processed with an in-house python script to remove identity information and reformat according to the VCF template provided by the PrecisionFDA.
**Variant normalization**:
- Break complex indels into smaller components: Yes
- Decompose multiallelic variants: Yes # "Yes" means multiple variants at the same location are put into separate lines, each line contains one biallelic variant.
- Left aligned: Yes # This is required.
- Parsimony/Trim: Yes # e.g., GAC>GC --trim--> GA>G. Trim is required.

[1] Author A., Author B. A great tool for processing raw reads. *Great Journal* 66, 77-88 (2020)
[2] Author C. Read mapping tools. https://github.com/My-Lab/MapTool/
[3] SomeTools. https://myinstitute.github.io/sometools/
[4] Author D., et al. Local alignment optimization will increate the accuracy of variant calling. *Another Great Journal* 199, 123-45 (2022)
[5] CallVar. https://www.ngs-tools.com/callvar/