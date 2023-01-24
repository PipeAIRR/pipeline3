**A pipeline for BCR repertoire libraries from the form of - Illumina MiSeq 2x250 BCR mRNA. that were produced in the same fashion as those in [Greiff et al. 2014](https://bmcimmunol.biomedcentral.com/articles/10.1186/s12865-014-0040-5)**

<u>Sequence library:</u>

The sequences were amplified using CPrimers and VPrimers. They were sequences with Illumins MiSeq 2x250. 

Each 250 base-pair read was sequenced from one end of the target cDNA, so that the two reads together cover the entire variable region of the Ig heavy chain. The V(D)J reading frame proceeds from the start of read 2 to the start of read 1. Read 1 is in the opposite orientation (reverse complement), and contains the C-region primer sequence. Both reads begin with a random sequences of 4 nucleotides.


<u>Input files:</u>

1. The read can downloaded from the EBI European Nucleotide Archive under accession ID: ERP003950 or from SRA accession: SRR1383456
2. The primers sequences available online at the table below.

<u>Output files:</u>

1. sample_collapse-unique.fastq
2. sample_atleast-2.fastq
3. log and log tab file for each step.

<u>Sequence processing:</u>

* Paired-end assembly

	1. AssemblePairs align
* Quality control and primer annotation

	2. FilterSeq quality
	2. MaskPrimer score
* Deduplication and filtering

	3. CollapseSeq
	4. SplitSeq group





**CPrimers**

| Header     | Primer |
| ----------- | ----------- |
| IGHG   | CARKGGATRRRCHGATGGGG       |




**VPrimers**

| Header     | Primer |
| ----------- | ----------- |
| VH-FW1   | GAKGTRMAGCTTCAGGAGTC    |
| VH-FW2   | GAGGTBCAGCTBCAGCAGTC    |
| VH-FW3   | CAGGTGCAGCTGAAGSASTC    |
| VH-FW4   | GAGGTCCARCTGCAACARTC    |
| VH-FW5   | CAGGTYCAGCTBCAGCARTC    |
| VH-FW6   | CAGGTYCARCTGCAGCAGTC    |
| VH-FW7   | CAGGTCCACGTGAAGCAGTC    |
| VH-FW8   | GAGGTGAASSTGGTGGAATC      |
| VH-FW9   | GAVGTGAWGYTGGTGGAGTC      |
| VH-FW10  | GAGGTGCAGSKGGTGGAGTC      |
| VH-FW11   | GAKGTGCAMCTGGTGGAGTC      |
| VH-FW12   | GAGGTGAAGCTGATGGARTC      |
| VH-FW13   | GAGGTGCARCTTGTTGAGTC      |
| VH-FW14   | GARGTRAAGCTTCTCGAGTC      |
| VH-FW15   | GAAGTGAARSTTGAGGAGTC      |
| VH-FW16   | CAGGTTACTCTRAAAGWGTSTG    |
| VH-FW17   | CAGGTCCAACTVCAGCARCC      |
| VH-FW18   | GATGTGAACTTGGAAGTGTC      |
| VH-FW19	| GAGGTGAAGGTCATCGAGTC		|


