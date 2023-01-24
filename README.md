
**A pipeline for BCR repertoire libraries from the form of - UMI Barcoded Illumina MiSeq 325+275 paired-end 5’RACE BCR mRNA. that were produced in the same fashion as those in [VanderHeiden et al. 2017](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in)**

<u>Sequence library:</u>

The sequences were amplified using  1. AbSeq R1 Human IG Primers.fasta , 2. AbSeq R2 TS.fasta 3. AbSeq Human IG InternalCRegion.fasta.
They were sequences with Illumins MiSeq 325+275. 

Each read was sequenced from one end of the target cDNA so that the two reads together cover the entire variable region of the Ig heavy chain. Sequencing was performed using an uneven number of cycles for read 1 and read 2 using a 2x300 kit. The V(D)J reading frame proceeds from the start of read 2 to the start of read 1. Read 1 is in the opposite orientation (reverse complement), contains a partial C-region, and is 325 nucleotides in length. Read 2 contains the 5’RACE template switch site with a 17 nucleotide UMI barcode preceding it, and is 275 nucleotides in length.


<u>Input files:</u>

1. The read can downloaded from the NCBI Sequence Read Archive under BioProject accession ID: PRJNA248475 or downloaded first 25,000 sequences of sample HD09_N_AB8KB (accession: SRR4026043) using fastq-dump from the SRA Toolkit
2. The primers sequences available at the table below.

<u>Output files:</u>

1. sample_collapse-unique.fastq
2. sample_atleast-2.fastq
3. log and log tab file for each step.

<u>Sequence processing:</u>

* Quality control, UMI annotation and primer masking
	1. FilterSeq quality
	2. MaskPrimer score
* Generation of UMI consensus sequences
	3. PairSeq
	4. BuildConsensus
* Paired-end assembly of UMI consensus sequences
	5. PairSeq	
	6. AssemblePairs sequential 
* Deduplication and filtering
	7. MaskPrimer align
	8. ParseHeaders collapse
	9. CollapseSeq
	10. SplitSeq group




**AbSeq R1 Human IG Primers**

| Header     | Primer |
| ----------- | ----------- |
| Human-IGHM   |	GAATTCTCACAGGAGACGAGG      |
| Human-IGHD   |	TGTCTGCACCCTGATATGATGG     |
| Human-IGHA   |	GGGTGCTGYMGAGGCTCA  	   |
| Human-IGHE   |	TTGCAGCAGCGGGTCAAGG 	   |
| Human-IGHG   |	CCAGGGGGAAGACSGATG  	   |
| Human-IGK    |	GACAGATGGTGCAGCCACAG       |
| Human-IGL    |	AGGGYGGGAACAGAGTGAC        |



**AbSeq R2 TS**

| Header     | Primer |
| ----------- | ----------- |
| TS-shift0   |	TACGGG      |
| TS-shift1   |	ATACGGG     |
| TS-shift2   |	TCTACGGG    |
| TS-shift3   |	CGATACGGG   |
| TS-shift4   |	GATCTACGGG  |





**AbSeq Human IG InternalCRegion**

| Header     | Primer |
| ----------- | ----------- |
| Human-IGHA      |	GGCTGGTCGGGGATGC       |
| Human-IGHD      |	GAGCCTTGGTGGGTGC       |
| Human-IGHE      |	GGCTCTGTGTGGAGGC  	   |
| Human-IGHG      |	GGCCCTTGGTGGARGC 	   |
| Human-IGHM      |	GGGCGGATGCACTCCC  	   |
| Human-IGKC      |	TTCGTTTRATHTCCAS       |
| Human-IGLC-1    |	TGGGGTTGGCCTTGGG       |
| Human-IGLC-2    |	AGGGGGCAGCCTTGGG  	   |
| Human-IGLC-3    |	YRGCCTTGGGCTGACC       |
| Human-IGLC-4    |	GCTGCCAAACATGTGC       |




