# Key Features for RSeQC
Comprehensively evaluate high throughput sequence data, especially RNA-seq data

1. [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py) checks the mapping statistics of reads that:
    * QC fail
	* Uniquely map
	* Splice map
	* Map in a proper pair

2. [geneBody_coverage.py](http://rseqc.sourceforge.net/#genebody-coverage-py) scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position
    * Generates a plot illustrating the coverage profile along the gene body  
![](https://cloud.githubusercontent.com/assets/11708268/15725008/bce81f50-2817-11e6-9b03-f205bde446f3.png)

3. [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)*to speculate the experimental design of the samples
    * Samples a subset of the reads from the BAM file and compares it to a reference gene model
		+ Determines if the RNA-seq has been sequenced with paired-end or single-end reads
		+ Also detects if sequencing is strand-specific and how the reads are stranded

4. [tin.py](http://rseqc.sourceforge.net/#tin-py) evalues RNA integrity at transcript level, this is analogous to RIN. 

## Not Yet Implemented

1. *inner_distance.py* estimates the inner distance distribution between paired reads
    * Estimated inner distance should be consistent with gel size selection
	* Detects structural variation or aberrant splicing in RNA-seq data

2. *read_distribution.py* calculates the fraction of reads mapped to the following regions:
    * Coding regions
	* 5' UTR exons
	* 3' URT exons
	* Introns
	* Intergenic regions *based on the gene model provided*
	    + Background noise levels can be estimated by adding custom gene models

3. *RPKM_saturation.py* determines the precision of estimated *Reads Per Kilobase of transcript per Million* (RPKM)
    * Estimation is performed at current sequencing depth by resampling the total mapped reads
	* Percent relative error **(100 * |RPKMp<sub>obs</sub> - RPKM<sub>real</sub> | / RPKM<sub>real</sub>)** is used to measure RPKM  
![](https://cloud.githubusercontent.com/assets/11708268/15725749/4ba206ae-281b-11e6-9065-8264178f0aad.png)

4. *junction_saturation.py* determines if the current sequencing depth is sufficient to perform alternative splicing analysis
    * Similar to **RPKM_saturation.py** where the splice junctions are detected for each resampled subset of reads
		+ Number of detected splice junctions will increase as the resample percentage increases 
		+ Eventually the resample percentage will reach a fixed value
	* Very important tool for alternative splicing analysis to ensure saturated sequencing depth  
![](https://cloud.githubusercontent.com/assets/11708268/15746372/49616932-28a4-11e6-9aee-725ee306ac28.png)

## Basic RSeQC Modules
* Sequence quality
* Nucleotide composition bias
* PCR bias
* GC bias

## RNA-seq Specific
* Sequencing saturation
* Mapped reads distribution
* Coverage uniformity
* Strand specificity
* Transcript level RNA itegrity

## FastQC vs RSeQC
* FastQC only focuses on raw sequence-related metrics
    + Not enough to ensure the usability of RNA-seq data
* RSeQC accesses the quality of RNA-seq experiments
