# SpliceTools
![Untitled-1](https://user-images.githubusercontent.com/23278950/176218942-46068321-f8ec-4137-9144-e30eeb76b248.jpg)
## Citation ##
Erik K Flemington, Samuel A Flemington, Tina M O’Grady, Melody Baddoo, Trang Nguyen, Yan Dong, Nathan A Ungerleider, SpliceTools, a suite of downstream RNA splicing analysis tools to investigate mechanisms and impact of alternative splicing, Nucleic Acids Research, 2023;, gkad111, https://doi.org/10.1093/nar/gkad111
## Requirements ##
* Python 3
* Perl 5
## Dependencies ##
Python libraries:
* Matplotlib
* Seaborn
* Pandas

Perl libraries:
* Parallel::ForkManager 
## Installation ##
```
git clone https://github.com/flemingtonlab/SpliceTools.git

cd SpliceTools/bin
```
## File formats ##
### Differential splicing files ###
| File  | Description |
| :-- | :-- |
| AS files  | Must be in rMATS[^1][^2][^3] output tsv format. SpliceTools utilizes columns 2-11, 20 and 23 of A3SS, A5SS, RI, and SE files and 2-13, 22, and 25 of rMATS MXE files (see [rMATS_file_formats.xlsx](https://github.com/flemingtonlab/SpliceTools/files/8995041/rMATS_file_formats.xlsx)) for those interested in converting other file formats for use in SpliceTools. |
| BED12 annotation  | Bed12 file format should be used with Gene ID column (column 4) containing the ENSEMBL ID and gene ID separated by an underscore (for example, “ENST00000456328_DDX11L1“).  |
| Genome fasta | Standard genome fasta file (can be wrapped or unwrapped). |
| Gene expression file | TSV file with *gene ID* in first column, *control TPM values* in following columns, *test TPM values* in the columns following control TPM values |

**Output files will be saved in the same directory as input AS files

## Usage ##
### RIFractionExpressed ###
Determines the fraction of expressed genes with statistically significant RI events at genes with a minimum input control or test TPM value.

Generates:
* Lists of statistically significant negative and positive IncDiff RI events in expressed genes (includes gene expression values for each)
	  
* File with list of genes without statistically significant change in RI.

##### Usage #####
```
perl RIFractionExpressed.pl \
	-r <retained intron file (rMATS JCEC)> \
	-e <expression file> \
	-TPM <min TPMs condition 1,min TPMs condition 2> \
	-SN <sample number condition 1,sample number condition 2> \
	-f <FDR>
```

##### Required #####
	-r <retained intron file>
	-e <expression file>
	-TPM <min TPMs condition 1,min TPMs condition 2>
		Note: if condition 1 or 2 to be not considered, enter "-"
		(e.g. TPM 3,- or TPM -,3)
	-SN <sample number condition 1,sample number condition 2>
	-f <FDR>
	
##### Additional #####
	-h help
	
##### Example #####
```
perl RIFractionExpressed.pl -r PATH/RIfile.txt -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05
```

##### Batch script #####
```perl RIFractionExpressedBatch.pl -r PATH/RI_files_directory -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05```



### RIIntronExonSizes ###
Generates:
* Lists of sizes for upstream exon, retained intron, and downstream exon for statistically significant positive and negative IncDiff retained intron events
* Summary file with average sizes

##### Usage #####
```
perl RIIntronExonSizes.pl \
	-r <retained intron file (rMATS JCEC)> \
	-a <bed12 annotation file> \
	-f <FDR>
```
##### Required #####
	-r <retained intron file>
	-a <bed12 annotation file>
	-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl RIIntronExonSizes.pl -r PATH/RIfile.txt -a PATH/bed12_annotation_file -f 0.05```

##### Batch script #####
```perl RIIntronExonSizesBatch.pl -r PATH/RI_input_files_directory -a PATH/bed12_annotation_file -f 0.05```


### RIMedley ###
Automatically runs the following on an rMATS RI JCEC file:
* RIFractionExpressed.pl
* RIIntronExonSizes.pl
* RISpliceSiteScoring.pl

##### Usage #####

```
perl RIMedley.pl \
	-r <retained intron file (rMATS JCEC)> \
	-a <bed12 annotation file> \
	-g <genome fasta file> \
	-e <expression file> \
	-TPM <min TPMs condition 1,min TPMs condition 2> \
	-SN <sample number condition 1,sample number condition2> \
	-f <FDR>
```

##### Required #####
	-r <retained intron file>
	-a <bed12 annotation file>
	-g <genome fasta file>
	-e <expression file>
	-TPM <min TPMs condition 1,min TPMs condition 2>
	-SN <sample number condition 1,sample number condition 2>
	-f <FDR>

##### Additional #####
	-h help
	
##### Example ##### 
```perl RIMedley.pl -r PATH/RIfile.txt -a PATH/bed12_annotation.bed -g PATH/genome.fa -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05```


### RISpliceSiteScoring ###
Generates:
* Lists of splice site scores (for plotting score distributions) for upstream donor and downstream acceptor sites for statistically significantly changed retained intron events
* Summary file with average scores

*Note: Inclusion of annotation file is optional but will generate data for annotated events for comparison.*

##### Usage #####
```
perl RISpliceSiteScoring.pl \
	-r <retained intron file (rMATS JCEC)> \
	-g <genome fasta file> \
	-f <FDR> \
	-a <bed12 annotation file>
```

##### Required #####
	-r <retained intron file>
	-g <genome fasta file>
	-f <FDR>

##### Additional #####
	-a <bed12 annotation file> (optional)
	-h help

##### Example #####
```perl RISpliceSiteScoring.pl -r PATH/RIfile.txt -g PATH/genome.fa -a PATH/bed12_annotation.bed -f 0.05```

##### Batch script #####
```perl RISpliceSiteScoring.pl -r PATH/RI_input_files_directory -g PATH/genome.fa -a PATH/bed12_annotation.bed -f 0.05```


### SEFractionExpressed ###
Determines the fraction of expressed genes with statistically significant SE events at genes with a minimum input control or test TPM value.

Also generates:
* Lists of statistically significant negative and positive IncDiff SE events in expressed genes (includes gene expression values for each)
* File with list of genes without statistically significant exon skipping.

##### Usage #####
```
perl SEFractionExpressed.pl \
	-s <skipped exon file (rMATS JCEC)> \
	-e <expression file> \
	-TPM <min TPMs condition 1,min TPMs condition 2> \
	-SN <sample number condition 1,sample number condition 2> \
	-f <FDR>
```
##### Required #####
	-s <skipped exon file>
	-e <expression file>
	-TPM <min TPMs condition 1,min TPMs condition 2>
		Note: if condition 1 or 2 to be not considered, enter "-"
		(e.g. TPM 3,- or TPM -,3)
	-SN <sample number condition 1,sample number condition 2>
	-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl SEFractionExpressed.pl -s PATH/SEfile.txt -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05```

##### Batch script #####
```perl SEFractionExpressedBatch.pl -s PATH/SE_input_files_directory -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05```


### SEIntronExonSizes ###
Generates:
* Lists of sizes for upstream exon, upstream intron, skipped exon, downstream intron, and downstream exon for statistically significant positive and negative IncDiff skipped exon events.
* Summary file with average sizes (for comparison, includes average sizes for potential SE events based on input annotation file).

##### Usage ##### 
```
perl SEIntronExonSizes.pl \
-s <skipped exon file (rMATS JCEC)> \
-a <bed12 annotation file> \
-f <FDR>
```
##### Required #####
	-s <skipped exon file>
	-a <bed12 annotation file>
	-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl SEIntronExonSizes.pl -s PATH/SEfile.txt -a PATH/bed12_annotation_file.bed -f 0.05```


##### Batch script #####
```perl SEIntronExonSizesBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -f 0.05```


### SEMedley ###
Automatically runs the following on an rMATS SE JCEC file:
* SEFractionExpressed.pl
* SEIntronExonSizes.pl
* SENumberSkipped.pl
* SESpliceSiteScoring.pl
* SETranslateNMD.pl
* SEUnannotated.pl

##### Usage ##### 
```
SEMedley.pl \
	-s <skipped exon file (rMATS JCEC)> \
	-a <bed12 annotation file> \
	-g <genome fasta file> \
	-e <expression file> \
	-TPM <min TPMs condition 1,min TPMs condition 2> \
	-SN <sample number condition 1,sample number condition2> \
	-f <FDR>
```
##### Required #####
	-s <skipped exon file>
	-a <bed12 annotation file>
	-g <genome fasta file>
	-e <expression file>
	-TPM <min TPMs condition 1,min TPMs condition 2>
	-SN <sample number condition 1,sample number condition 2>-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl SEMedley.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -g PATH/genome.fa -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05```


### SENumberSkipped ###
Generates:
* Lists of all predicted SE isoforms with number of exons skipped
* Lists of isoforms with maximum number of skipped exons
* Statistics for number of exons skipped

*Note: Uses annotated exon information to determine number of exons skipped
(0 exons skipped means that skipped exons were not in annotation file)*

##### Usage #####
```
SENumberSkipped.pl \
	-s <skipped exon file (rMATS JCEC)> \
	-a <bed12 annotation file> \
	-f <FDR>
```

##### Required #####
	-s <skipped exon file>
	-a <bed12 annotation file>
	-f <FDR>

##### Additional #####
	-h help

#### Example ####
```perl SENumberSkipped.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -f 0.05```

##### Batch script #####

```perl SENumberSkippedBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -f 0.05```


### SESpliceSiteScoring ###
Generates:
* Lists of splice site scores (for plotting score distributions) for upstream donor, skipped exon, and downstream acceptor sites for statistically significantly changed exon skipping events.
* Summary file with average scores.

*Note: Inclusion of annotation file is optional but will generate data for annotated events for comparison.*

##### Usage #####
```
perl SESpliceSiteScoring.pl \
	-s <skipped exon file (rMATS JCEC)> \
	-g <genome fasta file> \
	-f <FDR> \
	-a <bed12 annotation file>
```
##### Required #####
	-s <skipped exon file>
	-g <genome fasta file>
	-f <FDR>

##### Additional #####
	-a <bed12 annotation file> (optional)
	-h help

##### Example #####
```perl SESpliceSiteScoring.pl -s PATH/SEfile.txt -g PATH/genome.fa -a PATH/bed12_annotation.bed -f 0.05```

##### Batch script #####
```perl SESpliceSiteScoringBatch.pl -s PATH/SE_input_files_directory -g PATH/genome.fa -a PATH/bed12_annotation.bed -f 0.05 ```


### SETranslateNMD ###

Program will:
* Generate predicted RNA(DNA) isoform sequences based on annotated exon structures upstream and downstream from skipping event
* Generate translated isoform sequences for all such events
* Generate separate listings of frameshifted protein isoforms
* Generate candidate neopeptides produced by frameshifts
* Generate lists of skipped exon protein sequences for in-frame events for BATCH submission to NCBI conserved domain search (https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi)
* Identify events predicted to undergo Nonsense Mediated RNA Decay (NMD) and generate statistics on skipping events predicted to undergo and not undergo NMD

##### Usage #####
```
perl SETranslateNMD.pl \
	-s <skipped exon file (rMATS)> \
	-a <bed12 annotation file> \
	-g <genome fasta file> \
	-f <FDR>
```
##### Required #####
	-s <skipped exon file>
	-a <bed12 annotation file>
	-g <genome fasta file>
	-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl SETranslateNMD.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -g PATH/genome.fa -f 0.05```

##### Batch script #####
```perl SETranslateNMDBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -g PATH/genome.fa -f 0.05```


### SEUnannotated ###
Generates statistics on the number of skipping events that are present (annotated) and not present (unannotated) in a bed12 input annotation file. Also outputs lists of annotated and unannotated events (pos IncDiff and neg IncDiff).

##### Usage #####
```
perl SEUnannotated.pl \
	-s <skipped exon file (rMATS)> \
	-a <bed12 annotation file> \
	-f <FDR>
```

##### Required #####
	-s <skipped exon file>
	-a <bed12 annotation file>
	-f <FDR>

##### Additional #####
	-h help

##### Example #####
```perl SEUnannotated.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -f 0.05```

##### Batch script #####
```perl SEUnannotatedBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -f 0.05```


### SpliceCompare ###
Takes groups of rMATS differential splicing files and outputs lists of statistically significant events with positive IncDiff and negative IncDiff and summary files. It then compares differential splicing changes across experiments to assess functional relationships (hypergeometric test), outputting pval matrices, -log10(pval) matrices, and .svg graphical representation files displaying -log10(pval) of all comparisons.

Important notes:
 * Input file names must end in one or more of the following suffixes:
		A3SS.MATS.JCEC.txt
		A5SS.MATS.JCEC.txt
		MXE.MATS.JCEC.txt
		RI.MATS.JCEC.txt
		SE.MATS.JCEC.txt

* By default, SpliceCompare.pl generates a file of common splicing changes for each comparison. For large numbers of input files, this can lead to the generation of millions of comparison files. Use option -p 0 to suppress this output.

* By default, -log10(pval) maximum is 280. For experiments with less significant hypergeometric test values, max -log10(pval) can be decreased using -m option.*

##### Usage #####
```
perl SpliceCompare.pl \
	-i <input files directory (rMATS JCEC)> \
	-o <path to output directory> \
	-f <FDR>
```
##### Required #####
	-i <input files directory (rMATS JCEC)>
	-o <path to output directory>
	-f <FDR>

##### Optional #####
	-p <1=output overlap files (default), 0=no overlap files>
	-m <max cluster array value (default = 280)>

##### Additional #####
	-h help

##### Example #####
```perl SpliceCompare.pl -i PATH/input_files_directory -o PATH/output_files_directory -m 280 -p 1 -f 0.05```


## References ##

[^1]: Shen S., Park JW., Lu ZX., Lin L., Henry MD., Wu YN., Zhou Q., Xing Y. rMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data. PNAS, 111(51):E5593-601. doi: 10.1073/pnas.1419161111

[^2]: Park JW., Tokheim C., Shen S., Xing Y. Identifying differential alternative splicing events from RNA sequencing data using RNASeq-MATS. Methods in Molecular Biology: Deep Sequencing Data Analysis, 2013;1038:171-179 doi: 10.1007/978-1-62703-514-9_10

[^3]: Shen S., Park JW., Huang J., Dittmar KA., Lu ZX., Zhou Q., Carstens RP., Xing Y. MATS: A Bayesian Framework for Flexible Detection of Differential Alternative Splicing from RNA-Seq Data. Nucleic Acids Research, 2012;40(8):e61 doi: 10.1093/nar/gkr1291
