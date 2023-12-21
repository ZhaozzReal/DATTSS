# DATTSS


## Description

DATTSS is developed for Dynamic analyses of Alternative Tandem TSS events from standard RNA-seq data. Given the read coverage profiles of two or multiple RNA-seq samples, DATTSS first identifies the most distal TSS in each first exon region that has evidence of being used in the samples of interest. Then DATTSS infers the location of proximal tandem TSS based on “change point” model, where the usage of proximal tandem TSS would lead to a profound increase in RNA-seq read coverage in the first exon region. Once the proximal tandem TSSs are identified, library size-normalized expression levels and the relative usage of distal tandem TSS are calculated.


## Diagram illuminates the DATTSS algorithm

![image](https://github.com/ZhaozzReal/DATTSS/blob/main/diagram.png)


## Installation

DATTSS is built on Python, which requiring packages ```HTSeq```, ```numpy```, ```multiprocessing``` and ```argparse```.

Clone the lastest development version of DATTSS and change directory:

```
  git clone https://github.com/ZhaozzReal/DATTSS.git
  cd DATTSS
```

## Usage

```
python DATTSS_main.py -b /path/to/allbamfiles.txt -anno /path/to/hg38_first_exon_annotation_forDATTSS.txt -p 10 -o /path/to/DATTSS_output.txt
```

* ```/path/to/allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. 
The expected format is `,` to separate different files:
```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

* ```/path/to/hg38_first_exon_annotation_forDATTSS.txt``` contains first exon regions of all transcripts from protein-coding genes and annotated tandem TSSs within them, which is built on hg38 version.


## Output

In ```/path/to/DATTSS_output.txt```, each row corresponds to one tandem TSS event and its related features.

The explanation of each column is as follows:
 
 * genename: the HUGO gene symbol
 * first_exon_region：the genomic region of first exon region ranging from distal TSS to first exon end
 * Proximal_TSS：the position of proximal TSS inferred by DATTSS based on joint RNA-seq read coverage profile
 * MSE_ratio：the ratio of mean square error (MSE) calculated at the inferred proximal tandem TSS
 * sample1.bam：the distal TSS usage of given tandem TSS event in sample1 (```None``` means that the mean coverage of common 5'UTR region is lower than certain cutoff (default=20), which should be discarded from downstream analysis)
 * sample2.bam：the distal TSS usage of given tandem TSS event in sample2
 * sample3.bam：the distal TSS usage of given tandem TSS event in sample3
 * sampleN.bam：the distal TSS usage of given tandem TSS event in sampleN



# Compare alternative tandem TSS usage between conditions with replicates

**Command** 

```
python DATTSS_compare.py -b /path/to/allbamfiles.txt -anno /path/to/hg38_first_exon_annotation_forDATTSS.txt -p 10 -r /path/to/hg38_GEOCODE_ExonRegion_annotation.txt -d /path/to/exonCount/ -o /path/to/DATTSS_output.txt
```

allbamfiles.txt contains all filename of bamfile between two conditions, as shown below:

```
condition1=/path/to/ctrl1.bam,/path/to/ctrl2.bam 
condition2=/path/to/case1.bam,/path/to/case2.bam
```

Following counting reads mapped to GEOCODE-annotated exons and inferred alternative terminal exons, DATTSS outputs the results into directory ```/path/to/exonCount/```. 
```
/path/to/exonCount/
  |-- ctrl1_exoncount.txt
  |-- ctrl2_exoncount.txt
  |-- case1_exoncount.txt
  |-- case2_exoncount.txt
```

## Infer statistically differential tandem TSS usage between conditions

DATTSS utilizes DEXSeq, the model for differential exon usage analysis based on standard RNA-seq data, to detect differential usage of alternative 5' terminal exon. This statistical framework could account for biological variability between replicates and is robust to changes in isoform abundance between conditions.

**Command**

```
Rscript Infer_DU_tandemTSS.R -b /path/to/allbamfiles.txt -I /path/to/DATTSS_output.txt -d /path/to/exonCount -o /path/to/DATTSS_tandem_DU.txt
```

Final results will be saved in the file ```DATTSS_tandem_DU.txt```.




## Tips
DATTSS detects and quantifies dynamic tandem TSS usage based on prior TSS annotations, which could ensure the discovery of authentic and reliable tandem TSS events in diverse physiological and pathological processes by exploiting the huge amount of standard RNA-seq data. And the fast calculation speed of DATTSS is friendly to population-level analyses of dynamic tandem TSS usage, which could be applied in RNA-seq data from TCGA and GTEx project.
