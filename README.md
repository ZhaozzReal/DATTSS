# DATTSS


## Description

DATTSS is developed for Dynamic analyses of Alternative Tandem TSS events from standard RNA-seq data. A key aspect of DATTSS is that it takes advantage of previously reported tandem TSSs. DATTSS incorporates these known TSSs as prior knowledge and then identifies the used tandem TSSs based on ‘change point’ model, which forms a linear regression model to infer the location of proximal tandem TSSs that can best explain the localized changes of read coverage profiles in the first exon regions of transcripts.


## Diagram illuminates the DATTSS algorithm

![image](https://github.com/ZhaozzReal/DATTSS/blob/main/diagram.png)


## Installation

DATTSS is built on Python, which requires packages ```HTSeq```, ```numpy```, ```multiprocessing``` and ```argparse```.

Clone the lastest development version of DATTSS and change directory:

```
  git clone https://github.com/ZhaozzReal/DATTSS.git
  cd DATTSS
```

## Usage

```
python DATTSS_main.py -b /path/to/allbamfiles.txt -anno /path/to/hgXX_first_exon_annotation_forDATTSS.txt -p 10 -o /path/to/DATTSS_output.txt
```

* ```/path/to/allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. 
The expected format is `,` to separate different files:
```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

* ```/path/to/hgXX_first_exon_annotation_forDATTSS.txt``` contains first exon regions of all transcripts from protein-coding genes and annotated internal tandem TSSs within them, which is built on hg38 or hg19 genome version.


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



## Tips
DATTSS detects and quantifies dynamic tandem TSS usage based on prior TSS annotations, which could ensure the discovery of authentic and reliable tandem TSS events in diverse physiological and pathological processes by exploiting the huge amount of standard RNA-seq data. And the fast calculation speed of DATTSS is friendly to population-level analyses of dynamic tandem TSS usage, which could be applied in RNA-seq data from TCGA and GTEx project.


<br/>
<br/>




## Compare alternative tandem TSS usage between conditions

***Step1: Detect and quantify alternative tandem TSS events***


```
python DATTSS_compare.py -b /path/to/allbamfiles.txt -anno /path/to/hgXX_first_exon_annotation_forDATTSS.txt -p 10 -r /path/to/hgXX_GEOCODE_ExonRegion_annotation.txt -d /path/to/exonCount/ -o /path/to/DATTSS_output.txt
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

***Step2: Infer statistically differential tandem TSS usage between conditions***

DATTSS utilizes DEXSeq, the model for differential exon usage analysis based on standard RNA-seq data, to detect differential usage of alternative 5' terminal exon. This statistical framework could account for biological variability between replicates and is robust to changes in isoform abundance between conditions.


```
Rscript Infer_DU_tandemTSS.R -b /path/to/allbamfiles.txt -I /path/to/DATTSS_output.txt -d /path/to/exonCount -o /path/to/DATTSS_tandem_DU.txt
```

Final results will be saved in the file ```DATTSS_tandem_DU.txt```.


## Citation
Please cite the following articles if you use DATTSS in your research:

Zhao Z, Chen Y, Zou X, Lin L, Zhou X, Cheng X, Yang G, Xu Q, Gong L, Li L, Ni T. Pan-cancer transcriptome analysis reveals widespread regulation through alternative tandem transcription initiation. Sci Adv. 2024 Jul 12;10(28):eadl5606. PMID: 38985880.

