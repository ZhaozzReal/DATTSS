# DATTSS


## Description

DATTSS is developped for systematic analysis of alternative tandem TSS events from standard RNA-seq data. Given the read coverage profiles of two or multiple RNA-seq samples, DATTSS first identifies the most distal TSS in each first exon region that has evidence of being used in the samples of interest. Then DATTSS infers the location of proximal tandem TSS based on “change point” model, where the usage of internal tandem TSS would lead to a profound increase in RNA-seq read coverage in the first exon region of a given transcript. Once the proximal tandem TSSs are identified, library size-normalized expression levels and the relative usage of tandem TSS are calculated.


## Diagram illuminates the DATTSS algorithm

![image](https://github.com/ZhaozzReal/DATTSS/blob/main/diagram.png)


## Installation

DATTSS is built on Python, which requiring packages HTSeq, numpy, multiprocessing and argparse.

Clone the lastest development version of DATTSS and change directory:

```
  git clone https://github.com/ZhaozzReal/DATTSS.git
  cd DATTSS
```

## Usage

```
python DATTSS_main.py -b /path/to/allbamfiles.txt -anno /path/to/hg38_first_exon_annotation.txt -p 10 -o /path/to/DATTSS_output.txt
```

* ```/path/to/allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. 
The expected format is `,` to separate different files:
```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

* ```/path/to/hg38_first_exon_annotation.txt``` contains first exon regions of all transcripts and annotated tandem TSSs within them, which is built on hg38 version.


## Output

In ```/path/to/DATTSS_output.txt```, each row corresponds to one tandem TSS event and its related features.

The explanation of each column is as follows:
 
      * genename:
 * first_exon_region：
 * strand：
 * Annotated_TSSs：
 * Proximal_TSS：
 * MSE_ratio：
 * sample1.bam：
 * sample2.bam：
 * sample3.bam：
 * sampleN.bam：



## Tips
DATTSS detects and quantifies dynamic tandem TSS usage based on prior TSS annotations, which could ensure the discovery of authentic and reliable tandem TSS events in diverse physiological and pathological processes by exploiting the huge amount of standard RNA-seq data. And the fast calculation speed of DATTSS is friendly to population-level analyses of dynamic tandem TSS usage, which could be used in RNA-seq data from TCGA and GTEx project.
