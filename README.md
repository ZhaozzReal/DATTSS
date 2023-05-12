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

```allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. The expected format is `,` to separate different files.

***
