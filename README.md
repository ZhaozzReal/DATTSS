# DATTSS


## Description

DATTSS is developped for systematic analysis of alternative tandem TSS events from standard RNA-seq data. Given the read coverage profiles of two or multiple RNA-seq samples, DATTSS first identifies the most distal TSS in each first exon region that has evidence of being used in the samples of interest. Then DATTSS infers the location of proximal tandem TSS within the longest first exon region based on “change point” model, where the usage of internal tandem TSS would lead to a profound increase in RNA-seq read coverage in the first exon region of a given transcript. DATTSS progressively segments coverage profiles of the first exon regions at annotated TSSs, identifying the potential internal tandem TSS where the squared deviation decreases most from the mean coverage of the first exon when dividing the region into two segments compared with considering it as a single segment. Once the proximal tandem TSSs are identified, library size-normalized expression levels and the relative usage of tandem TSS are calculated.


## Diagram illuminates the DATTSS algorithm

![image](https://github.com/ZhaozzReal/DATTSS/blob/main/diagram.png)
