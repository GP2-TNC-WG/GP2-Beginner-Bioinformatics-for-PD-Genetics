## Module V. Variant annotation and prioritization using whole exome and whole genome sequencing
	
* Created by: GP2 Training and Networking

## Table of Contents

#### [0. Getting Started](#0)

#### [1. VCF file conversion](#1)

#### [2. Annotation](#2)

#### [3. Output explanation](#3)

---
<a id="0"></a>

## 0. Getting started

ANNOVAR is an efficient tool that uses updated information to functionally annotate genetic variants detected from different genomes (eg, human genome in versions hg18, hg19, hg38 as well as other species such as mouse, yeast, fly, ecc.).
This package can make the following annotations:

**Annotation at the gene level**

It allows identifying if a SNP or structural variant (eg CNV) causes changes at the protein level and which amino acids are affected. Users can use independent data sets for nominal genes (RefSeq, UCSC, ENSEMBL, GENCODE).

**Annotation at the region level**

It allows identifying variants in specific regions of the genome (eg regions conserved among 44 species, regions that encode for transcription factors, loci found from GWAS, or any other annotation based on genomic intervals).

**Annotation based on a specific filter**

It allows identifying variants documented in different databases (eg, whether or not a certain variant is reported in dbSNP, what is the allele frequency in 1000 genomes or ExAC, as well as calculating the scores from SIFT / PolyPhen / LRT / MutationTaster / MutationAssessor / FATHMM / MetaSVM / MetaLR, identify intergenic variants using GERP ++ score <2, among others).

### Set up your directory ..

```
cd /data/LNG/saraB/ANNO/
```
### Make a folder to store outputs ..

```
mkdir ANNOVAR_output
```
### Download ANNOVAR

```
http://download.openbioinformatics.org/cgi-bin/annovar_download.cgi
```

### Uncompress ANNOVAR
```
tar xvfz annovar.latest.tar.gz
```

### Download resources needed to perform annotation
```
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
```
---
<a id="1"></a>

## 1. Convert *.VCF into *.avinput
```
convert2annovar.pl -format vcf4 /data/LNG/saraB/ANNO/ANNOVAR_input/FILTERED.ataxia_chr1.WES.vcf > /data/LNG/saraB/ANNO/ANNOVAR_input/FILTERED.ataxia_chr1.WES.avinput
```
VCF is the gold standard format that most researchers use. Normally we start from a * .VCF file and make it a more manageable * .avinput input.

---
<a id="2"></a>

## 2. Annotation
```
table_annovar.pl /data/LNG/saraB/ANNO/ANNOVAR_input/FILTERED.ataxia_chr1.WES.avinput ANNOVAR_DATA/hg38 -buildver hg38 \
```
AMP-PD WGS data has been aligned with the human genome construct version GRCh38. We must map our variants in the same way.

The table_annovar.pl argument is a command used to annotate your input and generate a tab-delimited output that contains representative columns for each of the annotations.

### Specify directory and output name
```
-out /data/LNG/saraB/ANNO/ANNOVAR_output/FILTERED.ataxia_chr1.WES \
```
### Specify type of annotation and databases to be used
```
-remove -protocol refGene,ensGene,ljb26_all,gnomad211_genome,exac03,avsnp147,dbnsfp30a -operation g,g,f,f,f,f,f -nastring . -csvout
```

The -operation argument tells ANNOVAR what operations to use for each of the protocols.
The protocols can be "g" that indicates "annotation at the level of the gene" or "f" that means "annotation based on a specific filter", among others.
The -nastring argument will add a "." for those variants for which there is no information.
The -csvout argument will generate an output in csv format easily readable by excel.

### Annotate splicing variants
```
-arg '-splicing 15',,, \
```
---
<a id="3"></a>

## 3. Output explanation
The output file contains multiple columns.
Each of the columns corresponds to each of the protocols specified in your script.
The columns Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene contain annotations on how mutations affect the structure of the gene and protein.
The ExAC * columns represent the allele frequency in different sub-populations within the Exome Aggregation Consortium, while avsnp147 refers to the rs ID of each variant in version 147 of dbSNP.
The remaining columns contain pathogenicity predictors for non-synonymous variants using the SIFT tools, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, GERP ++ scores, CADD scores, DANN scores, PhyloP scores and SiPhy scores.

<a id="4"></a>
