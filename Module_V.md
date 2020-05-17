## Module V. Variant annotation and prioritization using whole exome and whole genome sequencing data
	
* Created by: GP2 Training and Networking

## Table of Contents

#### [0. Getting Started](#0)

#### [1. Background](#1)

#### [2. VCF file conversion](#2)

#### [3. Annotation](#3)

#### [4. Output explanation](#3)

---
<a id="0"></a>

## 0. Getting started

Files that you will need:

- **Variant calling format (VCF)** that you would like to annotate

- **ANNOVAR resources**

#### Set up your directory ..

```
cd Desktop/TEST/
```
#### Make a folder to store outputs ..

```
mkdir ANNOVAR_output
```
#### Download ANNOVAR

```
https://annovar.openbioinformatics.org/en/latest/user-guide/download/
```

#### Uncompress ANNOVAR
```
tar xvfz annovar.latest.tar.gz
```

#### Download resources needed to perform annotation
```
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ljb26_all humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20140902 humandb/

```
---
<a id="1"></a>

## 1. Background

**ANNOVAR** is an efficient tool that uses updated information to functionally annotate genetic variants detected from different genomes (eg, human genome in versions hg18, hg19, hg38 as well as other species such as mouse, yeast, fly, ecc.).
This package can make the following annotations:

**Annotation at the gene level**

It allows identifying if a SNP or structural variant (eg CNV) causes changes at the protein level and which amino acids are affected. Users can use independent datasets to call genes (RefSeq, UCSC, ENSEMBL, GENCODE).

**Annotation at the region level**

It allows identifying variants in specific regions of the genome (eg regions conserved among 44 species, regions that encode for transcription factors, loci found from GWAS, or any other annotation based on genomic intervals).

**Annotation based on a specific filter**

It allows identifying variants documented in different databases (eg, whether or not a certain variant is reported in dbSNP, what is the allele frequency in 1000 genomes or ExAC, as well as calculating the scores from SIFT / PolyPhen / LRT / MutationTaster / MutationAssessor / FATHMM / MetaSVM / MetaLR, identify intergenic variants using GERP ++ score <2, among others).

<a id="2"></a>

## 2. Convert *.VCF into *.avinput
```
perl convert2annovar.pl -format vcf4 TEST.vcf > TEST.avinput
```
VCF is the gold standard format that most researchers use. Normally we start from a * .VCF file and make it a more manageable * .avinput input.

---
<a id="3"></a>

## 3. Annotation
```
table_annovar.pl TEST.vcf humandb/ -buildver hg38 -out TEST.annovar -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20140902 -operation g,f,f,f -nastring . -vcfinput
```
 WGS data usually aligned with the human genome construct version GRCh38. We must map our variants in the same way.

The table_annovar.pl argument is a command used to annotate your input and generate a tab-delimited output that contains representative columns for each of the annotations.

The -operation argument tells ANNOVAR what operations to use for each of the protocols.
The protocols can be "g" that indicates "annotation at the level of the gene" or "f" that means "annotation based on a specific filter", among others.
The -nastring argument will add a "." for those variants for which there is no information.
The -csvout argument will generate an output in csv format easily readable by excel.


## 4. Output explanation
The output file contains multiple columns.
Each of the columns corresponds to each of the protocols specified in your script.
The columns Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene contain annotations on how mutations affect the structure of the gene and protein.
The ExAC * columns represent the allele frequency in different sub-populations within the Exome Aggregation Consortium, while avsnp147 refers to the rs ID of each variant in version 147 of dbSNP.
The remaining columns contain pathogenicity predictors for non-synonymous variants using the SIFT tools, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, GERP ++ scores, CADD scores, DANN scores, PhyloP scores and SiPhy scores.

<a id="5"></a>
