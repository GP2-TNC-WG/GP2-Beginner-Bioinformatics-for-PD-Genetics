## Module I. Genotyping Quality Control

#### Information

* Created by: GP2 Training and Networking

## Table of Contents

#### [0. Getting Started](#0)

#### [1. Sample QC - Genotyping call rates](#1)

#### [2. Sample QC - Heterozygosity](#2)

#### [3. Sample QC - Gender checking](#3)

#### [4. Sample QC - Ancestry](#4)

#### [5. Sample QC - Relatedness](#5)

#### [6. Variant QC - Missingness per variant](#6)

#### [7. Variant QC - Missingness by haplotype](#7)

#### [8. Variant QC - HWE](#8)

---
<a id="0"></a>
## 0. Getting Started

Genotyping quality control (QC) is crucial step to avoid spurious associations and bias. 
QC is usually performed at a sample and at a variant level. 

Files that you will need:

- **Raw PLINK Binary Files** containing cases and controls (.fam, .bim, .bed files)

This tutorial is written in R, and these are the packages you will need to install/load 

```
# Download the necessary packages 
if (!require(dplyr)) install.packages('dplyr')
if (!require(ggplot2)) install.packages('ggplot2')

# Load the necessary packages 
library(dplyr)
library(ggplot2)
```

Other programs you will need are:
- PLINK v1.9, GCTA

---
<a id="1"></a>

## 1. Sample QC - Genotyping call rates

```

NOTE: do all samples have a gender??
if not those samples that do not have a gender will be removed

NOTE: do all your samples have an affection status??
if not this will cause trouble at the end of this script
```
### Calculate genotyping call rates per sample

```
./plink --bfile RAW.test --missing --out call_rates
```
### Remove call rate outliers

```
./plink --bfile RAW.test --mind 0.05 --make-bed --out RAW.test_call_rate
```

```
mv raw.test_call_rate.irem CALL_RATE_OUTLIERS.txt
rm call_rates.lmiss
rm call_rates.log
rm call_rates.hh
rm call_rates.nosex
mv call_rates.imiss CALL_RATES_ALL_SAMPLES.txt
```

NOTES:
- All call rates outliers are in CALL_RATE_OUTLIERS.txt
- Open file call_rates.imiss and F_MISS - 1 = callrate
- All call rates saved in CALL_RATES_ALL_SAMPLES.txt
- Removed intermediate files

---
<a id="2"></a>

## 2. Sample QC - Heterozygosity

```
./plink --bfile RAW.test --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
./plink --bfile RAW.test --extract pruning.prune.in --make-bed --out pruned_data
./plink --bfile pruned_data --het --out prunedHet

awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt

cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt

mv prunedHet.het HETEROZYGOSITY_DATA.txt

./plink --bfile RAW.test_call_rate --remove all_outliers.txt --make-bed --out after_heterozyg
```

```
rm pruning.hh
rm pruning.nosex
rm pruned_data.bed
rm pruned_data.bim
rm pruned_data.fam
rm pruned_data.hh
rm pruned_data.log
rm pruned_data.nosex
rm prunedHet.hh
rm prunedHet.log
rm prunedHet.nosex
rm pruning.log
rm pruning.prune.in
rm pruning.prune.out
rm outliers1.txt
rm outliers2.txt
rm all_outliers.txt
```

NOTES:
- --het from LD pruned data > use F cut-off of -0.15 and <- 0.15 for inclusion
- outliers stored here -> all_outliers.txt
- all heterozygosity is stored here -> HETEROZYGOSITY_DATA.txt
- Removed intermediate files

---
<a id="3"></a>

## 3. Sample QC - Gender checking

```
./plink --bfile after_heterozyg_call_rate --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
./plink --bfile after_heterozyg_call_rate --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 
```

NOTES:
- Gender failures are stored in GENDER_FAILURES.txt
- Gender checks are stored in GENDER_CHECK1.txt and GENDER_CHECK2.txt

```
grep "PROBLEM" gender_check1.sexcheck > problems1.txt
grep "PROBLEM" gender_check2.sexcheck > problems2.txt
cat problems1.txt problems2.txt > GENDER_FAILURES.txt
cut -f 1,2 GENDER_FAILURES.txt > samples_to_remove.txt
```

```
./plink --bfile after_heterozyg_call_rate --remove samples_to_remove.txt --make-bed --out after_gender
```

```
rm gender_check2.hh
rm gender_check2.log
rm gender_check2.nosex
rm gender_check1.hh
rm gender_check1.log
rm gender_check1.nosex
rm after_heterozyg_call_rate.bed
rm after_heterozyg_call_rate.bim
rm after_heterozyg_call_rate.fam
rm after_heterozyg_call_rate.log
mv gender_check1.sexcheck GENDER_CHECK1.txt
mv gender_check2.sexcheck GENDER_CHECK2.txt
```

---
<a id="4"></a>

## 4. Sample QC - Ancestry

NOTES:
- This part is optional 
- No ancestry outliers -> based on Hapmap3 PCA plot, should be near combined CEU/TSI
- Keep in mind that this comparison with hapmap is based on the number of SNPs that overlap between your input dataset and hapmap

```
./plink --bfile after_gender --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
./plink --bfile after_gender --flip hapmap3_bin_snplis-merge.missnp --make-bed --out after_gender3
./plink --bfile after_gender3 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
./plink --bfile after_gender3 --exclude hapmap3_bin_snplis-merge.missnp --out after_gender4 --make-bed
./plink --bfile after_gender4 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
./plink --bfile hapmap3_bin_snplis --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
./plink --bfile hapmap3_bin_snplis --extract pruning.prune.in --make-bed --out pruned_data
./plink --bfile pruned_data --het --out prunedHet
./plink --bfile pruned_data --geno 0.01 --out pca --make-bed --pca 4
```

```
grep "EUROPE" pca.eigenvec > eur.txt
grep "ASIA" pca.eigenvec > asia.txt
grep "AFRICA" pca.eigenvec > afri.txt
grep -v -f eur.txt pca.eigenvec | grep -v -f asia.txt | grep -v -f afri.txt > new_samples.txt
cut -d " " -f 3 after_gender.fam > new_samples_add.txt
paste new_samples_add.txt new_samples.txt > new_samples2.txt
paste eur_add.txt eur.txt > euro.txt
paste asia_add.txt asia.txt > asiao.txt
paste afri_add.txt afri.txt > afrio.txt
cat new_samples2.txt euro.txt asiao.txt afrio.txt > pca.eigenvec2
```

R script for PCA plotting and filtering
```
R < PCA_in_R.R --no-save  
```

NOTES:
- Extract population of interest. In this case europeans, but can be any populations that is present in base comparison dataset
- This creates several plots and lists based on genetic ancestry

```
./plink --bfile after_gender --keep PCA_filtered_europeans.txt --make-bed --out after_gender_heterozyg_hapmap

cat PCA_filtered_asians.txt PCA_filtered_africans.txt PCA_filtered_mixed_race.txt > hapmap_outliers33.txt
```
---
<a id="5"></a>

## 5. Sample QC - Relatedness

NOTES: 
- This part is optional
- PIHAT threshold of 0.125 

```
./plink --bfile after_gender_heterozyg_hapmap --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
./plink --bfile after_gender_heterozyg_hapmap --extract pruning.prune.in --make-bed --out pruned_data
./plink --bfile pruned_data --het --out prunedHet

./gcta64 --bfile pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 
./gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
./plink --bfile after_gender_heterozyg_hapmap --keep GRM_matrix_0125.grm.id --make-bed --out after_gender_heterozyg_hapmap_pihat

cut -f 1,2 after_gender_heterozyg_hapmap.fam > IDs_before_relatedness_filter.txt
cut -f 1,2 after_gender_heterozyg_hapmap_pihat.fam > IDs_after_relatedness_filter.txt
```
---
<a id="6"></a>

## 6. Variant QC - Missingness per variant

```
./plink --bfile after_gender_heterozyg_hapmap_pihat --make-bed --out after_gender_heterozyg_pihat_mind --geno 0.05
grep "(--geno)" after_gender_heterozyg_pihat_mind.log > MISSINGNESS_SNPS.txt
```

NOTES:
- Missingness by case control P > 1E-4 # needs case control status

```
./plink --bfile after_gender_heterozyg_pihat_mind --test-missing --out missing_snps 
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
./plink --bfile after_gender_heterozyg_pihat_mind --exclude missing_snps_1E4.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing1
sort -u missing_snps_1E4.txt > VARIANT_TEST_MISSING_SNPS.txt
```
---
<a id="7"></a>

## 7. Variant QC - Missingness by haplotype

```
./plink --bfile after_gender_heterozyg_pihat_mind_missing1 --test-mishap --out missing_hap 
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\
/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt

sort -u missing_haps_1E4_final.txt > HAPLOTYPE_TEST_MISSING_SNPS.txt
./plink --bfile after_gender_heterozyg_pihat_mind_missing1 --exclude missing_haps_1E4_final.txt --make-bed --out after_gender_heterozyg_pihat_hapmap_mind_missing12
```

---
<a id="8"></a>

## 8. Variant QC - HWE

```
./plink --bfile after_gender_heterozyg_pihat_hapmap_mind_missing12 --filter-controls --hwe 1E-4 --write-snplist --out HWE_snps
./plink --bfile after_gender_heterozyg_pihat_hapmap_mind_missing12 --extract HWE_snps.snplist --make-bed --out after_gender_heterozyg_pihat_hapmap_mind_missing123
```

```
mv after_gender_heterozyg_pihat_hapmap_mind_missing123.bim FILTERED.test.bim
mv after_gender_heterozyg_pihat_hapmap_mind_missing123.bed FILTERED.test.bed
mv after_gender_heterozyg_pihat_hapmap_mind_missing123.fam FILTERED.test.fam
```
---
<a id="8"></a>

#######CLEAN FOLDER...

```
rm hapmap3_bin_snplis.bed
rm hapmap3_bin_snplis.bim
rm hapmap3_bin_snplis.fam
rm hapmap3_bin_snplis.log
rm hapmap3_bin_snplis.hh
rm after_gender4.bed
rm after_gender4.bim
rm after_gender4.fam
rm after_gender4.hh
rm after_gender4.log
rm hapmap3_bin_snplis-merge.missnp
rm after_gender3.bed
rm after_gender3.bim
rm after_gender3.fam
rm after_gender3.hh
rm after_gender3.log
rm after_gender.bed
rm after_gender.bim
rm after_gender.fam
rm after_gender.hh
rm after_gender.log
rm afri.txt
rm asia.txt
rm eur.txt
rm PCA.bed
rm PCA.bim
rm PCA.fam
rm PCA.log
rm PCA.eigenval
rm afrio.txt
rm asiao.txt
rm euro.txt
rm new_samples_add.txt
rm new_samples.txt
rm new_samples2.txt
rm after_gender_heterozyg_pihat_mind.hh
rm after_gender_heterozyg_hapmap_pihat.bed
rm after_gender_heterozyg_hapmap_pihat.bim
rm after_gender_heterozyg_hapmap_pihat.fam
rm GRM_matrix_0125.grm.gz
rm GRM_matrix_0125.grm.id
rm GRM_matrix.grm.gz
rm GRM_matrix.grm.id
rm after_gender_heterozyg_hapmap.bim
rm after_gender_heterozyg_hapmap.log
rm hapmap_outliers33.txt
rm after_gender_heterozyg_pihat_mind_missing1.bed
rm after_gender_heterozyg_pihat_mind_missing1.fam
rm after_gender_heterozyg_pihat_mind_missing1.hh
rm missing_snps_1E4.txt
rm after_gender_heterozyg_pihat_mind.bed
rm after_gender_heterozyg_pihat_mind.bim
rm after_gender_heterozyg_pihat_mind.fam
rm after_gender_heterozyg_pihat_mind.log
rm missing_snps.hh
rm missing_snps.log
rm missing_snps.missing
rm missing_haps_1E4_final.txt
rm missing_haps_1E4.txt
rm after_gender_heterozyg_pihat_mind_missing1.bim
rm after_gender_heterozyg_pihat_mind_missing1.log
rm missing_hap.hh
rm missing_hap.log
rm missing_hap.missing.hap
rm after_gender_heterozyg_hapmap.bed
rm after_gender_heterozyg_hapmap.fam
rm after_gender_heterozyg_hapmap.hh
rm plink.snplist
rm *.hh
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.bed
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.bim
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.fam
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.hh
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.log
rm plink.log
rm pruning.hh
rm pruning.nosex
rm pruned_data.bed
rm pruned_data.bim
rm pruned_data.fam
rm pruned_data.hh
rm pruned_data.log
rm pruned_data.nosex
rm prunedHet.hh
rm prunedHet.log
rm prunedHet.nosex
rm pruning.log
rm pruning.prune.in
rm pruning.prune.out
```

