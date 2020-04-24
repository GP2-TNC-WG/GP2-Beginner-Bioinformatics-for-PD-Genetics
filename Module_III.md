## Module III.  GWAS for Binary and Quantitative Traits 

### Information
- 


## Table of Contents

#### [0. Getting Started](#0)

#### [1. Generation of a Covariates File](#1)

#### [2. Running a GWAS in PLINK](#2)

#### [3. Generating Summary Statistics](#3)

#### [4. Prep for FUMA GWAS](#4)

#### [5. Data Visualization: Scree Plot](#5)

#### [6. Data Visualization: QQ Plot](#6)

#### [7. Data Visualization: Manhattan Plot](#7)

#### [8. Logistic versus Linear Regressions](#8)

#### [9. Meta-Analysis](#9) --> Might not be in this module?

---
<a id="0"></a>
## 0. Getting Started
Files that you will need:
- **QC’d Pre-Imputation PLINK Binary Files** containing cases and controls (.fam, .bim, .bed files)
	- These have gone through the quality control (QC) steps outlined in Module I

-   **QC’d Imputed PLINK Binary Files** containing cases and controls (.fam, .bim, .bed files)
	-   These have gone through Michigan Server imputation and the soft/hard-call QC steps outlined in Module II
   
-   **Covariates File** containing the covariates you would like to correct by (more on this in section 1)
	-  	This file at a minimum includes sample information (such as ID, SEX, PHENO) and principle components

This tutorial is written in R, and these are the packages you will need to install/load 
```R
# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(qqman)) install.packages('qqman')

# Load the necessary packages 
library(tidyverse)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(qqman)
```
Other programs you will need are:
- PLINK v1.9
- RVTests [optional]
- flashPCA [optional]

---
<a id="1"></a>
## 1. Generation of a Covariates File

### Background Information
-   What is a covariate file?
	-   A covariate file is a file that includes the sample information, principal components, and other information you would like to correct for
-   Why do we need it for a GWAS?
	-   Adding covariates is beneficial because non-confounding covariates can explain some of the variance → reducing noise
- What is included in a covariates file?
	- This file at a minimum includes sample information (such as ID, SEX, PHENO) and principle components
	- If you have additional covariates you'd like to use (such as AGE, FAMILY HISTORY, EDUCATION etc.) you can merge these into the covariate file 
- How do you generate the PCs?
	- This is done in PLINK

### Generating Principal Components in PLINK
```bash
# Make sure to use high-quality SNPs 
plink --bfile EXAMPLE_UNIMPUTED_QC --maf 0.01 --geno 0.05 --hwe 1E-6 --make-bed --out EXAMPLE_UNIMPUTED

# Prune out unnecessary SNPs (only need to do this to generate PCs)
plink --bfile EXAMPLE_UNIMPUTED --indep-pairwise 50 5 0.5 --out prune 

# Keep only pruned SNPs (only need to do this to generate PCs)
plink --bfile EXAMPLE_UNIMPUTED --extract prune.prune.in --make-bed --out prune 

# Generate PCs 
plink --bfile prune --pca --out EXAMPLE_UNIMPUTED.PCA
```


### Generating a Covariate File in R 
```R
# Generate a covariate file 

# Read in the PCA Eigenvalues and Eigenvectors 
print("Read in pca.eigenvec file from PLINK")
eigenvec <- read.delim("EXAMPLE_UNIMPUTED.eigenvec", sep ="\t", header = T, stringsAsFactors = F)

# Read in the .fam file
fam <- fread("EXAMPLE_UNIMPUTED.fam", header = F)

# Read in other covariates 
  # RANDOM AGES AND EDUCATION LEVELS !!
  # Format
    # FID IID AGE EDUCATION ETC
additional_cov <- fread("SAMPLE_COVARIATES.txt")

# Rename the columns 
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# Combine the covariate file with the fam file 
pheno_pcs <- left_join(fam, additional_cov, by=c("FID", "IID"))

# Now combine the additional covariates 
combined <- left_join(pheno_pcs, eigenvec, by=c("FID", "IID"))

# Save out the pheno_pcs file 
write.table(combined, file = "EXAMPLE_covariateFile_10PCs.txt", row.names=FALSE, na="", quote = FALSE, sep="\t")

```

Output Example :
```
FID	IID	PAT	MAT	SEX	PHENO	AGE	EDUCATION	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
001_10	001_10	0	0	1	1	40	2	0.00551348	-0.00146207	0.00468116	0.00435601	0.00743252	-0.00315913	-0.00430362	-0.00242458	-0.00223629	0.00381585
002_08	002_08	0	0	2	1	72	1	0.00405574	-0.000306384	-0.00226897	0.00280894	-0.00584511	-0.00364357	-0.000183857	-0.000944975
```

---
<a id="2"></a>
## 2. Running a GWAS in PLINK

### Background Information
- What is a GWAS?
	-   A GWAS is a popular method to test common genetic variants across the exomes/genotypes/genomes of many individuals to identify genotype-phenotype associations
    -   “If certain genetic variations are found to be significantly more frequent in people with the disease compared to people without disease, the variations are said to be ‘associated’ with the disease” [[Source]](https://www.genome.gov/about-genomics/fact-sheets/Genome-Wide-Association-Studies-Fact-Sheet)
-   What kind of information can we learn from this type of analysis?
    -   Loci that pass a genome-wide significance threshold for a large group of people that are more frequent in the disease phenotype than in controls

You can choose to run a GWAS in RVTests
A GWAS will only be run post-imputation and if it has been QC'd, as shown in Modules I and II.

### Running a GWAS in PLINK 
```bash
plink --bfile EXAMPLE_IMPUTED \
--logistic --ci 0.95 \
--hide-covar \
--covar EXAMPLE_covariateFile_10PCs.txt \
--covar-name AGE,EDUCATION,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out IMPUTED.test

# The --hide-covar option removes each individual test from the output
```

Output Example:
```
 CHR           SNP         BP   A1       TEST    NMISS         OR       SE      L95      U95         STAT            P
   1      1:753538     753538    A        ADD     2000      1.263    1.416  0.07868    20.28        0.165       0.8689
   1      1:754396     754396    C        ADD     2000      2.212    1.231   0.1979    24.71       0.6446       0.5192
   1      1:758591     758591    C        ADD     2000      2.212    1.231   0.1979    24.71       0.6446       0.5192
   1      1:761531     761531    A        ADD     2000      1.263    1.416  0.07868    20.28        0.165       0.8689
```

### Running a GWAS in RVTests
```

							Under development -MM.

```

-   What are the limitations with a GWAS?
    -   Only looks at common variants (MAF >0.05), and cannot provide information on rare variants
    - Only nominates loci and the nearest protein-coding gene (though, functionally, this might not be the causal gene)
    - Explain a number of small-effect associations
    - Literature shows that loci nominated in one population might not translate to other populations

---
<a id="3"></a>
## 3. Generating Summary Statistics

### Background Information
-   What are summary statistics? What information can be found?
	- Summary statistics have aggregate p-values and association data for every variant analyzed in a study
    -   In addition to SNP and P-value, information such as effect allele, directionality, frequency, standard error, etc. can also be found in summary statistics
 - Can summary statistics be made public?
	 - Yes! In fact, you should strive to make all your summary statistics public
    

-   Where can I find more summary statistics?
	- NHGRI-EBI have a website where they curate GWAS summary statistics, available for public download with their accompanying manuscripts at [GWAS Catalog](https://www.ebi.ac.uk/gwas/)


### Generating Summary Statistics
```

							Under development -MM.

```

---
<a id="4"></a>
## 4. Prep for FUMA GWAS

### Background Information

From their [website](https://fuma.ctglab.nl/):
-   “FUMA is a platform that can be used to annotate, prioritize, visualize and interpret GWAS results”  
-   “The ***SNP2GENE*** function takes GWAS summary statistics as an input, and provides extensive functional annotation for all SNPs in genomic areas identified by lead SNPs”

-   How do I get started with FUMA?
	-   Register with your email address
    -   Summary statistics ready with rsIDs (can merge with HRC)
    
-   What kind of information can I get?
	-  GWAS-based and gene-based Manhattan and QQ plots
	-  Gene set analysis
	-  Tissue expression analysis

SOMETHING HERE ABOUT LOOK AT THE SLIDES TO GET STARTED? SAVE IMAGES IN THE README?
```

							Under development -MM.

```


---
<a id="5"></a>
## 5. Data Visualization: Scree Plot

### Background Information
-   What is a Scree plot?
    -   Method to visualize how much each PC individually contributes to the amount of variance we see in the data
-   Why do we do this?
	-   Good way to check that the first few PCs you will use later for your GWAS do indeed capture most of the variance we see  
	-   If you do not see this, this indicates issues in the data or the processing of the data upstream
-   How do you read a Scree plot?
	-   The percent of variance is on the Y-axis, with the number of PCs at the bottom
	-  The higher the point, the more % variance can be explained by that PC
	- Typically for a GWAS, the plateau happens at PC4 or PC5

### Making a Scree Plot in R 
```R
# Making a Scree plot 

# Read in the PCA Eigenvalues and Eigenvectors 
print("Read in pca.eigenval files from PLINK")
eigenval <- read.delim("IPDGC_all_to_include_plink2.eigenval", sep ="\t", header = F, stringsAsFactors = F)

# Update column names
colnames(eigenval)[1] <- "Eigenvalues"
eigenval$PC <- as.numeric(rownames(eigenval))
eigenval$VarianceExplained <- eigenval$Eigenvalues/sum(eigenval$Eigenvalues)*100

# Keeping only the first 10 PCs
eigenval2 <- head(eigenval,10)

# Generating the plot 
scree <- ggplot(data = eigenval2, aes(x = PC, y = VarianceExplained)) +
  geom_line() + 
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_x_continuous(name = "Principal Components", breaks = seq(0,10,1), limits = c(NA,10)) +
  scale_y_continuous(name = "Percent (%) Variance Explained", breaks = seq(0,50,5), limits = c(0,50)) +
  ggtitle("Scree Plot: \n IPDGC Samples \n (2,1478 Cases; 2,4388 Controls)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Saving out the plot as PDF and JPEG
ggsave("ScreePlot_IPDGC_unrelated.pdf", scree, width = 5, height = 3.5, units = "in")
ggsave("ScreePlot_IPDGC_unrelated.jpg", scree, width = 5, height = 3.5, units = "in")
```
---
<a id="6"></a>
## 6. Data Visualization: QQ Plot

### Background Information
-   What is a QQ plot?
	-   A Q–Q (quantile-quantile) plot is a probability plot that compares two probability distributions by plotting their quantiles against each other
-   Why do we do this?
	-   “To search for evidence of systematic bias (from unrecognized population structure, analytical approach, genotyping artifacts, etc.)” [[Source](https://www.broadinstitute.org/diabetes-genetics-initiative/plotting-genome-wide-association-results)] otherwise known as the lambda (λ), or the genomic inflation factor
	-   λ scales with sample size, so also important to report the λ1000  (the inflation factor of a study of 1000 cases and 1000 controls)
	-  Visually represents the extent that the observed distribution of the test statistic follows the expected (null) distribution
- How do you read a QQ Plot?
	- The expected -log (a normal distribution) is drawn as a line, and you want to see how many of your points follow this distribution
    -   If the λ1000 is…
	    - <1 : Study is underpowered
	    - ~1: Data observed is close to expected values
	    - \>1: Study/data has inconsistencies
	- You can generate a QQ plot pre- or post-imputation
	-  Post-imputation will be closer to the normal distribution and will look better

### Making a QQ Plot in R [UNIMPUTED DATA]

```R
# Making a QQ Plot
  # You can make a QQ plot pre- or post- imputation
  # This script shows you how to plot it pre-imputation 
  # FUMA GWAS will generate a QQ Plot post imputation for you :) 

# Read in the .fam file
fam.file  <- fread("UMIMPUTED_PD_september_2018_no_duplicates.fam")

# Separate out the cases
case <- sum(fam.file[,6] == 2)

# Separate out the controls 
control <-  sum(fam.file[,6] == 1)

# Read in the .assoc file  
assoc.df <- fread("UMIMPUTED_PD_september_2018_no_duplicates.assoc.assoc")

# Non-factorization is important
assoc.df = assoc.df[complete.cases(assoc.df), ]

# Remove any NA lines
assoc.df$CHISQ <- qchisq(assoc.df$P, 1, lower.tail=FALSE)

# qchisq(assoc.df$P,1,lower.tail=FALSE) can convert p-value (even < 5.5e-17) to chisq
    # while qchisq(1-assoc.df$P,1) fails to convert small p-value (Inf in this case) 

# Calculate lambda 
lambda <- median(assoc.df$CHISQ) / qchisq(0.5, 1)

# Calculate the lambda 1000 
lambda1000 <- 1 + (lambda - 1) * (1 / case + 1 / control) * 500
 
# qchisq(0.5,1) will give 0.4549364
# some people directly / 0.454 which may cause slightly different result
# chisq = z**2
# z <- qnorm(p-value/2)
 
# Make the QQ Plot
qqplot <- qq(assoc.df$P, main = "Q-Q plot of Pre-Imputed GWAS p-values")

# Print out the lambda statistics 
lambda1 <- paste('The lambda value is: ', lambda, sep="")
lambda10001 <- paste('The lambda 1000 value is: ', lambda1000, sep="")

# Write out the lambda statistics to an external .txt file 
lambda_stats<-file("lambda_stats.txt")
writeLines(c(paste(lambda1, lambda10001)), lambda_stats)
close(lambda_stats)

# Save the QQ Plot as PDF and JPEG
ggsave("QQPlot_UNIMPUTED.pdf", qqplot, width = 20, height = 20, units = "in")
ggsave("QQPlot_UNIMPUTED.jpg", qqplot, width = 5, height = 5, units = "in")
```

---
<a id="7"></a>
## 7. Data Visualization: Manhattan Plot

### Background Information 

```

							Under development -MM.

```
