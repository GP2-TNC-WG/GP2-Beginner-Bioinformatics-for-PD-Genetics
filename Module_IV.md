## Module IV. Polygenic risk score in Parkinson Disease

#### Information

* Created by: GP2 Training and Networking

## Table of Contents

#### [0. Getting Started](#0)

#### [1. GRS versus disease status](#1)

#### [2. GRS versus age at onset](#2)

#### [3. Data visualization: Violin Plots](#3)

#### [4. Data visualization: Quantile plots](#4)

#### [5. ROC calculation ](#5)

#### [6. Data Visualization: ROC plots](#6)

#### [7. Data Visualization: Density plots](#7)

---
<a id="0"></a>
## 0. Getting Started

Genetic risk scores (AKA: polygenic scores, polygenic risk scores, or genome-wide scores) is a summary measure of a set of risk-associated genetic variants and can be easily calculated using PLINK

Files that you will need:

**QC’d Imputed PLINK Binary Files** containing cases and controls (.fam, .bim, .bed files)

- These have gone through Michigan Server imputation and the soft/hard-call QC steps outlined in Module II

**Covariates File** containing the covariates you would like to correct by

- This file at a minimum includes sample information (such as ID, SEX, PHENO) and principal components

**SCORE file** containing  the 90 risk loci (variant identifier, A1, Beta value from Nalls et al.,  2019.)

- See META5_GRS_chr_bp.txt or META5_GRS_RSid.txt depending on your format

This tutorial is written in R, and these are the packages you will need to install/load 
```R
# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyr')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(caret)) install.packages('caret')


# Load the necessary packages 
library(tidyr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(plotROC)
library(caret)
```
Other programs you will need are:
- PLINK v1.9

---
<a id="1"></a>

## GRS versus disease status (Nalls et al., 2019)

### Calculate Score (profile) in PLINK

```
cd Desktop/TEST/
./plink --bfile test --score META5_GRS_chr_bp.txt --out GRS_PD_test
```

Where:
```
test = standard binary file prefix (will point to .bed, .bim, and .fam files)
GRS_PD_test.profile = whatever you want it to be, the output will have the extension .profile
META5_GRS_chr_bp.txt = file with variant-name, allele and score-value
```

### Read PLINK output, merge with covariate file and recode CASE (1) and CONTROL (0)

```
setwd("~/Desktop/TEST/")
temp_data <- read.table("GRS_PD_test.profile", header = T) 
temp_covs <- read.table("test_covs.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO.x - 1
```

### Normalize Score to Z-Score 

```
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
```

### Perform logistic regression adjusted by covariates

```
grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AAO, family="binomial", data = data)
summary(grsTests)
```
---
<a id="2"></a>

## GRS versus age at onset

### Subset ONLY cases perform linear regression adjusted by covariates

```
cases <- subset(data, CASE == 1)
meanPop <- mean(cases$SCORE)
sdPop <- sd(cases$SCORE)
cases$zSCORE <- (cases$SCORE - meanPop)/sdPop
grsTests <- lm(AAO ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cases)
summary(grsTests)
```
---
<a id="3"></a>

## Data visualization - Violin plots

```
data$PHENO[data$CASE ==0] <- "Controls"
data$PHENO[data$CASE ==1] <- "PD"

p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("PD GRS (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("PD_GRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

```
---
<a id="4"></a>

## Data visualization -  Quantile plots

* Make quantiles

```
data$quantile1 <- ifelse(data$zSCORE <= quantile(data$zSCORE)[2], 1, 0)
data$quantile2 <- ifelse(data$zSCORE > quantile(data$zSCORE)[2] & data$zSCORE <= quantile(data$zSCORE)[3], 1, 0)
data$quantile3 <- ifelse(data$zSCORE > quantile(data$zSCORE)[3] & data$zSCORE <= quantile(data$zSCORE)[4], 1, 0)
data$quantile4 <- ifelse(data$zSCORE > quantile(data$zSCORE)[4], 1, 0)
data$quantiles <- 1
data$quantiles[data$quantile2 == 1] <- 2
data$quantiles[data$quantile3 == 1] <- 3
data$quantiles[data$quantile4 == 1] <- 4
quintileTests <- glm(CASE ~ as.factor(data$quantiles) + AAO + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)

```
* Summarize the regression and export a table

```
summary(quintileTests)
summary_stats <- data.frame(summary(quintileTests)$coef[2:4,1:2])
names(summary_stats) <- c("BETA","SE")
summary_stats$QUANTILE <- c("2nd","3rd","4th")
summary_stats[4,] <- c(0,0,"1st")
summary_stats_sorted <- summary_stats[order(summary_stats$QUANTILE),]
write.table(summary_stats_sorted, "quantile_table.csv", quote = F, row.names = F, sep = ",")
```

* Make quantile plot

```
to_plot <- read.table("quantile_table.csv", header = T, sep = ",")
to_plot$low <- to_plot$BETA - (1.96*to_plot$SE)
to_plot$high <- to_plot$BETA + (1.96*to_plot$SE)
plotted <- ggplot(to_plot, aes(QUANTILE, BETA)) + geom_pointrange(aes(ymin = low, ymax = high))
ggsave(plot = plotted, filename = "plotQuantile.png", width = 4, height = 4, units = "in", dpi = 300)
```
---
<a id="5"></a>

### ROC calculation

* Load additional packages

```
packageList <- c("caret","ggplot2","data.table","plotROC")
lapply(packageList, library, character.only = TRUE)
```
* Run regression model 

```
Model <- glm(CASE ~ SCORE, data = data, family = 'binomial')
```

* Make predictions

```
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")
```
---
<a id="6"></a>

### Data visualization - ROC plots

```
overlayedRocs <- ggplot(data, aes(d = CASE, m = probDisease)) + geom_roc(labels = FALSE) + geom_rocci() + style_roc(theme = theme_gray) + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "plotRoc.png", width = 8, height = 5, units = "in", dpi = 300)

```

* Show the confusion matrix (specificity and sensitivity)

```
confMat <- confusionMatrix(data = as.factor(data$predicted), reference = as.factor(data$reported), positive = "DISEASE")
confMat
```

---
<a id="7"></a>

### Data visualization - Density plots

```
densPlot <- ggplot(data, aes(probDisease, fill = CASE, color = CASE)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity.png", width = 8, height = 5, units = "in", dpi = 300)

```
---
<a id="8"></a>
