# Calculating a Genetic Risk Score

- **Author(s):** Sara Bandres-Ciga
- **Date Last Updated:** April. 2020

## Introduction

Genetic risk scores (AKA: polygenic scores, polygenic risk scores, or genome-wide scores) is a summary measure of a set of risk-associated genetic variants and can be easily calculated using PLINK

## GRS versus disease status (Nalls et al., 2019)

### Calculate profile in PLINK

```
plink --bfile test --score META5_GRS_chr_bp.txt --out GRS_PD_test.profile
```

Where:
test = standard binary file prefix (will point to .bed, .bim, and .fam files)
GRS_PD_test.profile = whatever you want it to be, the output will have the extension .profile
META5_GRS_chr_bp.txt = file with variant-name, allele and score-value

### Load R libraries

```
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
```

### Calculate Score in cases versus controls 

```
temp_data <- read.table("GRS_PD_test.profile", header = T) 
temp_covs <- read.table("test_covs.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
```

### Normalize Score to Z-Score 

```
data$zSCORE <- (data$SCORE - meanControls)/sdControls
```

### Normalize Score to Z-Score and perform logistic regression adjusted by covariates

```
grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CONSENSUS_AGE, family="binomial", data = data)
summary(grsTests)
```

## GRS versus age at onset (Blauwendraat et al., 2019)

### Perform linear regression adjusted by covariates

```
cases <- subset(data, PHENO == 1)
meanPop <- mean(cases$SCORE)
sdPop <- sd(cases$SCORE)
cases$zSCORE <- (cases$SCORE - meanPop)/sdPop
grsTests <- lm(AGE_onset ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cases)
summary(grsTests)
```

## Data visualization

### Violin plots

```
data$PHENO[data$PHENO ==1] <- "Controls"
data$PHENO[data$PHENO ==2] <- "PD"

p <- ggplot(all_data, aes(x= reorder(as.factor(PHENO), zSCORE), y=zSCORE, fill=as.factor(PHENO))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("PD GRS (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("PD_GRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

```

### Quantile plots

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
quintileTests <- glm(CASE ~ as.factor(data$quantiles) + as.factor(data$DATASET) + AGE + SEX_COV + PC1 + PC2 + PC3 + PC4 + PC5, family="binomial", data = data)

```
* Sumarize the regression and export a table

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

### ROC calculation

* Load additional packages

```
packageList <- c("caret","ggplot2","data.table","plotROC")
lapply(packageList, library, character.only = TRUE)
```
* Build your basic model as best you can

```
bestModel <- glm(PHENO ~ SCORE, data = data, family = 'binomial')
```

* Make predictions

```
data$predicted <- predict(bestModel, data)
data$probDisease <- predict(bestModel, data, type = "prob")[2]
```

* Make ROC plots

```
overlayedRocs <- ggplot(data, aes(d = PHENO, m = probDisease)) + geom_roc(labels = FALSE) + geom_rocci() + style_roc(theme = theme_gray) + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "plotRoc.png", width = 8, height = 5, units = "in", dpi = 300)
```

* Show the confusion matrix (specificity and sensitivity)

```
confMat <- confusionMatrix(data = as.factor(trained$predicted), reference = as.factor(trained$PHENO), positive = "DISEASE")
confMat
```

### Density plots

```
densPlot <- ggplot(data, aes(probDisease, fill = PHENO, color = PHENO)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity.png", width = 8, height = 5, units = "in", dpi = 300)

```
