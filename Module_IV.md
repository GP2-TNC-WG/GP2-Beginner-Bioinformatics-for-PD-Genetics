# Calculating a Genetic Risk Score

- **Author(s):** Sara Bandres-Ciga
- **Date Last Updated:** April. 2020

## GRS versus disease status (Nalls et al., 2019)

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

### Normalize Score to Z-Score and perform linear regression adjusted by covariates

```
grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CONSENSUS_AGE, family="binomial", data = data)
summary(grsTests)
```

## GRS versus age at onset (Blauwendraat et al., 2019)

### Perform logistic regression adjusted by covariates

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
### Density plots

### Quantile plots

### ROC calculation

```
library("pROC")
rocAuc <- roc(data$label,data$predictedLabel)
auc(rocAuc)
ci(rocAuc, of="auc")
coords(rocAuc, "best")
thresh <- coords(rocAuc, "best")[1]
```

```
png("rocCurve.png", width = 5.5, height = 4, units = "in", res = 300)
rocobj <- plot.roc(data$label, data$predictedLabel,  main="ROC curve", percent=FALSE,  ci=TRUE, print.auc=TRUE, col = "darkred")
ciobj <- ci.se(rocobj, specificities=seq(0, 1, 0.05))
plot(ciobj, type="shape", col=alpha("dark red", 1/2))
dev.off()
```
