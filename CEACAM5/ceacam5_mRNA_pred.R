# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# Pull data from cbioportal
CEACAM5_AMP <- read.delim(file = "CEACAM5/ceacam5_alterations_across_samples.tsv")# 0=living 1=deceased
#CEACAM5_mRNA_rel <- read.delim(file = "CEACAM5/ceacam5_mRNA expression z-Scores relative to all samples (log microarray).txt")# 0=living 1=deceased
CEACAM5_META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
CEACAM5_MICROARRAY_METABRIC <- read.delim(file = "CEACAM5/ceacam5_mRNA expression (microarray).txt")
# Merge data
names(CEACAM5_MICROARRAY_METABRIC)[names(CEACAM5_MICROARRAY_METABRIC)=="SAMPLE_ID"] <- "Sample.ID"
colnames(CEACAM5_MICROARRAY_METABRIC)
colnames(CEACAM5_AMP)

names(CEACAM5_MICROARRAY_METABRIC)[names(CEACAM5_MICROARRAY_METABRIC)=="STUDY_ID"] <- "Study.ID"
MERGE_CEACAM5 <- merge(CEACAM5_AMP, CEACAM5_MICROARRAY_METABRIC, by = "Sample.ID")
colnames(MERGE_CEACAM5)
colnames(CEACAM5_META_METABRIC)
MERGE_CEACAM5_2 <- merge(MERGE_CEACAM5, CEACAM5_META_METABRIC, by = "Patient.ID")
colnames(MERGE_CEACAM5_2)
head(MERGE_CEACAM5_2)
MERGE_CEACAM5_2$Overall.Survival.Status <- gsub(":.*","", MERGE_CEACAM5_2$Overall.Survival.Status) # remove everything after :
# Subset merged data
colnames(MERGE_CEACAM5_2)
head(MERGE_CEACAM5_2)
subset_ceacam5 <- MERGE_CEACAM5_2[c(1,4,13,38:39)] # ceacam5.y is expression
colnames(subset_ceacam5)

# Change sructure
str(subset_ceacam5)
subset_ceacam5$Overall.Survival.Status <- as.numeric(as.character(subset_ceacam5$Overall.Survival.Status))
str(subset_ceacam5)
subset_ceacam5$Altered <- as.factor(as.integer(subset_ceacam5$Altered))


colnames(subset_ceacam5)
# Add col classifying AKR1C3 exp as +/-
subset_ceacam5$CEACAM5.y_HL <- ifelse(subset_ceacam5$CEACAM5.y <= 5.55, 'Low',
                           ifelse(subset_ceacam5$CEACAM5.y >=5.55 & subset_ceacam5$CEACAM5.y <=13, 'up',
                                  ifelse(subset_ceacam5$CEACAM5.y >=13, 'High', 'something else')
                           ))
# Change sructure of new col
str(subset_ceacam5)
subset_ceacam5$CEACAM5.y_HL <- as.factor(as.character(subset_ceacam5$CEACAM5.y_HL))
str(subset_ceacam5)
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CEACAM5.y_HL, data = subset_ceacam5)

#ceacam5_plot_alt<-coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Altered, data = subset_ceacam5)
#colnames(subset)

## plot ----
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~CEACAM5.y_HL, data=subset_ceacam5)
plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Low", "High"), legend.title="CEACAM5",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for CEACAM5 mRNA expression (METABRIC Breast Cancer)", 
           #risk.table.height=.15,
           data = subset_ceacam5)

## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     mRNA exression relative to all samples      ^^^^^^^^^^^^^^^^^^^^^^^^^^


# Pull data from cbioportal
CEACAM5_AMP <- read.delim(file = "CEACAM5/ceacam5_alterations_across_samples.tsv")# 0=living 1=deceased
CEACAM5_mRNA_rel <- read.delim(file = "CEACAM5/ceacam5_mRNA expression z-Scores relative to all samples (log microarray).txt")# 0=living 1=deceased
CEACAM5_META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
#CEACAM5_MICROARRAY_METABRIC <- read.delim(file = "CEACAM5/ceacam5_mRNA expression (microarray).txt")
# Merge data
names(CEACAM5_mRNA_rel)[names(CEACAM5_mRNA_rel)=="SAMPLE_ID"] <- "Sample.ID"
colnames(CEACAM5_mRNA_rel)
colnames(CEACAM5_AMP)

names(CEACAM5_mRNA_rel)[names(CEACAM5_mRNA_rel)=="STUDY_ID"] <- "Study.ID"
MERGE_CEACAM5 <- merge(CEACAM5_AMP, CEACAM5_mRNA_rel, by = "Sample.ID")
colnames(MERGE_CEACAM5)
colnames(CEACAM5_META_METABRIC)
MERGE_CEACAM5_2 <- merge(MERGE_CEACAM5, CEACAM5_META_METABRIC, by = "Patient.ID")
colnames(MERGE_CEACAM5_2)
head(MERGE_CEACAM5_2)
MERGE_CEACAM5_2$Overall.Survival.Status <- gsub(":.*","", MERGE_CEACAM5_2$Overall.Survival.Status) # remove everything after :
# Subset merged data
colnames(MERGE_CEACAM5_2)
head(MERGE_CEACAM5_2)
subset_ceacam5 <- MERGE_CEACAM5_2[c(1,4,13,38:39)] # ceacam5.y is expression
colnames(subset_ceacam5)

# Change sructure
str(subset_ceacam5)
subset_ceacam5$Overall.Survival.Status <- as.numeric(as.character(subset_ceacam5$Overall.Survival.Status))
str(subset_ceacam5)
subset_ceacam5$Altered <- as.factor(as.integer(subset_ceacam5$Altered))


colnames(subset_ceacam5)
# Add col classifying AKR1C3 exp as +/-
subset_ceacam5$CEACAM5.y_HL <- ifelse(subset_ceacam5$CEACAM5.y <= 0, 'Negative',
                                      ifelse(subset_ceacam5$CEACAM5.y >=0 & subset_ceacam5$CEACAM5.y <=13, 'Positive',
                                             ifelse(subset_ceacam5$CEACAM5.y >=13, 'High', 'something else')
                                      ))
# Change sructure of new col
str(subset_ceacam5)
subset_ceacam5$CEACAM5.y_HL <- as.factor(as.character(subset_ceacam5$CEACAM5.y_HL))
str(subset_ceacam5)
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CEACAM5.y_HL, data = subset_ceacam5)

#ceacam5_plot_alt<-coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Altered, data = subset_ceacam5)
#colnames(subset)

## plot ----
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~CEACAM5.y_HL, data=subset_ceacam5)
plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Negative", "Positive"), legend.title="CEACAM5",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for CEACAM5 mRNA expression (METABRIC Breast Cancer)\nRelative to all samples", 
           #risk.table.height=.15,
           data = subset_ceacam5)
