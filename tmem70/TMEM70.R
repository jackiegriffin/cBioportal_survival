# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# Pull data from cbioportal
ANGPT2_AMP <- read.delim(file = "ANGPT2/ANGPT2_alterations_across_samples.tsv")# 0=living 1=deceased
ANGPT2_META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
  ANGPT2_MICROARRAY_METABRIC <- read.delim(file = "ANGPT2/ANGPT2_mRNA expression (microarray).txt")
# Merge data
names(ANGPT2_MICROARRAY_METABRIC)[names(ANGPT2_MICROARRAY_METABRIC)=="SAMPLE_ID"] <- "Sample.ID"
colnames(ANGPT2_MICROARRAY_METABRIC)
colnames(ANGPT2_AMP)

names(ANGPT2_MICROARRAY_METABRIC)[names(ANGPT2_MICROARRAY_METABRIC)=="STUDY_ID"] <- "Study.ID"
MERGE_ANGPT2 <- merge(ANGPT2_AMP, ANGPT2_MICROARRAY_METABRIC, by = "Sample.ID")
colnames(MERGE_ANGPT2)
colnames(ANGPT2_META_METABRIC)
MERGE_ANGPT2_2 <- merge(MERGE_ANGPT2, ANGPT2_META_METABRIC, by = "Patient.ID")
colnames(MERGE_ANGPT2_2)
head(MERGE_ANGPT2_2)
MERGE_ANGPT2_2$Overall.Survival.Status <- gsub(":.*","", MERGE_ANGPT2_2$Overall.Survival.Status) # remove everything after :
# Subset merged data
colnames(MERGE_ANGPT2_2)
head(MERGE_ANGPT2_2)
subset_ANGPT2 <- MERGE_ANGPT2_2[c(1,4,13,38:39)] # ANGPT2.y is expression
colnames(subset_ANGPT2)

# Change sructure
str(subset_ANGPT2)
subset_ANGPT2$Overall.Survival.Status <- as.numeric(as.character(subset_ANGPT2$Overall.Survival.Status))
str(subset_ANGPT2)
subset_ANGPT2$Altered <- as.factor(as.integer(subset_ANGPT2$Altered))


colnames(subset_ANGPT2)
# Add col classifying AKR1C3 exp as +/-
#visualize expression data distribution
colnames(subset_ANGPT2)
hist(subset_ANGPT2[,3])

subset_ANGPT2$ANGPT2.y_HL <- ifelse(subset_ANGPT2$ANGPT2.y <= 6.8, 'Negative',
                                      ifelse(subset_ANGPT2$ANGPT2.y >=6.8 & subset_ANGPT2$ANGPT2.y <=13, 'Positive',
                                             ifelse(subset_ANGPT2$ANGPT2.y >=13, 'High', 'something else')
                                      ))
# Change sructure of new col
str(subset_ANGPT2)
subset_ANGPT2$ANGPT2.y_HL <- as.factor(as.character(subset_ANGPT2$ANGPT2.y_HL))
str(subset_ANGPT2)
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ ANGPT2.y_HL, data = subset_ANGPT2)

#ANGPT2_plot_alt<-coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Altered, data = subset_ANGPT2)
#colnames(subset)

## plot ----
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~ANGPT2.y_HL, data=subset_ANGPT2)
#plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Low<6.8", "high>6.8"), legend.title="ANGPT2",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for ANGPT2 mRNA expression (METABRIC Breast Cancer)", 
           #risk.table.height=.15,
           data = subset_ANGPT2)
