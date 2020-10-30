# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# Pull data from cbioportal
ANGPTL4_AMP <- read.delim(file = "ANGPTL4/ANGPTL4_alterations_across_samples.tsv")# 0=living 1=deceased
ANGPTL4_META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
  ANGPTL4_MICROARRAY_METABRIC <- read.delim(file = "ANGPTL4/ANGPTL4_mRNA expression z-Scores relative to all samples (log microarray).txt")
# Merge data
names(ANGPTL4_MICROARRAY_METABRIC)[names(ANGPTL4_MICROARRAY_METABRIC)=="SAMPLE_ID"] <- "Sample.ID"
colnames(ANGPTL4_MICROARRAY_METABRIC)
colnames(ANGPTL4_AMP)

names(ANGPTL4_MICROARRAY_METABRIC)[names(ANGPTL4_MICROARRAY_METABRIC)=="STUDY_ID"] <- "Study.ID"
MERGE_ANGPTL4 <- merge(ANGPTL4_AMP, ANGPTL4_MICROARRAY_METABRIC, by = "Sample.ID")
colnames(MERGE_ANGPTL4)
colnames(ANGPTL4_META_METABRIC)
MERGE_ANGPTL4_2 <- merge(MERGE_ANGPTL4, ANGPTL4_META_METABRIC, by = "Patient.ID")
colnames(MERGE_ANGPTL4_2)
head(MERGE_ANGPTL4_2)
MERGE_ANGPTL4_2$Overall.Survival.Status <- gsub(":.*","", MERGE_ANGPTL4_2$Overall.Survival.Status) # remove everything after :
# Subset merged data
colnames(MERGE_ANGPTL4_2)
head(MERGE_ANGPTL4_2)
subset_ANGPTL4 <- MERGE_ANGPTL4_2[c(1,4,13,38:39)] # ANGPTL4.y is expression
colnames(subset_ANGPTL4)

# Change sructure
str(subset_ANGPTL4)
subset_ANGPTL4$Overall.Survival.Status <- as.numeric(as.character(subset_ANGPTL4$Overall.Survival.Status))
str(subset_ANGPTL4)
subset_ANGPTL4$Altered <- as.factor(as.integer(subset_ANGPTL4$Altered))


colnames(subset_ANGPTL4)
# Add col classifying AKR1C3 exp as +/-
#visualize expression data distribution
colnames(subset_ANGPTL4)
hist(subset_ANGPTL4[,3])

subset_ANGPTL4$ANGPTL4.y_HL <- ifelse(subset_ANGPTL4$ANGPTL4.y <= 0, 'Negative',
                                      ifelse(subset_ANGPTL4$ANGPTL4.y >=0 & subset_ANGPTL4$ANGPTL4.y <=13, 'Positive',
                                             ifelse(subset_ANGPTL4$ANGPTL4.y >=13, 'High', 'something else')
                                      ))
# Change sructure of new col
str(subset_ANGPTL4)
subset_ANGPTL4$ANGPTL4.y_HL <- as.factor(as.character(subset_ANGPTL4$ANGPTL4.y_HL))
str(subset_ANGPTL4)
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ ANGPTL4.y_HL, data = subset_ANGPTL4)

#ANGPTL4_plot_alt<-coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Altered, data = subset_ANGPTL4)
#colnames(subset)

## plot ----
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~ANGPTL4.y_HL, data=subset_ANGPTL4)
#plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Low<6.8", "high>6.8"), legend.title="ANGPTL4",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for ANGPTL4 mRNA expression (METABRIC Breast Cancer)", 
           #risk.table.height=.15,
           data = subset_ANGPTL4)
