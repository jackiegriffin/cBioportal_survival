# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# Pull data from cbioportal
AKR1C3_AMP <- read.delim(file = "alterations_across_samples_meta.tsv")# 0=living 1=deceased
META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
AKR1C3_MICROARRAY_METABRIC <- read.delim(file = "mRNA expression z-Scores relative to all samples (log microarray).txt")
# Merge data
names(AKR1C3_MICROARRAY_METABRIC)[names(AKR1C3_MICROARRAY_METABRIC)=="SAMPLE_ID"] <- "Sample.ID"
names(AKR1C3_MICROARRAY_METABRIC)[names(AKR1C3_MICROARRAY_METABRIC)=="STUDY_ID"] <- "Study.ID"
MERGE <- merge(AKR1C3_AMP, AKR1C3_MICROARRAY_METABRIC, by = "Sample.ID")
MERGE_2 <- merge(MERGE, META_METABRIC, by = "Patient.ID")
MERGE_2$Overall.Survival.Status <- gsub(":.*","", MERGE_2$Overall.Survival.Status) # remove everything after :
# Subset merged data
colnames(MERGE_2)
subset <- MERGE_2[c(1,5,11,17:19,22:23,27:28,32,36:37)]
colnames(subset)

# Change sructure
str(subset)
subset$Overall.Survival.Status <- as.numeric(as.character(subset$Overall.Survival.Status))
str(subset)
# Flip factor levels, WANT: no_alteration/negative = 1 and amp/positive = 2 
subset$AKR1C3.x <- fct_rev(subset$AKR1C3.x)
str(subset)
# Remove missing surivival data
subset<-na.omit(subset, cols = c("Overall.Survival..Months.", "Overall.Survival.Status"))# remove rows without OS data

# cox model for gene expression
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ AKR1C3.y, data = subset)
colnames(subset)
hist(subset[,3])

# Add col classifying AKR1C3 exp as +/-
subset$AKR1C3_HL <- ifelse(subset$AKR1C3.y <= -0, 'Negative',
                           ifelse(subset$AKR1C3.y >=0 & subset$AKR1C3.y <=10, 'positive',
                                  ifelse(subset$AKR1C3.y >=10, 'High', 'something else')
                           ))


# Change sructure of new col
str(subset)
subset$AKR1C3_HL <- as.factor(as.character(subset$AKR1C3_HL))
str(subset)

metabric_plot<-coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ AKR1C3_HL, data = subset)
colnames(subset)
## plot 
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~AKR1C3_HL
                , data=subset)
plot(sfit)
ggsurvplot(sfit, data = subset)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Negative", "positive"), legend.title="AKR1C3 expression",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for AKR1C3 espression (METABRIC Breast Cancer)", 
           risk.table.height=.15, data = subset)

# cox model with nodes
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ AKR1C3 + Lymph.nodes.examined.positive, data = merge)



