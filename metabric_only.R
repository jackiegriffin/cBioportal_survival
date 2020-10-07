library("survival")
library("survminer")
library("plyr")
library("tibble")
install.packages("tibble")
# METABRIC ONLY - entire dataset
# microarray of AKR1C3 expression data 
#AKR1C3_AMP <- read.delim(file = "alterations_across_samples_meta.tsv")# 0=living 1=deceased

META_METABRIC <- read.delim(file = "brca_metabric_clinical_data.tsv")# 0=living 1=deceased
META_METABRIC$Overall.Survival.Status <- gsub(":.*","", META_METABRIC$Overall.Survival.Status) # remove everything after :



AKR1C3_MICROARRAY <- read.delim(file = "mRNA expression z-Scores relative to all samples (log microarray).txt")
colnames(META_METABRIC)
colnames(AKR1C3_MICROARRAY)

names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="SAMPLE_ID"] <- "Sample.ID"
names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="STUDY_ID"] <- "Study.ID"

merge <- merge(META_METABRIC, AKR1C3_MICROARRAY, by = "Sample.ID")
merge$Overall.Survival.Status <- as.numeric(as.character(merge$Overall.Survival.Status))

colnames(merge)
str(merge)



# cox model
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ AKR1C3_HL, data = merge)


# cox model with nodes
coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ AKR1C3 + Lymph.nodes.examined.positive, data = merge)



###add new column of data classifying AKR1C3 expression as high and low
merge$AKR1C3_HL <- ifelse(merge$AKR1C3 <= -0, 'Negative',
                  ifelse(merge$AKR1C3 >=0 & merge$AKR1C3 <=10, 'positive',
                         ifelse(merge$AKR1C3 >=10, 'High', 'something else')
                  ))
hist(merge[,38])

## plot 
sfit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status)~AKR1C3_HL, data=merge)
plot(sfit)
ggsurvplot(sfit, data = merge)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Negative", "positive"), legend.title="AKR1C3 expression",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for AKR1C3 espression (METABRIC Breast Cancer)", 
           risk.table.height=.15, data = merge)


