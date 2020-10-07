library("survival")
library("survminer")
library("plyr")
library("tibble")
install.packages("tibble")
# METABRIC ONLY - entire dataset
# microarray of AKR1C3 expression data 
#AKR1C3_AMP <- read.delim(file = "alterations_across_samples_meta.tsv")# 0=living 1=deceased

META <- read.delim(file = "meta/combined_study_clinical_data.tsv")# 0=living 1=deceased
colnames(META)
MERGE <- merge(META, AKR1C3_MICROARRAY, by = "Sample.ID")
colnames(OS)
MERGE_2 <- merge(MERGE, OS, by ="Patient.ID")
colnames(MERGE_2)

subset <- MERGE_2[c(1,294,296:297)]

#META_METABRIC$Overall.Survival.Status <- gsub(":.*","", META_METABRIC$Overall.Survival.Status) # remove everything after :
OS <- read.delim(file = "KM_Plot__Overall_Survival_.txt")
AKR1C3_MICROARRAY <- read.delim(file = "meta/mRNA expression z-Scores relative to all samples (log microarray).txt")
colnames(OS)
colnames(AKR1C3_MICROARRAY)
names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="SAMPLE_ID"] <- "Sample.ID"
names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="STUDY_ID"] <- "Study.ID"

merge <- merge(META_METABRIC, AKR1C3_MICROARRAY, by = "Sample.ID")

merge <- merge(META_METABRIC, AKR1C3_MICROARRAY, by = "Sample.ID")
merge$Overall.Survival.Status <- as.numeric(as.character(merge$Overall.Survival.Status))

colnames(merge)
str(merge)



# cox model
colnames(subset)
head(subset)
str(subset)
subset$OS_STATUS <- gsub(":.*","", subset$OS_STATUS) # remove everything after :
subset$OS_STATUS <- as.numeric(as.character(subset$OS_STATUS))

coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3_HL, data = subset)

###add new column of data classifying AKR1C3 expression as high and low
subset$AKR1C3_HL <- ifelse(subset$AKR1C3 <= 0, 'Negative',
                          ifelse(subset$AKR1C3 >=0 & subset$AKR1C3 <=10, 'positive',
                                 ifelse(subset$AKR1C3 >=10, 'High', 'something else')
                          ))
hist(subset[,2])

## plot 
sfit <- survfit(Surv(OS_MONTHS, OS_STATUS)~AKR1C3_HL, data=subset)
plot(sfit)
ggsurvplot(sfit, data = subset)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Negative", "positive"), legend.title="AKR1C3 expression",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Overall survival for AKR1C3 espression (METABRIC Breast Cancer)", 
           risk.table.height=.15, data = subset)


