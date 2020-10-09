R.Version()


# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# Pull data from cbioportal
META <- read.delim(file = "meta/combined_study_clinical_data.tsv")# 0=living 1=deceased
AKR1C3_MICROARRAY <- read.delim(file = "meta/mRNA expression z-Scores relative to all samples (log microarray).txt")
OS <- read.delim(file = "KM_Plot__Overall_Survival_.txt")
# Merge data by Patient.ID
names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="SAMPLE_ID"] <- "Sample.ID"
names(AKR1C3_MICROARRAY)[names(AKR1C3_MICROARRAY)=="STUDY_ID"] <- "Study.ID"
MERGE <- merge(META, AKR1C3_MICROARRAY, by = "Sample.ID")
MERGE_2 <- merge(MERGE, OS, by ="Patient.ID")
# Subset merged data
colnames(MERGE_2)
subset <- MERGE_2[c(1,294,296:297)]
colnames(subset)
subset$OS_STATUS <- gsub(":.*","", subset$OS_STATUS) # remove everything after :
# Change sructure
str(subset)
subset$OS_STATUS <- as.numeric(as.character(subset$OS_STATUS))
str(subset)
# Remove unwanted data

subset<-na.omit(subset, cols = c("OS_STATUS", "OS_MONTHS"))# remove rows without OS data

#visualize expression data distribution
colnames(subset)
hist(subset[,2])
str(subset)
aka_cox_expression <- coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = subset)
coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = subset)
beta_coef_i_want <- coef(summary(aka_cox_expression))[1]
hazard_ratio_i_want <- coef(summary(aka_cox_expression))[2]
p_value_i_want <- coef(summary(aka_cox_expression))[5]
summary(aka_cox_expression)
#95% CI 1.01-1.70





# Add col classifying AKR1C3 exp as +/-
subset$AKR1C3_HL <- ifelse(subset$AKR1C3 <= -0, 'Negative',
                          ifelse(subset$AKR1C3 >=0 & subset$AKR1C3 <=10, 'positive',
                                 ifelse(subset$AKR1C3 >=10, 'High', 'something else')
                          ))
# Change sructure
str(subset)
subset$AKR1C3_HL <- as.factor(as.character(subset$AKR1C3_HL))
str(subset)


# Flip factor levels, WANT: no_alteration/negative = 1 and amp/positive = 2 
#subset$AKR1C3_HL <- fct_rev(subset$AKR1C3_HL)

# cox model
#hazard ratio - positive value indicated greater risk (opposite of linear model)
aka_cox_exp <- coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3_HL, data = subset)
coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3_HL, data = subset)
beta_coef_i_want <- coef(summary(aka_cox_exp))[1]
hazard_ratio_i_want <- coef(summary(aka_cox_exp))[2]
p_value_i_want <- coef(summary(aka_cox_exp))[5]
summary(aka_cox_exp)
#95% CI 1.01-1.70

## plot 
sfit <- survfit(Surv(OS_MONTHS, OS_STATUS)~AKR1C3_HL, data=subset)
plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Negative", "Positive"), legend.title="AKR1C3",  
           palette=c("orchid2", "dodgerblue2"), 
           title="AKR1C3 Expression Risk in Lobular Breast Cancer", 
           #risk.table.height=.15, 
           data = subset)

head(sub_merge_AKR1C3_OS$AKR1C3)