# Load libraries
library(survival)
library(survminer)
library(plyr)
library(tibble)
library(forcats)

# pull data from cbioportal
cbio_OS_data <- read.delim(file = "KM_Plot__Overall_Survival_.txt")# 0=living 1=deceased
cbio_AKR1C3amp_data <- read.delim(file = "alterations_across_samples.tsv")
colnames(cbio_OS_data)
colnames(cbio_AKR1C3amp_data)
# merge data by patient ID
merge_AKR1C3_OS <- merge(cbio_OS_data, cbio_AKR1C3amp_data, by = "Patient.ID")
merge_AKR1C3_OS$OS_STATUS <- gsub(":.*","", merge_AKR1C3_OS$OS_STATUS) # remove everything after :
# subset merged data
colnames(merge_AKR1C3_OS)
sub_merge_AKR1C3_OS <-merge_AKR1C3_OS[c(1,3:4,6:8)]
colnames(sub_merge_AKR1C3_OS)

# change sructure
str(sub_merge_AKR1C3_OS)
sub_merge_AKR1C3_OS$OS_STATUS <- as.numeric(as.character(sub_merge_AKR1C3_OS$OS_STATUS))
sub_merge_AKR1C3_OS$AKR1C3 <- as.character(as.factor(sub_merge_AKR1C3_OS$AKR1C3))
str(sub_merge_AKR1C3_OS)
# remove unwanted data
sub_merge_AKR1C3_OS<-sub_merge_AKR1C3_OS[!(sub_merge_AKR1C3_OS$AKR1C3=="not profiled"),] # remove rows with AKR1C3 'not profiled'
sub_merge_AKR1C3_OS<-na.omit(sub_merge_AKR1C3_OS, cols = c("OS_STATUS", "OS_MONTHS"))# remove rows without OS data
sub_merge_AKR1C3_OS$AKR1C3 <- revalue(sub_merge_AKR1C3_OS$AKR1C3, c("C154F"="AMP"))
# change AKR1C3 back to factor; want 2 factors
str(sub_merge_AKR1C3_OS)
sub_merge_AKR1C3_OS$AKR1C3 <- as.factor(as.character(sub_merge_AKR1C3_OS$AKR1C3))
str(sub_merge_AKR1C3_OS)
#flip levels, we want amp as 2 and no alteration as 1 
sub_merge_AKR1C3_OS$AKR1C3 <- fct_rev(sub_merge_AKR1C3_OS$AKR1C3)


# cox model
#hazard ratio - positive value indicated greater risk (opposite of linear model)
aka_cox <- coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = sub_merge_AKR1C3_OS)
coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = sub_merge_AKR1C3_OS)

beta_coef_i_want <- coef(summary(aka_cox))[1]
hazard_ratio_i_want <- coef(summary(aka_cox))[2]
p_value_i_want <- coef(summary(aka_cox))[5]
summary(aka_cox)
#95% CI 1.58-6.50

## plot 
sfit <- survfit(Surv(OS_MONTHS, OS_STATUS)~AKR1C3, data=sub_merge_AKR1C3_OS)
plot(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("No alteration", "Amplified"), legend.title="AKR1C3",  
           palette=c("orchid2", "dodgerblue2"), 
           title="AKR1C3 Amplification Risk in Lobular Breast Cancer", 
           #risk.table.height=.15, 
           data = sub_merge_AKR1C3_OS)

head(sub_merge_AKR1C3_OS$AKR1C3)
