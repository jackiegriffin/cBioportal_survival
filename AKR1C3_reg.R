
library("survival")
library("survminer")
library("plyr")

# pull data from cbioportal
cbio_OS_data <- read.delim(file = "KM_Plot__Overall_Survival_.txt")# 0=living 1=deceased
cbio_AKR1C3amp_data <- read.delim(file = "alterations_across_samples.tsv")

cnv <- read.delim(file = "cna.txt")


names(cnv)[names(cnv)=="SAMPLE_ID"] <- "Sample.ID"
colnames(cnv)
colnames(merge_AKR1C3_OS)
colnames(cbio_OS_data)
colnames(cbio_AKR1C3amp_data)

# merge data by patient ID
merge_AKR1C3_OS <- merge(cbio_OS_data, cbio_AKR1C3amp_data, by = "Patient.ID")

merge_cnv <-merge(merge_AKR1C3_OS, cnv, by="Sample.ID")

merge_cnv$OS_STATUS <- gsub(":.*","", merge_cnv$OS_STATUS) # remove everything after :
sub_merge_cnv<-sub_merge_cnv[!(sub_merge_cnv$AKR1C3.x=="not profiled"),] # remove rows with AKR1C3 'not profiled'
write.csv(sub_merge_cnv_omit, file = "s.csv")
colnames(merge_cnv)

# subset data
sub_merge_cnv <-merge_cnv[c(1:8,13:14)]
new <- read.csv(file = "s.csv", stringsAsFactors = FALSE, row.names = 1)
colnames(sub_AKR1C3_OS)

# need to remove data with NAs
sub_merge_cnv_omit<-na.omit(sub_merge_cnv, cols = c("OS_STATUS", "OS_MONTHS"))

str(sub_merge_cnv_omit)

sub_merge_cnv_omit$OS_STATUS <- as.numeric(as.character(sub_merge_cnv_omit$OS_STATUS))
sub_merge_cnv_omit$AKR1C3.y <- as.numeric(as.integer(sub_merge_cnv_omit$AKR1C3.y))


## model ----- using complete
relation <- lm(Complete$OS_MONTHS ~ Complete$AKR1C3.x)

print(relation)


# fit <- survfit(Surv(Days, status) ~ Tx, data = zr75_3b)

Complete$AKR1C3 <- revalue(Complete$AKR1C3, c("C154F"="AMP"))
Complete["1880", "AKR1C3"] <- AMP

coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3.x, data = new)

# make vector of patient_ids
colnames(AKR1C3)


patient_id <- reg_sub$Patient.ID
AKR1C3 <- reg_sub[c(2:3)]

