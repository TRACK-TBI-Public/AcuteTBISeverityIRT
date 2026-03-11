# TITLE: Accompanies
#	Nelson, L.D. et al. (2025). Beyond mild, moderate, and severe traumatic 
#   brain injury: Modelling severity from clinical, neuroimaging, and blood-based 
#    indicators. eBioMedicine, 121, 106001.
# Script 2/2 - performs analysis of CENTER-TBI sample, scoring it using the
# IRT model established in the TRACK-TBI sample

install.packages("haven")
install.packages("mirt")
install.packages("plotly")
install.packages("reshape2")
install.packages("ggridges")
install.packages("fmsb")
install.packages("DescTools")
install.packages("gridExtra")
install.packages("psych")
install.packages("ggplot2")
install.packages("MplusAutomation")
install.packages("dplyr")

library(haven)
library(mirt)
library(plotly)
library(reshape2)
library(ggridges)
library(fmsb)
library(DescTools)
library(gridExtra)
library(psych)
library(ggplot2)
library(MplusAutomation)
library(dplyr)

#Set the working directory, e.g.,
setwd("I:/FolderName/")

#Open the CENTER-TBI dataset
dat <- read.csv("CENTER_MergedData_Wide.csv")

# Table 2

describe(dat$Age)
quantile(dat$Age, probs = c(0.25, 0.75), na.rm = TRUE)

dat$Sex <- factor(dat$Sex, levels = c(1, 2), labels = c("Male", "Female"))
table(dat$Sex)
prop.table(table(dat$Sex))

## Race variable was not harmonized to TRACK - need to check category labels
#dat$Race[is.na(dat$Race)] <- 88
#dat$Race <- factor(dat$Race, levels = c(1, 2, 3, 4, 5, 6, 7, 88),
#                   labels = c("Indian",
#                              "Alaska Native/Inuit",
#                              "Asian",
#                              "Black",
#                              "Native Hawaiian/Pacific Islander",
#                              "White",
#                              "Mixed race",
#                              "Unknown"))
# Combine race categories
#dat$Race_r <- factor(dat$Race, labels = c("Other/Unknown",
#                              "Other/Unknown",
#                              "Asian",
#                              "Black",
#                              "Other/Unknown",
#                              "White",
#                              "Other/Unknown",
#                              "Other/Unknown"))

#table(dat$Race_r)
#prop.table(table(dat$Race_r))

## No ethnicity variable in CENTER
#dat$Ethnicity <- factor(dat$Ethnicity, levels = c(1, 2, 88),
#                        labels = c("Hispanic", "Non-Hispanic", "Unknown"))
#table(dat$Ethnicity)
#prop.table(table(dat$Ethnicity))

dat$EduYearsOfEducation <- replace(dat$EduYearsOfEducation, dat$EduYearsOfEducation == 88, NA)

describe(dat$EduYearsOfEducation)
quantile(dat$EduYearsOfEducation, probs = c(0.25, 0.75), na.rm = T)

## N/A in CENTER
#dat$HealthInsurance <- factor(dat$HealthInsurance, levels = c(1, 2, 3, 4, 5, 6, 7, 88, 99),
#                              labels = c("self-pay (uninsured)",
#                              "Insurance through a current or former employer (incl. thru family member)",
#                              "Insurance purchased directly from an insurance company or on the health insurance",
#                              "Medicare, for people 65 and older, or people with certain disabilities",
#                              "Medicaid, Medical Assistance, “the State” or any kind of government-assistance plan for low income/disability",
#                              "Medicaid Pending",
#                              "TRICARE, VA or other military health care",
#                              "Unknown",
#                              "Any other type of health insurance or health coverage plan"))
#dat$HealthInsurance_r <- factor(dat$HealthInsurance,
#                              labels = c("Medicaid/Uninsured",
#                                         "Other Insurance",
#                                         "Other Insurance",
#                                         "Other Insurance",
#                                         "Medicaid/Uninsured",
#                                         "Medicaid/Uninsured",
#                                         "Other Insurance",
#                                         "Unknown",
#                                         "Other Insurance"))
#table(dat$HealthInsurance_r)
#prop.table(table(dat$HealthInsurance_r))

describe(dat$gcser)
quantile(dat$gcser, probs = c(0.25, 0.75), na.rm = TRUE)

dat$gcser_r <- ifelse(dat$gcser >= 3 & dat$gcser <= 8, "3-8",
                      ifelse(dat$gcser >= 9 & dat$gcser <= 12, "9-12", "13-15"))
  
table(dat$gcser_r)
prop.table(table(dat$gcser_r))

table(dat$ctpos)
prop.table(table(dat$ctpos))

table(dat$NumberNonreactingPupils)
prop.table(table(dat$NumberNonreactingPupils))

dat$InjCause <- factor(dat$InjCause, levels = c(1, 2, 3, 4, 5, 6, 99),
                       labels = c("Road traffic incident",
                                  "Incidental fall",
                                  "Other non-intentional injury",
                                  "Violence/assault",
                                  "Act of mass violence",
                                  "Suicide attempt",
                                  "Other"))

dat$InjCause_r <- factor(dat$InjCause,
                       labels = c("Motor Vehicle/Traffic Crash",
                                  "Fall",
                                  "Other/Unknown",
                                  "Assault/Violence",
                                  "Other/Unknown",
                                  "Other/Unknown",
                                  "Other/Unknown"))

table(dat$InjCause_r)
prop.table(table(dat$InjCause_r))

dat$levcare <- factor(dat$levcare, levels = c(1, 2, 3),
                      labels = c("Emergency Department", "Inpatient Floor", "ICU"))

table(dat$levcare)
prop.table(table(dat$levcare))

dat$LOCDuration_r <- factor(dat$LOCDuration_R, levels = c(1, 2, 3, 4, 5, 6, 7),
                            labels = c("None",
                                       "< 1 Min",
                                       "1-29 Min",
                                       "30-59 Min",
                                       "1-24 Hours",
                                       "24 Hours - 7 Days",
                                       "> 7 Days"))

table(dat$LOCDuration_r)
prop.table(table(dat$LOCDuration_r))

dat$LOCPTADuration_r <- factor(dat$LOCPTADuration_R, levels = c(1, 2, 3, 4, 5, 6, 7),
                            labels = c("None",
                                       "< 1 Min",
                                       "1-29 Min",
                                       "30-59 Min",
                                       "1-24 Hours",
                                       "24 Hours - 7 Days",
                                       "> 7 Days"))

table(dat$LOCPTADuration_r)
prop.table(table(dat$LOCPTADuration_r))

#

table(dat$GCSEDArrVerbal_R)
table(dat$GCSEDArrMotor_R)
table(dat$GCSEDArrEyes_R)
#Proportion for Table 2
prop.table(table(dat$GCSEDArrVerbal_R))
prop.table(table(dat$GCSEDArrMotor_R))
prop.table(table(dat$GCSEDArrEyes_R))

#Recode letter categories of GCS to missing so IRT model can treat as ordinal
dat$GCSEDArrVerbal_R <- as.numeric(as.character(dat$GCSEDArrVerbal_R))
dat$GCSEDArrVerbal_R[dat$GCSEDArrVerbal_R == "T"] <- NA
#dat$GCSEDArrVerbal_R <- 6 - dat$GCSEDArrVerbal_R

dat$GCSEDArrMotor_R <- as.numeric(as.character(dat$GCSEDArrMotor_R))
dat$GCSEDArrMotor_R[dat$GCSEDArrMotor_R == "P"] <- NA
#dat$GCSEDArrMotor_R <- 7 - dat$GCSEDArrMotor_R

dat$GCSEDArrEyes_R <- as.numeric(as.character(dat$GCSEDArrEyes_R))
dat$GCSEDArrEyes_R[dat$GCSEDArrEyes_R == "S"] <- NA
#dat$GCSEDArrEyes_R <- 5 - dat$GCSEDArrEyes_R

table(dat$GCSEDArrVerbal_R)
table(dat$GCSEDArrMotor_R)
table(dat$GCSEDArrEyes_R)
prop.table(table(dat$GCSEDArrVerbal_R))


# Missing data

sum(is.na(dat$CT_SkullFx_R))/length(dat$CT_SkullFx_R)
sum(is.na(dat$CT_Contusion_R))/length(dat$CT_Contusion_R)
sum(is.na(dat$CT_Shear_R))/length(dat$CT_Shear_R)
sum(is.na(dat$CT_ExtraaxHematoma_R))/length(dat$CT_ExtraaxHematoma_R)
sum(is.na(dat$CT_EDH_R))/length(dat$CT_EDH_R)
sum(is.na(dat$CT_SDH_R))/length(dat$CT_SDH_R)
sum(is.na(dat$CT_SAH_R))/length(dat$CT_SAH_R)
sum(is.na(dat$CT_IVH_R))/length(dat$CT_IVH_R)
#sum(is.na(dat$CT_Edema_R))/length(dat$CT_Edema_R)
#sum(is.na(dat$CT_DownwardHerniation_R))/length(dat$CT_DownwardHerniation_R)
#sum(is.na(dat$CT_UpwardCerebellarHerniation_R))/length(dat$CT_UpwardCerebellarHerniation_R)
#sum(is.na(dat$CT_DuretHemorrhage_R))/length(dat$CT_DuretHemorrhage_R)
sum(is.na(dat$CT_MidlineShift_R))/length(dat$CT_MidlineShift_R)

sum(is.na(dat$NumberNonreactingPupils))/length(dat$NumberNonreactingPupils)

sum(is.na(dat$LOCDuration_R))/length(dat$LOCDuration_R)
sum(is.na(dat$LOCPTADuration_R))/length(dat$LOCPTADuration_R)

sum(is.na(dat$GCSEDArrVerbal_R_R))/length(dat$GCSEDArrVerbal_R_R)
sum(is.na(dat$GCSEDArrMotor_R))/length(dat$GCSEDArrMotor_R)
sum(is.na(dat$GCSEDArrEyes_R))/length(dat$GCSEDArrEyes_R)

sum(is.na(dat$gfap_cat))/length(dat$gfap_cat)
sum(is.na(dat$uchl_cat))/length(dat$uchl_cat)
sum(is.na(dat$nse_cat))/length(dat$nse_cat)
sum(is.na(dat$s100_cat))/length(dat$s100_cat)
#sum(is.na(dat$crp_cat))/length(dat$crp_cat)

# Variables for the CENTER-TBI IRT model

allvar <- c("Subject_ID", "CT_SkullFx_R", "CT_Contusion_R", "CT_Shear_R", "CT_ExtraaxHematoma_R", "CT_EDH_R", "CT_SDH_R",
            "CT_SAH_R", "CT_IVH_R", "CT_MidlineShift_R", 
            "GCSEDArrVerbal_R", "GCSEDArrMotor_R", "GCSEDArrEyes_R",
            "gfap_cat", "uchl_cat", "nse_cat", "s100_cat", "NumberNonreactingPupils")

allvar <- dat[allvar] # Create dataset with just IRT variables

## Explore missingness patterns - Added by Lin for CENTER data
install.packages("naniar")
library(naniar)
vis_miss(allvar)
miss_var_summary(allvar)
print(miss_case_summary(allvar),n=30)
#There are 7 missing 94.7% of data - they may be causing an issue with mirt so will 
# exclude them below.

##Lin Added missing_prop to remove the 7 people with 94.7% missingness
missing_prop <- rowMeans(is.na(allvar))
allvardata <- allvar[missing_prop <= 0.94,] # Remove cases with 94% missing
  #Compare original to new data - confirmed went from 4507 to 4500
str(allvar)
str(allvardata)


# Score CENTER TBI cases on the IRT model established in the TRACK dataset###

#Open TRACK data
datTRACK <- read.csv("dataset_1168_1773155061.csv")

table(datTRACK$GCSEDArrVerbal_R)
table(datTRACK$GCSEDArrMotor_R)
table(datTRACK$GCSEDArrEyes_R)

datTRACK$GCSEDArrVerbal_R <- as.numeric(as.character(datTRACK$GCSEDArrVerbal_R))
datTRACK$GCSEDArrVerbal_R[datTRACK$GCSEDArrVerbal_R == "T"] <- NA
datTRACK$GCSEDArrVerbal_R <- 6 - datTRACK$GCSEDArrVerbal_R

datTRACK$GCSEDArrMotor_R <- as.numeric(as.character(datTRACK$GCSEDArrMotor_R))
datTRACK$GCSEDArrMotor_R[datTRACK$GCSEDArrMotor_R == "P"] <- NA
datTRACK$GCSEDArrMotor_R <- 7 - datTRACK$GCSEDArrMotor_R

datTRACK$GCSEDArrEyes_R <- as.numeric(as.character(datTRACK$GCSEDArrEyes_R))
datTRACK$GCSEDArrEyes_R[datTRACK$GCSEDArrEyes_R == "S"] <- NA
datTRACK$GCSEDArrEyes_R <- 5 - datTRACK$GCSEDArrEyes_R

table(datTRACK$GCSEDArrVerbal_R)
table(datTRACK$GCSEDArrMotor_R)
table(datTRACK$GCSEDArrEyes_R)

# Fit TRACK model in the TRACK dataset

allvarTRACK <- c("Subject_ID", "CT_SkullFx_R", "CT_Contusion_R", "CT_Shear_R", "CT_ExtraaxHematoma_R", "CT_EDH_R", "CT_SDH_R",
                 "CT_SAH_R", "CT_IVH_R", "CT_Edema_R", "CT_DownwardHerniation_R", "CT_UpwardCerebellarHerniation_R",
                 "CT_DuretHemorrhage_R", "CT_MidlineShift_R", "LOCDuration_R", "LOCPTADuration_R",
                 "GCSEDArrVerbal_R", "GCSEDArrMotor_R", "GCSEDArrEyes_R",
                 "gfap_cat", "uchl_cat", "nse_cat", "s100_cat", "crp_cat", "NumberNonreactingPupils")

allvarTRACK <- datTRACK[allvarTRACK] # Create dataset with just IRT variables
allvarTRACKdata <- allvarTRACK[rowSums(is.na(allvarTRACK)) != ncol(allvarTRACK), ] # Remove cases with all missing

allvarTRACK.uni <- mirt(allvarTRACKdata[ , 2:25], itemtype = c(rep("graded", 24)), model = 1, SE = T)

coef(allvarTRACK.uni, IRTpars = T, simplify = T)

# Score CENTER-TBI subjects on the TRACK-TBI model

# First, create the missing variables in the CENTER dataset so response patterning will match
#  Get item names from training set
#train_items <- allvarTRACK
# Find missing items in new dataset vs. training TRACK set
missing_items <- setdiff(colnames(allvarTRACK), colnames(allvardata))
#Add missing item columns filled with NA
allvardata[, missing_items] <- NA
#Reorder columns to match training data
allvardataT2C <- allvardata[, colnames(allvarTRACK)]

#score Center cases on TRACK parameters, I called it T2C b/c we used TRACK-TBI model 2 score CENTER-TBI data
allvar.uni.scores.T2C <- fscores(allvarTRACK.uni, response.pattern = allvardataT2C[,2:25],
                                 method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)

allvar.uni.scores.T2C <- data.frame(allvar.uni.scores.T2C)
allvardataT2C$EAPT2C <- allvar.uni.scores.T2C$F1
allvardataT2C$SET2C <- allvar.uni.scores.T2C$SE_F1

# Merge dataset with IRT scores with original dataset (final dataset called "alldat")

scores <- allvardataT2C %>% select(Subject_ID, EAPT2C, SET2C)
alldatT2C <- merge(dat, scores, by = "Subject_ID")

## Plots for CENTER scoring of TRACK IRT model

#histogram of Acute TBI Severity IRT scores 
ggplot(alldatT2C, aes(x = EAPT2C)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "white") +
  labs(title = "IRT Scores CENTER data [TRACK model]", x = "TBI Severity", y = "Frequency") +
  # scale_x_continuous(limits = c(-4, 4)) +
  theme_minimal()

# Number unique IRT scores (CENTER data, TRACK model)
unique_counts <- alldatT2C %>%
  group_by(TBISev_GCS, ctpos) %>%
  dplyr::summarize(unique_count = n_distinct(EAPT2C))
print(unique_counts)

unique_counts <- alldatT2C %>%
  group_by(ctpos) %>%
  dplyr::summarize(unique_count = n_distinct(EAPT2C))
print(unique_counts)

unique_counts <- alldatT2C %>%
  dplyr::summarize(unique_count = n_distinct(EAPT2C))
print(unique_counts)

# Scatterplot GCS-based classification vs. IRT scores (CENTER data, TRACK model)
alldatT2C <- alldatT2C %>% filter( !is.na(TBISev_GCS))

alldatT2C$TBISev_GCS <- factor(alldatT2C$TBISev_GCS,
                               levels = c(1, 2, 3),
                               labels = c("13-15", "9-12", "3-8"))

g1 <- ggplot(alldatT2C, aes(x = EAPT2C, y = TBISev_GCS)) +
  geom_density_ridges(
    jittered_points = TRUE, position = "raincloud",
    alpha = 0.7, scale = 0.9) +
  labs(title = "GCS-Based TBI Classification", x = "TBI Severity (CENTER data, TRACK model)", y = "") +
  theme_bw() +
  scale_x_continuous(name="TBI Severity", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

g1

# Predictive validity analyses

#alldatT2C$goseALL.T1 <- ifelse(alldat$goseALL.T1 < 10, alldat$goseALL.T1, NA)
alldatT2C$goseALL.T3 <- ifelse(alldatT2C$goseALL.T3 < 10, alldatT2C$goseALL.T3, NA)

# Create categorical outcome variables

alldatT2C$dead <- ifelse(alldatT2C$goseALL.T3 == 1, 1, 0)
alldatT2C$unfavorable <- ifelse(alldatT2C$goseALL.T3 < 4, 1, 0)
alldatT2C$incomplete <- ifelse(alldatT2C$goseALL.T3 < 8, 1, 0)

alldatT2C$TBISev_GCS <- factor(alldatT2C$TBISev_GCS)
alldatT2C$TBISev_GCS <- relevel(alldatT2C$TBISev_GCS, ref = "13-15")

#Predictive models vs GCS 

logit.death <- glm(dead ~ TBISev_GCS, data = alldatT2C, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAPT2C <- glm(dead ~ TBISev_GCS + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.death.EAPT2C)
anova(logit.death.EAPT2C, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAPT2C, test = "Chisq")
NagelkerkeR2(logit.death.EAPT2C)
PseudoR2(logit.death.EAPT2C, which = "CoxSnell")
exp(logit.death.EAPT2C$coefficients)
exp(cbind(OR = coef(logit.death.EAPT2C), confint(logit.death.EAPT2C)))
with(logit.death.EAPT2C, null.deviance - deviance)
with(logit.death.EAPT2C, df.null - df.residual)
with(logit.death.EAPT2C, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

logit.unfavorableT2C <- glm(unfavorable ~ TBISev_GCS, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C)
anova(logit.unfavorableT2C, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C)
PseudoR2(logit.unfavorableT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorableT2C), confint(logit.unfavorableT2C)))

logit.unfavorableT2C.EAP <- glm(unfavorable ~ TBISev_GCS + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C.EAP)
anova(logit.unfavorableT2C.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorableT2C, logit.unfavorableT2C.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C.EAP)
PseudoR2(logit.unfavorableT2C.EAP, which = "CoxSnell")
exp(logit.unfavorableT2C.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorableT2C.EAP), confint(logit.unfavorableT2C.EAP)))
with(logit.unfavorableT2C.EAP, null.deviance - deviance)
with(logit.unfavorableT2C.EAP, df.null - df.residual)
with(logit.unfavorableT2C.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

logit.incompleteT2C <- glm(incomplete ~ TBISev_GCS, data = alldatT2C, family = "binomial")
summary(logit.incompleteT2C)
anova(logit.incompleteT2C, test = "Chisq")
NagelkerkeR2(logit.incompleteT2C)
PseudoR2(logit.incompleteT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.incompleteT2C), confint(logit.incompleteT2C)))

logit.incompleteT2C.EAP <- glm(incomplete ~ TBISev_GCS + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.incompleteT2C.EAP)
anova(logit.incompleteT2C.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.incompleteT2C, logit.incompleteT2C.EAP, test = "Chisq")
NagelkerkeR2(logit.incompleteT2C.EAP)
PseudoR2(logit.incompleteT2C.EAP, which = "CoxSnell")
exp(logit.incompleteT2C.EAP$coefficients)
exp(cbind(OR = coef(logit.incompleteT2C.EAP), confint(logit.incompleteT2C.EAP)))
with(logit.incompleteT2C.EAP, null.deviance - deviance)
with(logit.incompleteT2C.EAP, df.null - df.residual)
with(logit.incompleteT2C.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

#Predictive models vs IMPACT

logit.death <- glm(dead ~ IMPACTCoreLPMort, data = alldatT2C, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAPT2C <- glm(dead ~ IMPACTCoreLPMort + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.death.EAPT2C)
anova(logit.death.EAPT2C, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAPT2C, test = "Chisq")
NagelkerkeR2(logit.death.EAPT2C)
PseudoR2(logit.death.EAPT2C, which = "CoxSnell")
exp(logit.death.EAPT2C$coefficients)
exp(cbind(OR = coef(logit.death.EAPT2C), confint(logit.death.EAPT2C)))
with(logit.death.EAPT2C, null.deviance - deviance)
with(logit.death.EAPT2C, df.null - df.residual)
with(logit.death.EAPT2C, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

logit.unfavorableT2C <- glm(unfavorable ~ IMPACTCoreLPUnfav, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C)
anova(logit.unfavorableT2C, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C)
PseudoR2(logit.unfavorableT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorableT2C), confint(logit.unfavorableT2C)))

logit.unfavorableT2C.EAP <- glm(unfavorable ~ IMPACTCoreLPUnfav + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C.EAP)
anova(logit.unfavorableT2C.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorableT2C, logit.unfavorableT2C.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C.EAP)
PseudoR2(logit.unfavorableT2C.EAP, which = "CoxSnell")
exp(logit.unfavorableT2C.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorableT2C.EAP), confint(logit.unfavorableT2C.EAP)))
with(logit.unfavorableT2C.EAP, null.deviance - deviance)
with(logit.unfavorableT2C.EAP, df.null - df.residual)
with(logit.unfavorableT2C.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death <- glm(dead ~ IMPACTExtendedLPMort, data = alldatT2C, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAPT2C <- glm(dead ~ IMPACTExtendedLPMort  + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.death.EAPT2C)
anova(logit.death.EAPT2C, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAPT2C, test = "Chisq")
NagelkerkeR2(logit.death.EAPT2C)
PseudoR2(logit.death.EAPT2C, which = "CoxSnell")
exp(logit.death.EAPT2C$coefficients)
exp(cbind(OR = coef(logit.death.EAPT2C), confint(logit.death.EAPT2C)))
with(logit.death.EAPT2C, null.deviance - deviance)
with(logit.death.EAPT2C, df.null - df.residual)
with(logit.death.EAPT2C, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

logit.unfavorableT2C <- glm(unfavorable ~ IMPACTExtendedLPUnfav, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C)
anova(logit.unfavorableT2C, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C)
PseudoR2(logit.unfavorableT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorableT2C), confint(logit.unfavorableT2C)))

logit.unfavorableT2C.EAP <- glm(unfavorable ~ IMPACTExtendedLPUnfav + EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C.EAP)
anova(logit.unfavorableT2C.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorableT2C, logit.unfavorableT2C.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C.EAP)
PseudoR2(logit.unfavorableT2C.EAP, which = "CoxSnell")
exp(logit.unfavorableT2C.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorableT2C.EAP), confint(logit.unfavorableT2C.EAP)))
with(logit.unfavorableT2C.EAP, null.deviance - deviance)
with(logit.unfavorableT2C.EAP, df.null - df.residual)
with(logit.unfavorableT2C.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

#IRT scores only

logit.death <- glm(dead ~ EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.unfavorableT2C <- glm(unfavorable ~ EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.unfavorableT2C)
anova(logit.unfavorableT2C, test = "Chisq")
NagelkerkeR2(logit.unfavorableT2C)
PseudoR2(logit.unfavorableT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorableT2C), confint(logit.unfavorableT2C)))

logit.incompleteT2C <- glm(incomplete ~ EAPT2C, data = alldatT2C, family = "binomial")
summary(logit.incompleteT2C)
anova(logit.incompleteT2C, test = "Chisq")
NagelkerkeR2(logit.incompleteT2C)
PseudoR2(logit.incompleteT2C, which = "CoxSnell")
exp(cbind(OR = coef(logit.incompleteT2C), confint(logit.incompleteT2C)))



## START OF CODE TO FIT CENTER MODEL IN CENTER DATA NOT USING THIS FOR REVISION

allvarC <- c("Subject_ID", "CT_SkullFx_R", "CT_Contusion_R", "CT_Shear_R", "CT_ExtraaxHematoma_R", "CT_EDH_R", "CT_SDH_R",
            "CT_SAH_R", "CT_IVH_R", "CT_MidlineShift_R", "LOCDuration_R",
            "GCSEDArrVerbal_R", "GCSEDArrMotor_R", "GCSEDArrEyes_R",
            "gfap_cat", "uchl_cat", "nse_cat", "s100_cat", "NumberNonreactingPupils")

allvarC <- dat[allvarC] # Create dataset with just IRT variables

missing_prop <- rowMeans(is.na(allvarC))
allvarC_R <- allvarC[missing_prop <= 0.94,] # Remove cases with 94% missing
#Compare original to new data - confirmed went from 4507 to 4500
str(allvarC_R)

# Write Mplus data file
##For CENTER THIS currently COMES FROM AN SPSS DATA PREP SCRIPT

#prepareMplusData(
#  allvar,
#  filename = "ShareDataMplus.dat",
#  inpfile = TRUE)

# Unidimensional IRT model

allvarC.uni <- mirt(allvarC_R[ , 2:19], itemtype = c(rep("graded", 18)), model = 1, SE = T)

coef(allvarC.uni, IRTpars = T, simplify = T)

# Estimate IRT scores for full model

allvarC.uni.scores <- fscores(allvarC.uni, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
allvarC.uni.scores <- data.frame(allvarC.uni.scores)
allvarC_R$EAP <- allvarC.uni.scores$F1
allvarC_R$SE <- allvarC.uni.scores$SE_F1

# Merge dataset with IRT scores with original dataset (final dataset called "alldat")

scores <- allvarC_R %>% select(Subject_ID, EAP, SE)
alldatC <- merge(alldatT2C, scores, by = "Subject_ID")

#Check correlation between EAP scores from Center-derived model and TRACK-derived model
# r = .98
cor(alldatC$EAPT2C, alldatC$EAP, use="complete.obs",method="pearson")

## Plots

#Scatterplot of CENTER IRT scores using the CENTER vs. TRACK model
g2 <- ggplot(alldatC, aes(x = EAP, y = EAPT2C)) +
  geom_point(alpha =  0.5) + 
  labs(title = "Acute TBI Severity IRT Scores (CENTER-TBI Data)", x = "CENTER Model", y = "TRACK Model (External Validation)") +
  theme_bw() +
  scale_x_continuous(name="CENTER Model", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

g2

#Put histogram of the two EAP scores from CENTER side-by-side
library(patchwork) 

hist1 <- ggplot(alldatC, aes(x = EAP)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "white") +
  labs(title = "IRT Scores (CENTER data, CENTER model)", x = "TBI Severity", y = "Frequency") +
 # scale_x_continuous(limits = c(-4, 4)) +
  theme_minimal()
hist2 <- ggplot(alldatC, aes(x = EAPT2C)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "white") +
  labs(title = "IRT Scores (CENTER data, TRACK model)", x = "TBI Severity", y = "Frequency") +
  # scale_x_continuous(limits = c(-4, 4)) +
  theme_minimal()
grid.arrange(hist1,hist2,ncol=2)


###LIN LEFT OFF HERE- NEED TO REFERENCE THE RIGHT DATA BELOW 

# Extract information for each domain and entire test

HeadCT <- 0
Clinical.All <- 0
Clinical.GCS <- 0
Clinical.LOC.Pupil <- 0
Biomarkers <- 0

Theta <- seq(-3, 3, 0.1)

for(i in 1:9){
  HeadCT <- HeadCT + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in c(10:13,18)){
  Clinical.All <- Clinical.All + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in c(10,18)){
  Clinical.PTA.LOC.Pupil <- Clinical.PTA.LOC.Pupil + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in 11:13){
  Clinical.GCS <- Clinical.GCS + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in 14:17){
  Biomarkers <- Biomarkers + iteminfo(extract.item(allvar.uni, i), Theta)
}

Test <- testinfo(allvar.uni, Theta)

# Information plot (no test information)

InfoData <- cbind(Theta, HeadCT, Clinical.PTA.LOC.Pupil, Clinical.GCS, Biomarkers)
InfoData <- data.frame(InfoData)
names(InfoData) <- c("Theta", "Head CT", "Clinical-LOC/Pupil", "Clinical-GCS", "Biomarkers")

InfoDataLong <- melt(InfoData, id.vars=c("Theta"))
names(InfoDataLong) <- c("Theta", "Domain", "Information")

thinned <- floor(seq(from = 1, to = dim(InfoDataLong)[1], length = 70))

ggplot(data = InfoDataLong, aes(x = Theta, y = Information, color = Domain, shape = Domain)) +
  geom_line(aes(color = Domain), linewidth = 0.75) +
  geom_point(data = InfoDataLong[thinned,], aes(shape = Domain), size = 2.5) + 
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Information by Domain") +
  scale_color_manual(name = "Domain",
                     values = c("#0000CC", "#006633", "#33FF99", "#CC0000")) +
  theme_bw()

# Information plot (include test information)

InfoData <- cbind(Theta, HeadCT, Clinical.PTA.LOC.Pupil, Clinical.GCS, Biomarkers, Test)
InfoData <- data.frame(InfoData)
names(InfoData) <- c("Theta", "Head CT", "Clinical-LOCPupil", "Clinical-GCS", "Biomarkers", "Test")

InfoDataLong <- melt(InfoData, id.vars=c("Theta"))
names(InfoDataLong) <- c("Theta", "Domain", "Information")

thinned <- floor(seq(from = 1, to = dim(InfoDataLong)[1], length = 70))

ggplot(data = InfoDataLong, aes(x = Theta, y = Information, color = Domain, shape = Domain)) +
  geom_line(aes(color = Domain), linewidth = 0.75) +
  geom_point(data = InfoDataLong[thinned,], aes(shape = Domain), size = 2.5) + 
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Information by Domain") +
  scale_color_manual(name = "Domain",
                     values = c("#0000CC", "#006633", "#33FF99", "#CC0000", "black")) +
  theme_bw()

# Add domains one by one

TestGCS <- Clinical.GCS
TestGCSCT <- Clinical.GCS + HeadCT
TestGCSCTLOCPupil <- Clinical.GCS + HeadCT + Clinical.PTA.LOC.Pupil
TestGCSCTLOCPupilBio <- Clinical.GCS + HeadCT + Clinical.PTA.LOC.Pupil + Biomarkers

SE.TestGCS <- 1/sqrt(TestGCS)
SE.TestGCSCT <- 1/sqrt(TestGCSCT)
SE.TestGCSCTLOCPupil <- 1/sqrt(TestGCSCTLOCPupil)
SE.TestGCSCTLOCBioPupil <- 1/sqrt(TestGCSCTLOCPupilBio)

InfoDataIncremental <- cbind(Theta, TestGCS, TestGCSCT, TestGCSCTLOCPupil, TestGCSCTLOCPupilBio)
InfoDataIncremental <- data.frame(InfoDataIncremental)
names(InfoDataIncremental) <- c("Theta", "GCS", "GCS + Head CT", "GCS + Head CT + LOC/Pupil", "GCS + Head CT + LOC/Pupil + Biomarkers")
InfoDataIncremental <- melt(InfoDataIncremental, id.vars=c("Theta"))
names(InfoDataIncremental) <- c("Theta", "Domain", "Information")

SEDataIncremental <- cbind(Theta, SE.TestGCS, SE.TestGCSCT, SE.TestGCSCTLOCPupil, SE.TestGCSCTLOCBioPupil)
SEDataIncremental <- data.frame(SEDataIncremental)
names(SEDataIncremental) <- c("Theta", "GCS", "GCS + Head CT", "GCS + Head CT + LOC/Pupil", "GCS + Head CT + LOC/Pupil + Biomarkers")
SEDataIncremental <- melt(SEDataIncremental, id.vars=c("Theta"))
names(SEDataIncremental) <- c("Theta", "Domain", "SE")

thinned1 <- floor(seq(from = 1, to = dim(InfoDataIncremental)[1], length = 70))

ggplot(data = InfoDataIncremental, aes(x = Theta, y = Information, color = Domain)) +
  geom_line(aes(color = Domain), linewidth = 0.75) +
  geom_point(data = InfoDataIncremental[thinned1,], aes(), size = 2.5) + 
  scale_x_continuous(name = "TBI Severity", limits = c(-4, 4), n.breaks = 10) +
  labs(title = "Incremental Information by Domain") +
  scale_color_manual(name = "Domain",
                     values = c("#0000CC", "#006633", "#33FF99", "#CC0000")) +
  theme_bw() 

ggplot(data = InfoDataIncremental, aes(x = Theta, y = Information, color = Domain)) +
  geom_line(aes(color = Domain), linewidth = 0.75) +
  geom_line(data = SEDataIncremental, aes(x = Theta, y = SE*0.85, color = Domain), linewidth = 0.75, linetype = "dashed") + # Secondary y-axis scaled
  scale_y_continuous(
    name = "Information",
    sec.axis = sec_axis(~ ./0.85, name = "Standard Error")
  ) +
 # geom_point(data = InfoDataIncremental[thinned1,], aes(), size = 2.5) + 
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Incremental Information by Domain") +
  scale_color_manual(name = "Domain",
                     values = c("#0000CC", "#006633", "#33FF99", "#CC0000")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title=element_blank())

InfoData2 <- cbind(Theta, HeadCT, Clinical.All, Biomarkers)
InfoData2 <- data.frame(InfoData2)
names(InfoData2) <- c("Theta", "Head CT", "Clinical", "Biomarkers")

InfoDataLong2 <- melt(InfoData2, id.vars=c("Theta"))
names(InfoDataLong2) <- c("Theta", "Domain", "Information")

thinned2 <- floor(seq(from = 1, to = dim(InfoDataLong2)[1], length = 70))

ggplot(data = InfoDataLong2, aes(x = Theta, y = Information, color = Domain)) +
  geom_line(aes(color = Domain), linewidth = 0.75) +
  #geom_point(data = InfoDataLong2[thinned2,], aes(), size = 2.5) + 
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Information by Domain") +
  scale_color_manual(name = "Domain",
                     values = c("#0000CC", "#999900", "#CC0000")) +
  theme_bw() 

# Extract information for each item

CT_SkullFx <- iteminfo(extract.item(allvar.uni, 1), Theta)
CT_Contusion <- iteminfo(extract.item(allvar.uni, 2), Theta)
CT_Shear <- iteminfo(extract.item(allvar.uni, 3), Theta)
CT_ExtraaxHematoma <- iteminfo(extract.item(allvar.uni, 4), Theta)
CT_EDH <- iteminfo(extract.item(allvar.uni, 5), Theta)
CT_SDH <- iteminfo(extract.item(allvar.uni, 6), Theta)
CT_SAH <- iteminfo(extract.item(allvar.uni, 7), Theta)
CT_IVH <- iteminfo(extract.item(allvar.uni, 8), Theta)
#CT_Edema <- iteminfo(extract.item(allvar.uni, 9), Theta)
#CT_DownwardHerniation <- iteminfo(extract.item(allvar.uni, 10), Theta)
#CT_UpwardCereHerniation <- iteminfo(extract.item(allvar.uni, 11), Theta)
#CT_DuretHemorrhage <- iteminfo(extract.item(allvar.uni, 12), Theta)
CT_MidlineShift <- iteminfo(extract.item(allvar.uni, 9), Theta)
LOCDuration <- iteminfo(extract.item(allvar.uni, 10), Theta)
#LOCPTADuration <- iteminfo(extract.item(allvar.uni, 15), Theta)
GCSEDArrVerbal_R <- iteminfo(extract.item(allvar.uni, 11), Theta)
GCSEDArrMotor_R <- iteminfo(extract.item(allvar.uni, 12), Theta)
GCSEDArrEyes_R <- iteminfo(extract.item(allvar.uni, 13), Theta)
gfap_cat <- iteminfo(extract.item(allvar.uni, 14), Theta)
uchl_cat <- iteminfo(extract.item(allvar.uni, 15), Theta)
nse_cat <- iteminfo(extract.item(allvar.uni, 16), Theta)
s100_cat <- iteminfo(extract.item(allvar.uni, 17), Theta)
#crp_cat <- iteminfo(extract.item(allvar.uni, 23), Theta)
NumberNonreactingPupils <- iteminfo(extract.item(allvar.uni, 18), Theta)
Test <- testinfo(allvar.uni, Theta)

InfoData <- cbind(Theta, CT_SkullFx, CT_Contusion, CT_Shear, CT_ExtraaxHematoma,
                  CT_EDH, CT_SDH, CT_SAH, CT_IVH,
                  CT_MidlineShift, LOCDuration, GCSEDArrVerbal_R, GCSEDArrMotor_R,
                  GCSEDArrEyes_R, NumberNonreactingPupils, gfap_cat, uchl_cat, nse_cat, s100_cat)
InfoData <- data.frame(InfoData)
names(InfoData) <- c("Theta", "Skull Fracture", "Contusion", "Shear", "Extraaxial Hematoma", "EDH", "SDH",
                     "SAH", "IVH", "Midline Shift", "LOC Duration",
                     "GCS-Verbal", "GCS-Motor", "GCS-Eyes", "Pupil Reactivity",
                     "GFAP", "UCH-L1", "NSE", "S100B")

InfoDataLong <- melt(InfoData, id.vars=c("Theta"))
names(InfoDataLong) <- c("Theta", "Item", "Information")

thinned <- floor(seq(from = 1, to = dim(InfoDataLong)[1], length = 70))

ggplot(data = InfoDataLong, aes(x = Theta, y = Information, color = Item)) +
  geom_line(aes(linetype = Item), linewidth = 0.75) +
  scale_linetype_manual(values=c("solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
#                                 "dotted",
#                                 "dotdash",
#                                 "longdash",
#                                 "twodash",
                                 "solid",
                                 "dashed",
#                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
#                                 "longdash",
                                 "twodash")) +
  scale_color_manual(values = c("#006666",
                                "#009999",
                                "#00CCCC",
                                "#00FFFF",
                                "#99FFFF",
                                "#003366",
                                "#004C99",
                                "#0066CC",
#                                "#0080FF",
#                                "#66B2FF",
#                                "#0000FF",
#                                "#6666FF",
                                "#CCCCFF",
                                "#666600",
#                                "#999900",
                                "#CCCC00",
                                "#DAA520",
                                "#B8860B",
                                "#FFD700",
                                "#660000",
                                "#990000",
                                "#CC0000",
#                                "#FF0000",
                                "#FF3333")) +
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Item Information", size = 2) +
  theme_bw() +
  theme(legend.position = "bottom")

thinned2 <- floor(seq(from = 1, to = dim(InfoDataLong)[1], length = 500))

ggplot(data = InfoDataLong, aes(x = Theta, y = Information, color = Item)) +
  geom_line(aes(linetype = Item), linewidth = 0.75) +
  geom_point(data = InfoDataLong[thinned2,], aes(shape = Item), size = 2.5) + 
  scale_linetype_manual(values=c("solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
#                                 "dotted",
#                                 "dotdash",
#                                 "longdash",
#                                 "twodash",
                                 "solid",
                                 "dashed",
#                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
#                                 "dotdash",
                                 "solid")) +
  scale_color_manual(values = c("#006666",
                                "#009999",
                                "#00CCCC",
                                "#00FFFF",
                                "#99FFFF",
                                "#003366",
                                "#004C99",
                                "#0066CC",
#                                "#0080FF",
#                                "#66B2FF",
#                                "#0000FF",
#                                "#6666FF",
                                "#CCCCFF",
                                "#666600",
#                                "#999900",
                                "#CCCC00",
                                "#FFFF33",
                                "#FFFF99",
                                "#660000",
                                "#990000",
                                "#CC0000",
                                "#FF0000",
                                "#FF9999",
#                                "darkgrey",
                                "black")) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4,
                                5, 6, 7, 8, 9,
                                10, 11, 12, 13, 14,
                                15, 16, 17, 18, 19,
                                20, 21, 22, 23, 24)) +
  scale_x_continuous(name = "TBI Severity", limits = c(-3, 3), n.breaks = 10) +
  labs(title = "Item Information") +
  theme_bw() +
  theme(legend.position = "bottom")

# Predictive validity analyses

#alldat$TBISev_Civilian4grp <- factor(alldat$TBISev_Civilian4grp,
#                                     levels = c(1, 2, 3, 4),
#                                     labels = c("u-MTBI", "c-MTBI", "Mod.", "Severe"))

alldat$TBISev_GCS <- factor(alldat$TBISev_GCS,
                                     levels = c(1, 2, 3),
                                     labels = c("13-15", "9-12", "3-8"))

#alldat$TBISev_VA3grp <- factor(alldat$TBISev_VA3grp,
#                            levels = c(1, 2, 3),
#                            labels = c("Mild", "Mod.", "Severe"))


alldat$goseALL.T1 <- ifelse(alldat$goseALL.T1 < 10, alldat$goseALL.T1, NA)
alldat$goseALL.T3 <- ifelse(alldat$goseALL.T3 < 10, alldat$goseALL.T3, NA)

alldat$goseTBI.T1 <- ifelse(alldat$goseTBI.T1 < 10, alldat$goseTBI.T1, NA)
alldat$goseTBI.T3 <- ifelse(alldat$goseTBI.T3 < 10, alldat$goseTBI.T3, NA)

# Create categorical outcome variables

alldat$dead <- ifelse(alldat$goseALL.T3 == 1, 1, 0)
alldat$unfavorable <- ifelse(alldat$goseALL.T3 < 4, 1, 0)
alldat$incomplete <- ifelse(alldat$goseALL.T3 < 8, 1, 0)

# Filter out NA values
filtered_data <- alldat %>% 
  filter(!is.na(goseALL.T3))

scatter <- ggplot(filtered_data, aes(x = EAP, y = as.factor(goseALL.T3))) +
  geom_density_ridges(
    jittered_points = TRUE, position = "raincloud",
    alpha = 0.7, scale = 0.9) +
  labs(title = "IRT Scores at Each Ordinal GOSE Score", x = "TBI Severity IRT Score", y = "Ordinal GOSE") +
  theme_bw() +
  scale_x_continuous(name="TBI Severity IRT Score", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

scatter


#alldat <- alldat %>% filter( !is.na(TBISev_Civilian4grp))
alldat <- alldat %>% filter( !is.na(TBISev_GCS))
#alldat <- alldat %>% filter( !is.na(TBISev_VA3grp))


unique_counts <- alldat %>%
  group_by(TBISev_GCS, ctpos) %>%
  dplyr::summarize(unique_count = n_distinct(EAP))
print(unique_counts)

unique_counts <- alldat %>%
  group_by(ctpos) %>%
  dplyr::summarize(unique_count = n_distinct(EAP))
print(unique_counts)

unique_counts <- alldat %>%
  dplyr::summarize(unique_count = n_distinct(EAP))
print(unique_counts)

g1 <- ggplot(alldat, aes(x = EAP, y = TBISev_GCS)) +
  geom_density_ridges(
    jittered_points = TRUE, position = "raincloud",
    alpha = 0.7, scale = 0.9) +
  labs(title = "GCS-Based TBI Classification", x = "TBI Severity", y = "") +
  theme_bw() +
  scale_x_continuous(name="TBI Severity", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

g1

## Logistic regressions predicting dichotomized GOSE as death, unfavorable outcome, complete recovery
#  Models run hierchically with single-predictors first, followed by the addition of TBI Severity IRT scores

alldat$TBISev_GCS <- factor(alldat$TBISev_GCS)
alldat$TBISev_GCS <- relevel(alldat$TBISev_GCS, ref = "13-15")

logit.death <- glm(dead ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAP <- glm(dead ~ TBISev_GCS + EAP, data = alldat, family = "binomial")
summary(logit.death.EAP)
anova(logit.death.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAP, test = "Chisq")
NagelkerkeR2(logit.death.EAP)
PseudoR2(logit.death.EAP, which = "CoxSnell")
exp(logit.death.EAP$coefficients)
exp(cbind(OR = coef(logit.death.EAP), confint(logit.death.EAP)))
with(logit.death.EAP, null.deviance - deviance)
with(logit.death.EAP, df.null - df.residual)
with(logit.death.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death <- glm(dead ~ IMPACTCoreLPMort, data = alldat, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAP <- glm(dead ~ IMPACTCoreLPMort + EAP, data = alldat, family = "binomial")
summary(logit.death.EAP)
anova(logit.death.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAP, test = "Chisq")
NagelkerkeR2(logit.death.EAP)
PseudoR2(logit.death.EAP, which = "CoxSnell")
exp(logit.death.EAP$coefficients)
exp(cbind(OR = coef(logit.death.EAP), confint(logit.death.EAP)))
with(logit.death.EAP, null.deviance - deviance)
with(logit.death.EAP, df.null - df.residual)
with(logit.death.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death <- glm(dead ~ IMPACTExtendedLPMort, data = alldat, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))

logit.death.EAP <- glm(dead ~ IMPACTExtendedLPMort + EAP, data = alldat, family = "binomial")
summary(logit.death.EAP)
anova(logit.death.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.death, logit.death.EAP, test = "Chisq")
NagelkerkeR2(logit.death.EAP)
PseudoR2(logit.death.EAP, which = "CoxSnell")
exp(logit.death.EAP$coefficients)
exp(cbind(OR = coef(logit.death.EAP), confint(logit.death.EAP)))
with(logit.death.EAP, null.deviance - deviance)
with(logit.death.EAP, df.null - df.residual)
with(logit.death.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death <- glm(dead ~ EAP, data = alldat, family = "binomial")
summary(logit.death)
anova(logit.death, test = "Chisq")
NagelkerkeR2(logit.death)
PseudoR2(logit.death, which = "CoxSnell")
exp(cbind(OR = coef(logit.death), confint(logit.death)))


logit.unfavorable <- glm(unfavorable ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.unfavorable)
anova(logit.unfavorable, test = "Chisq")
NagelkerkeR2(logit.unfavorable)
PseudoR2(logit.unfavorable, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable), confint(logit.unfavorable)))

logit.unfavorable.EAP <- glm(unfavorable ~ TBISev_GCS + EAP, data = alldat, family = "binomial")
summary(logit.unfavorable.EAP)
anova(logit.unfavorable.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorable, logit.unfavorable.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.EAP)
PseudoR2(logit.unfavorable.EAP, which = "CoxSnell")
exp(logit.unfavorable.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.EAP), confint(logit.unfavorable.EAP)))
with(logit.unfavorable.EAP, null.deviance - deviance)
with(logit.unfavorable.EAP, df.null - df.residual)
with(logit.unfavorable.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable <- glm(unfavorable ~ IMPACTCoreLPUnfav, data = alldat, family = "binomial")
summary(logit.unfavorable)
anova(logit.unfavorable, test = "Chisq")
NagelkerkeR2(logit.unfavorable)
PseudoR2(logit.unfavorable, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable), confint(logit.unfavorable)))

logit.unfavorable.EAP <- glm(unfavorable ~ IMPACTCoreLPUnfav + EAP, data = alldat, family = "binomial")
summary(logit.unfavorable.EAP)
anova(logit.unfavorable.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorable, logit.unfavorable.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.EAP)
PseudoR2(logit.unfavorable.EAP, which = "CoxSnell")
exp(logit.unfavorable.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.EAP), confint(logit.unfavorable.EAP)))
with(logit.unfavorable.EAP, null.deviance - deviance)
with(logit.unfavorable.EAP, df.null - df.residual)
with(logit.unfavorable.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable <- glm(unfavorable ~ IMPACTExtendedLPUnfav, data = alldat, family = "binomial")
summary(logit.unfavorable)
anova(logit.unfavorable, test = "Chisq")
NagelkerkeR2(logit.unfavorable)
PseudoR2(logit.unfavorable, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable), confint(logit.unfavorable)))

logit.unfavorable.EAP <- glm(unfavorable ~ IMPACTExtendedLPUnfav + EAP, data = alldat, family = "binomial")
summary(logit.unfavorable.EAP)
anova(logit.unfavorable.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.unfavorable, logit.unfavorable.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.EAP)
PseudoR2(logit.unfavorable.EAP, which = "CoxSnell")
exp(logit.unfavorable.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.EAP), confint(logit.unfavorable.EAP)))
with(logit.unfavorable.EAP, null.deviance - deviance)
with(logit.unfavorable.EAP, df.null - df.residual)
with(logit.unfavorable.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable <- glm(unfavorable ~ EAP, data = alldat, family = "binomial")
summary(logit.unfavorable)
anova(logit.unfavorable, test = "Chisq")
NagelkerkeR2(logit.unfavorable)
PseudoR2(logit.unfavorable, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable), confint(logit.unfavorable)))

logit.incomplete <- glm(incomplete ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.incomplete)
anova(logit.incomplete, test = "Chisq")
NagelkerkeR2(logit.incomplete)
PseudoR2(logit.incomplete, which = "CoxSnell")
exp(cbind(OR = coef(logit.incomplete), confint(logit.incomplete)))

logit.incomplete.EAP <- glm(incomplete ~ TBISev_GCS + EAP, data = alldat, family = "binomial")
summary(logit.incomplete.EAP)
anova(logit.incomplete.EAP, test = "Chisq")
#Likelihood ratio test
anova(logit.incomplete, logit.incomplete.EAP, test = "Chisq")
NagelkerkeR2(logit.incomplete.EAP)
PseudoR2(logit.incomplete.EAP, which = "CoxSnell")
exp(logit.incomplete.EAP$coefficients)
exp(cbind(OR = coef(logit.incomplete.EAP), confint(logit.incomplete.EAP)))
with(logit.incomplete.EAP, null.deviance - deviance)
with(logit.incomplete.EAP, df.null - df.residual)
with(logit.incomplete.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

logit.incomplete <- glm(incomplete ~ EAP, data = alldat, family = "binomial")
summary(logit.incomplete)
anova(logit.incomplete, test = "Chisq")
NagelkerkeR2(logit.incomplete)
PseudoR2(logit.incomplete, which = "CoxSnell")
exp(cbind(OR = coef(logit.incomplete), confint(logit.incomplete)))

## END OF CODE TO FIT CENTER MODEL IN CENTER DATA

##Prediction of missingness START

#Create binary variables reflecting missingness for missingness analysis
dat$CompleteDataCT <- ifelse(complete.cases(dat[, c(
  "CT_SkullFx_R",
  "CT_Contusion_R",
  "CT_Shear_R",
  "CT_ExtraaxHematoma_R",
  "CT_EDH_R",
  "CT_SDH_R",
  "CT_SAH_R",
  "CT_IVH_R",
#  "CT_Edema_R",
#  "CT_DownwardHerniation_R",
#  "CT_UpwardCerebellarHerniation_R",
#  "CT_DuretHemorrhage_R",
  "CT_MidlineShift_R"
)]), 1, 0)
table(dat$CompleteDataCT)

dat$CompleteDataGCS <- ifelse(complete.cases(dat[, c(
  "GCSEDArrVerbal_R",
  "GCSEDArrMotor_R",
  "GCSEDArrEyes_R"
)]), 1, 0)
table(dat$CompleteDataGCS)

dat$CompleteDataGFAPUCH <- ifelse(complete.cases(dat[, c(
  "gfap_cat",
  "uchl_cat"
)]), 1, 0)
table(dat$CompleteDataGFAPUCH)

#Set categorical predictors as factors
categorical_vars <- c(
  "Sex", "ctpos", "InjCause_r", "levcare")

dat[categorical_vars] <- lapply(dat[categorical_vars], factor)

ModelMissingCT <- glm(
  CompleteDataCT ~ Age + Sex + EduYearsOfEducation +
    gcser + InjCause_r + levcare ,
  data = dat,
  family = binomial(link = "logit")
)
summary(ModelMissingCT)
anova(ModelMissingCT, test = "Chisq")
NagelkerkeR2(ModelMissingCT)
exp(cbind(OR = coef(ModelMissingCT), confint(ModelMissingCT)))

ModelMissingGCS <- glm(
  CompleteDataGCS ~ Age + Sex + EduYearsOfEducation +
    ctpos + InjCause_r + levcare,
  data = dat,
  family = binomial(link = "logit")
)
summary(ModelMissingGCS)
anova(ModelMissingGCS, test = "Chisq")
NagelkerkeR2(ModelMissingGCS)
exp(cbind(OR = coef(ModelMissingGCS), confint(ModelMissingGCS)))

ModelMissingGFAPUCH <- glm(
  CompleteDataGFAPUCH ~ Age + Sex + EduYearsOfEducation +
    gcser + ctpos + InjCause_r + levcare,
  data = dat,
  family = binomial(link = "logit")
)
summary(ModelMissingGFAPUCH)
anova(ModelMissingGFAPUCH, test = "Chisq")
NagelkerkeR2(ModelMissingGFAPUCH)
exp(cbind(OR = coef(ModelMissingGFAPUCH), confint(ModelMissingGFAPUCH)))
#^ Model gives warning - with only 123 missing this probably too many predictors. 

##Prediction of missingness END

