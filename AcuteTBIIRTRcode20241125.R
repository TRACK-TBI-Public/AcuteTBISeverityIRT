# TITLE: Nelson et al. - A model of traumatic brain injury severity
# from clinical, neuroimaging, and blood-based biomarkers - EFA 
# Script 1/2 - performs descriptive, IRT, and regression modeling
# Saves data file for use in Mplus factor modeling

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
setwd("I:/Neurosurgery Research -Clinical/Nelson, Lin/Manuscripts/R01_IPCOM/TRACK_AcuteIRT/Data/06_SharedVersion/20241125_Version")

dat <- read.csv("AcuteTBIIRTData20241125.csv")

# Table 2

describe(dat$Age)
quantile(dat$Age, probs = c(0.25, 0.75), na.rm = TRUE)

dat$Sex <- factor(dat$Sex, levels = c(1, 2), labels = c("Male", "Female"))
table(dat$Sex)
prop.table(table(dat$Sex))

dat$Race[is.na(dat$Race)] <- 88

dat$Race <- factor(dat$Race, levels = c(1, 2, 3, 4, 5, 6, 7, 88),
                   labels = c("Indian",
                              "Alaska Native/Inuit",
                              "Asian",
                              "Black",
                              "Native Hawaiian/Pacific Islander",
                              "White",
                              "Mixed race",
                              "Unknown"))


# Combine race categories

dat$Race_r <- factor(dat$Race, labels = c("Other/Unknown",
                              "Other/Unknown",
                              "Asian",
                              "Black",
                              "Other/Unknown",
                              "White",
                              "Other/Unknown",
                              "Other/Unknown"))

table(dat$Race_r)
prop.table(table(dat$Race_r))

dat$Ethnicity <- factor(dat$Ethnicity, levels = c(1, 2, 88),
                        labels = c("Hispanic", "Non-Hispanic", "Unknown"))
table(dat$Ethnicity)
prop.table(table(dat$Ethnicity))

dat$EduYearsOfEducation <- replace(dat$EduYearsOfEducation, dat$EduYearsOfEducation == 88, NA)

describe(dat$EduYearsOfEducation)
quantile(dat$EduYearsOfEducation, probs = c(0.25, 0.75), na.rm = T)

dat$HealthInsurance <- factor(dat$HealthInsurance, levels = c(1, 2, 3, 4, 5, 6, 7, 88, 99),
                              labels = c("self-pay (uninsured)",
                              "Insurance through a current or former employer (incl. thru family member)",
                              "Insurance purchased directly from an insurance company or on the health insurance",
                              "Medicare, for people 65 and older, or people with certain disabilities",
                              "Medicaid, Medical Assistance, “the State” or any kind of government-assistance plan for low income/disability",
                              "Medicaid Pending",
                              "TRICARE, VA or other military health care",
                              "Unknown",
                              "Any other type of health insurance or health coverage plan"))

dat$HealthInsurance_r <- factor(dat$HealthInsurance,
                              labels = c("Medicaid/Uninsured",
                                         "Other Insurance",
                                         "Other Insurance",
                                         "Other Insurance",
                                         "Medicaid/Uninsured",
                                         "Medicaid/Uninsured",
                                         "Other Insurance",
                                         "Unknown",
                                         "Other Insurance"))

table(dat$HealthInsurance_r)
prop.table(table(dat$HealthInsurance_r))

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

dat$GCSEDArrVerbal_R <- as.numeric(as.character(dat$GCSEDArrVerbal_R))
dat$GCSEDArrVerbal_R[dat$GCSEDArrVerbal_R == "T"] <- NA
dat$GCSEDArrVerbal_R <- 6 - dat$GCSEDArrVerbal_R

dat$GCSEDArrMotor_R <- as.numeric(as.character(dat$GCSEDArrMotor_R))
dat$GCSEDArrMotor_R[dat$GCSEDArrMotor_R == "P"] <- NA
dat$GCSEDArrMotor_R <- 7 - dat$GCSEDArrMotor_R

dat$GCSEDArrEyes_R <- as.numeric(as.character(dat$GCSEDArrEyes_R))
dat$GCSEDArrEyes_R[dat$GCSEDArrEyes_R == "S"] <- NA
dat$GCSEDArrEyes_R <- 5 - dat$GCSEDArrEyes_R

table(dat$GCSEDArrVerbal_R)


# Missing data

sum(is.na(dat$CT_SkullFx_R))/length(dat$CT_SkullFx_R)
sum(is.na(dat$CT_Contusion_R))/length(dat$CT_Contusion_R)
sum(is.na(dat$CT_Shear_R))/length(dat$CT_Shear_R)
sum(is.na(dat$CT_ExtraaxHematoma_R))/length(dat$CT_ExtraaxHematoma_R)
sum(is.na(dat$CT_EDH_R))/length(dat$CT_EDH_R)
sum(is.na(dat$CT_SDH_R))/length(dat$CT_SDH_R)
sum(is.na(dat$CT_SAH_R))/length(dat$CT_SAH_R)
sum(is.na(dat$CT_IVH_R))/length(dat$CT_IVH_R)
sum(is.na(dat$CT_Edema_R))/length(dat$CT_Edema_R)
sum(is.na(dat$CT_DownwardHerniation_R))/length(dat$CT_DownwardHerniation_R)
sum(is.na(dat$CT_UpwardCerebellarHerniation_R))/length(dat$CT_UpwardCerebellarHerniation_R)
sum(is.na(dat$CT_DuretHemorrhage_R))/length(dat$CT_DuretHemorrhage_R)
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
sum(is.na(dat$crp_cat))/length(dat$crp_cat)

# Variables for the IRT model

allvar <- c("Subject_ID", "CT_SkullFx_R", "CT_Contusion_R", "CT_Shear_R", "CT_ExtraaxHematoma_R", "CT_EDH_R", "CT_SDH_R",
            "CT_SAH_R", "CT_IVH_R", "CT_Edema_R", "CT_DownwardHerniation_R", "CT_UpwardCerebellarHerniation_R",
            "CT_DuretHemorrhage_R", "CT_MidlineShift_R", "LOCDuration_R", "LOCPTADuration_R",
            "GCSEDArrVerbal_R", "GCSEDArrMotor_R", "GCSEDArrEyes_R",
            "gfap_cat", "uchl_cat", "nse_cat", "s100_cat", "crp_cat", "NumberNonreactingPupils")


# Variables for the IRT model (minus biomarkers)

allvar_nobiomarkers <- c("Subject_ID", "CT_SkullFx_R", "CT_Contusion_R", "CT_Shear_R", "CT_ExtraaxHematoma_R", "CT_EDH_R", "CT_SDH_R",
                         "CT_SAH_R", "CT_IVH_R", "CT_Edema_R", "CT_DownwardHerniation_R", "CT_UpwardCerebellarHerniation_R",
                         "CT_DuretHemorrhage_R", "CT_MidlineShift_R", "LOCDuration_R", "LOCPTADuration_R",
                         "GCSEDArrVerbal_R", "GCSEDArrMotor_R", "GCSEDArrEyes_R", "NumberNonreactingPupils")

allvar <- dat[allvar] # Create dataset with just IRT variables
allvar_nobiomarkers <- dat[allvar_nobiomarkers] # Create dataset with just IRT variables (minus biomarkers)

allvardata <- allvar[rowSums(is.na(allvar)) != ncol(allvar), ] # Remove cases with all missing
allvar_nobiomarkers <- allvar_nobiomarkers[rowSums(is.na(allvar_nobiomarkers)) != ncol(allvar_nobiomarkers), ] # Remove cases with all missing

# Write Mplus data file

prepareMplusData(
  allvar,
  filename = "ShareDataMplus.dat",
  inpfile = TRUE)

# Unidimensional IRT model

allvar.uni <- mirt(allvardata[ , 2:25], itemtype = c(rep("graded", 24)), model = 1, SE = T)

coef(allvar.uni, IRTpars = T, simplify = T)


# No Biomarkers

no.biomarkers.uni <- mirt(allvar_nobiomarkers[ , 2:19], itemtype = c("graded", "graded", "graded", "graded", "graded",
                                                     "graded", "graded", "graded", "graded", "graded",
                                                     "graded", "graded", "graded", "graded", "graded",
                                                     "graded", "graded", "graded"), model = 1, SE = T)
coef(no.biomarkers.uni, IRTpars = T, simplify = T)

# Estimate IRT scores for full model

allvar.uni.scores <- fscores(allvar.uni, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
allvar.uni.scores <- data.frame(allvar.uni.scores)
allvardata$EAP <- allvar.uni.scores$F1
allvardata$SE <- allvar.uni.scores$SE_F1

# Estimate IRT scores using same item parameters but excluding biomarkers

allvar.uni.scores.nobiomarkers.sep.model <- fscores(no.biomarkers.uni, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
allvar.uni.scores.nobiomarkers.sep.model <- data.frame(allvar.uni.scores.nobiomarkers.sep.model)
allvardata$EAP.no.bio.sep <- allvar.uni.scores.nobiomarkers.sep.model$F1

# Estimate IRT scores for separate model without biomarkers

allvar.uni.scores.nobiomarkers.subset <- fscores(allvar.uni, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE,
                                                 item_weights = c(rep(1, 19), rep(0, 5)))
allvar.uni.scores.nobiomarkers.subset <- data.frame(allvar.uni.scores.nobiomarkers.subset)
allvardata$EAP.no.bio.subset <- allvar.uni.scores.nobiomarkers.subset$F1

# Correlations of IRT scores

cor(allvardata$EAP.no.bio.sep, allvardata$EAP.no.bio.subset,
    use = "complete.obs")

cor(allvardata$EAP, allvardata$EAP.no.bio.subset,
    use = "complete.obs")

# Merge dataset with IRT scores with original dataset (final dataset called "alldat")

scores <- allvardata %>% select(Subject_ID, EAP, SE, EAP.no.bio.sep, EAP.no.bio.subset)
alldat <- merge(dat, scores, by = "Subject_ID")

## k-fold cross validation START

all_data <- allvardata

# Define the number of folds for k-fold internal cross-validation
nreps <- 10  # Adjust as needed

# Define function to fit the IRT model, calculate scores, and replace non-test scores with NA
fit_irt_and_replace_with_na <- function(train_data, all_data, test_data, all_items) {
  model <- mirt(train_data[, all_items], itemtype = "graded", model = 1, SE = TRUE,
                verbose = FALSE, technical = list(message = FALSE))
  item_params <- coef(model, IRTpars = TRUE, simplify = TRUE)
  all_scores <- fscores(model, response.pattern = all_data[, all_items], method = "EAP", full.scores = TRUE, na.rm = FALSE)
  all_scores_df <- data.frame(Subject_ID = all_data$Subject_ID, F1 = all_scores[, "F1"]) %>%
    mutate(F1 = if_else(Subject_ID %in% test_data$Subject_ID, F1, NA_real_))
  list(item_params = item_params, test_scores = all_scores_df)
}

# Prepare data
all_items <- 2:25  # Adjust item column indices as needed
set.seed(123)  # For reproducibility
subject_ids <- allvardata$Subject_ID
all_scores_df <- data.frame(Subject_ID = subject_ids)
item_params_list <- vector("list", nreps)

# Run nreps replications
for (i in seq_len(nreps)) {
  cat("\nRunning replication:", i, "\n")
  
  # Split data into training and test sets
  train_indices <- sample(seq_len(nrow(allvardata)), size = 0.9 * nrow(allvardata))
  train_data <- allvardata[train_indices, ]
  test_data <- allvardata[-train_indices, ]
  
  result <- fit_irt_and_replace_with_na(train_data, allvardata, test_data, all_items)
  
  # Store item parameters and add scores to all_scores_df
  item_params_list[[i]] <- result$item_params
  all_scores_df <- all_scores_df %>%
    left_join(result$test_scores, by = "Subject_ID") %>%
    rename(!!paste0("Score_Sample_", i) := F1)
}

# Output results
cat("\nFinal inspection of all_scores_df after all replications:\n")
print(head(all_scores_df))

for (i in seq_len(nreps)) {
  cat("\nItem parameters for replication", i, ":\n")
  print(item_params_list[[i]])
  cat("\n---------------------------------\n")
}

# Merge `all_scores_df` with the original `all_data` by `Subject_ID`
final_data <- all_data %>%
  left_join(all_scores_df, by = "Subject_ID")

# Display the first few rows of the merged dataframe to confirm
cat("Merged dataframe (first few rows):\n")
print(head(final_data))

# Identify columns in `dat` that are also in `final_data`, excluding `Subject_ID`
common_columns <- intersect(names(final_data), names(dat))
common_columns <- setdiff(common_columns, "Subject_ID")  # Exclude `Subject_ID` from common columns

# Remove common columns from `dat` before merging
merged_data <- final_data %>%
  left_join(select(dat, -one_of(common_columns)), by = "Subject_ID")

# Display the first few rows of the merged dataframe
cat("Merged dataframe with unique columns only (first few rows):\n")
#print(head(merged_data))
cor.test(merged_data$EAP, merged_data$Score_Sample_1,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_2,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_3,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_4,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_5,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_6,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_7,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_8,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_9,
         use = "complete.obs",
         conf.level = 0.95)
cor.test(merged_data$EAP, merged_data$Score_Sample_10,
         use = "complete.obs",
         conf.level = 0.95)
# Define the variables for the Score_Sample_X
variables <- paste0("Score_Sample_", 1:10)

# Initialize a vector to store the correlations
correlations <- numeric()

# Loop over the variables and compute the correlations
for (var in variables) {
  test <- try(cor.test(merged_data$EAP, merged_data[[var]], 
                       use = "complete.obs", conf.level = 0.95), silent = TRUE)
  if (inherits(test, "htest")) {
    correlations <- c(correlations, test$estimate)  # Extract correlation
  } else {
    correlations <- c(correlations, NA)  # Handle errors with NA
  }
}

# Create a data frame for the results
results_df <- data.frame(Sample = variables, Correlation = round(correlations, 3))

# Display the results
print(results_df)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_1,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_2,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_3,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_4,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_5,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_6,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_7,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_8,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_9,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$Score_Sample_10,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
cor.test(merged_data$goseALL.T3, merged_data$EAP,
         use = "complete.obs", method = "spearman",
         conf.level = 0.95)
# Define the variables to compute correlations
variables <- c(paste0("Score_Sample_", 1:10), "EAP")
# Initialize an empty vector to store correlations
correlations <- numeric()
# Loop over the variables to compute correlations
for (var in variables) {
  test <- try(cor.test(merged_data$goseALL.T3, merged_data[[var]], 
                       use = "complete.obs", method = "spearman"), silent = TRUE)
  if (inherits(test, "htest")) {
    correlations <- c(correlations, test$estimate)  # Extract correlation
  } else {
    correlations <- c(correlations, NA)  # Handle errors by adding NA
  }
}
# Create a data frame with the results
results_df <- data.frame(Sample = variables, Correlation = round(correlations, 3))
# Display the results
print(results_df)

## k-fold cross validation END


## Plots

ggplot(alldat, aes(x= EAP.no.bio.subset, y = EAP)) +
  geom_point(size = 2, shape = 23) +
  labs(title = "IRT Scores with and without Biomarkers", x = "IRT Scores without Biomarkers", y = "IRT Scores with Biomarkers") +
  scale_x_continuous(limits = c(-2, 2.5)) +
  scale_y_continuous(limits = c(-2, 2.5)) + 
  theme_minimal() +
  annotate("text", x = 2, y = -1.5, label = "r = .97", color = "black", size = 4)

ggplot(alldat, aes(x = EAP)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "white") +
  labs(title = "IRT Scores", x = "TBI Severity", y = "Frequency") +
 # scale_x_continuous(limits = c(-4, 4)) +
  theme_minimal()

# Extract information for each domain and entire test

HeadCT <- 0
Clinical.All <- 0
Clinical.GCS <- 0
Clinical.PTA.LOC.Pupil <- 0
Biomarkers <- 0

Theta <- seq(-3, 3, 0.1)

for(i in 1:13){
  HeadCT <- HeadCT + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in c(14:18,24)){
  Clinical.All <- Clinical.All + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in c(14:15,24)){
  Clinical.PTA.LOC.Pupil <- Clinical.PTA.LOC.Pupil + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in 16:18){
  Clinical.GCS <- Clinical.GCS + iteminfo(extract.item(allvar.uni, i), Theta)
}

for(i in 19:23){
  Biomarkers <- Biomarkers + iteminfo(extract.item(allvar.uni, i), Theta)
}

Test <- testinfo(allvar.uni, Theta)

# Information plot (no test information)

InfoData <- cbind(Theta, HeadCT, Clinical.PTA.LOC.Pupil, Clinical.GCS, Biomarkers)
InfoData <- data.frame(InfoData)
names(InfoData) <- c("Theta", "Head CT", "Clinical-LOC/PTA/Pupil", "Clinical-GCS", "Biomarkers")

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
names(InfoData) <- c("Theta", "Head CT", "Clinical-LOC/PTA/Pupil", "Clinical-GCS", "Biomarkers", "Test")

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
names(InfoDataIncremental) <- c("Theta", "GCS", "GCS + Head CT", "GCS + Head CT + LOC/PTA/Pupil", "GCS + Head CT + LOC/PTA/Pupil + Biomarkers")
InfoDataIncremental <- melt(InfoDataIncremental, id.vars=c("Theta"))
names(InfoDataIncremental) <- c("Theta", "Domain", "Information")

SEDataIncremental <- cbind(Theta, SE.TestGCS, SE.TestGCSCT, SE.TestGCSCTLOCPupil, SE.TestGCSCTLOCBioPupil)
SEDataIncremental <- data.frame(SEDataIncremental)
names(SEDataIncremental) <- c("Theta", "GCS", "GCS + Head CT", "GCS + Head CT + LOC/PTA/Pupil", "GCS + Head CT + LOC/PTA/Pupil + Biomarkers")
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
CT_Edema <- iteminfo(extract.item(allvar.uni, 9), Theta)
CT_DownwardHerniation <- iteminfo(extract.item(allvar.uni, 10), Theta)
CT_UpwardCereHerniation <- iteminfo(extract.item(allvar.uni, 11), Theta)
CT_DuretHemorrhage <- iteminfo(extract.item(allvar.uni, 12), Theta)
CT_MidlineShift <- iteminfo(extract.item(allvar.uni, 13), Theta)
LOCDuration <- iteminfo(extract.item(allvar.uni, 14), Theta)
LOCPTADuration <- iteminfo(extract.item(allvar.uni, 15), Theta)
GCSEDArrVerbal_R <- iteminfo(extract.item(allvar.uni, 16), Theta)
GCSEDArrMotor_R <- iteminfo(extract.item(allvar.uni, 17), Theta)
GCSEDArrEyes_R <- iteminfo(extract.item(allvar.uni, 18), Theta)
gfap_cat <- iteminfo(extract.item(allvar.uni, 19), Theta)
uchl_cat <- iteminfo(extract.item(allvar.uni, 20), Theta)
nse_cat <- iteminfo(extract.item(allvar.uni, 21), Theta)
s100_cat <- iteminfo(extract.item(allvar.uni, 22), Theta)
crp_cat <- iteminfo(extract.item(allvar.uni, 23), Theta)
NumberNonreactingPupils <- iteminfo(extract.item(allvar.uni, 24), Theta)
Test <- testinfo(allvar.uni, Theta)

InfoData <- cbind(Theta, CT_SkullFx, CT_Contusion, CT_Shear, CT_ExtraaxHematoma,
                  CT_EDH, CT_SDH, CT_SAH, CT_IVH,
                  CT_Edema, CT_DownwardHerniation, CT_UpwardCereHerniation, CT_DuretHemorrhage,
                  CT_MidlineShift, LOCDuration, LOCPTADuration, GCSEDArrVerbal_R, GCSEDArrMotor_R,
                  GCSEDArrEyes_R, NumberNonreactingPupils, gfap_cat, uchl_cat, nse_cat, s100_cat, crp_cat)
InfoData <- data.frame(InfoData)
names(InfoData) <- c("Theta", "Skull Fracture", "Contusion", "Shear", "Extraaxial Hematoma", "EDH", "SDH",
                     "SAH", "IVH", "Edema", "Downward Herniation", "Upward Cerebellar Herniation",
                     "Duret Hemorrhage", "Midline Shift", "LOC Duration", "PTA Duration",
                     "GCS-Verbal", "GCS-Motor", "GCS-Eyes", "Pupil Reactivity",
                     "GFAP", "UCH-L1", "NSE", "S100B", "hsCRP")

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
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash")) +
  scale_color_manual(values = c("#006666",
                                "#009999",
                                "#00CCCC",
                                "#00FFFF",
                                "#99FFFF",
                                "#003366",
                                "#004C99",
                                "#0066CC",
                                "#0080FF",
                                "#66B2FF",
                                "#0000FF",
                                "#6666FF",
                                "#CCCCFF",
                                "#666600",
                                "#999900",
                                "#CCCC00",
                                "#DAA520",
                                "#B8860B",
                                "#FFD700",
                                "#660000",
                                "#990000",
                                "#CC0000",
                                "#FF0000",
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
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "twodash",
                                 "solid",
                                 "dashed",
                                 "dotted",
                                 "dotdash",
                                 "longdash",
                                 "dotdash",
                                 "solid")) +
  scale_color_manual(values = c("#006666",
                                "#009999",
                                "#00CCCC",
                                "#00FFFF",
                                "#99FFFF",
                                "#003366",
                                "#004C99",
                                "#0066CC",
                                "#0080FF",
                                "#66B2FF",
                                "#0000FF",
                                "#6666FF",
                                "#CCCCFF",
                                "#666600",
                                "#999900",
                                "#CCCC00",
                                "#FFFF33",
                                "#FFFF99",
                                "#660000",
                                "#990000",
                                "#CC0000",
                                "#FF0000",
                                "#FF9999",
                                "darkgrey",
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

alldat$TBISev_Civilian4grp <- factor(alldat$TBISev_Civilian4grp,
                                     levels = c(1, 2, 3, 4),
                                     labels = c("u-MTBI", "c-MTBI", "Mod.", "Severe"))

alldat$TBISev_GCS <- factor(alldat$TBISev_GCS,
                                     levels = c(1, 2, 3),
                                     labels = c("13-15", "9-12", "3-8"))

alldat$TBISev_VA3grp <- factor(alldat$TBISev_VA3grp,
                            levels = c(1, 2, 3),
                            labels = c("Mild", "Mod.", "Severe"))


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


alldat <- alldat %>% filter( !is.na(TBISev_Civilian4grp))
alldat <- alldat %>% filter( !is.na(TBISev_GCS))
alldat <- alldat %>% filter( !is.na(TBISev_VA3grp))


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

g2 <- ggplot(alldat, aes(x = EAP, y = TBISev_VA3grp)) +
  geom_density_ridges(
    jittered_points = TRUE, position = "raincloud",
    alpha = 0.7, scale = 0.9) +
  labs(title = "U.S. Veteran's Affairs TBI Classification", x = "TBI Severity", y = "") +
  theme_bw() +
  xlim(-4, 4) +
  scale_x_continuous(name="TBI Severity", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

g2

g3 <- ggplot(alldat, aes(x = EAP, y = TBISev_Civilian4grp)) +
  geom_density_ridges(
    jittered_points = TRUE, position = "raincloud",
    alpha = 0.7, scale = 0.9) +
  labs(title = "Civilian 4-Group TBI Classification", x = "TBI Severity", y = "") +
  theme_bw() +
  scale_x_continuous(name="TBI Severity", limits=c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(plot.title = element_text(hjust = 0.5))

g3


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


# Sensitivity analysis using IRT scores without biomarkers

logit.death.EAP.nobio <- glm(dead ~ EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.death.EAP.nobio)
anova(logit.death.EAP.nobio, test = "Chisq")
NagelkerkeR2(logit.death.EAP.nobio)
PseudoR2(logit.death.EAP.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.death.EAP.nobio), confint(logit.death.EAP.nobio)))

logit.unfavorable.EAP.nobio <- glm(unfavorable ~ EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.unfavorable.EAP.nobio)
anova(logit.unfavorable.EAP.nobio, test = "Chisq")
NagelkerkeR2(logit.unfavorable.EAP.nobio)
PseudoR2(logit.unfavorable.EAP.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable.EAP.nobio), confint(logit.unfavorable.EAP.nobio)))

logit.incomplete.EAP <- glm(incomplete ~ EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.incomplete.EAP)
anova(logit.incomplete.EAP, test = "Chisq")
NagelkerkeR2(logit.incomplete.EAP)
PseudoR2(logit.incomplete.EAP, which = "CoxSnell")
exp(cbind(OR = coef(logit.incomplete.EAP), confint(logit.incomplete.EAP)))


logit.death.nobio <- glm(dead ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.death.nobio)
anova(logit.death.nobio, test = "Chisq")
NagelkerkeR2(logit.death.nobio)
PseudoR2(logit.death.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.death.nobio), confint(logit.death.nobio)))

logit.death.EAP.nobio <- glm(dead ~ TBISev_GCS + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.death.EAP.nobio)
anova(logit.death.EAP.nobio, test = "Chisq")
anova(logit.death.nobio, logit.death.EAP.nobio, test = "Chisq")
NagelkerkeR2(logit.death.EAP.nobio)
PseudoR2(logit.death.EAP.nobio, which = "CoxSnell")
exp(logit.death.EAP.nobio$coefficients)
exp(cbind(OR = coef(logit.death.EAP.nobio), confint(logit.death.EAP.nobio)))
with(logit.death.EAP.nobio, null.deviance - deviance)
with(logit.death.EAP.nobio, df.null - df.residual)
with(logit.death.EAP.nobio, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death.nobio <- glm(dead ~ IMPACTCoreLPUnfav, data = alldat, family = "binomial")
summary(logit.death.nobio)
anova(logit.death.nobio, test = "Chisq")
NagelkerkeR2(logit.death.nobio)
PseudoR2(logit.death.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.death.nobio), confint(logit.death.nobio)))

logit.death.EAP.nobio <- glm(dead ~ IMPACTCoreLPUnfav + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.death.EAP.nobio)
anova(logit.death.EAP.nobio, test = "Chisq")
anova(logit.death.nobio, logit.death.EAP.nobio, test = "Chisq")
NagelkerkeR2(logit.death.EAP.nobio)
PseudoR2(logit.death.EAP.nobio, which = "CoxSnell")
exp(logit.death.EAP.nobio$coefficients)
exp(cbind(OR = coef(logit.death.EAP.nobio), confint(logit.death.EAP.nobio)))
with(logit.death.EAP.nobio, null.deviance - deviance)
with(logit.death.EAP.nobio, df.null - df.residual)
with(logit.death.EAP.nobio, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.death.nobio <- glm(dead ~ IMPACTExtendedLPUnfav, data = alldat, family = "binomial")
summary(logit.death.nobio)
anova(logit.death.nobio, test = "Chisq")
NagelkerkeR2(logit.death.nobio)
PseudoR2(logit.death.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.death.nobio), confint(logit.death.nobio)))

logit.death.EAP.nobio <- glm(dead ~ IMPACTExtendedLPUnfav + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.death.EAP.nobio)
anova(logit.death.EAP.nobio, test = "Chisq")
anova(logit.death.nobio, logit.death.EAP.nobio, test = "Chisq")
NagelkerkeR2(logit.death.EAP.nobio)
PseudoR2(logit.death.EAP.nobio, which = "CoxSnell")
exp(logit.death.EAP.nobio$coefficients)
exp(cbind(OR = coef(logit.death.EAP.nobio), confint(logit.death.EAP.nobio)))
with(logit.death.EAP.nobio, null.deviance - deviance)
with(logit.death.EAP.nobio, df.null - df.residual)
with(logit.death.EAP.nobio, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable.nobio <- glm(unfavorable ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio)
anova(logit.unfavorable.nobio, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio)
PseudoR2(logit.unfavorable.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable.nobio), confint(logit.unfavorable.nobio)))

logit.unfavorable.nobio.EAP <- glm(unfavorable ~ TBISev_GCS + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio.EAP)
anova(logit.unfavorable.nobio.EAP, test = "Chisq")
anova(logit.unfavorable.nobio, logit.unfavorable.nobio.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio.EAP)
PseudoR2(logit.unfavorable.nobio.EAP, which = "CoxSnell")
exp(logit.unfavorable.nobio.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.nobio.EAP), confint(logit.unfavorable.nobio.EAP)))
with(logit.unfavorable.nobio.EAP, null.deviance - deviance)
with(logit.unfavorable.nobio.EAP, df.null - df.residual)
with(logit.unfavorable.nobio.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable.nobio <- glm(unfavorable ~ IMPACTCoreLPUnfav, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio)
anova(logit.unfavorable.nobio, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio)
PseudoR2(logit.unfavorable.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable.nobio), confint(logit.unfavorable.nobio)))

logit.unfavorable.nobio.EAP <- glm(unfavorable ~ IMPACTCoreLPUnfav + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio.EAP)
anova(logit.unfavorable.nobio.EAP, test = "Chisq")
anova(logit.unfavorable.nobio, logit.unfavorable.nobio.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio.EAP)
PseudoR2(logit.unfavorable.nobio.EAP, which = "CoxSnell")
exp(logit.unfavorable.nobio.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.nobio.EAP), confint(logit.unfavorable.nobio.EAP)))
with(logit.unfavorable.nobio.EAP, null.deviance - deviance)
with(logit.unfavorable.nobio.EAP, df.null - df.residual)
with(logit.unfavorable.nobio.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.unfavorable.nobio <- glm(unfavorable ~ IMPACTExtendedLPUnfav, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio)
anova(logit.unfavorable.nobio, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio)
PseudoR2(logit.unfavorable.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.unfavorable.nobio), confint(logit.unfavorable.nobio)))

logit.unfavorable.nobio.EAP <- glm(unfavorable ~ IMPACTExtendedLPUnfav + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.unfavorable.nobio.EAP)
anova(logit.unfavorable.nobio.EAP, test = "Chisq")
anova(logit.unfavorable.nobio, logit.unfavorable.nobio.EAP, test = "Chisq")
NagelkerkeR2(logit.unfavorable.nobio.EAP)
PseudoR2(logit.unfavorable.nobio.EAP, which = "CoxSnell")
exp(logit.unfavorable.nobio.EAP$coefficients)
exp(cbind(OR = coef(logit.unfavorable.nobio.EAP), confint(logit.unfavorable.nobio.EAP)))
with(logit.unfavorable.nobio.EAP, null.deviance - deviance)
with(logit.unfavorable.nobio.EAP, df.null - df.residual)
with(logit.unfavorable.nobio.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))


logit.incomplete.nobio <- glm(incomplete ~ TBISev_GCS, data = alldat, family = "binomial")
summary(logit.incomplete.nobio)
anova(logit.incomplete.nobio, test = "Chisq")
NagelkerkeR2(logit.incomplete.nobio)
PseudoR2(logit.incomplete.nobio, which = "CoxSnell")
exp(cbind(OR = coef(logit.incomplete.nobio), confint(logit.incomplete.nobio)))

logit.incomplete.nobio.EAP <- glm(incomplete ~ TBISev_GCS + EAP.no.bio.subset, data = alldat, family = "binomial")
summary(logit.incomplete.nobio.EAP)
anova(logit.incomplete.nobio.EAP, test = "Chisq")
anova(logit.incomplete.nobio, logit.incomplete.nobio.EAP, test = "Chisq")
NagelkerkeR2(logit.incomplete.nobio.EAP)
PseudoR2(logit.incomplete.nobio.EAP, which = "CoxSnell")
exp(logit.incomplete.nobio.EAP$coefficients)
exp(cbind(OR = coef(logit.incomplete.nobio.EAP), confint(logit.incomplete.nobio.EAP)))
with(logit.incomplete.nobio.EAP, null.deviance - deviance)
with(logit.incomplete.nobio.EAP, df.null - df.residual)
with(logit.incomplete.nobio.EAP, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = F))

## Sensitivity analyses to treat GOSE ordinally 

# Spearman Correlations between IRT scores and GOSE, stratified by GCS-based severity levels

GCS1315 <- filter(alldat, TBISev_GCS == "13-15") 
GCS912 <- filter(alldat, TBISev_GCS == "9-12")
GCS38 <- filter(alldat, TBISev_GCS == "3-8")

cor.test(GCS1315$EAP, GCS1315$goseALL.T1, method = "spearman", use = "complete.obs")
cor.test(GCS1315$EAP, GCS1315$goseTBI.T1, method = "spearman", use = "complete.obs")

cor.test(GCS1315$EAP, GCS1315$goseALL.T3, method = "spearman", use = "complete.obs")
cor.test(GCS1315$EAP, GCS1315$goseTBI.T3, method = "spearman", use = "complete.obs")

cor.test(GCS912$EAP, GCS912$goseALL.T1, method = "spearman", use = "complete.obs")
cor.test(GCS912$EAP, GCS912$goseTBI.T1, method = "spearman", use = "complete.obs")

cor.test(GCS912$EAP, GCS912$goseALL.T3, method = "spearman", use = "complete.obs")
cor.test(GCS912$EAP, GCS912$goseTBI.T3, method = "spearman", use = "complete.obs")

cor.test(GCS38$EAP, GCS38$goseALL.T1, method = "spearman", use = "complete.obs")
cor.test(GCS38$EAP, GCS38$goseTBI.T1, method = "spearman", use = "complete.obs")

cor.test(GCS38$EAP, GCS38$goseALL.T3, method = "spearman", use = "complete.obs")
cor.test(GCS38$EAP, GCS38$goseTBI.T3, method = "spearman", use = "complete.obs")

#Rank Regressions 

# Convert 6M GOSE to ranks
alldat$GOSE_rank <- rank(as.numeric(alldat$goseALL.T3))

# Hierarchical rank regression with GCS + TBI Severity IRT scores and F change
rank1 <- lm(GOSE_rank ~ TBISev_GCS, data = alldat)
summary(rank1)

rank2 <- lm(GOSE_rank ~ TBISev_GCS + EAP, data = alldat)
summary(rank2)

anova(rank1,rank2)

# Hierarchical rank regression with IMPACT Core + TBI Severity IRT scores and F change

rank3 <- lm(GOSE_rank ~ IMPACTCoreLPUnfav, data = alldat)
summary(rank3)

rank4 <- lm(GOSE_rank ~ IMPACTCoreLPUnfav + EAP, data = alldat)
summary(rank4)

anova(rank3,rank4)

# Hierarchical rank regression with IMPACT Extended + TBI Severity IRT scores and F change

rank5 <- lm(GOSE_rank ~ IMPACTExtendedLPUnfav, data = alldat)
summary(rank5)

rank6 <- lm(GOSE_rank ~ IMPACTExtendedLPUnfav + EAP, data = alldat)
summary(rank6)

anova(rank5,rank6)

