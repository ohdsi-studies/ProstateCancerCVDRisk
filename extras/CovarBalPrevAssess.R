# Install and load the sqldf package
install.packages('sqldf')
library(sqldf)

shinyDataFolder <- "D:/StudyResults/EPI720_2/shinyData"
studyFolder <- "D:/StudyResults/EPI720_2"

databasefiles <- list.files(shinyDataFolder, pattern = "database_")
databases <- sub(".*database_*(.*?) *.rds", "\\1", databasefiles)

getComparisonSummaries <- function(database) {
  comparisonSummary <- readRDS(file.path(shinyDataFolder, paste0("comparison_summary_", database, ".rds")))
  comparisonSummary <- comparisonSummary[c("target_id", "comparator_id", "database_id")]
  return(comparisonSummary)
}
comparisonSummaries <- lapply(databases, getComparisonSummaries)
comparisonSummaries <- do.call(rbind, comparisonSummaries)
comparisonSummaries <- comparisonSummaries[!duplicated(comparisonSummaries), ]


getCovariateSummaries <- function(database) {
  covariateSummary <- readRDS(file.path(shinyDataFolder, paste0("covariate_", database, ".rds")))
  return(covariateSummary)
}
covariateSummaries <- lapply(databases, getCovariateSummaries)
covariateSummaries <- do.call(rbind, covariateSummaries)
covariateSummaries <- covariateSummaries[c("covariate_id", "covariate_name")]
covariateSummaries <- covariateSummaries[!duplicated(covariateSummaries), ]

exposures <- readRDS(file.path(shinyDataFolder, paste0("exposure_of_interest_Optum.rds")))[, 1:2] # hardcoded, fix
analyses <- readRDS(file.path(shinyDataFolder, paste0("cohort_method_analysis_Optum.rds")))[, 1:2] # hardcoded, fix  (MITCH ADDED)

getCovarBalSummaries <- function(database, comparisonSummaries) { # database <- databases[1]
  comparisonSummary <- comparisonSummaries[comparisonSummaries$database_id == database, ]
  fullBal <- data.frame()
  for (i in 1:nrow(comparisonSummary)) { # i = 1
    targetId <- comparisonSummary$target_id[i]
    comparatorId <- comparisonSummary$comparator_id[i]
    database <- comparisonSummary$database_id[i]
    bal <- readRDS(file.path(shinyDataFolder, paste0("covariate_balance_t", targetId, "_c", comparatorId, "_", database, ".rds")))
    bal$covariate_analysis_id <- as.numeric(substr(bal$covariate_id, nchar(bal$covariate_id)-2, nchar(bal$covariate_id)))
    bal <- bal[!bal$covariate_analysis_id %in% c(901, 902, 903, 904, 920), ] # revmove charlson, dcsi, chads2vasc, chads
    fullBal <- rbind(fullBal, bal)
  }
  return(fullBal)
}

# MC NOTE: Here is where I apply individual blinding rules 

covarBalSummaries <- lapply(databases, getCovarBalSummaries, comparisonSummaries) # databases[c(2,8)]
covarBalSummaries <- do.call(rbind, covarBalSummaries)

covarBalSummaries$blindOriginal <- ifelse(abs(covarBalSummaries$std_diff_after) > 0.1 & 
                                    abs(pmin(covarBalSummaries$target_mean_after, 1-covarBalSummaries$target_mean_after) - pmin(covarBalSummaries$comparator_mean_after, 1-covarBalSummaries$comparator_mean_after)) > 0.05, 1, 0)
covarBalSummaries$blindStdzdiff <- ifelse(abs(covarBalSummaries$std_diff_after) > 0.1,1,0)
covarBalSummaries$blindAbsdiff <- ifelse(abs(pmin(covarBalSummaries$target_mean_after, 1-covarBalSummaries$target_mean_after) - pmin(covarBalSummaries$comparator_mean_after, 1-covarBalSummaries$comparator_mean_after)) > 0.05, 1, 0)
covarBalSummaries$blind <- ifelse(covarBalSummaries$blindStdzdiff==1 & covarBalSummaries$blindAbsdiff==1,1,0)

covarBalSummaries <- merge(covariateSummaries, covarBalSummaries)
covarBalSummaries <- merge(exposures, covarBalSummaries, by.x = "exposure_id", by.y = "comparator_id")
names(covarBalSummaries)[names(covarBalSummaries) == "exposure_name"] <- "comparator_name"
names(covarBalSummaries)[names(covarBalSummaries) == "exposure_id"] <- "comparator_id"
covarBalSummaries <- merge(exposures, covarBalSummaries, by.x = "exposure_id", by.y = "target_id")
names(covarBalSummaries)[names(covarBalSummaries) == "exposure_name"] <- "target_name"
names(covarBalSummaries)[names(covarBalSummaries) == "exposure_id"] <- "target_id"
covarBalSummaries <- merge(analyses, covarBalSummaries, by.x = "analysis_id", by.y = "analysis_id")
names(covarBalSummaries)[names(covarBalSummaries) == "description"] <- "analysis_name"

DataFolder <- "D:/StudyResults/EPI720_2/"

#covBalSummaries$comparator <- ifelse(covBalSummaries$comparator_id=13316, "[EPI_720 - C01] New users of enzalutamide with CRPC (controlled CVD in prior 6-mo)",
#                                ifelse(covBalSummaries$comparator_id=13319, "[EPI_720 - C01] New users of enzalutamide with CRPC (uncontrolled CVD in prior 6-mo)",
#                                  ifelse(covBalSummaries$comparator_id=13247,	"[EPI_720 - C01] New users of enzalutamide with CRPC (no CVD in prior 6-mo)",
#                                      ifelse(covBalSummaries$comparator_id=13246, "[EPI_720 - C01] New users of enzalutamide with CRPC (any CVD in prior 6-mo)",
#                                        ifelse(covBalSummaries$comparator_id=13245, "[EPI_720 - C01] New users of enzalutamide with CRPC", "ERROR")))))

#covarBalSummaries <- merge(analyses, covarBalSummaries, by.x = "")

# drop useless columns
#write.csv(covarBalSummaries, file.path(studyFolder, "covarBalSummaries.csv"), row.names = FALSE)

#toBlind <- covarBalSummaries[covarBalSummaries$blind == 1, ]
#blindingDecisions <- unique(covarBalSummaries[, c("database_id", "target_id", "target_name", "comparator_id", "analysis_id", "analysis_name", "comparator_name", "blind", "blindStdzdiff", "blindAbsdiff")])


covarBalSummaries$printImbalance <- ifelse(covarBalSummaries$database_id == "Optum" & 
                                           covarBalSummaries$target_id == 13186 &
                                           covarBalSummaries$comparator_id == 13245 &
                                           covarBalSummaries$analysis_id <= 8 &
                                           abs(covarBalSummaries$std_diff_after) > 0.1, 1, 0)
printImbalanceSummary <- covarBalSummaries[covarBalSummaries$printImbalance == 1, ]

keepVars <- c("database_id", "covariate_name", "analysis_name", "target_mean_before", "comparator_mean_before", "std_diff_before", "target_mean_after", "comparator_mean_after", "std_diff_after")
printImbalanceDetail <- printImbalanceSummary[keepVars]
write.csv(printImbalanceDetail, file.path(studyFolder, "printImbalanceDetail.csv"), row.names = FALSE)

printImbalanceSummary <- unique(printImbalanceSummary[, c("database_id", "covariate_name")])
write.csv(printImbalanceSummary, file.path(studyFolder, "printImbalanceSummary.csv"), row.names = FALSE)

# MC NOTE: There are 5 strata, 6 outcomes (not counting Alzheimers / dementia / cognitive decline), 2 databases, and 10 analyses
# MC Note: outcome_id = 13327 refers to the analysis with the Alzheimer's / cogn decline outcome and is excluded from the covariate inspection below.

blindingDecisions <- sqldf('select database_id, target_id, target_name, comparator_id, comparator_name, outcome_id, analysis_id, analysis_name,
                            max(blind) as blind, max(blindStdzdiff) as blindStdzdiff, max(blindAbsdiff) as blindAbsdiff
                            from (select * from covarBalSummaries where outcome_id != 13327)
                            group by database_id, target_id, target_name, comparator_id, comparator_name, outcome_id, analysis_id, analysis_name')

toBlind <- covarBalSummaries[covarBalSummaries$blind == 1 & covarBalSummaries$outcome_id != 13327, ]
toBlind <- unique(toBlind[, c("database_id", "target_id", "target_name", "comparator_id", "outcome_id", "analysis_id", "comparator_name", "analysis_name", "blind")])

# The code below generates a list of covariates that have StdzDiff > .1 but an absolute prevalence difference of < 0.05
# We plan to manually inspect the list below for pluasible confounders.

# First, we merge so that we do not consider any analyses that failed outright on the basis of having a covariate that failed both blinding criteria

toUnblind <- blindingDecisions[blindingDecisions$blind != 1, ]

# MC Note: The code below restricts to analyses that we plan to unblind and then pulls the covariate-specific information for manual inspection
toUnblindCheck <- sqldf('select a.database_id, a.target_id, a.target_name, a.comparator_id, a.comparator_name, a.outcome_id, a.analysis_id, a.analysis_name,
                                a.blind, a.blindStdzdiff, a.blindAbsdiff, 
                                b.covariate_name, b.blind as cov_blind, b.blindStdzdiff as cov_blindStdzdiff, b.blindAbsdiff as cov_blindAbsdiff
                            from (select * from toUnblind where outcome_id != 13327) as a left join (select * from covarBalSummaries where outcome_id != 13327) as b
                            on (a.database_id=b.database_id and
                                a.target_id=b.target_id and
                                a.target_name=b.target_name and
                                a.comparator_id=b.comparator_id and
                                a.comparator_name=b.comparator_name and
                                a.outcome_id=b.outcome_id and
                                a.analysis_id=b.analysis_id and
                                a.analysis_name=b.analysis_name)')

sqldf('select distinct target_id from blindingDecisions where blind != 1')

unique(toBlind$target_id)

unique(blindingDecisions$target_id)

unique(toUnblind$target_id)

# MC Note: These are just confirmations that my blinding variables have the proper ranges based on the above restrictions (everything looks good)
# min(toUnblindCheck$blind) 
# max(toUnblindCheck$blind)
# min(toUnblindCheck$blindStdzdiff)
# max(toUnblindCheck$blindStdzdiff)
# min(toUnblindCheck$blindAbsdiff)
# max(toUnblindCheck$blindAbsdiff)
# min(covarBalSummaries$blind) 
# max(covarBalSummaries$blind)
# min(covarBalSummaries$blindStdzdiff)
# max(covarBalSummaries$blindStdzdiff)
# min(covarBalSummaries$blindAbsdiff)
# max(covarBalSummaries$blindAbsdiff)

# In this step I subset to the variables that failed the standardized difference criterion but passed the absolute prevalence criterion.
# Among those covariates, I eliminate the duplicates and trim everything such that the same covariate isn't appear multiple times in the 
#  list for multiple different look-back periods
unBlindCheck <- toUnblindCheck[toUnblindCheck$cov_blindStdzdiff==1 & toUnblindCheck$cov_blindAbsdiff == 0, ]
unBlindCheck$covariateNameclipped <-gsub(".*:","",unBlindCheck$covariate_name)
unBlindCheck <- unBlindCheck[order(unBlindCheck$covariateNameclipped),]
unBlind <- unique(unBlindCheck[, c("covariateNameclipped")])
write.csv(unBlind, file.path(studyFolder, "unBlind.csv"), row.names = FALSE)

#unBlindTEST <- covarBalSummaries
#unBlindTEST$covariateNameclipped <- gsub(".*:","",covarBalSummaries$covariate_name)
#unBlindTEST <- unBlindTEST[order(unBlindTEST$covariateNameclipped),]
#unBlindTEST <- unique(unBlindTEST[, c("covariateNameclipped")])

saveRDS(toBlind, file.path(shinyDataFolder, "to_blind.rds"))
write.csv(toBlind, file.path(studyFolder, "toBlind.csv"), row.names = FALSE)

write.csv(blindingDecisions, file.path(studyFolder, "blindingDecisions.csv"), row.names = FALSE)


# Importing the results of the manual covariate review by Jamie, Mitch, Dina, Gowtham and Patrick
install.packages('XLConnect')
library(XLConnect)
wb_manual_cov_review = loadWorkbook("C:/Users/admin_mconove1/Documents/EPI720_ManualCovReview_PotentialConfounders (GROUP RESPONSES).xlsx")
manual_cov_review = readWorksheet(wb_manual_cov_review, sheet = 'manual_cov_review', header = TRUE)
manual_cov_review <- manual_cov_review[c("blind_GE2","blind_GE1","blind_total","blind_mitch","blind_jamie","blind_gowtham","blind_dina","blind_patrick","Covariate")]
  # I manually edited the index month and index year covariate names in the excel doc - just switching them back.
manual_cov_review$Covariate <- ifelse(startsWith(manual_cov_review$Covariate, "index"),
                                      gsub(".*:","",manual_cov_review$Covariate),manual_cov_review$Covariate)

# MC Note: Merging the results of the manual covariate review with the dataset containing blinding decisions at the covariate-/analysis-/target-/outcome-specific level
unBlindCheckManual <- sqldf('select a.*, 
                                    b.blind_GE2, b.blind_GE1, b.blind_total, 
                                    b.blind_mitch, b.blind_jamie, b.blind_gowtham, b.blind_dina, b.blind_patrick 
                            from unBlindCheck as a LEFT JOIN manual_cov_review as b
                            on a.covariateNameclipped=b.Covariate')

## MC Note: Summarizing the blinding decisions at the analysis-/target-/outcome-specific level
unBlindCheckManualCollapse <- sqldf('select database_id, target_id, target_name, comparator_id, outcome_id, analysis_id, comparator_name, analysis_name,
                                     max(blind_GE1) AS blind_GE1, max(blind_GE2) AS blind_GE2
                                     from unBlindCheckManual
                                     group by database_id, target_id, target_name, comparator_id, outcome_id, analysis_id, comparator_name, analysis_name')
# MC Note: Merging back with the master set of listed analyses.
# MC Note: There are 5 strata, 6 outcomes (not counting Alzheimers / dementia / cognitive decline), 2 databases, and 10 analyses
# MC Note: outcome_id = 13327 refers to the analysis with the Alzheimer's / cogn decline outcome and is excluded from the covariate inspection below.

finalUnblindSummary <- sqldf('select a.*, b.blind_GE2, b.blind_GE1,
                              case when a.outcome_id = 13320 then "Composite: ischemic stroke, hemorrhagic stroke, heart failure, acute myocardial infarction or sudden cardiac death)"
                                   when a.outcome_id = 13321 then "Ischemic stroke"
                                   when a.outcome_id = 13323 then "Hemorrhagic stroke"
                                   when a.outcome_id = 13324 then "Heart failure"
                                   when a.outcome_id = 13325 then "Acute myocardial infarction (AMI)"
                                   when a.outcome_id = 13326 then "Sudden cardiac death"
                              end as outcome_name 
                              from blindingDecisions as a LEFT JOIN unBlindCheckManualCollapse as b
                              on a.database_id=b.database_id
                                 and a.target_id = b.target_id
                                 and a.target_name = b.target_name
                                 and a.comparator_id = b.comparator_id
                                 and a.outcome_id = b.outcome_id
                                 and a.analysis_id = b.analysis_id
                                 and a.comparator_name = b.comparator_name
                                 and a.analysis_name = b.analysis_name
                              where a.outcome_id != 13327')

# MC Note: change NAs to 0s
finalUnblindSummary$blind_GE1<-ifelse(is.na(finalUnblindSummary$blind_GE1),0,finalUnblindSummary$blind_GE1)
finalUnblindSummary$blind_GE2<-ifelse(is.na(finalUnblindSummary$blind_GE2),0,finalUnblindSummary$blind_GE2)

# Final decision on blinding
  # If absolute AND standardized differences are both past threshold then blind
  # Otherwise, only blind if one of its covariates failed the imbalance inspection for at least 1 of the manual reviewers
finalUnblindSummary$finalBlind_GE1 <- ifelse(finalUnblindSummary$blindStdzdiff == 1 & finalUnblindSummary$blindAbsdiff == 1,1,
                                     ifelse(finalUnblindSummary$blind_GE1 == 1, 1, 0))

  # If absolute AND standardized differences are both past threshold then blind
  # Otherwise, only blind if one of its covariates failed the imbalance inspection for at least 2 of the manual reviewers
finalUnblindSummary$finalBlind_GE2 <- ifelse(finalUnblindSufmmary$blindStdzdiff == 1 & finalUnblindSummary$blindAbsdiff == 1,1,
                                     ifelse(finalUnblindSummary$blind_GE2 == 1, 1, 0))

FinalUnblind <- sqldf('select * from finalUnblindSummary where finalBlind_GE1 = 0')
write.csv(FinalUnblind, file.path(studyFolder, "FinalUnblind_GE1_031820.csv"), row.names = FALSE)
write.csv(finalUnblindSummary, file.path(studyFolder, "FinalBlindingSummary_031820.csv"), row.names = FALSE)
