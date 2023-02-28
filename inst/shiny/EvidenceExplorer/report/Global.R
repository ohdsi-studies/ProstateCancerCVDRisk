source("DataPulls.R")
source("PlotsAndTables.R")

shinySettings <- list(dataFolder = "D:/StudyResults/EPI720_3/shinyDataAll", blind = TRUE)
dataFolder <- shinySettings$dataFolder
blind <- shinySettings$blind
connection <- NULL
positiveControlOutcome <- NULL

splittableTables <- c("covariate_balance", "preference_score_dist", "kaplan_meier_dist")

files <- list.files(dataFolder, pattern = ".rds")

# Find part to remove from all file names (usually databaseId):
databaseFileName <- files[grepl("^database", files)]
removePart <- paste0(gsub("database", "", databaseFileName), "$")

# Remove data already in global environment:
tableNames <- gsub("_t[0-9]+_c[0-9]+$", "", gsub(removePart, "", files)) 
camelCaseNames <- SqlRender::snakeCaseToCamelCase(tableNames)
camelCaseNames <- unique(camelCaseNames)
camelCaseNames <- camelCaseNames[!(camelCaseNames %in% SqlRender::snakeCaseToCamelCase(splittableTables))]
suppressWarnings(
  rm(list = camelCaseNames)
)
# Load data from data folder. R data objects will get names derived from the filename:
loadFile <- function(file) {
  tableName <- gsub("_t[0-9]+_c[0-9]+$", "", gsub(removePart, "", file)) 
  camelCaseName <- SqlRender::snakeCaseToCamelCase(tableName)
  if (!(tableName %in% splittableTables)) {
    newData <- readRDS(file.path(dataFolder, file))
    colnames(newData) <- SqlRender::snakeCaseToCamelCase(colnames(newData))
    if (exists(camelCaseName, envir = .GlobalEnv)) {
      existingData <- get(camelCaseName, envir = .GlobalEnv)
      newData <- rbind(existingData, newData)
    }
    assign(camelCaseName, newData, envir = .GlobalEnv)
  }
  invisible(NULL)
}
lapply(files, loadFile)

# create filter variables
primaryCohortIds <- c(13186, 13245)
priorCVDCohortIds <- c(13248, 13247)
noPriorCVDCohortIds <- c(13185, 13246)
severeCVDCohortIds <- c(13191, 13319)
lessSevereCVDCohortIds <- c(13187, 13316)
# afibCohortIds <- c(11400, 11403, 11401, 11402, 13033, 13234, 13199, 13200)
# vteCohortIds <- c(11404, 11407, 11405, 11406, 13203, 13235, 13204, 13205)
# cadPadCohortIds <- c(11408, 11411, 11409, 11410, 13207, 13236, 13208, 13209)
# thrTkrCohortIds <- c(11412, 11415, 11413, 11414, 13211, 13237, 13212, 13213)

primaryTarCohortIds <- primaryCohortIds
priorCVDTarCohortIds <- priorCVDCohortIds
noPriorCVDTarCohortIds <- noPriorCVDCohortIds
severeCVDTarCohortIds <- severeCVDCohortIds
lessSevereTarCVDCohortIds <- lessSevereCVDCohortIds
sensitivityTarTarCohortIds <- c(priorCVDCohortIds, noPriorCVDCohortIds, severeCVDCohortIds, lessSevereCVDCohortIds)

# primaryTarCohortIds <- c(11400, 11403, 11401, 11402, 11404, 11407, 11405, 11406,
#                          11408, 11411, 11409, 11410, 11412, 11415, 11413, 11414)
# sensitivityTarCohortIds <- c(13033, 13234, 13199, 13200, 13203, 13235, 13204, 13205,
#                              13207, 13236, 13208, 13209, 13211, 13237, 13212, 13213)

tcos <- unique(cohortMethodResult[, c("targetId", "comparatorId", "outcomeId")])
tcos <- tcos[tcos$outcomeId %in% outcomeOfInterest$outcomeId, ]
tcos$indicationId <- NA
tcos$indicationId[tcos$targetId %in% primaryCohortIds & tcos$comparatorId %in% primaryCohortIds] <- 1
tcos$indicationId[tcos$targetId %in% priorCVDCohortIds & tcos$comparatorId %in% priorCVDCohortIds] <- 2
tcos$indicationId[tcos$targetId %in% noPriorCVDCohortIds & tcos$comparatorId %in% noPriorCVDCohortIds] <- 3
tcos$indicationId[tcos$targetId %in% severeCVDCohortIds & tcos$comparatorId %in% severeCVDCohortIds] <- 4
tcos$indicationId[tcos$targetId %in% lessSevereCVDCohortIds & tcos$comparatorId %in% lessSevereCVDCohortIds] <- 5

# indications
exposureIndications <- data.frame(indicationId = c(1:5),
                                  indicationName = c("Primary Comparison",
                                                     "CVD in prior 180 days",
                                                     "No CVD in prior 180 days",
                                                     "Severe CVD (>=3) in prior 180 days",
                                                     "Less Severe CVD (<3) in prior 180 days"),
                                  stringsAsFactors = FALSE)

exposureOfInterest$definition <- NULL
exposureOfInterest$indicationId <- NA
exposureOfInterest$indicationId[exposureOfInterest$exposureId %in% primaryCohortIds] <- 1
exposureOfInterest$indicationId[exposureOfInterest$exposureId %in% priorCVDCohortIds] <- 2
exposureOfInterest$indicationId[exposureOfInterest$exposureId %in% noPriorCVDCohortIds] <- 3
exposureOfInterest$indicationId[exposureOfInterest$exposureId %in% severeCVDCohortIds] <- 4
exposureOfInterest$indicationId[exposureOfInterest$exposureId %in% lessSevereCVDCohortIds] <- 5

exposureOfInterest <- merge(exposureOfInterest, exposureIndications)

cohortMethodResult$indicationId <- NA
cohortMethodResult$indicationId[cohortMethodResult$targetId %in% primaryCohortIds & cohortMethodResult$comparatorId %in% primaryCohortIds] <- 1
cohortMethodResult$indicationId[cohortMethodResult$targetId %in% priorCVDCohortIds & cohortMethodResult$comparatorId %in% priorCVDCohortIds] <- 2
cohortMethodResult$indicationId[cohortMethodResult$targetId %in% noPriorCVDCohortIds & cohortMethodResult$comparatorId %in% noPriorCVDCohortIds] <- 3
cohortMethodResult$indicationId[cohortMethodResult$targetId %in% severeCVDCohortIds & cohortMethodResult$comparatorId %in% severeCVDCohortIds] <- 4
cohortMethodResult$indicationId[cohortMethodResult$targetId %in% lessSevereCVDCohortIds & cohortMethodResult$comparatorId %in% lessSevereCVDCohortIds] <- 5

primaryTar <- cohortMethodResult$targetId %in% primaryTarCohortIds & cohortMethodResult$comparatorId %in% primaryTarCohortIds & cohortMethodResult$analysisId %in% c(5)
  # MC Note <- I'm only going to produce reports on the analysis variant that was unblinded (5). Otherwise I could add additional analysisIds to the line below
#sensitivityTar <- cohortMethodResult$targetId %in% sensitivityTarCohortIds & cohortMethodResult$comparatorId %in% sensitivityTarCohortIds & cohortMethodResult$analysisId %in% c(5)
  # MC Note - I dropped the lines below - I don't think they're necessary since I don't have any redundant analyses
#ittTar <- cohortMethodResult$targetId %in% primaryTarCohortIds & cohortMethodResult$comparatorId %in% primaryTarCohortIds & cohortMethodResult$analysisId %in% c(2, 4)
#cohortMethodResult <- cohortMethodResult[primaryTar | sensitivityTar | ittTar, ] # this drops 1:1 SENS ITT amd 1:100 SENS ITT (since the same as 1:1 ITT and 1:100 ITT)

#cohortMethodAnalysis$definition <- NULL
#typeof(cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == 5])
#cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == 5]
#cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == 5] <- "As-treated (keep all) / 1:100 variable ratio matching"
#newAnalyses <- data.frame(analysisId = c(5, 6), 
#                          description = c("1:1 PS match / on-treatment (long gap days)",
#                                          "1:100 PS match / on-treatment (long gap days)"))
#cohortMethodAnalysis <- rbind(cohortMethodAnalysis, newAnalyses)
cohortMethodAnalysis$order <- match(cohortMethodAnalysis$analysisId, c(5))
cohortMethodAnalysis <- cohortMethodAnalysis[order(cohortMethodAnalysis$order), ]
cohortMethodAnalysis$order <- NULL

#analysisId5cm <- cohortMethodResult$targetId %in% sensitivityTarCohortIds & cohortMethodResult$comparatorId %in% sensitivityTarCohortIds & cohortMethodResult$analysisId == 1 & cohortMethodResult$databaseId != "Meta-analysis"
#analysisId6cm <- cohortMethodResult$targetId %in% sensitivityTarCohortIds & cohortMethodResult$comparatorId %in% sensitivityTarCohortIds & cohortMethodResult$analysisId == 3 & cohortMethodResult$databaseId != "Meta-analysis"
#cohortMethodResult$analysisId[analysisId5cm] <- 5
#cohortMethodResult$analysisId[analysisId6cm] <- 6
# 
# analysisId5fu <- cmFollowUpDist$targetId %in% sensitivityTarCohortIds & cmFollowUpDist$comparatorId %in% sensitivityTarCohortIds & cmFollowUpDist$analysisId == 1
# analysisId6fu <- cmFollowUpDist$targetId %in% sensitivityTarCohortIds & cmFollowUpDist$comparatorId %in% sensitivityTarCohortIds & cmFollowUpDist$analysisId == 3
# cmFollowUpDist$analysisId[analysisId5fu] <- 5
# cmFollowUpDist$analysisId[analysisId6fu] <- 6
# 
# analysisId5atr <- attrition$targetId %in% sensitivityTarCohortIds & attrition$comparatorId %in% sensitivityTarCohortIds & attrition$analysisId == 1
# analysisId6atr <- attrition$targetId %in% sensitivityTarCohortIds & attrition$comparatorId %in% sensitivityTarCohortIds & attrition$analysisId == 3
# attrition$analysisId[analysisId5atr] <- 5
# attrition$analysisId[analysisId6atr] <- 6
# 
# covariate56 <- covariate[covariate$analysisId %in% c(1, 3), ]
# covariate56$analysisId[covariate56$analysisId == 1] <- 5
# covariate56$analysisId[covariate56$analysisId == 3] <- 6
# covariate <- rbind(covariate, covariate56)
# rm(covariate56)
# 
# analysisId5pm <- propensityModel$targetId %in% sensitivityTarCohortIds & propensityModel$comparatorId %in% sensitivityTarCohortIds & propensityModel$analysisId == 1
# analysisId6pm <- propensityModel$targetId %in% sensitivityTarCohortIds & propensityModel$comparatorId %in% sensitivityTarCohortIds & propensityModel$analysisId == 3
# propensityModel$analysisId[analysisId5pm] <- 5
# propensityModel$analysisId[analysisId6pm] <- 6

  # MC Note: I removed the lines below since we are not completing a meta-analysis
#metaAnalysisId <- data.frame(databaseId = "Meta-analysis",
#                             databaseName = "Meta-analysis",
#                             description = "Meta-analysis",
#                             isMetaAnalysis = 1)
#database <- rbind(database, metaAnalysisId)

# renaming
exposureOfInterest$exposureName <- gsub("[EPI_720 - T01] ", "", exposureOfInterest$exposureName, fixed = TRUE)
exposureOfInterest$exposureName <- gsub("[EPI_720 - C01] ", "", exposureOfInterest$exposureName, fixed = TRUE)
#levels(exposureOfInterest$exposureName) <- gsub("[680] ", "", levels(exposureOfInterest$exposureName), fixed = TRUE)

outcomeOfInterest$definition <- NULL
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O01] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O02] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O03] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O04] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O05] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O06] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)
outcomeOfInterest$outcomeName <- gsub("[EPI_720 - O07] ", "", outcomeOfInterest$outcomeName, fixed = TRUE)

#levels(outcomeOfInterest$outcomeName) <- gsub("[680] ", "", levels(outcomeOfInterest$outcomeName), fixed = TRUE)