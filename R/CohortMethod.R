# Copyright 2019 Observational Health Data Sciences and Informatics
#
# This file is part of EPI720
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Run CohortMethod package
#'
#' @details
#' Run the CohortMethod package, which implements the comparative cohort design.
#'
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cohortDatabaseSchema Schema name where intermediate data can be stored. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the exposure and outcome cohorts used in this
#'                             study.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param outputFolder         Name of local folder where the results were generated; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param maxCores             How many parallel cores should be used? If more cores are made available
#'                             this can speed up the analyses.
#'
#' @export
runCohortMethod <- function(connectionDetails,
                            cdmDatabaseSchema,
                            cohortDatabaseSchema,
                            cohortTable,
                            oracleTempSchema,
                            outputFolder,
                            timeAtRisk,
                            maxCores) {
  cmOutputFolder <- file.path(outputFolder, "cmOutput")
  if (!file.exists(cmOutputFolder)) {
    dir.create(cmOutputFolder)
  }
  cmAnalysisListFile <- system.file("settings",
                                    "cmAnalysisList.json",
                                    package = "EPI720")
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  
  
  if (timeAtRisk == "AsTreated") {
    cmAnalysisList <- cmAnalysisList[c(2,4,5,6,9,10)]
    # cmAnalysisList[2] "As-treated (keep first) / 1:1 matching"
    # cmAnalysisList[4] "As-treated (keep first) / 1:100 variable ratio matching"
    # cmAnalysisList[[5]] "As-treated (keep all) /  1:100 variable ratio matching"
    # cmAnalysisList[[6]] "As-treated (keep all) / 1:1 matching
    # cmAnalysisList[[9]] "As-treated (remove all) / 1:100 variable ratio matching"
    # cmAnalysisList[[10]] "As-treated (remove all) / 1:1 matching"

  }
  if (timeAtRisk == "IntentToTreat") {
    cmAnalysisList <- cmAnalysisList[c(1,3,7,8)]
    # cmAnalysisList[1] "Intent-to-treat (keep first) / 1:1 matching"
    # cmAnalysisList[3] "Intent-to-treat (keep first) / 1:100 variable ratio matching"
    # cmAnalysisList[[7]] "Intent-to-treat (remove all) / 1:100 variable ratio matching"
    # cmAnalysisList[[8]]  "Intent-to-treat (remove all)  / 1:1 matching"
  }
  
  
  tcosList <- createTcos(outputFolder = outputFolder)
  
  #  tcosList1 <- lapply(1:nrow(tcs), createTco)        
  
  #The code below was added by Jamie on Dec-12 2019
  
   pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "EPI720")
   tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE)
   exposureIds <- unique(c(tcosOfInterest$targetId, tcosOfInterest$comparatorId))
  
   pathToCohortCounts <- file.path(studyFolder,databaseName,"CohortCounts.csv") 
   CohortCounts <- read.csv(pathToCohortCounts, stringsAsFactors = FALSE)
   smallExposureCountId <- CohortCounts$cohortDefinitionId[CohortCounts$cohortDefinitionId %in% exposureIds & CohortCounts$personCount < 50]

  #The code below was added by Jamie on Dec-12 2019
#   newTcosList <- list()
#   for (i in 1:length(tcosList)) { # i = 1
#     if (tcosList[[i]]$targetId %in% smallExposureCountId | tcosList[[i]]$comparatorId %in% smallExposureCountId) {
#       tcosList[[i]] <- NULL
#     }
#     
#   }
  #The code above was added by Jamie on Dec-12 2019

   #The code below was added by Mitch on Dec-12 2019
      numDrops <- 0
      for (i in 1:length(tcosList)) { # i = 1
        if (tcosList[[i-numDrops]]$targetId %in% smallExposureCountId | tcosList[[i-numDrops]]$comparatorId %in% smallExposureCountId) {
          tcosList[[i-numDrops]] <- NULL
          numDrops <- numDrops+1
        }
      }
      numDrops <- 0
  #The code above was added by Mitch on Dec-12 2019
      
  #  tcosList2 <- merge(x=tcosList1,y=CohortCountsMC[,c("cohortDefinitionId","cohortCount")],by.x="targetId", by.y="cohortDefinitionId",all.x=TRUE)
  #  tcosList3 <- merge(x=tcosList2,y=CohortCountsMC[,c("cohortDefinitionId","cohortCount")],by.x="comparatorId", by.y="cohortDefinitionId",all.x=TRUE,
  #                     suffixes=c("Target","Comparator"))
  #  tcosList <- tcosList3[ which(tcosList3$cohortCountTarget >= 50 & tcosList3$cohortCountComparator >= 50 ), ]
  #  tcosList$cohortCountTarget <- tcosList$cohortCountComparator <- NULL
  
  outcomesOfInterest <- getOutcomesOfInterest()
  results <- CohortMethod::runCmAnalyses(connectionDetails = connectionDetails,
                                         cdmDatabaseSchema = cdmDatabaseSchema,
                                         exposureDatabaseSchema = cohortDatabaseSchema,
                                         exposureTable = cohortTable,
                                         outcomeDatabaseSchema = cohortDatabaseSchema,
                                         outcomeTable = cohortTable,
                                         outputFolder = cmOutputFolder,
                                         oracleTempSchema = oracleTempSchema,
                                         cmAnalysisList = cmAnalysisList,
                                         targetComparatorOutcomesList = tcosList,
                                         getDbCohortMethodDataThreads = min(3, maxCores),
                                         createStudyPopThreads = min(3, maxCores),
                                         createPsThreads = max(1, round(maxCores/10)),
                                         psCvThreads = min(10, maxCores),
                                         trimMatchStratifyThreads = min(10, maxCores),
                                         fitOutcomeModelThreads = max(1, round(maxCores/4)),
                                         outcomeCvThreads = min(4, maxCores),
                                         refitPsForEveryOutcome = FALSE,
                                         outcomeIdsOfInterest = outcomesOfInterest)
  
  ParallelLogger::logInfo("Summarizing results")
  analysisSummary <- CohortMethod::summarizeAnalyses(referenceTable = results, 
                                                     outputFolder = cmOutputFolder)
  analysisSummary <- addCohortNames(analysisSummary, "targetId", "targetName")
  analysisSummary <- addCohortNames(analysisSummary, "comparatorId", "comparatorName")
  analysisSummary <- addCohortNames(analysisSummary, "outcomeId", "outcomeName")
  analysisSummary <- addAnalysisDescription(analysisSummary, "analysisId", "analysisDescription")
  write.csv(analysisSummary, file.path(outputFolder, "analysisSummary.csv"), row.names = FALSE)
  
  ParallelLogger::logInfo("Computing covariate balance") 
  balanceFolder <- file.path(outputFolder, "balance")
  if (!file.exists(balanceFolder)) {
    dir.create(balanceFolder)
  }
  subset <- results[results$outcomeId %in% outcomesOfInterest,]
  subset <- subset[subset$strataFile != "", ]
  if (nrow(subset) > 0) {
    subset <- split(subset, seq(nrow(subset)))
    cluster <- ParallelLogger::makeCluster(min(3, maxCores))
    ParallelLogger::clusterApply(cluster, subset, computeCovariateBalance, cmOutputFolder = cmOutputFolder, balanceFolder = balanceFolder)
    ParallelLogger::stopCluster(cluster)
  }
}

computeCovariateBalance <- function(row, cmOutputFolder, balanceFolder) {
  outputFileName <- file.path(balanceFolder,
                              sprintf("bal_t%s_c%s_o%s_a%s.rds", row$targetId, row$comparatorId, row$outcomeId, row$analysisId))
  if (!file.exists(outputFileName)) {
    ParallelLogger::logTrace("Creating covariate balance file ", outputFileName)
    cohortMethodDataFolder <- file.path(cmOutputFolder, row$cohortMethodDataFolder)
    cohortMethodData <- CohortMethod::loadCohortMethodData(cohortMethodDataFolder)
    strataFile <- file.path(cmOutputFolder, row$strataFile)
    strata <- readRDS(strataFile)
    balance <- CohortMethod::computeCovariateBalance(population = strata, cohortMethodData = cohortMethodData)
    saveRDS(balance, outputFileName)
  }
}

addAnalysisDescription <- function(data, IdColumnName = "analysisId", nameColumnName = "analysisDescription") {
  cmAnalysisListFile <- system.file("settings",
                                    "cmAnalysisList.json",
                                    package = "EPI720")
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  idToName <- lapply(cmAnalysisList, function(x) data.frame(analysisId = x$analysisId, description = as.character(x$description)))
  idToName <- do.call("rbind", idToName)
  names(idToName)[1] <- IdColumnName
  names(idToName)[2] <- nameColumnName
  data <- merge(data, idToName, all.x = TRUE)
  # Change order of columns:
  idCol <- which(colnames(data) == IdColumnName)
  if (idCol < ncol(data) - 1) {
    data <- data[, c(1:idCol, ncol(data) , (idCol+1):(ncol(data)-1))]
  }
  return(data)
}

createTcos <- function(outputFolder) {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "EPI720")
  tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE)
  allControls <- getAllControls(outputFolder)
  tcs <- unique(rbind(tcosOfInterest[, c("targetId", "comparatorId")],
                      allControls[, c("targetId", "comparatorId")]))
  createTco <- function(i) {
    targetId <- tcs$targetId[i]
    comparatorId <- tcs$comparatorId[i]
    outcomeIds <- as.character(tcosOfInterest$outcomeIds[tcosOfInterest$targetId == targetId & tcosOfInterest$comparatorId == comparatorId])
    outcomeIds <- as.numeric(strsplit(outcomeIds, split = ";")[[1]])
    outcomeIds <- c(outcomeIds, allControls$outcomeId[allControls$targetId == targetId & allControls$comparatorId == comparatorId])
    excludeConceptIds <- as.character(tcosOfInterest$excludedCovariateConceptIds[tcosOfInterest$targetId == targetId & tcosOfInterest$comparatorId == comparatorId])
    if (length(excludeConceptIds) == 1 && is.na(excludeConceptIds)) {
      excludeConceptIds <- c()
    } else if (length(excludeConceptIds) > 0) {
      excludeConceptIds <- as.numeric(strsplit(excludeConceptIds, split = ";")[[1]])
    }
    includeConceptIds <- as.character(tcosOfInterest$includedCovariateConceptIds[tcosOfInterest$targetId == targetId & tcosOfInterest$comparatorId == comparatorId])
    if (length(includeConceptIds) == 1 && is.na(includeConceptIds)) {
      includeConceptIds <- c()
    } else if (length(includeConceptIds) > 0) {
      includeConceptIds <- as.numeric(strsplit(excludeConceptIds, split = ";")[[1]])
    }
    tco <- CohortMethod::createTargetComparatorOutcomes(targetId = targetId,
                                                        comparatorId = comparatorId,
                                                        outcomeIds = outcomeIds,
                                                        excludedCovariateConceptIds = excludeConceptIds,
                                                        includedCovariateConceptIds = includeConceptIds)
    return(tco)
  }
tcosList <- lapply(1:nrow(tcs), createTco)
#Mitch added the lines below on Dec-12, 2019
#  tcosList1 <- lapply(1:nrow(tcs), createTco)        
#  pathToCohortCounts <- file.path(studyFolder,databaseName,"CohortCounts.csv") 
#  CohortCountsMC <- read.csv(pathToCohortCounts, stringsAsFactors = FALSE)
#  tcosList2 <- merge(x=tcosList1,y=CohortCountsMC[,c("cohortDefinitionId","cohortCount")],by.x="targetId", by.y="cohortDefinitionId",all.x=TRUE)
#  tcosList3 <- merge(x=tcosList2,y=CohortCountsMC[,c("cohortDefinitionId","cohortCount")],by.x="comparatorId", by.y="cohortDefinitionId",all.x=TRUE,
#                     suffixes=c("Target","Comparator"))
#  tcosList <- tcosList3[ which(tcosList3$cohortCountTarget >= 50 & tcosList3$cohortCountComparator >= 50 ), ]
#  tcosList$cohortCountTarget <- tcosList$cohortCountComparator <- NULL
#Mitch added the lines above on Dec-12, 2019
  return(tcosList)
}

getOutcomesOfInterest <- function() {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "EPI720")
  tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE) 
  outcomeIds <- as.character(tcosOfInterest$outcomeIds)
  outcomeIds <- do.call("c", (strsplit(outcomeIds, split = ";")))
  outcomeIds <- unique(as.numeric(outcomeIds))
  return(outcomeIds)
}

getAllControls <- function(outputFolder) {
  allControlsFile <- file.path(outputFolder, "AllControls.csv")
  if (file.exists(allControlsFile)) {
    # Positive controls must have been synthesized. Include both positive and negative controls.
    allControls <- read.csv(allControlsFile)
  } else {
    # Include only negative controls
    pathToCsv <- system.file("settings", "NegativeControls.csv", package = "EPI720")
    allControls <- read.csv(pathToCsv)
    allControls$oldOutcomeId <- allControls$outcomeId
    allControls$targetEffectSize <- rep(1, nrow(allControls))
  }
  return(allControls)
}
