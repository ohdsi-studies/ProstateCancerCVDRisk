#source("S:/MiscCode/SetEnvironmentVariables.R")
source("C:/Users/admin_mconove1/Documents/MiscCode/ConnectionDetails.R")

connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "pdw",
                                                                server = Sys.getenv("server"),
                                                                user = NULL,
                                                                password = NULL,
                                                                port = Sys.getenv("port"))

options(fftempdir = "D:/FFTemp")
oracleTempSchema = NULL
baseUrl <- Sys.getenv("baseUrl")
studyResults <- "D:/StudyResults/EPI720_3/CohortDiagnostics_Outcomes"

writeToCsv <- function(data, fileName, incremental = FALSE, ...) {
  colnames(data) <- SqlRender::camelCaseToSnakeCase(colnames(data))
  if (incremental) {
    params <- list(...)
    names(params) <- SqlRender::camelCaseToSnakeCase(names(params))
    params$data = data
    params$fileName = fileName
    do.call(saveIncremental, params)
  } else {
    readr::write_csv(data, fileName)
  }
}


# Optum DOD settings ----------------------------------------------------------- done
databaseId <- "OptumDOD"
exportFolder <- file.path(studyResults, databaseId)
inclusionStatisticsFolder <- file.path(exportFolder, "incStats")
cdmDatabaseSchema <- "CDM_OPTUM_EXTENDED_DOD_V1194.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "stratPopAbirateroneAll"


# databases
database <- tibble::tibble(databaseId = databaseId,
                           databaseName = databaseId,
                           description = databaseId,
                           isMetaAnalysis = 0)
writeToCsv(database, file.path(exportFolder, "database.csv"))

# cohorts
cohort1 <- data.frame(cohort_full_name = "Abiraterone with HF outcome",
                      cohort_id = 1,
                      cohort_name = "Abiraterone with HF outcome",
                      json = "",
                      sql = "")
cohort2 <- data.frame(cohort_full_name = "Abiraterone with AMI outcome",
                      cohort_id = 2,
                      cohort_name = "Abiraterone with AMI outcome",
                      json = "",
                      sql = "")
cohort3 <- data.frame(cohort_full_name = "Abiraterone with ischemic stroke outcome",
                      cohort_id = 3,
                      cohort_name = "Abiraterone with ischemic stroke  outcome",
                      json = "",
                      sql = "")
cohorts <- rbind(cohort1, cohort2, cohort3)
writeToCsv(cohorts, file.path(exportFolder, "cohort.csv"))

# cohort counts
counts1 <- CohortDiagnostics::getCohortCounts(connectionDetails = connectionDetails,
                                              cohortDatabaseSchema = cohortDatabaseSchema,
                                              cohortTable = cohortTable,
                                              cohortIds = 1) 

counts2 <- CohortDiagnostics::getCohortCounts(connectionDetails = connectionDetails,
                                              cohortDatabaseSchema = cohortDatabaseSchema,
                                              cohortTable = cohortTable,
                                              cohortIds = 2) 

counts3 <- CohortDiagnostics::getCohortCounts(connectionDetails = connectionDetails,
                                              cohortDatabaseSchema = cohortDatabaseSchema,
                                              cohortTable = cohortTable,
                                              cohortIds = 3) 
counts <- cbind(rbind(counts1, counts2, counts3),databaseId)

writeToCsv(counts, file.path(exportFolder, "cohort_count.csv"), incremental = FALSE, cohortId = counts$cohortId)

covariateSettings <- FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE,
                                                                 useDemographicsAge = TRUE,
                                                                 useDemographicsAgeGroup = TRUE,
                                                                 useDemographicsRace = TRUE,
                                                                 useDemographicsEthnicity = TRUE,
                                                                 useDemographicsIndexYear = TRUE,
                                                                 useDemographicsIndexMonth = TRUE,
                                                                 useDemographicsPriorObservationTime = TRUE,
                                                                 useDemographicsPostObservationTime = TRUE,
                                                                 useDemographicsTimeInCohort = TRUE,
                                                                 useConditionOccurrenceAnyTimePrior = TRUE,
                                                                 useConditionOccurrenceLongTerm = TRUE,
                                                                 useConditionOccurrenceMediumTerm = TRUE,
                                                                 useConditionOccurrenceShortTerm = TRUE,
                                                                 useConditionEraAnyTimePrior = TRUE,
                                                                 useConditionEraLongTerm = TRUE,
                                                                 useConditionEraMediumTerm = TRUE,
                                                                 useConditionEraShortTerm = TRUE,
                                                                 useConditionEraOverlapping = TRUE,
                                                                 useConditionEraStartLongTerm = TRUE,
                                                                 useConditionEraStartMediumTerm = TRUE,
                                                                 useConditionEraStartShortTerm = TRUE,
                                                                 useConditionGroupEraAnyTimePrior = TRUE,
                                                                 useConditionGroupEraLongTerm = TRUE,
                                                                 useConditionGroupEraMediumTerm = TRUE,
                                                                 useConditionGroupEraShortTerm = TRUE,
                                                                 useConditionGroupEraOverlapping = TRUE,
                                                                 useConditionGroupEraStartLongTerm = TRUE,
                                                                 useConditionGroupEraStartMediumTerm = TRUE,
                                                                 useConditionGroupEraStartShortTerm = TRUE,
                                                                 useConditionOccurrencePrimaryInpatientLongTerm = TRUE,
                                                                 useDrugExposureAnyTimePrior = TRUE,
                                                                 useDrugExposureLongTerm = TRUE,
                                                                 useDrugExposureMediumTerm = TRUE,
                                                                 useDrugExposureShortTerm = TRUE,
                                                                 useDrugEraAnyTimePrior = TRUE,
                                                                 useDrugEraLongTerm = TRUE,
                                                                 useDrugEraMediumTerm = TRUE,
                                                                 useDrugEraShortTerm = TRUE,
                                                                 useDrugEraOverlapping = TRUE,
                                                                 useDrugEraStartLongTerm = TRUE,
                                                                 useDrugEraStartMediumTerm = TRUE,
                                                                 useDrugEraStartShortTerm = TRUE,
                                                                 useDrugGroupEraAnyTimePrior = TRUE,
                                                                 useDrugGroupEraLongTerm = TRUE,
                                                                 useDrugGroupEraMediumTerm = TRUE,
                                                                 useDrugGroupEraShortTerm = TRUE,
                                                                 useDrugGroupEraOverlapping = TRUE,
                                                                 useDrugGroupEraStartLongTerm = TRUE,
                                                                 useDrugGroupEraStartMediumTerm = TRUE,
                                                                 useDrugGroupEraStartShortTerm = TRUE,
                                                                 useProcedureOccurrenceAnyTimePrior = TRUE,
                                                                 useProcedureOccurrenceLongTerm = TRUE,
                                                                 useProcedureOccurrenceMediumTerm = TRUE,
                                                                 useProcedureOccurrenceShortTerm = TRUE,
                                                                 useDeviceExposureAnyTimePrior = TRUE,
                                                                 useDeviceExposureLongTerm = TRUE,
                                                                 useDeviceExposureMediumTerm = TRUE,
                                                                 useDeviceExposureShortTerm = TRUE,
                                                                 useMeasurementAnyTimePrior = TRUE,
                                                                 useMeasurementLongTerm = TRUE,
                                                                 useMeasurementMediumTerm = TRUE,
                                                                 useMeasurementShortTerm = TRUE,
                                                                 useMeasurementValueAnyTimePrior = TRUE,
                                                                 useMeasurementValueLongTerm = TRUE,
                                                                 useMeasurementValueMediumTerm = TRUE,
                                                                 useMeasurementValueShortTerm = TRUE,
                                                                 useMeasurementRangeGroupAnyTimePrior = TRUE,
                                                                 useMeasurementRangeGroupLongTerm = TRUE,
                                                                 useMeasurementRangeGroupMediumTerm = TRUE,
                                                                 useMeasurementRangeGroupShortTerm = TRUE,
                                                                 useObservationAnyTimePrior = TRUE,
                                                                 useObservationLongTerm = TRUE,
                                                                 useObservationMediumTerm = TRUE,
                                                                 useObservationShortTerm = TRUE,
                                                                 useCharlsonIndex = TRUE,
                                                                 useDcsi = TRUE,
                                                                 useChads2 = TRUE,
                                                                 useChads2Vasc = TRUE,
                                                                 useDistinctConditionCountLongTerm = TRUE,
                                                                 useDistinctConditionCountMediumTerm = TRUE,
                                                                 useDistinctConditionCountShortTerm = TRUE,
                                                                 useDistinctIngredientCountLongTerm = TRUE,
                                                                 useDistinctIngredientCountMediumTerm = TRUE,
                                                                 useDistinctIngredientCountShortTerm = TRUE,
                                                                 useDistinctProcedureCountLongTerm = TRUE,
                                                                 useDistinctProcedureCountMediumTerm = TRUE,
                                                                 useDistinctProcedureCountShortTerm = TRUE,
                                                                 useDistinctMeasurementCountLongTerm = TRUE,
                                                                 useDistinctMeasurementCountMediumTerm = TRUE,
                                                                 useDistinctMeasurementCountShortTerm = TRUE,
                                                                 useVisitCountLongTerm = TRUE,
                                                                 useVisitCountMediumTerm = TRUE,
                                                                 useVisitCountShortTerm = TRUE,
                                                                 longTermStartDays = -365,
                                                                 mediumTermStartDays = -180,
                                                                 shortTermStartDays = -30,
                                                                 endDays = 0,
                                                                 includedCovariateConceptIds = c(),
                                                                 addDescendantsToInclude = FALSE,
                                                                 excludedCovariateConceptIds = c(),
                                                                 addDescendantsToExclude = FALSE,
                                                                 includedCovariateIds = c())

# cohort characteristics
# HF
data1 <- CohortDiagnostics::getCohortCharacteristics(connectionDetails = connectionDetails,
                                                    oracleTempSchema = oracleTempSchema,
                                                    cdmDatabaseSchema = cdmDatabaseSchema,
                                                    cohortDatabaseSchema = cohortDatabaseSchema,
                                                    cohortTable = cohortTable,
                                                    cohortId = 1,
                                                    covariateSettings = covariateSettings)
data1 <- data1[round(data1$mean, 3) != 0, ]

covariates1 <- unique(data1[, c("covariateId", "covariateName", "analysisId")])
colnames(covariates1)[[3]] <- "covariateAnalysisId"
#writeToCsv(covariates, file.path(exportFolder, "covariate.csv"), incremental = FALSE, covariateId = covariates$covariateId)

data1$covariateName <- NULL
data1$analysisId <- NULL
if (nrow(data1) > 0) {
  data1$databaseId <- databaseId
  data1 <- merge(data1, counts1[, c("cohortId", "cohortEntries")])
  #data1 <- enforceMinCellValue(data1, "mean", minCellCount/data1$cohortEntries)
  data1$sd[data1$mean < 0] <- NA
  data1$cohortEntries <- NULL
  data1$mean <- round(data1$mean, 3)
  data1$sd <- round(data1$sd, 3)
}

# AMI
data2 <- CohortDiagnostics::getCohortCharacteristics(connectionDetails = connectionDetails,
                                                     oracleTempSchema = oracleTempSchema,
                                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                                     cohortTable = cohortTable,
                                                     cohortId = 2,
                                                     covariateSettings = covariateSettings)
data2 <- data2[round(data2$mean, 3) != 0, ]

# pile up the covariates files and dedup them (since they just refer to the variable labels and not the actual data) then save deduped file as covariates

covariates2 <- unique(data2[, c("covariateId", "covariateName", "analysisId")])
colnames(covariates2)[[3]] <- "covariateAnalysisId"
#writeToCsv(covariates, file.path(exportFolder, "covariate.csv"), incremental = FALSE, covariateId = covariates$covariateId)

data2$covariateName <- NULL
data2$analysisId <- NULL
if (nrow(data2) > 0) {
  data2$databaseId <- databaseId
  data2 <- merge(data2, counts2[, c("cohortId", "cohortEntries")])
  #data2 <- enforceMinCellValue(data2, "mean", minCellCount/data2$cohortEntries)
  data2$sd[data2$mean < 0] <- NA
  data2$cohortEntries <- NULL
  data2$mean <- round(data2$mean, 3)
  data2$sd <- round(data2$sd, 3)
}

# Stroke
data3 <- CohortDiagnostics::getCohortCharacteristics(connectionDetails = connectionDetails,
                                                     oracleTempSchema = oracleTempSchema,
                                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                                     cohortTable = cohortTable,
                                                     cohortId = 3,
                                                     covariateSettings = covariateSettings)

#data3store <- data3

data3 <- data3[round(data3$mean, 3) != 0, ]

covariates3 <- unique(data3[, c("covariateId", "covariateName", "analysisId")])
colnames(covariates3)[[3]] <- "covariateAnalysisId"
#writeToCsv(covariates, file.path(exportFolder, "covariate.csv"), incremental = FALSE, covariateId = covariates$covariateId)

data3$covariateName <- NULL
data3$analysisId <- NULL
if (nrow(data3) > 0) {
  data3$databaseId <- databaseId
  data3 <- merge(data3, counts3[, c("cohortId", "cohortEntries")])
  #data3 <- enforceMinCellValue(data3, "mean", minCellCount/data3$cohortEntries)
  data3$sd[data3$mean < 0] <- NA
  data3$cohortEntries <- NULL
  data3$mean <- round(data3$mean, 3)
  data3$sd <- round(data3$sd, 3)
}

covariates <- unique(rbind(covariates1, covariates2, covariates3))
writeToCsv(covariates, file.path(exportFolder, "covariate.csv"), incremental = FALSE, covariateId = covariates$covariateId)

data <- rbind(data1, data2, data3)
data <- unique(data[which(!is.na(data$covariateId)), ])
writeToCsv(data, file.path(exportFolder, "covariate_value.csv"), incremental = FALSE, cohortId = subset$cohortId)

#data[which(data$covariateId==1007 & data$cohortId==1),]
#data1[which(data1$covariateId==1007),]
#data2[which(data2$covariateId==1007),]
#data3[which(data3$covariateId==1007),]


zipName <- file.path(exportFolder, paste0("Results_", databaseId, ".zip"))
files <- list.files(exportFolder, pattern = ".*\\.csv$")
oldWd <- setwd(exportFolder)
on.exit(setwd(oldWd), add = TRUE)
DatabaseConnector::createZipFile(zipFile = zipName, files = files)


CohortDiagnostics::preMergeDiagnosticsFiles(exportFolder)
CohortDiagnostics::launchDiagnosticsExplorer(exportFolder)


#read.csv("D:/StudyResults/EPI720_3/CohortDiagnostics_Outcomes/OptumDOD/cohort_count.csv")
#read.csv("D:/StudyResults/EPI720_3/CohortDiagnostics_Outcomes/OptumDOD/cohort.csv")
#read.csv("D:/StudyResults/EPI720_3/CohortDiagnostics_Outcomes/OptumDOD/database.csv")

#unique(read.csv("D:/StudyResults/EPI720_3/CohortDiagnostics_Outcomes/OptumDOD/covariate.csv")
       
     #  covariates[which(covariates$covariateId==30753101),]
     #  data[which(data$covariateId==30753101),]
       
