source("C:/Users/admin_mconove1/Documents/MiscCode/ConnectionDetails.R")

setwd("D:/StudyResults/EPI720_3/Optum/AsTreated/cmOutput/")
getwd()

# server connection:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "pdw",
                                                                server = Sys.getenv("server"),
                                                                user = NULL,
                                                                password = NULL,
                                                                port = Sys.getenv("port"))


# The code below was used to confirm / identify the correct StratPop_*.rds file that corresponds to the heart failure analysis that we are interested in analyzing further

# Primary comparison for HF outcome
# T=13186
# C=13245
# O=13324

#13186,[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC,13186,EPI_720__T01_New_Users_of_abiraterone__prednisone_with_CRPC
#13245,[EPI_720 - C01] New users of enzalutamide with CRPC,13245,EPI_720__C01_New_users_of_enzalutamide_with_CRPC
#13321,[EPI_720 - O02] Cardiovascular morbidity - ischemic stroke,13321,EPI_720__O02_Cardiovascular_morbidity__ischemic_stroke
#13324,[EPI_720 - O04] Cardiovascular morbidity - heart failure,13324,EPI_720__O04_Cardiovascular_morbidity__heart_failure
#13325,[EPI_720 - O05] Cardiovascular morbidity - acute myocardial infarction (AMI),13325,EPI_720__O05_Cardiovascular_morbidity__acute_myocardial_infarction_AMI
# 
# test_l1_s1 <- readRDS("StratPop_l1_s1_p1_t13186_c13245_s1_o13324.rds", refhook = NULL)
# test_l2_s1 <- readRDS("StratPop_l2_s1_p1_t13186_c13245_s1_o13324.rds", refhook = NULL)
# test_l3_s1 <- readRDS("StratPop_l3_s1_p1_t13186_c13245_s1_o13324.rds", refhook = NULL)
# 
# test_l1_s2 <- readRDS("StratPop_l1_s1_p1_t13186_c13245_s2_o13324.rds", refhook = NULL)
# test_l2_s2 <- readRDS("StratPop_l2_s1_p1_t13186_c13245_s2_o13324.rds", refhook = NULL)
# test_l3_s2 <- readRDS("StratPop_l3_s1_p1_t13186_c13245_s2_o13324.rds", refhook = NULL)
# 
# test_l1_s1_E1 <- test_l1_s1[which(test_l1_s1$treatment==1),]
# test_l1_s1_E0 <- test_l1_s1[which(test_l1_s1$treatment==0),]
#   test_l2_s1_E1 <- test_l2_s1[which(test_l2_s1$treatment==1),]
#   test_l2_s1_E0 <- test_l2_s1[which(test_l2_s1$treatment==0),]
# test_l3_s1_E1 <- test_l3_s1[which(test_l3_s1$treatment==1),]
# test_l3_s1_E0 <- test_l3_s1[which(test_l3_s1$treatment==0),]
# 
# test_l1_s2_E1 <- test_l1_s2[which(test_l1_s2$treatment==1),]
# test_l1_s2_E0 <- test_l1_s2[which(test_l1_s2$treatment==0),]
#   test_l2_s2_E1 <- test_l2_s2[which(test_l2_s2$treatment==1),]
#   test_l2_s2_E0 <- test_l2_s2[which(test_l2_s2$treatment==0),]
# test_l3_s2_E1 <- test_l3_s2[which(test_l3_s2$treatment==1),]
# test_l3_s2_E0 <- test_l3_s2[which(test_l3_s2$treatment==0),]
# 
# length((test_l1_s1_E1$subjectId))
# length((test_l1_s1_E0$subjectId))
# 
# length((test_l1_s2_E1$subjectId))
# length((test_l1_s2_E0$subjectId))
# 
# length((test_l2_s1_E1$subjectId))
# length((test_l2_s1_E0$subjectId))
# 
# length(unique(test_l2_s2_E1$subjectId))
# length(unique(test_l2_s2_E0$subjectId))
# 
# length((test_l3_s1_E1$subjectId))
# length((test_l3_s1_E0$subjectId))
# 
# length((test_l3_s2_E1$subjectId))
# length((test_l3_s2_E0$subjectId))


# Load the stratPop_*.rds files and generate Kaplan Meier plots

# HF
  stratPopHF <- readRDS("StratPop_l2_s1_p1_t13186_c13245_s2_o13324.rds", refhook = NULL)
  # Subset to the treated (abiraterone) group
  stratPopHFAbiraterone <- stratPopHF[which(stratPopHF$treatment==1),]
  # Subset to those who experienced outcomes during follow-up
  stratPopAbirateroneHF <- stratPopHF[which(stratPopHF$treatment==1 & stratPopHF$outcomeCount != 0),]
  stratPopEnzalutamideHF <- stratPopHF[which(stratPopHF$treatment==0 & stratPopHF$outcomeCount != 0),]
  
  
  kmAll <- readRDS("D:/StudyResults/EPI720_3/shinyDataAll/kaplan_meier_dist_t13186_c13245.rds", refhook = NULL)
  kmAll_HF <- kmAll[kmAll$analysis_id==5 & kmAll$database_id == "Optum" & kmAll$outcome_id==13324,]
  names(kmAll_HF) <- SqlRender::snakeCaseToCamelCase(names(kmAll_HF))
  plotKaplanMeierShiny(kmAll_HF, targetName="Abiraterone", comparatorName="Enzalutamide")
  
  #table(kmAll$analysis_id, kmAll$outcome_id, kmAll$database_id)

# AMI
  stratPopAMI <- readRDS("StratPop_l2_s1_p1_t13186_c13245_s2_o13325.rds", refhook = NULL)
  # Subset to the treated (abiraterone) group
  stratPopAMIAbiraterone <- stratPopAMI[which(stratPopAMI$treatment==1),]
  # Subset to those who experienced outcomes during follow-up
  stratPopAbirateroneAMI <- stratPopAMI[which(stratPopAMI$treatment==1 & stratPopAMI$outcomeCount != 0),]
  stratPopEnzalutamideAMI <- stratPopAMI[which(stratPopAMI$treatment==0 & stratPopAMI$outcomeCount != 0),]
  
  kmAll <- readRDS("D:/StudyResults/EPI720_3/shinyDataAll/kaplan_meier_dist_t13186_c13245.rds", refhook = NULL)
  kmAll_AMI <- kmAll[kmAll$analysis_id==5 & kmAll$database_id == "Optum" & kmAll$outcome_id==13325,]
  names(kmAll_AMI) <- SqlRender::snakeCaseToCamelCase(names(kmAll_AMI))
  plotKaplanMeierShiny(kmAll_AMI, targetName="Abiraterone", comparatorName="Enzalutamide")


# Ischemic stroke
  stratPopStroke <- readRDS("StratPop_l2_s1_p1_t13186_c13245_s2_o13321.rds", refhook = NULL)
  # Subset to the treated (abiraterone) group
  stratPopStrokeAbiraterone <- stratPopStroke[which(stratPopStroke$treatment==1),]
  # Subset to those who experienced outcomes during follow-up
  stratPopAbirateroneStroke <- stratPopStroke[which(stratPopStroke$treatment==1 & stratPopStroke$outcomeCount != 0),]
  stratPopEnzalutamideStroke <- stratPopStroke[which(stratPopStroke$treatment==0 & stratPopStroke$outcomeCount != 0),]
  
  kmAll <- readRDS("D:/StudyResults/EPI720_3/shinyDataAll/kaplan_meier_dist_t13186_c13245.rds", refhook = NULL)
  kmAll_Stroke <- kmAll[kmAll$analysis_id==5 & kmAll$database_id == "Optum" & kmAll$outcome_id==13321,]
  names(kmAll_Stroke) <- SqlRender::snakeCaseToCamelCase(names(kmAll_Stroke))
  plotKaplanMeierShiny(kmAll_Stroke, targetName="Abiraterone", comparatorName="Enzalutamide")
  

install.packages("rlang")
library(rlang)

# Generating the median event times for the different treatment/comparator/outcomes of interest 

median(stratPopAbirateroneHF$survivalTime)
median(stratPopEnzalutamideHF$survivalTime)

median(stratPopAbirateroneAMI$survivalTime)
median(stratPopEnzalutamideAMI$survivalTime)

median(stratPopAbirateroneStroke$survivalTime)
median(stratPopEnzalutamideStroke$survivalTime)


quantile(stratPopAbirateroneHF$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))
quantile(stratPopEnzalutamideHF$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))

quantile(stratPopAbirateroneAMI$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))
quantile(stratPopEnzalutamideAMI$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))

quantile(stratPopAbirateroneStroke$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))
quantile(stratPopEnzalutamideStroke$survivalTime, c(0, .01, .05, .1, .25, .5, .75, .90, .95, .99, 1))



# Load the insertDbPopulation function so that we can insert the cohorts into the database

#' Insert a population into a database
#'
#' @details
#' Inserts a population table into a database. The table in the database will have the same structure
#' as the COHORT table in the Common Data Model.
#'
#' @param population             Either an object of type [CohortMethodData] or a population
#'                               object generated by functions like [createStudyPopulation()].
#' @param cohortIds              The IDs to be used for the target and comparator cohort,
#'                               respectively.
#' @param connectionDetails      An R object of type `connectionDetails` created using the
#'                               [DatabaseConnector::createConnectionDetails()] function.
#' @param cohortDatabaseSchema   The name of the database schema where the data will be written.
#'                               Requires write permissions to this database. On SQL Server, this
#'                               should specify both the database and the schema, so for example
#'                               'cdm_instance.dbo'.
#' @param cohortTable            The name of the table in the database schema where the data will be
#'                               written.
#' @param createTable            Should a new table be created? If not, the data will be inserted into
#'                               an existing table.
#' @param dropTableIfExists      If `createTable = TRUE` and the table already exists it will be
#'                               overwritten.
#' @param cdmVersion             Define the OMOP CDM version used: currently supports "5".
#'
#' @export
insertDbPopulation <- function(population,
                               cohortIds = c(1, 0),
                               connectionDetails,
                               cohortDatabaseSchema,
                               cohortTable = "cohort",
                               createTable = FALSE,
                               dropTableIfExists = TRUE,
                               cdmVersion = "5") {
  if (is(population, "cohortMethodData")) {
    population <- population$cohorts %>%
      collect()
  }
  newCohortIds <- plyr::mapvalues(population$treatment, c(1, 0), cohortIds)
  population <- population[, c("subjectId", "cohortStartDate")]
  population$cohortDefinitionId <- newCohortIds
  population$cohortEndDate <- NA
  colnames(population) <- SqlRender::camelCaseToSnakeCase(colnames(population))
  connection <- DatabaseConnector::connect(connectionDetails)
  ParallelLogger::logInfo("Writing ",
                          nrow(population),
                          " rows to ",
                          cohortDatabaseSchema, ".",
                          cohortTable)
  start <- Sys.time()
  if (!createTable) {
    if (cdmVersion == "4") {
      sql <- "DELETE FROM @table WHERE cohort_concept_id IN (@cohort_ids);"
    } else {
      sql <- "DELETE FROM @table WHERE cohort_definition_id IN (@cohort_ids);"
    }
    sql <- SqlRender::render(sql,
                             table = paste(cohortDatabaseSchema, cohortTable, sep = "."),
                             cohort_ids = cohortIds)
    sql <- SqlRender::translate(sql, targetDialect = connectionDetails$dbms)
    DatabaseConnector::executeSql(connection = connection,
                                  sql = sql,
                                  progressBar = FALSE,
                                  reportOverallTime = FALSE)
  }
  DatabaseConnector::insertTable(connection = connection,
                                 tableName = paste(cohortDatabaseSchema, cohortTable, sep = "."),
                                 data = population,
                                 dropTableIfExists = dropTableIfExists,
                                 createTable = createTable,
                                 tempTable = FALSE,
                                 oracleTempSchema = NULL)
  DatabaseConnector::disconnect(connection)
  delta <- Sys.time() - start
  ParallelLogger::logInfo(paste("Inserting rows took", signif(delta, 3), attr(delta, "units")))
  invisible(TRUE)
}

# Insert the populations into the scratch.dbo database

insertDbPopulation(population = stratPopAbirateroneHF,
                 #  cohortIds = c(101,100),
                   connectionDetails = connectionDetails,
                   cohortDatabaseSchema = "scratch.dbo",
                   cohortTable = "stratPopabirateroneHF",
                   createTable = TRUE,
                   dropTableIfExists = TRUE,
                   cdmVersion = "5")


insertDbPopulation(population = stratPopAbirateroneAMI,
                   #  cohortIds = c(101,100),
                   connectionDetails = connectionDetails,
                   cohortDatabaseSchema = "scratch.dbo",
                   cohortTable = "stratPopabirateroneAMI",
                   createTable = TRUE,
                   dropTableIfExists = TRUE,
                   cdmVersion = "5")


insertDbPopulation(population = stratPopAbirateroneStroke,
                   #  cohortIds = c(101,100),
                   connectionDetails = connectionDetails,
                   cohortDatabaseSchema = "scratch.dbo",
                   cohortTable = "stratPopabirateroneStroke",
                   createTable = TRUE,
                   dropTableIfExists = TRUE,
                   cdmVersion = "5")

# Double-checking that the populations were properly inserted into the scratch.dbo database
library(CohortMethod)
library(SqlRender)
connection <- DatabaseConnector::connect(connectionDetails)
DatabaseConnector::querySql(connection, "SELECT * FROM scratch.dbo.stratPopabirateroneHF")
DatabaseConnector::querySql(connection, "SELECT * FROM scratch.dbo.stratPopabirateroneAMI")
DatabaseConnector::querySql(connection, "SELECT * FROM scratch.dbo.stratPopabirateroneStroke")
