source("C:/Users/admin_mconove1/Documents/MiscCode/ConnectionDetails.R")

library(EPI720)

options(fftempdir = "D:/FFTemp")
maxCores <- parallel::detectCores()
studyFolder <- "D:/StudyResults/EPI720_3"

# server connection:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "pdw",
                                                                server = Sys.getenv("server"),
                                                                user = NULL,
                                                                password = NULL,
                                                                port = Sys.getenv("port"))

#connection <- DatabaseConnector::connect(connectionDetails)

mailSettings <- list(from = Sys.getenv("emailAddress"),
                     to = c(Sys.getenv("emailAddress")),
                     smtp = list(host.name = Sys.getenv("emailHost"), port = 25,
                                 user.name = Sys.getenv("emailAddress"),
                                 passwd = Sys.getenv("emailPassword"), ssl = FALSE),
                     authenticate = FALSE,
                     send = TRUE)

# MDCR settings ----------------------------------------------------------------
databaseId <- "MDCR"
databaseName <- "MDCR"
databaseDescription <- "MDCR"
cdmDatabaseSchema = "CDM_IBM_MDCR_V1062.dbo"
outputFolder <- file.path(studyFolder, databaseName)
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "EPI720_MDCR"

# CCAE settings ----------------------------------------------------------------
#databaseId <- "CCAE"
#databaseName <- "CCAE"
#databaseDescription <- "CCAE"
#cdmDatabaseSchema = "CDM_IBM_CCAE_V1061.dbo"
#outputFolder <- file.path(studyFolder, databaseName)
#cohortDatabaseSchema <- "scratch.dbo"
#cohortTable <- "EPI720_CCAE"

# OPTUM_EXTENDED_DOD settings ----------------------------------------------------------------
databaseId <- "Optum"
databaseName <- "Optum"
databaseDescription <- "Optum DOD"
cdmDatabaseSchema = "CDM_OPTUM_EXTENDED_DOD_V1064.dbo"
outputFolder <- file.path(studyFolder, databaseName)
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "EPI720_Optum"

# Run --------------------------------------------------------------------------
OhdsiRTools::runAndNotify(expression = {
        execute(connectionDetails = connectionDetails,
                cdmDatabaseSchema = cdmDatabaseSchema,
                cohortDatabaseSchema = cohortDatabaseSchema,
                cohortTable = cohortTable,
                oracleTempSchema = NULL,
                outputFolder = outputFolder,
                databaseId = databaseId,
                databaseName = databaseName,
                databaseDescription = databaseDescription,
                createCohorts = FALSE,
                synthesizePositiveControls = FALSE,
                runAnalyses = FALSE,
                runDiagnostics = FALSE,
                packageResults = TRUE,
                maxCores = maxCores)
}, mailSettings = mailSettings, label = paste0("EPI720 ", databaseId), stopOnWarning = FALSE)

resultsZipFile <- file.path(outputFolder, "exportFull", paste0("Results", databaseId, ".zip"))
dataFolder <- file.path(outputFolder, "shinyData")
prepareForEvidenceExplorer(resultsZipFile = resultsZipFile, dataFolder = dataFolder)
launchEvidenceExplorer(dataFolder = dataFolder, blind = TRUE, launch.browser = FALSE)


compileShinyData(studyFolder = studyFolder, databases = c("MDCR", "Optum"))  
