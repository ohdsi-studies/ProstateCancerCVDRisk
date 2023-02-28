
# This code was written during a 1-hour troubleshooting session where Mitch & Jamie tried to figure out a (possible) error where 20% of the T/C cohorts appear to have
# hematologic neoplasms. We tried manually running feature extractoin to pull that variable but it produces 0 rows.  There may be some issues with what is rolling up
# to the hematologic neoplasms concept.  For now, we are deleting the hematologic neoplasms row from all of the Table 1s included in the topline report.
# we will troubleshoot this issue during the extension of the EPI_720 project analyzing prior/subsequent exposure to Chemo and ADT, which is currently being submitted to Wrike

Sys.getenv("PATH")

Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                        "C:\\Users\\admin_mconove1\\Documents\\R\\rtools40",
                        "C:\\Users\\admin_mconove1\\Documents\\R\\rtools40\\mingw64\\bin", 
                        sep = ";"))
# Did it work (look at the end)?
Sys.getenv("PATH")

test <- readRDS("D:/StudyResults/EPI720_3/Optum/AsTreated/cmOutput/outcomeModelReference.rds")

testOmFile <- "D:/StudyResults/EPI720_3/Optum/AsTreated/cmOutput/CmData_l2_t13186_c13245/"

test2 <- CohortMethod::loadCohortMethodData(file.path(testOmFile))


test2<- readRDS(testOmFile)
covRetest <- as.data.frame(test2$covariateRef)



remotes::install_github("ohdsi/CohortMethod", ref = "v3.0.2", dependencies=TRUE)

# Confirming 20% prevalence of hematologic malignancies
settings <- FeatureExtraction::createCovariateSettings(#includedCovariateIds=4044013210,
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
                                                       includedCovariateConceptIds=4044013)

covariateData <- FeatureExtraction::getDbCovariateData(connectionDetails = connectionDetails,
                                    cdmDatabaseSchema = cdmDatabaseSchema,
                                    cohortDatabaseSchema = cohortDatabaseSchema,
                                    cohortTable = cohortTable,
                                    cohortId = 13186,
                                    rowIdField = "subject_id",
                                    covariateSettings = settings,
                                    aggregated = TRUE)

as.data.frame(covariateData$covariateRef)
as.data.frame(covariateData$covariate)

