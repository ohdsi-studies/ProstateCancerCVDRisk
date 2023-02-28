setwd("C:/Users/admin_mconove1/Documents/epi_720/EPI720_3/inst/shiny/EvidenceExplorer/")
#reportFolder <- "./report"
reportFolder <- "C:/Users/admin_mconove1/Documents/epi_720/EPI720_3/inst/shiny/EvidenceExplorer/report"
#source("global.R")
source("C:/Users/admin_mconove1/Documents/epi_720/EPI720_3/inst/shiny/EvidenceExplorer/report/global.R")
source("C:/Users/admin_mconove1/Documents/epi_720/EPI720_3/inst/shiny/EvidenceExplorer/report/ReportPlotsAndTables.R")
library("magrittr")
library("flextable")
library("CohortMethod")

source("C:/Users/admin_mconove1/Documents/MiscCode/ConnectionDetails.R")

library(EPI720)


# server connection:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "pdw",
                                                                server = Sys.getenv("server"),
                                                                user = NULL,
                                                                password = NULL,
                                                                port = Sys.getenv("port"))


# Global parameters ------------------------------------------------------------
#13186	[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC
#13245	[EPI_720 - C01] New users of enzalutamide with CRPC
#
# 13248	[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC (no CVD in prior 6-mo)
# 13247	[EPI_720 - C01] New users of enzalutamide with CRPC (no CVD in prior 6-mo)
# 
# 13185	[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC (any CVD in prior 6-mo)
# 13246	[EPI_720 - C01] New users of enzalutamide with CRPC (any CVD in prior 6-mo)
# 
# 13191	[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC (uncontrolled CVD in prior 6-mo)
# 13319	[EPI_720 - C01] New users of enzalutamide with CRPC (uncontrolled CVD in prior 6-mo)
# 
# 13187	[EPI_720 - T01] New Users of abiraterone + prednisone with CRPC (controlled CVD in prior 6-mo)
# 13316	[EPI_720 - C01] New users of enzalutamide with CRPC (controlled CVD in prior 6-mo)

databaseIds <- database$databaseId
outcomeId <- 13324 # <- c(13324, 13325, 13321)
outcomeLabels <- c("HF", "AMI", "Ischemic Stroke")
cbind(c(13324, 13325, 13321),c("Heart Failure", "Acute Myocardial Infarction", "Ischemic Stroke"))


outcomeLabels <-
data.frame("outcomeId" = c(13320, 13321,13323,13324,13325,13326,13327), 
           "outcomeName"= c("Cardiovascular morbidity (composite)",
                             "Ischemic stroke","Hemorrhagic stroke","Heart failure","Acute myocardial infarction","Sudden cardiac death","Alzheimers/dementia/cognative decline"))

cohortMethodAnalysis$shortDescription <- 
  ifelse(cohortMethodAnalysis$description == "As-treated (keep all) /  1:100 variable ratio matching", "AT/keep all/1:100",
         ifelse(cohortMethodAnalysis$description == "As-treated (keep first) / 1:1 matching", "AT/keep 1st/1:1",
                ifelse(cohortMethodAnalysis$description == "As-treated (keep first) / 1:100 variable ratio matching", "AT/keep 1st/1:100",
                       ifelse(cohortMethodAnalysis$description == "As-treated (keep all) / 1:1 matching", "AT/keep all/1:1",
                              ifelse(cohortMethodAnalysis$description == "As-treated (remove all) / 1:100 variable ratio matching", "AT/remove all/1:100",
                                     ifelse(cohortMethodAnalysis$description == "As-treated (remove all) / 1:1 matching", "AT/remove all/1:1",
                                            ifelse(cohortMethodAnalysis$description == "Intent-to-treat (keep first) / 1:1 matching", "ITT/keep 1st/1:1",
                                                   ifelse(cohortMethodAnalysis$description == "Intent-to-treat (keep first) / 1:100 variable ratio matching", "ITT/keep 1st/1:100",
                                                          ifelse(cohortMethodAnalysis$description == "Intent-to-treat (remove all) / 1:100 variable ratio matching", "ITT/remove all/1:100",
                                                                 ifelse(cohortMethodAnalysis$description == "Intent-to-treat (remove all)  / 1:1 matching", "ITT/remove all/1:1","ERROR"))))))))))



# 13320,"Cardiovascular morbidity (ischemic stroke, hemorrhagic stroke, heart failure, acute myocardial infarction or sudden cardiac death)"
# 13321,"Ischemic stroke"
# 13323,"Hemorrhagic stroke"
# 13324,"Heart failure"
# 13325,"Acute myocardial infarction"
# 13326,"Sudden cardiac death"
# 13327,"Alzheimers/dementia/cognative decline"


indicationLabels <-c("PrimaryComparison","PriorCVD","NoPriorCVD","SevereCVD","LessSevereCVD")
primaryCohortIds <- c(13186, 13245)
priorCVDCohortIds <- c(13248, 13247)
noPriorCVDCohortIds <- c(13185, 13246)
severeCVDCohortIds <- c(13191, 13319)
lessSevereCVDCohortIds <- c(13187, 13316)

#outcomeGroups <- list(list(outcome = "HF", cohortId = 13324), 
#                      list(outcome = "AMI", cohortId = 13325),
#                      list(outcome = "IschemStroke", cohortId = 13321))
indicationGroups <- list(list(indication = "PrimaryComparison", cohortIDs = primaryCohortIds),
                         list(indication = "PriorCVD", cohortIDs = priorCVDCohortIds),
                         list(indication = "NoPriorCVD", cohortIDs = noPriorCVDCohortIds),
                         list(indication = "SevereCVD", cohortIDs = severeCVDCohortIds),
                         list(indication = "LessSevereCVD", cohortIDs = lessSevereCVDCohortIds))

# databaseIds <- database$databaseId
# outcomeId <- 14347
# indicationLabels <- c("NVAF", "VTE", "TKR/THR")
# nvafCohortIds <- c(11400, 11403, 11401, 11402, 13033, 13234, 13199, 13200)
# vteCohortIds <- c(11404, 11407, 11405, 11406, 13203, 13235, 13204, 13205)
# thrTkrCohortIds <- c(11412, 11415, 11413, 11414, 13211, 13237, 13212, 13213)
# indicationGroups <- list(list(indication = "NVAF", cohortIds = nvafCohortIds),
#                          list(indication = "VTE", cohortIds = vteCohortIds),
#                          list(indication = "THR/TKR", cohortIds = thrTkrCohortIds))

blank <- ""
blankPlot <- ggplot2::ggplot() + ggplot2::theme_void()

headingFormat <- officer::fp_text(font.size = 12, bold = TRUE, font.family = "Calibri")
titleFormat <- officer::fp_text(font.size = 10, bold = FALSE, font.family = "Calibri")

indicationTitles <- list()
indicationTitles[[1]] <- officer::fpar(officer::ftext("PrimaryComparison", prop = titleFormat))
indicationTitles[[7]] <- officer::fpar(officer::ftext("PriorCVD", prop = titleFormat))
indicationTitles[[13]] <- officer::fpar(officer::ftext("NoPriorCVD", prop = titleFormat))
indicationTitles[[19]] <- officer::fpar(officer::ftext("SevereCVD", prop = titleFormat))
indicationTitles[[25]] <- officer::fpar(officer::ftext("LessSevereCVD", prop = titleFormat))


outcomeIds <- outcomeOfInterest$outcomeId
databaseIds <- database$databaseId
exposureOfInterest$shortName <- exposureOfInterest$exposureName
exposureOfInterest$shortName <- gsub("New users of ", "", gsub(" with CRPC", "", gsub("New Users of ", "", exposureOfInterest$exposureName, fixed = TRUE), fixed = TRUE), fixed = TRUE)


# 13327,"Alzheimers/dementia/cognative decline"
# 13320,"Cardiovascular morbidity (ischemic stroke, hemorrhagic stroke, heart failure, acute myocardial infarction or sudden cardiac death)"
# 13323,"Hemorrhagic stroke"
# 13326,"Sudden cardiac death"
# 13321,"Ischemic stroke"
# 13324,"Heart failure"
# 13325,"Acute myocardial infarction"



# indicationTitles <- list()
# indicationTitles[[1]] <- officer::fpar(officer::ftext("NVAF", prop = titleFormat))
# indicationTitles[[7]] <- officer::fpar(officer::ftext("VTE", prop = titleFormat))
# indicationTitles[[13]] <- officer::fpar(officer::ftext("TKR/THR", prop = titleFormat))

#outcomeTitles <- list()
#outcomeTitles[[1]] <- officer::fpar(officer::ftext("HF", prop = titleFormat))
#outcomeTitles[[7]] <- officer::fpar(officer::ftext("AMI", prop = titleFormat))
#outcomeTitles[[13]] <- officer::fpar(officer::ftext("IschemStroke", prop = titleFormat))

# Exposure cohort counts--------------------------------------------------------

section1 <- "Study population counts"
section1 <- officer::fpar(officer::ftext(section1, prop = headingFormat))

counts <- getPatientCounts(attrition, "Original cohorts") # 1st exposure, 1st cohort only, restrict to common period, 365 days of obs prior, no prior outcome
counts$subjects <- formatC(counts$subjects, big.mark = ",", format = "d")
counts$exposureSubjects <- formatC(counts$exposureSubjects, big.mark = ",", format = "d")
counts$databaseSubjects <- formatC(counts$databaseSubjects, big.mark = ",", format = "d")


# by exposure
exposureCountTitle <- "Counts by exposure cohort"
exposureCountTitle <- officer::fpar(officer::ftext(exposureCountTitle, prop = titleFormat))
countsByExposure <- counts
countsByExposure$indicationName <- sub('\\:.*', "", countsByExposure$exposureName)
countsByExposure$exposureName <-  capitalize(sub(".*? ", "", countsByExposure$exposureName))
countsByExposure$totalExposureCount <- sprintf("%s (n = %s)", countsByExposure$exposureName, counts$exposureSubjects)
#MC Note: I added the line below to incorporate outcomes in the output
countsByExposure <- merge(countsByExposure, outcomeLabels, by="outcomeId")

countsByExposure <- countsByExposure[, c("outcomeName", "totalExposureCount", "databaseId", "subjects", "exposurePercent")]
countsByExposure$subjects <- gsub("^-", "<", countsByExposure$subjects)
countsByExposure$exposurePercent <- gsub("^-", "<", countsByExposure$exposurePercent)
header <- c("Outcome", "Exposure", "Database", "Patients", "Percent")
# MC Note: I had to revise the line below to make it work (original line commented out below)
names(countsByExposure) <- header
# MC Note: I had to add the line below to make it work
countsByExposure <- countsByExposure[order(countsByExposure$Outcome, countsByExposure$Exposure, countsByExposure$Database),]
#countsByExposure <- rbind(header, countsByExposure)
countsByExposure <- createCountsFlextable(countsByExposure)

# by database
countsByDatabase <- counts
databaseCountTitle <- "Counts by database"
databaseCountTitle <- officer::fpar(officer::ftext(databaseCountTitle, prop = titleFormat))
countsByDatabase$databaseOrder <- match(countsByDatabase$databaseId, countsByDatabase$databaseId)
countsByDatabase$exposureOrder <- match(countsByDatabase$exposureName, countsByDatabase$exposureName)
countsByDatabase <- countsByDatabase[order(countsByDatabase$databaseOrder, countsByDatabase$exposureOrder), ]
countsByDatabase[, c("databaseOrder", "exposureOrder")] <- NULL
countsByDatabase$totalDatabaseCount <- sprintf("%s (n = %s)", countsByDatabase$databaseId, countsByDatabase$databaseSubjects)
#MC Note: I added the line below to incorporate outcomes in the output
countsByDatabase <- merge(countsByDatabase, outcomeLabels, by="outcomeId")
countsByDatabase$indicationName <- sub('\\:.*', "", countsByDatabase$exposureName)
countsByDatabase$exposureName <-  capitalize(sub(".*? ", "", countsByDatabase$exposureName))
countsByDatabase <- countsByDatabase[, c("outcomeName", "totalDatabaseCount", "exposureName", "subjects", "databasePercent")]
countsByDatabase$subjects <- gsub("^-", "<", countsByDatabase$subjects)
countsByDatabase$databasePercent <- gsub("^-", "<", countsByDatabase$databasePercent)
header <- c("Outcome", "Database", "Exposure", "Patients", "Percent")
# MC Note: I had to revise the line below to make it work (original line commented out below)
names(countsByDatabase) <- header
# MC Note: I had to add the line below to make it work
countsByDatabase <- countsByDatabase[order(countsByDatabase$Outcome, countsByDatabase$Database, countsByDatabase$Exposure),]
#countsByDatabase <- rbind(header, countsByDatabase)
countsByDatabase <- createCountsFlextable(countsByDatabase)

doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section1, style = "heading 1") %>%
  officer::body_add_fpar(exposureCountTitle, style = "heading 2") %>%
  flextable::body_add_flextable(countsByExposure) %>%
  officer::body_add_break() %>%
  officer::body_add_fpar(databaseCountTitle, style = "heading 2") %>%
  flextable::body_add_flextable(countsByDatabase) %>% # Having an issue with this line - in WordPad the second table doesn't display
  print(target = file.path(reportFolder, "studyPopCounts.docx"))



# Baseline characteristics -----------------------------------------------------

section2 <- "Patient baseline characteristics"
section2 <- officer::fpar(officer::ftext(section2, prop = headingFormat))

table1s <- list()
table1Titles <- list()

for (analysisId in c(5)) {  # analysisId <- 5
  analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
  analysisName <- sub(".*? ", "", analysisName)
  tcs <- tcos[tcos$targetId %in% primaryTarCohortIds & tcos$comparatorId %in% primaryTarCohortIds & tcos$outcomeId %in% c(13321,13324,13325), ]
  
  for (databaseId in c("Optum")){ #databaseIds[databaseIds != "Meta-analysis"]) { # databaseId <- "MDCD"
    for (i in 1:nrow(tcs)) { # i=14
      # MC Note: I had to revise the line below to make it work (original line commented out below)
      targetLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcs$targetId[i]]
      #targetLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$targetId[i]]
      targetLabel <- capitalize(sub(".*? ", "", targetLabel))
      # MC Note: I had to revise the line below to make it work (original line commented out below)
      comparatorLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcs$comparatorId[i]]
      #comparatorLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$comparatorId[i]]
      comparatorLabel <- capitalize(sub(".*? ", "", comparatorLabel))
      # MC Note: I had to add the two lines below to provide outcome labels in my output
      outcomeLabel <- outcomeLabels$outcomeName[outcomeLabels$outcomeId == tcs$outcomeId[i]]
      outcomeLabel <- capitalize(sub(".*? ", "", outcomeLabel))
      
      balance <- getCovariateBalance(targetId = tcs$targetId[i],
                                     comparatorId = tcs$comparatorId[i],
                                     databaseId = databaseId,
                                     analysisId = analysisId,
     # MC Note: I had to revise the line below to make it work (original line commented out below)
                                     outcomeId = tcs$outcomeId[i]) 
                                     #outcomeId = outcomeId) 
      if (!is.null(balance)) {
        table1 <- prepareTable1(balance,
                                targetLabel = sub(".*? ", "", targetLabel),
                                comparatorLabel = comparatorLabel)
        facs <- sapply(table1, is.factor)
        table1[facs] <- lapply(table1[facs], as.character)
        rm(facs)
        table1 <- rbind(c(blank, "Before matching", blank, blank, "After matching", blank, blank), table1)
        colnames(table1) <- letters[1:length(colnames(table1))]
        
        table1 <- flextable::qflextable(table1)
        table1 <- flextable::delete_part(table1, part = "header")
        table1 <- flextable::fontsize(table1, part = "all", size = 6)
        table1 <- flextable::align(table1, j = 1, align = "left", part = "all")
        table1 <- flextable::merge_h(table1, i = 1, part = "body")
        table1 <- flextable::autofit(table1, add_w = 0, add_h = 0)
        table1 <- flextable::padding(table1, padding = 0, part = "all")
        border <- officer::fp_border(color = "black", width = 1)
        table1 <- flextable::border_inner(table1, border = border, part = "all")
        table1 <- flextable::border_outer(table1, border = border, part = "all")
      } else {
        table1 <- data.frame(a = NA, b = NA, c = NA)
        table1 <- flextable::qflextable(table1)
      }
      table1s[[length(table1s) + 1]] <- table1
      title <- sprintf("%s vs. %s", targetLabel, comparatorLabel)
      #title <- sprintf("%s (%s), %s vs. %s", databaseId, analysisName, targetLabel, comparatorLabel)
      title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
      table1Titles[[length(table1Titles) + 1]] <- title
    }
  }
}
table1Pairs <- list(table1Titles, table1s)
doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section2, style = "heading 1")
for(i in 1:length(table1s)) { #i=1
#   if (i %in% c(1, 7, 13)) {
    doc <- doc %>% 
      officer::body_add_fpar( officer::fpar(officer::ftext(paste(databaseId,analysisName), 
                                                           prop = titleFormat)), style = "heading 2") %>%
      officer::body_add_fpar(table1Pairs[[1]][[i]], style = "heading 3") %>%
      #MC Note: below added
      officer::body_add_fpar( officer::fpar(officer::ftext(as.character(outcomeLabels$outcomeName[outcomeLabels$outcomeId==tcs$outcomeId[i]]), 
                                                           prop = titleFormat)), style = "heading 3") %>%
      #MC Note: above added
      flextable::body_add_flextable(table1Pairs[[2]][[i]]) %>%
      officer::body_add_break()
    
#  } else {
#    doc <- doc %>%
#      officer::body_add_fpar(table1Pairs[[1]][[i]], style = "heading 3") %>%
#      flextable::body_add_flextable(table1Pairs[[2]][[i]]) %>%
#      officer::body_add_break()
#  }
 }
print(doc, target = file.path(reportFolder, "table1s.docx"))


# IR tables after matching ------------------------------------------------------------------

section4 <- "Incidence rates"
section4 <- officer::fpar(officer::ftext(section4, prop = headingFormat))

eventTables <- list()
eventTableTitles <- list()

#exposureOfInterest$shortName[exposureOfInterest$exposureId == 2] <- "HCQ"
#exposureOfInterest$shortName[exposureOfInterest$exposureId == 28] <- "SSZ"
#To produce estimates 
#getdbcohortmethoddata
#createstudypopulation
#rather than stratapop - studypop
# studypop objects have removed prior outcomes

# The code below pulls post-matching exposure/event counts and IR estimates
# Action 1 - confirm dina wants crude tables


for (analysisId in c(5)) { # analysisId=5
  analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
  # if you want crude tables - incorporate getIRIP function here and run on corresponding studypop file.
  mainResults <- getMainResults(targetIds = 13186, 
                                comparatorIds = exposureOfInterest$exposureId, 
                                outcomeIds = c(13321,13324,13325), 
                                databaseIds = "Optum",
                                analysisIds = analysisId) 
  mainResults <- mainResults[!is.na(mainResults$seLogRr), ]
  #irTable <- prepareReportIrTable(mainResults, outcomeOfInterest[which(outcomeOfInterest$outcomeId %in% c("13321","13324","13325")),])
  irTable <- prepareReportIrTable(mainResults, outcomeOfInterest)
  irTable$outcomeOrder <- match(irTable$outcomeId,outcomeIds)
  irTable$dbOrder <- match(irTable$databaseId, databaseIds)
  irTable <- irTable[order(irTable$outcomeOrder, irTable$dbOrder, irTable$indicationName),]
  irTable[, c("outcomeId", "targetId", "outcomeOrder", "dbOrder")] <- NULL
 # if (analysisId == 1) {
#    irTable <- irTable[!(irTable$outcomeName == "Hospitalization for psychosis" & irTable$databaseId == "Meta-analysis"), ]  
#  }
  header1 <- c(blank, blank, blank, "Patients", blank, "TAR", blank, "Events", blank, "IR", blank, blank)
  header2 <- c("Outcome", "Prior CVD Stratification", "Database", rep(c("T", "C"), 4), "MDRR")    
  irTable <- rbind(header1, header2, irTable)
  irTable <- createIrFlexTable(irTable)
  eventTables[[length(eventTables) + 1]] <- irTable
  title <- paste0("Optum: ",analysisName)
  title <- officer::fpar(officer::ftext(title, prop = titleFormat))
  eventTableTitles[[length(eventTableTitles) + 1]] <- title 
}
section3 <- "Post-matching: Event counts and incidence rates"
section3 <- officer::fpar(officer::ftext(section3, prop = headingFormat))
eventTablePairs <- list(eventTableTitles, eventTables)
 doc <- officer::read_docx() %>% 
   officer::body_add_fpar(section3, style = "heading 1")
 for(i in 1:length(eventTables)) { #i=1
   doc <- doc %>% 
     officer::body_add_fpar(eventTablePairs[[1]][[i]], style = "heading 2") %>%
     officer::body_add_break() %>%
     flextable::body_add_flextable(eventTablePairs[[2]][[i]]) %>%
     officer::body_add_break()
 }
 print(doc, target = file.path(reportFolder, "irTables.docx"))
 

 
# Crude IRs -------------------------------------------------------------------------

crudeIrFolder <- file.path(reportFolder, "crudeIrFolder")
if (!file.exists(crudeIrFolder)) {
  dir.create(crudeIrFolder)
}

section3 <- "Unadjusted incidence rates"
section3 <- officer::fpar(officer::ftext(section3, prop = headingFormat))

#source("S:/MiscCode/SetEnvironmentVariables.R")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "pdw",
                                                                server = Sys.getenv("server"),
                                                                user = NULL,
                                                                password = NULL,
                                                                port = Sys.getenv("port"))

covariateSettings <- FeatureExtraction::createCovariateSettings(useDemographicsAge = TRUE,
                                                                useDemographicsIndexYear = TRUE)

dbs <- data.frame(name = c("MDCR", "Optum"), 
                  cdmDatabaseSchema = c("CDM_IBM_MDCR_V1062.dbo", "CDM_OPTUM_EXTENDED_DOD_V1064.dbo"),
                  #cohortDatabaseSchema = c("CDM_IBM_MDCR_V1062.dbo", "CDM_OPTUM_EXTENDED_DOD_V1064.dbo"), # 
                  cohortDatabaseSchema = rep("scratch.dbo", 2),
                  cohortTable = c("EPI720_MDCR", "EPI720_OPTUM"), # Ask Jamie - do I need to revise this?
                  stringsAsFactors = FALSE)

irResults <- data.frame()
attritionResults <- data.frame()

for (i in 1:nrow(dbs)) { # i=1
  cdmDatabaseSchema <- dbs$cdmDatabaseSchema[i]
  cohortDatabaseSchema <- dbs$cohortDatabaseSchema[i]
  cohortTable <- dbs$cohortTable[i] # ASK JAMIE - Do I need to revise this
  dbName <- dbs$name[i]
  
  for (j in 1:nrow(exposureOfInterest)) { # j=1
    if (exposureOfInterest$exposureId[j] %in% primaryTarCohortIds) {
      strataName <- "All patients"
    } else if (exposureOfInterest$exposureId[j] %in% priorCVDTarCohortIds) {
      strataName <- "Prior CVD"
    } else if (exposureOfInterest$exposureId[j] %in% noPriorCVDTarCohortIds) {
      strataName <- "No Prior CVD"
    } else if (exposureOfInterest$exposureId[j] %in% severeCVDTarCohortIds) {
      strataName <- "Severe CVD"
    } else if (exposureOfInterest$exposureId[j] %in% lessSevereTarCVDCohortIds) {
      strataName <- "Less Severe CVD"
    }
    
    for (h in 1:nrow(outcomeLabels)) { # h=4
      exposureId <- exposureOfInterest$exposureId[j]
      exposureName <- exposureOfInterest$exposureName[j]
      #indication <- sub('\\:.*', "", exposureName)
      exposureName <- capitalize(sub(".*? ", "", exposureName))
      #MC Note: I added the lines below
      outcomeId <- outcomeLabels$outcomeId[h]
      outcomeName <- as.character(outcomeLabels$outcomeName[h])
  
      studyPopFile <- sprintf("studyPop_e%s_o%s_%s.rds", exposureId, outcomeId, dbName)
      studyPopFile <- file.path(crudeIrFolder, studyPopFile)
      if (!file.exists(studyPopFile)) {
        cmData <- CohortMethod::getDbCohortMethodData(connectionDetails = connectionDetails,
                                                      studyStartDate = "20120831",
                                                      studyEndDate = "",
                                                      excludeDrugsFromCovariates = FALSE,
                                                      cdmDatabaseSchema = cdmDatabaseSchema,
                                                      targetId = exposureId,
                                                      comparatorId = 0,
                                                      outcomeIds = outcomeId,
                                                      exposureDatabaseSchema = cohortDatabaseSchema,
                                                      exposureTable = cohortTable,
                                                      outcomeDatabaseSchema = cohortDatabaseSchema,
                                                      outcomeTable = cohortTable,
                                                      covariateSettings = covariateSettings,
                                                      firstExposureOnly = FALSE,
                                                      restrictToCommonPeriod = FALSE,
                                                      removeDuplicateSubjects = "keep all")
        studyPop <- CohortMethod::createStudyPopulation(cohortMethodData = cmData, 
                                                        outcomeId = outcomeId,
                                                        removeSubjectsWithPriorOutcome = TRUE,
                                                        minDaysAtRisk = 1,
                                                        startAnchor = "cohort start",
                                                        riskWindowStart = 30,
                                                        endAnchor = "cohort end",
                                                        riskWindowEnd = 0) # <- set to zero to produce the as-treated analysis numbers
        saveRDS(studyPop, studyPopFile)
      } else {
        studyPop <- readRDS(studyPopFile) 
      }
      irIpTable <- getIrIp(studyPop)
      irIpTable <- cbind(strataName, outcomeName, exposureName, dbName, irIpTable)
      irResults <- rbind(irResults, irIpTable)
      
      studyPopAttrition <- attr(studyPop, "metaData")$attrition
      studyPopAttrition <- studyPopAttrition[, c(1,2)]
      studyPopAttrition <- t(studyPopAttrition)
      colNames <- studyPopAttrition[1, ]
      studyPopAttrition <- as.data.frame(studyPopAttrition, row.names = FALSE)
      studyPopAttrition <- studyPopAttrition[-1, ]
      names(studyPopAttrition) <- colNames
      studyPopAttrition <- cbind(strataName, outcomeName, exposureName, dbName, studyPopAttrition)
      attritionResults <- rbind(attritionResults, studyPopAttrition)
      
      # if (exposureOfInterest$exposureId[j] %in% primaryTarCohortIds) {
      #   analysisName <- "ITT" # ASK JAMIE - WHAT IS THIS?
      #   tarId <- 3
      #   studyPopFile <- sprintf("studyPop_e%s_tar%s_%s.rds", exposureId, tarId, dbName)
      #   studyPopFile <- file.path(crudeIrFolder, studyPopFile)
      #   if (!file.exists(studyPopFile)) {
      #     studyPop <- CohortMethod::createStudyPopulation(cohortMethodData = cmData, 
      #                                                     outcomeId = outcomeId[i],
      #                                                     removeSubjectsWithPriorOutcome = TRUE,
      #                                                     minDaysAtRisk = 1,
      #                                                     startAnchor = "cohort start",
      #                                                     riskWindowStart = 30,
      #                                                     endAnchor = "cohort start",
      #                                                     riskWindowEnd = 9999)
      #     saveRDS(studyPop, studyPopFile)
      #   } else {
      #     studyPop <- readRDS(studyPopFile)
      #   }
      #   irIpTable <- getIrIp(studyPop)
      #   irIpTable <- cbind(analysisName, outcomeName, exposureName, dbName, irIpTable)
      #   irResults <- rbind(irResults, irIpTable)
      #   
      #   studyPopAttrition <- attr(studyPop, "metaData")$attrition
      #   studyPopAttrition <- studyPopAttrition[, c(1,2)]
      #   studyPopAttrition <- t(studyPopAttrition)
      #   colNames <- studyPopAttrition[1, ]
      #   studyPopAttrition <- as.data.frame(studyPopAttrition, row.names = FALSE)
      #   studyPopAttrition <- studyPopAttrition[-1, ]
      #   names(studyPopAttrition) <- colNames
      #   studyPopAttrition <- cbind(analysisName, outcomeName, exposureName, dbName, studyPopAttrition)
      #   attritionResults <- rbind(attritionResults, studyPopAttrition)
      #   }
    }
  }
}
crudeIrTable <- irResults
crudeIrTable$exposureName <- gsub("Users", "New users", crudeIrTable$exposureName, fixed = TRUE)
crudeIrTable$analysisOrder <- match(crudeIrTable$strataName, c("All patients", "Prior CVD","No Prior CVD","Severe CVD","Less Severe CVD"))
#crudeIrTable$indicationOrder <- match(crudeIrTable$indication, c("NVAF", "VTE", "TKR/THR"))
#crudeIrTable$exposureOrder <- match(crudeIrTable$exposureName, c("Rivaroxaban", "Apixaban", "Dabigatran", "Warfarin"))
crudeIrTable$dbOrder <- match(crudeIrTable$dbName, databaseIds)
crudeIrTable <- crudeIrTable[order(crudeIrTable$analysisOrder, crudeIrTable$dbOrder), ]
crudeIrTable[, c("analysisOrder", "dbOrder")] <- NULL
keeps <- c("events", "pys", "ir1k", "ip1k")
#onTreatment <- crudeIrTable[crudeIrTable$analysisName == "On-treatment", -1]
#onTreatment30 <- crudeIrTable[crudeIrTable$analysisName == "On-treatment (30-day gap)", keeps]
#itt <- crudeIrTable[crudeIrTable$analysisName == "ITT", keeps]
#crudeIrTable <- cbind(onTreatment, onTreatment30, itt)
facs <- sapply(crudeIrTable, is.factor)
crudeIrTable[facs] <- lapply(crudeIrTable[facs], as.character)
rm(facs)
#header1 <- c(blank, blank, blank, blank, rep("On-treatment", 4), rep("On-treatment (30-day gap)", 4), rep("Intent-to-treat", 4))
header2 <- c("Strata", "Outcome", "Exposure", "Database", "Patients", rep(c("Events", "PYs", "IR/1,000 PYs", "IP/1,1000 P"), 1))
crudeIrTable <- rbind(header2, crudeIrTable)
names(crudeIrTable) <- letters[1:ncol(crudeIrTable)]

crudeIrTable <- flextable::qflextable(crudeIrTable)
crudeIrTable <- flextable::delete_part(crudeIrTable, part = "header")
crudeIrTable <- flextable::fontsize(crudeIrTable, part = "all", size = 6)
crudeIrTable <- flextable::merge_v(crudeIrTable, j = 1:3, part = "body")
crudeIrTable <- flextable::merge_h(crudeIrTable, i = 1, part = "body")
border <- officer::fp_border(color = "black", width = 1)
crudeIrTable <- flextable::border_inner(crudeIrTable, border = border, part = "all")
crudeIrTable <- flextable::border_outer(crudeIrTable, border = border, part = "all")
crudeIrTable <- flextable::align(crudeIrTable, j = 1:3, align = "left", part = "all")
crudeIrTable <- flextable::align(crudeIrTable, i = 1:2, align = "left", part = "all")
crudeIrTable <- flextable::autofit(crudeIrTable, add_w = 0.1, add_h = 0.1)
crudeIrTable <- flextable::padding(crudeIrTable, padding = 0, part = "all")

doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section3, style = "heading 1") %>%
  flextable::body_add_flextable(crudeIrTable) %>%
  officer::body_add_break() %>%
  print(target = file.path(reportFolder, "crudeIrTable.docx"))

write.csv(attritionResults, file.path(reportFolder, "crudeIrAttrition.csv"), row.names = FALSE)
# 
# # After matching IRs ----------------------------------------------------------
# 
# section4 <- "After matching incidence rates"
# section4 <- officer::fpar(officer::ftext(section4, prop = headingFormat))
# 
# eventTables <- list()
# eventTableTitles <- list()
# 
# header1 <- c(blank, blank, blank, "Patients", blank, "TAR", blank, "Events", blank, "IR", blank, blank)
# header2 <- c("Indication", "T vs C", "Database", rep(c("T", "C"), 4), "MDRR")
# 
# for (analysisId in c(1, 3, 8, 9, 2, 4)) { # analysisId <- 8
#   if (analysisId == 8) { # restrict to primary comparisons for short gap follow-up and ITT follow-up
#     titleSuffix <- " (30-day gap)"
#     analysisId <- 1
#     tarCohortIds <- sensitivityTarCohortIds
#   } else if (analysisId == 9) {
#     titleSuffix <- " (30-day gap)"
#     analysisId <- 3
#     tarCohortIds <- sensitivityTarCohortIds
#   } else {
#     titleSuffix <- ""
#     tarCohortIds <- primaryTarCohortIds
#   }
#   analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
#   analysisName <- paste0(analysisName, titleSuffix)
#   
#   mainResults <- getMainResults(targetIds = primaryTarCohortIds, 
#                                 comparatorIds = primaryTarCohortIds, 
#                                 outcomeIds = outcomeOfInterest$outcomeId, 
#                                 databaseIds = databaseIds,
#                                 analysisIds = analysisId)
#   #mainResults <- mainResults[!is.na(mainResults$seLogRr), ]
#   irTable <- prepareReportIrTable(mainResults, exposureOfInterest)
#   irTable <- reorderTable(irTable, indicationLabels)
#   irTable <- rbind(header1, header2, irTable)
#   irTable <- createFlextable(irTable)
#   title <- analysisName
#   title <- officer::fpar(officer::ftext(title, prop = titleFormat))
#   eventTables[[length(eventTables) + 1]] <- irTable
#   eventTableTitles[[length(eventTableTitles) + 1]] <- title 
# }
# eventTablePairs <- list(eventTableTitles, eventTables)
# doc <- officer::read_docx() %>% 
#   officer::body_add_fpar(section4, style = "heading 1")
# for(i in 1:length(eventTables)) { #i=1
#   doc <- doc %>% 
#     officer::body_add_fpar(eventTablePairs[[1]][[i]], style = "heading 2") %>%
#     flextable::body_add_flextable(eventTablePairs[[2]][[i]]) %>%
#     officer::body_add_break()
# }
# print(doc, target = file.path(reportFolder, "irTables.docx"))




# diagnostics for main paper (includes only unblinded diagnostics)) ------------------------------------------------------------------

section3 <- "Study dianosticis"
section3 <- officer::fpar(officer::ftext(section3, prop = headingFormat))

diagnosticsTitles <- list()
diagnosticsFileNames <- list()
compSummary <- comparisonSummary[which(comparisonSummary$targetId == 13186),]
comparisonSummaries <- rbind(compSummary[which(compSummary$databaseId == "MDCR"), ], compSummary[which(compSummary$databaseId == "Optum"), ])

#comparisonSummaries <- list(comparisonSummary[1:5, ], comparisonSummary[7:11, ])
#comparisonSummaries <- list(comparisonSummary[c(1,2,3,4,5,7,8,9,10,11), ])

#plotNumber <- 0
for (i in 1:nrow(comparisonSummaries)) { # i=2
  # print("print cs")
  # print(cs)
  # cs <- comparisonSummaries[2] # DELETE THIS LINE
  # plotNumber <- plotNumber + 1
  psPlots <- list()
  balPlots <- list()
  calPlots <- list()
  
  fileName <- file.path(reportFolder, sprintf("diagnostics_%s.png", comparisonSummaries$databaseId[i]))
  diagnosticsFileNames[[length(diagnosticsFileNames) + 1]] <- fileName
  title <- sprintf("Diagnostics plot: %s", comparisonSummaries$databaseId[i])
  title <- officer::fpar(officer::ftext(title, prop = titleFormat))
  diagnosticsTitles[[length(diagnosticsTitles) + 1]] <- title
  
  #for (i in 1:nrow(cs)) { # i = 1  should this be as.data.frame(cs)?
  targetId <- comparisonSummaries$targetId[i]
  targetLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$targetId[i]]
  comparatorId <- comparisonSummaries$comparatorId[i]
  comparatorLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$comparatorId[i]]
  databaseId <- comparisonSummaries$databaseId[i]
  for (analysisId in c(5)) { # analysisId = 5
    ps <- getPs(connection, targetId, comparatorId, analysisId = analysisId, databaseId) #Doesn't this need to be outcome specific?  The problem here appears to be that in the shinyDataAll folder the file preference_score_dist_t13191_c13319.rds doesn't exist (model only includes intercept term)
    psPlot <- plotPs2(ps, targetLabel, comparatorLabel)
    psPlots[[length(psPlots) + 1]] <- psPlot
    bal <- getCovariateBalance(connection, targetId, comparatorId, databaseId, analysisId = analysisId, outcomeIds[4])
    balPlot <- plotCovariateBalanceScatterPlot2(bal,beforeLabel = "Before Matching",afterLabel = "After Matching")
    balPlots[[length(balPlots) + 1]] <- balPlot
    
    # for (analysisId in cohortMethodAnalysis$analysisId) { # analysisId = 2
    controlResults <- getControlResults(connection = connection,
                                        targetId = targetId,
                                        comparatorId = comparatorId,
                                        analysisId = analysisId,
                                        databaseId = databaseId)
    controlEstimates <- controlResults[controlResults$effectSize == 1, ]
    calPlot <- plotLargeScatter2(controlEstimates, xLabel = "HR")
    calPlots[[length(calPlots) + 1]] <- calPlot
  }
  # }
  col0 <- grid::textGrob("")
  col1 <- grid::textGrob("PS distribution", gp = grid::gpar(fontsize = 12))
  col2 <- grid::textGrob("Covariate balanace", gp = grid::gpar(fontsize = 12))
  col3 <- grid::textGrob("Null dist.", gp = grid::gpar(fontsize = 12))
  row1 <- grid::textGrob(cohortMethodAnalysis$shortDescription[5], rot = 90, gp = grid::gpar(fontsize = 8))
  plotGrid <- gridExtra::grid.arrange(col0, col1, col2, col3, #col4,
                                      row1, psPlots[[1]], balPlots[[1]], calPlots[[1]],
                                      nrow = 2,
                                      heights = c(0.5, 4),
                                      widths = c(0.25, 3.5, 3.5, 3.5))
  ggplot2::ggsave(fileName, plotGrid, width = 8, height = 4, dpi = 400)
}
diagnosticsPlotPairs <- list(diagnosticsTitles, diagnosticsFileNames)
doc <- officer::read_docx() %>%
  officer::body_add_fpar(section3, style = "heading 1")
for(i in 1:length(diagnosticsFileNames)) { #i=1
  doc <- doc %>%
    officer::body_add_fpar(diagnosticsTitles[[i]], style = "heading 2") %>%
    officer::body_add_img(diagnosticsPlotPairs[[2]][[i]], width = 5, height = 2) %>%
    officer::body_add_break()
}
print(doc, target = file.path(reportFolder, "diagnosticPlotsUnblinded.docx"))




# diagnostics for supplemental appendix (includes all analyses) ------------------------------------------------------------------

section3 <- "Study dianosticis"
section3 <- officer::fpar(officer::ftext(section3, prop = headingFormat))

diagnosticsTitles <- list()
diagnosticsFileNames <- list()
compSummary <- comparisonSummary[which(comparisonSummary$targetId == 13186),]
comparisonSummaries <- rbind(compSummary[which(compSummary$databaseId == "MDCR"), ], compSummary[which(compSummary$databaseId == "Optum"), ])

#comparisonSummaries <- list(comparisonSummary[1:5, ], comparisonSummary[7:11, ])
#comparisonSummaries <- list(comparisonSummary[c(1,2,3,4,5,7,8,9,10,11), ])

#plotNumber <- 0
for (i in 1:nrow(comparisonSummaries)) { # i=2

  psPlots <- list()
  balPlots <- list()
  calPlots <- list()
  
  fileName <- file.path(reportFolder, sprintf("diagnostics_%s.png", comparisonSummaries$databaseId[i]))
  diagnosticsFileNames[[length(diagnosticsFileNames) + 1]] <- fileName
  title <- sprintf("Diagnostics plot: %s", comparisonSummaries$databaseId[i])
  title <- officer::fpar(officer::ftext(title, prop = titleFormat))
  diagnosticsTitles[[length(diagnosticsTitles) + 1]] <- title

  targetId <- comparisonSummaries$targetId[i]
  targetLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$targetId[i]]
  comparatorId <- comparisonSummaries$comparatorId[i]
  comparatorLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$comparatorId[i]]
  databaseId <- comparisonSummaries$databaseId[i]

  
  for (analysisId in cohortMethodAnalysis$analysisId) { # analysisId = 5
    ps <- getPs(connection, targetId, comparatorId, analysisId = analysisId, databaseId) #Doesn't this need to be outcome specific?  The problem here appears to be that in the shinyDataAll folder the file preference_score_dist_t13191_c13319.rds doesn't exist (model only includes intercept term)
    psPlot <- plotPs2(ps, targetLabel, comparatorLabel)
    psPlots[[length(psPlots) + 1]] <- psPlot
    bal <- getCovariateBalance(connection, targetId, comparatorId, databaseId, analysisId = analysisId, outcomeIds[4])
    balPlot <- plotCovariateBalanceScatterPlot2(bal,beforeLabel = "Before Matching",afterLabel = "After Matching")
    balPlots[[length(balPlots) + 1]] <- balPlot
    
    controlResults <- getControlResults(connection = connection,
                                        targetId = targetId,
                                        comparatorId = comparatorId,
                                        analysisId = analysisId,
                                        databaseId = databaseId)
    controlEstimates <- controlResults[controlResults$effectSize == 1, ]
    calPlot <- plotLargeScatter2(controlEstimates, xLabel = "HR")
    calPlots[[length(calPlots) + 1]] <- calPlot
    }
 # }
  col0 <- grid::textGrob("")
  col1 <- grid::textGrob("PS distribution", gp = grid::gpar(fontsize = 12))
  col2 <- grid::textGrob("Covariate balanace", gp = grid::gpar(fontsize = 12))
  col3 <- grid::textGrob("Null dist.", gp = grid::gpar(fontsize = 12))
 # col4 <- grid::textGrob("Null dist. (on-treatment)", gp = grid::gpar(fontsize = 12))
  row1 <- grid::textGrob(cohortMethodAnalysis$shortDescription[1], rot = 90, gp = grid::gpar(fontsize = 8))
  row2 <- grid::textGrob(cohortMethodAnalysis$shortDescription[2], rot = 90, gp = grid::gpar(fontsize = 8))
  row3 <- grid::textGrob(cohortMethodAnalysis$shortDescription[3], rot = 90, gp = grid::gpar(fontsize = 8))
  row4 <- grid::textGrob(cohortMethodAnalysis$shortDescription[4], rot = 90, gp = grid::gpar(fontsize = 8))
  row5 <- grid::textGrob(cohortMethodAnalysis$shortDescription[5], rot = 90, gp = grid::gpar(fontsize = 8))
  row6 <- grid::textGrob(cohortMethodAnalysis$shortDescription[6], rot = 90, gp = grid::gpar(fontsize = 8))
  row7 <- grid::textGrob(cohortMethodAnalysis$shortDescription[7], rot = 90, gp = grid::gpar(fontsize = 8))
  row8 <- grid::textGrob(cohortMethodAnalysis$shortDescription[8], rot = 90, gp = grid::gpar(fontsize = 8))
  row9 <- grid::textGrob(cohortMethodAnalysis$shortDescription[9], rot = 90, gp = grid::gpar(fontsize = 8))
  row10 <- grid::textGrob(cohortMethodAnalysis$shortDescription[10], rot = 90, gp = grid::gpar(fontsize = 8))
  #row3 <- grid::textGrob(cs$databaseId[3], rot = 90, gp = grid::gpar(fontsize = 12))
  #row4 <- grid::textGrob(cs$databaseId[4], rot = 90, gp = grid::gpar(fontsize = 12))
  #row5 <- grid::textGrob(cs$databaseId[5], rot = 90, gp = grid::gpar(fontsize = 12))
  plotGrid <- gridExtra::grid.arrange(col0, col1, col2, col3, #col4,
                                      row1, psPlots[[1]], balPlots[[1]], calPlots[[1]], 
                                      row2, psPlots[[2]], balPlots[[2]], calPlots[[2]], 
                                      row3, psPlots[[3]], balPlots[[3]], calPlots[[3]], 
                                      row4, psPlots[[4]], balPlots[[4]], calPlots[[4]],
                                      row5, psPlots[[5]], balPlots[[5]], calPlots[[5]],
                                      row6, psPlots[[6]], balPlots[[6]], calPlots[[6]],
                                      row7, psPlots[[7]], balPlots[[7]], calPlots[[7]],
                                      row8, psPlots[[8]], balPlots[[8]], calPlots[[8]],
                                      row9, psPlots[[9]], balPlots[[9]], calPlots[[9]],
                                      row10, psPlots[[10]], balPlots[[10]], calPlots[[10]],
                                      nrow = 11,
                                      heights = c(0.5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
                                      widths = c(0.25, 3.5, 3.5, 3.5))
  ggplot2::ggsave(fileName, plotGrid, width = 8, height = 12, dpi = 400)
}
diagnosticsPlotPairs <- list(diagnosticsTitles, diagnosticsFileNames)
doc <- officer::read_docx() %>%
  officer::body_add_fpar(section3, style = "heading 1")
for(i in 1:length(diagnosticsFileNames)) { #i=1
  doc <- doc %>%
    officer::body_add_fpar(diagnosticsTitles[[i]], style = "heading 2") %>%
    officer::body_add_img(diagnosticsPlotPairs[[2]][[i]], width = 5, height = 7) %>%
    officer::body_add_break()
}
print(doc, target = file.path(reportFolder, "diagnosticPlots.docx"))



# hrs and forest plots ---------------------------------------------------------

section5 <- "Calibrated hazard ratios"
section5 <- officer::fpar(officer::ftext(section5, prop = headingFormat))

hrTitles <- list()
hrOutcomeTitles <- list()
hrFileNames <- list()
stackResults <- list()
outcomeUnblind <- outcomeOfInterest[which(outcomeOfInterest$outcomeId %in% c(13324, 13325, 13321)),]


for (i in 1:nrow(outcomeUnblind)) { # i=1
  mainResults <- getMainResults(targetIds = c(13186),
                                comparatorIds = c(13245),
                                outcomeIds = outcomeUnblind$outcomeId[i],
                                databaseIds = c("Optum"),
                                analysisIds = c(5))
  mainResults <- merge(mainResults, cohortMethodAnalysis[, c("analysisId", "description")])
  mainResults$analysisOrder <- match(mainResults$analysisId, cohortMethodAnalysis$analysisId)
  mainResults$dbOrder <- match(mainResults$databaseId, databaseIds)
  mainResults <- mainResults[order(mainResults$analysisOrder, mainResults$dbOrder), ]
  names(mainResults)[names(mainResults) == "outcomeId"] <- "outcomeId"
  mainResults <- merge(mainResults, outcomeUnblind)
  mainResults$outcomeName <- outcomeLabels$outcomeName[which(outcomeLabels$outcomeId == outcomeUnblind$outcomeId[i])]
  mainResults$outcomeName <- ifelse(mainResults$outcomeId == 13325,"AMI",mainResults$outcomeName)
  stackResults <- rbind(stackResults,mainResults)
}
  forestPlot <- plotReportForest(stackResults, exposureOfInterest$shortName[1], exposureOfInterest$shortName[2])
  
  fileName <- file.path(reportFolder, "hr_plots.png")
  ggplot2::ggsave(fileName, forestPlot, width = 8.5, height = 4, dpi = 400)
  hrFileNames[[length(hrFileNames) + 1]] <- fileName
  title <- "Optum: (keep all) / 1:100 variable ratio matching" #outcomeUnblind$outcomeName[i] # NEED TO FIX THIS LINE AND THEN INCORPORATE OUTCOME TITLES INTO THE OFFICER DOC BELOW
  title <- officer::fpar(officer::ftext(title, prop = titleFormat))
  hrTitles[[length(hrTitles) + 1]] <- title

doc <- officer::read_docx() %>%
  officer::body_add_fpar(section5, style = "heading 1") %>%
 # officer::body_add_fpar(hrTitles[[2]], style = "heading 2") %>%
 # officer::body_add_img(hrFileNames[[2]], width = 5, height = 4) %>%
 # officer::body_add_fpar(hrTitles[[3]], style = "heading 2") %>%
 # officer::body_add_img(hrFileNames[[3]], width = 5, height = 4) %>%
  officer::body_add_fpar(hrTitles[[1]], style = "heading 2") %>%
  officer::body_add_img(hrFileNames[[1]], width = 5, height = 4) %>%
  print(doc, target = file.path(reportFolder, "hrPlots.docx"))




outcomeIds <- c(13324, 13325, 13321)
targetId <- c(13186)
comparatorId <- c(13245)
databaseId <- "Optum"
analysisId <- c(5)

for(i in 1:length(outcomeIds)) { i <- 1
  outcomeId <- outcomeIds[i]
  getKaplanMeier(connection, targetId, comparatorId, outcomeId, databaseId, analysisId)
}




# Primary comparison diagnostics -----------------------------------------------
# 
# primaryFolder <- file.path(reportFolder, "primaryFolder")
# if (!file.exists(primaryFolder)) {
#   dir.create(primaryFolder)
# }
# 
# section5 <- "Primary analysis dianosticis"
# section5 <- officer::fpar(officer::ftext(section5, prop = headingFormat))
# 
# # restrict to primary comparisons for one analysisId
# tcds <- comparisonSummary[comparisonSummary$targetId %in% primaryTarCohortIds & comparisonSummary$comparatorId %in% primaryTarCohortIds, ]
# primaryAnalysisId <- 1
# 
# primaryTitles <- list()
# primaryFileNames <- list()
# 
# for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
#   indication <- indicationGroup$indication
#   indicationCohortIds <- indicationGroup$cohortIds
#   tcs <- unique(tcds[tcds$targetId %in% indicationCohortIds & tcds$comparatorId %in% indicationCohortIds, c("targetId", "comparatorId")])
#   
#   for (i in 1:nrow(tcs)) { # i = 1
#     targetId <- tcs$targetId[i]
#     targetLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$targetId[i]]
#     targetLabel <- capitalize(sub(".*? ", "", targetLabel))
#     comparatorId <- tcs$comparatorId[i]
#     comparatorLabel <- exposureOfInterest$exposureName[exposureOfInterest$exposureId == tcos$comparatorId[i]]
#     comparatorLabel <- capitalize(sub(".*? ", "", comparatorLabel))
#     
#     fileName <- file.path(primaryFolder, sprintf("primary_diagnostics_%s_%s_%s.png", indication, targetLabel, comparatorLabel))
#     fileName <- sub("THR/TKR", "THRTKR", fileName)
#     primaryFileNames[[length(primaryFileNames) + 1]] <- fileName
#     
#     primaryTitle <- sprintf("%s vs %s", targetLabel, comparatorLabel)
#     #primaryTitle <- sprintf("%s: %s vs %s", indication, targetLabel, comparatorLabel)
#     primaryTitle <- officer::fpar(officer::ftext(primaryTitle, prop = titleFormat))
#     primaryTitles[[length(primaryTitles) + 1]] <- primaryTitle 
#     
#     if (!file.exists(fileName)) {
#       psPlots <- list()
#       balPlots <- list()
#       calPlots <- list()
#       
#       databases <- databaseIds[databaseIds != "Meta-analysis"]
#       for (databaseId in databases) { # databaseId <- "CCAE"
#         
#         ps <- getPs(connection, targetId, comparatorId, primaryAnalysisId, databaseId)
#         if (!is.null(ps)) {
#           tcSizes <- getAttrition(connection, targetId, comparatorId, outcomeId, primaryAnalysisId, databaseId)
#           targetSize <- tcSizes$targetPersons[tcSizes$description == "Have at least 1 days at risk"]
#           comparatorSize <- tcSizes$comparatorPersons[tcSizes$description == "Have at least 1 days at risk"]
#           psPlot <- plotReportPs(ps, targetSize = targetSize, comparatorSize = comparatorSize)
#         } else {
#           psPlot <- blankPlot
#         }
#         psPlots[[length(psPlots) + 1]] <- psPlot
#         
#         bal <- getCovariateBalance(connection, targetId, comparatorId, databaseId, primaryAnalysisId, outcomeId)
#         if (!is.null(bal)) {
#           balPlot <- plotReportCovariateBalanceScatterPlot(bal)
#         } else {
#           balPlot <- blankPlot
#         }
#         balPlots[[length(balPlots) + 1]] <- balPlot
#         
#         controlResults <- getControlResults(connection = connection,
#                                             targetId = targetId,
#                                             comparatorId = comparatorId,
#                                             analysisId = primaryAnalysisId,
#                                             databaseId = databaseId)
#         if (nrow(controlResults > 0)) {
#           controlEstimates <- controlResults[controlResults$effectSize == 1, ]
#           calPlot <- plotReportScatter(controlEstimates)
#         } else {
#           calPlot <- blankPlot
#         }
#         calPlots[[length(calPlots) + 1]] <- calPlot
#         
#         # add KM plots?
#       }
#       col0 <- grid::textGrob("")
#       col1 <- grid::textGrob("PS distribution", gp = grid::gpar(fontsize = 16))
#       col2 <- grid::textGrob("Covariate balanace", gp = grid::gpar(fontsize = 16))
#       col3 <- grid::textGrob("Empirical null distribution", gp = grid::gpar(fontsize = 16))
#       row1 <- grid::textGrob(databaseIds[1], rot = 90, gp = grid::gpar(fontsize = 16))
#       row2 <- grid::textGrob(databaseIds[2], rot = 90, gp = grid::gpar(fontsize = 16))
#       row3 <- grid::textGrob(databaseIds[3], rot = 90, gp = grid::gpar(fontsize = 16))
#       row4 <- grid::textGrob(databaseIds[4], rot = 90, gp = grid::gpar(fontsize = 16))
#       plotGrid <- gridExtra::grid.arrange(col0, col1, col2, col3,
#                                           row1, psPlots[[1]], balPlots[[1]], calPlots[[1]],
#                                           row2, psPlots[[2]], balPlots[[2]], calPlots[[2]],
#                                           row3, psPlots[[3]], balPlots[[3]], calPlots[[3]],
#                                           row4, psPlots[[4]], balPlots[[4]], calPlots[[4]],
#                                           nrow = 5,
#                                           heights = c(0.25, 4, 4, 4, 4),
#                                           widths = c(0.25, 4, 4, 4))
#       ggplot2::ggsave(fileName, plotGrid, width = 12.5, height = 12.5, dpi = 400)
#     }
#   }
# }
# primaryPlotPairs <- list(primaryTitles, primaryFileNames)
# doc <- officer::read_docx() %>%
#   officer::body_add_fpar(section5, style = "heading 1")
# for(i in 1:length(primaryFileNames)) { #i=1
#   if (i %in% c(1, 7, 13)) {
#     doc <- doc %>%
#       officer::body_add_fpar(indicationTitles[[i]], style = "heading 2") %>%
#       officer::body_add_fpar(primaryPlotPairs[[1]][[i]], style = "heading 3") %>%
#       officer::body_add_img(primaryPlotPairs[[2]][[i]], width = 6, height = 6) %>%
#       officer::body_add_break()
#   } else {
#     doc <- doc %>%
#       officer::body_add_fpar(primaryPlotPairs[[1]][[i]], style = "heading 3") %>%
#       officer::body_add_img(primaryPlotPairs[[2]][[i]], width = 6, height = 6) %>%
#       officer::body_add_break()
#   }
# }
# print(doc, target = file.path(reportFolder, "primaryDiagnosticPlots.docx"))

# HR forest plots long ---------------------------------------------------------

hrFolder <- file.path(reportFolder, "hrFolder")
if (!file.exists(hrFolder)) {
  dir.create(hrFolder)
}

section6 <- "Calibrated hazard ratios"
section6 <- officer::fpar(officer::ftext(section6, prop = headingFormat))

hrTitles <- list()
hrFileNames <- list()

for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
  indication <- indicationGroup$indication
  indicationCohortIds <- indicationGroup$cohortIds
  indicationComparisons <- comparisonSummary[comparisonSummary$targetId %in% indicationCohortIds & comparisonSummary$comparatorId %in% indicationCohortIds, ]
  
  for (analysisId in c(1, 3, 8, 9, 2, 4)) { # analysisId <- 3
    if (analysisId == 8) { # restrict to primary comparisons for short gap follow-up and ITT follow-up
      fileNameAnalysisId <- 8
      titleSuffix <- " (30-day gap)"
      analysisId <- 1
      tarCohortIds <- sensitivityTarCohortIds
    } else if (analysisId == 9) {
      fileNameAnalysisId <- 9
      titleSuffix <- " (30-day gap)"
      analysisId <- 3
      tarCohortIds <- sensitivityTarCohortIds
    } else {
      fileNameAnalysisId <- analysisId
      titleSuffix <- ""
      tarCohortIds <- primaryTarCohortIds
    }
    indicationAnalysis <- indicationComparisons[indicationComparisons$targetId %in% tarCohortIds & indicationComparisons$comparatorId %in% tarCohortIds, ]
    targetIds <- unique(indicationAnalysis$targetId)
    comparatorIds <- unique(indicationAnalysis$comparatorId)
    
    analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
    analysisName <- paste0(analysisName, titleSuffix)
    fileName <- file.path(hrFolder, sprintf("hr_plots_%s_a%s.png", indication, fileNameAnalysisId))
    fileName <- sub("THR/TKR", "THRTKR", fileName)
    
    mainResults <- getMainResults(targetIds = targetIds, 
                                  comparatorIds = comparatorIds, 
                                  outcomeIds = outcomeId, 
                                  databaseIds = databaseIds,
                                  analysisIds = analysisId)
    mainResults <- merge(mainResults, exposureOfInterest, by.x = "targetId", by.y = "exposureId", all.x = TRUE)
    names(mainResults)[names(mainResults) == "exposureName"] <- "targetName"
    mainResults$targetName <- capitalize(sub(".*? ", "", mainResults$targetName))
    mainResults$targetName <- sub(" \\(30d gap)", "", mainResults$targetName)
    mainResults <- merge(mainResults, exposureOfInterest, by.x = "comparatorId", by.y = "exposureId", all.x = TRUE)
    names(mainResults)[names(mainResults) == "exposureName"] <- "comparatorName"
    mainResults$comparatorName <- capitalize(sub(".*? ", "", mainResults$comparatorName))
    mainResults$comparatorName <- sub(" \\(30d gap)", "", mainResults$comparatorName)
    mainResults$comparison <- paste(mainResults$targetName, mainResults$comparatorName, sep = " vs ")
    mainResults$targetOrder <- match(mainResults$targetId, targetIds)
    mainResults$comparatorOrder <- match(mainResults$comparatorId, comparatorIds)
    mainResults$dbOrder <- match(mainResults$databaseId, databaseIds)
    mainResults <- mainResults[order(mainResults$targetOrder, mainResults$comparatorOrder, mainResults$dbOrder), ]
    mainResults[, c("targetOrder", "comparatorOrder", "dbOrder")] <- NULL
    
    forestPlot <- plotReportMultipleForest(mainResults)
    ggplot2::ggsave(fileName, forestPlot, width = 12.5, height = 9, dpi = 400)
    hrFileNames[[length(hrFileNames) + 1]] <- fileName
    title <- analysisName
    #title <- sprintf("%s: %s", indication, analysisName)
    title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
    hrTitles[[length(hrTitles) + 1]] <- title
  }
}
hrPlotPairs <- list(hrTitles, hrFileNames)
doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section6, style = "heading 1")
for(i in 1:length(hrFileNames)) { #i=1
  if (i %in% c(1, 7, 13)) {
    doc <- doc %>% 
      officer::body_add_fpar(indicationTitles[[i]], style = "heading 2") %>%
      officer::body_add_fpar(hrPlotPairs[[1]][[i]], style = "heading 3") %>%
      officer::body_add_img(hrPlotPairs[[2]][[i]], width = 6, height = 5) %>%
      officer::body_add_break()
  } else {
    doc <- doc %>% 
      officer::body_add_fpar(hrPlotPairs[[1]][[i]], style = "heading 3") %>%
      officer::body_add_img(hrPlotPairs[[2]][[i]], width = 6, height = 5) %>%
      officer::body_add_break()
  }
}
print(doc, target = file.path(reportFolder, "hrPlots.docx"))

source("report/buildReport.R")


# HR forest plots grid ---------------------------------------------------------

hrGridFolder <- file.path(reportFolder, "hrGridFolder")
if (!file.exists(hrGridFolder)) {
  dir.create(hrGridFolder)
}

section8 <- "Calibrated hazard ratios and confidence interval with forest plots comparing rivaroxabab, apixaban, dabigatran, and warfarin"
section8 <- officer::fpar(officer::ftext(section8, prop = headingFormat))

hrGridTitles <- list()
hrGridFileNames <- list()

for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
  indication <- indicationGroup$indication
  indicationCohortIds <- indicationGroup$cohortIds
  indicationComparisons <- comparisonSummary[comparisonSummary$targetId %in% indicationCohortIds & comparisonSummary$comparatorId %in% indicationCohortIds, ]
  
  for (analysisId in c(1, 3, 8, 9, 2, 4)) { # analysisId <- 1
    
    if (analysisId == 8) { # restrict to primary comparisons for short gap follow-up and ITT follow-up
      fileNameAnalysisId <- 8
      titleSuffix <- " (30-day gap)"
      analysisId <- 1
      tarCohortIds <- sensitivityTarCohortIds
    } else if (analysisId == 9) {
      fileNameAnalysisId <- 9
      titleSuffix <- " (30-day gap)"
      analysisId <- 3
      tarCohortIds <- sensitivityTarCohortIds
    } else {
      fileNameAnalysisId <- analysisId
      titleSuffix <- ""
      tarCohortIds <- primaryTarCohortIds
    }
    indicationAnalysis <- indicationComparisons[indicationComparisons$targetId %in% tarCohortIds & indicationComparisons$comparatorId %in% tarCohortIds, ]
    
    analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
    analysisName <- paste0(analysisName, titleSuffix)
    fileName <- file.path(hrGridFolder, sprintf("hr_plots_grid_%s_a%s.png", indication, fileNameAnalysisId))
    fileName <- sub("THR/TKR", "THRTKR", fileName)
    
    hrGridPlots <- list()
    
    tcs <- unique(indicationAnalysis[, c("targetId", "comparatorId")])
    
    for (i in 1:nrow(tcs)) { # i=1
      targetId <- tcs$targetId[i]
      comparatorId <- tcs$comparatorId[i]
      
      mainResults <- getMainResults(targetIds = targetId, 
                                    comparatorIds = comparatorId, 
                                    outcomeIds = outcomeId, 
                                    databaseIds = databaseIds,
                                    analysisIds = analysisId)
      
      mainResults$dbOrder <- match(mainResults$databaseId, databaseIds)
      mainResults <- mainResults[order(mainResults$dbOrder), ]
      mainResults$dbOrder <- NULL
      
      forestPlot <- plotReportSingleForest(mainResults)
      hrGridPlots[[length(hrGridPlots) + 1]] <- forestPlot
    }
    row1 <- grid::textGrob("Rivaroxaban", rot = 90, gp = grid::gpar(fontsize = 16))
    row2 <- grid::textGrob("Apixaban", rot = 90, gp = grid::gpar(fontsize = 16))
    row3 <- grid::textGrob("Dabigatran", rot = 90, gp = grid::gpar(fontsize = 16))
    row4 <- grid::textGrob("")
    col1 <- grid::textGrob("Apixaban", gp = grid::gpar(fontsize = 16))
    col2 <- grid::textGrob("Dabigatran", gp = grid::gpar(fontsize = 16))
    col3 <- grid::textGrob("Warfarin", gp = grid::gpar(fontsize = 16))
    plotGrid <- gridExtra::grid.arrange(row1, hrGridPlots[[1]], hrGridPlots[[2]], hrGridPlots[[3]],
                                        row2, blankPlot,        hrGridPlots[[4]], hrGridPlots[[5]],
                                        row3, blankPlot,        blankPlot,        hrGridPlots[[6]],
                                        row4, col1, col2, col3,
                                        nrow = 4,
                                        heights = c(4, 4, 4, 0.25),
                                        widths = c(0.25, 4, 4, 4))
    ggplot2::ggsave(fileName, plotGrid, width = 12.5, height = 12.5, dpi = 400)
    
    
    
    
    calFileNames[[length(calFileNames) + 1]] <- fileName
    title <- sprintf("%s: %s, %s", indication, databaseId, analysisName)
    title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
    calTitles[[length(calTitles) + 1]] <- title
  }
}


# Balance plots ----------------------------------------------------------------

balFolder <- file.path(reportFolder, "balFolder")
if (!file.exists(balFolder)) {
  dir.create(balFolder)
}

section6 <- "Covariate balance sensitivity analysis"
section6 <- officer::fpar(officer::ftext(section6, prop = headingFormat))

# restrict to primary comparisons
tcds <- comparisonSummary[comparisonSummary$targetId %in% primaryTarCohortIds & comparisonSummary$comparatorId %in% primaryTarCohortIds, ]

balTitles <- list()
balFileNames <- list()

for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
  indication <- indicationGroup$indication
  indicationCohortIds <- indicationGroup$cohortIds
  
  for (databaseId in databaseIds[databaseIds != "Meta-analysis"]) { # databaseId <- "CCAE"
    
    for (analysisId in c(5)) { # analysisId <- 1
      analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
      analysisName <- sub("On-treatment, ", "", analysisName)
      fileName <- file.path(balFolder, sprintf("bal_plots_%s_%s_a%s.png", indication, databaseId, analysisId))
      fileName <- sub("THR/TKR", "THRTKR", fileName)
      
      if (!file.exists(fileName)) {
        tcs <- tcds[tcds$databaseId == databaseId & tcds$targetId %in% indicationCohortIds & tcds$comparatorId %in% indicationCohortIds, ]
        balPlots <- list()
        
        for (i in 1:nrow(tcs)) { # i=1
          targetId <- tcs$targetId[i]
          comparatorId <- tcs$comparatorId[i]
          bal <- getCovariateBalance(connection, targetId, comparatorId, databaseId, analysisId, outcomeId)
          if (!is.null(bal)) {
            balPlot <- plotReportCovariateBalanceScatterPlot(bal)
          } else {
            balPlot <- blankPlot
          }
          balPlots[[length(balPlots) + 1]] <- balPlot
        }
        row1 <- grid::textGrob("Rivaroxaban", rot = 90, gp = grid::gpar(fontsize = 16))
        row2 <- grid::textGrob("Apixaban", rot = 90, gp = grid::gpar(fontsize = 16))
        row3 <- grid::textGrob("Dabigatran", rot = 90, gp = grid::gpar(fontsize = 16))
        row4 <- grid::textGrob("")
        col1 <- grid::textGrob("Apixaban", gp = grid::gpar(fontsize = 16))
        col2 <- grid::textGrob("Dabigatran", gp = grid::gpar(fontsize = 16))
        col3 <- grid::textGrob("Warfarin", gp = grid::gpar(fontsize = 16))
        plotGrid <- gridExtra::grid.arrange(row1, balPlots[[1]], balPlots[[2]], balPlots[[3]],
                                            row2, blankPlot,     balPlots[[4]], balPlots[[5]],
                                            row3, blankPlot,     blankPlot,     balPlots[[6]],
                                            row4, col1, col2, col3,
                                            nrow = 4,
                                            heights = c(4, 4, 4, 0.25),
                                            widths = c(0.25, 4, 4, 4))
        ggplot2::ggsave(fileName, plotGrid, width = 12.5, height = 12.5, dpi = 400)
      }
      balFileNames[[length(balFileNames) + 1]] <- fileName
      title <- sprintf("%s: %s, %s", indication, databaseId, analysisName)
      title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
      balTitles[[length(balTitles) + 1]] <- title
    }
  }
}
balPlotPairs <- list(balTitles, balFileNames)
doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section6, style = "heading 1")
for(i in 1:length(balFileNames)) { #i=1
  doc <- doc %>% 
    officer::body_add_fpar(balPlotPairs[[1]][[i]], style = "heading 2") %>%
    officer::body_add_img(balPlotPairs[[2]][[i]], width = 6, height = 6) %>%
    officer::body_add_break()
}
print(doc, target = file.path(reportFolder, "balPlots.docx"))

# Calibration plots ------------------------------------------------------------

calFolder <- file.path(reportFolder, "calFolder")
if (!file.exists(calFolder)) {
  dir.create(calFolder)
}

section6 <- "Empirical null distributions from negative control outomces comparing rivaroxabab, apixaban, dabigatran, and warfarin"
section6 <- officer::fpar(officer::ftext(section6, prop = headingFormat))

calTitles <- list()
calFileNames <- list()

for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
  indication <- indicationGroup$indication
  indicationCohortIds <- indicationGroup$cohortIds
  
  for (databaseId in databaseIds[databaseIds != "Meta-analysis"]) { # databaseId <- "CCAE"
    
    for (analysisId in c(3, 8, 9, 2, 4)) { # analysisId <- 3
      
      # restrict to primary comparisons for short gap follow-up and ITT follow-up
      if (analysisId == 8) {
        fileNameAnalysisId <- 8
        titleSuffix <- " (30-day gap)"
        analysisId <- 1
        tarCohortIds <- sensitivityTarCohortIds
      } else if (analysisId == 9) {
        fileNameAnalysisId <- 9
        titleSuffix <- " (30-day gap)"
        analysisId <- 3
        tarCohortIds <- sensitivityTarCohortIds
      } else {
        fileNameAnalysisId <- analysisId
        titleSuffix <- ""
        tarCohortIds <- primaryTarCohortIds
      }
      tcds <- comparisonSummary[comparisonSummary$targetId %in% tarCohortIds & comparisonSummary$comparatorId %in% tarCohortIds, ]
      
      analysisName <- cohortMethodAnalysis$description[cohortMethodAnalysis$analysisId == analysisId]
      analysisName <- paste0(analysisName, titleSuffix)
      fileName <- file.path(calFolder, sprintf("cal_plots_%s_%s_a%s.png", indication, databaseId, fileNameAnalysisId))
      fileName <- sub("THR/TKR", "THRTKR", fileName)
      
      if (!file.exists(fileName)) {
        tcs <- tcds[tcds$databaseId == databaseId & tcds$targetId %in% indicationCohortIds & tcds$comparatorId %in% indicationCohortIds, ]
        calPlots <- list()
        
        for (i in 1:nrow(tcs)) { # i=1
          targetId <- tcs$targetId[i]
          comparatorId <- tcs$comparatorId[i]
          controlResults <- getControlResults(connection = connection,
                                              targetId = targetId,
                                              comparatorId = comparatorId,
                                              analysisId = analysisId,
                                              databaseId = databaseId)
          
          if (nrow(controlResults > 0)) {
            controlEstimates <- controlResults[controlResults$effectSize == 1, ]
            calPlot <- plotReportScatter(controlEstimates)
          } else {
            calPlot <- blankPlot
          }
          calPlots[[length(calPlots) + 1]] <- calPlot
        }
        row1 <- grid::textGrob("Rivaroxaban", rot = 90, gp = grid::gpar(fontsize = 16))
        row2 <- grid::textGrob("Apixaban", rot = 90, gp = grid::gpar(fontsize = 16))
        row3 <- grid::textGrob("Dabigatran", rot = 90, gp = grid::gpar(fontsize = 16))
        row4 <- grid::textGrob("")
        col1 <- grid::textGrob("Apixaban", gp = grid::gpar(fontsize = 16))
        col2 <- grid::textGrob("Dabigatran", gp = grid::gpar(fontsize = 16))
        col3 <- grid::textGrob("Warfarin", gp = grid::gpar(fontsize = 16))
        plotGrid <- gridExtra::grid.arrange(row1, calPlots[[1]], calPlots[[2]], calPlots[[3]],
                                            row2, blankPlot,     calPlots[[4]], calPlots[[5]],
                                            row3, blankPlot,     blankPlot,     calPlots[[6]],
                                            row4, col1, col2, col3,
                                            nrow = 4,
                                            heights = c(4, 4, 4, 0.25),
                                            widths = c(0.25, 4, 4, 4))
        ggplot2::ggsave(fileName, plotGrid, width = 12.5, height = 12.5, dpi = 400)
      }
      calFileNames[[length(calFileNames) + 1]] <- fileName
      title <- sprintf("%s: %s, %s", indication, databaseId, analysisName)
      title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
      calTitles[[length(calTitles) + 1]] <- title
    }
  }
}
calPlotPairs <- list(calTitles, calFileNames)
doc <- officer::read_docx() %>% 
  officer::body_add_fpar(section6, style = "heading 1")
for(i in 1:length(calFileNames)) { #i=1
  doc <- doc %>% 
    officer::body_add_fpar(calPlotPairs[[1]][[i]], style = "heading 2") %>%
    officer::body_add_img(calPlotPairs[[2]][[i]], width = 6, height = 6) %>%
    officer::body_add_break()
}
print(doc, target = file.path(reportFolder, "calPlots.docx"))


# KM plots (TODO??) ---------------------------------------------------------------------
# PS plots (depreciated) -------------------------------------------------------

# psFolder <- file.path(reportFolder, "psFolder")
# if (!file.exists(psFolder)) {
#   dir.create(psFolder)
# }
# 
# section4 <- "Preference score distributions comparing rivaroxabab, apixaban, dabigatran, and warfarin"
# section4 <- officer::fpar(officer::ftext(section4, prop = headingFormat))
# 
# # restrict to primary comparisons for one analysisId
# tcds <- comparisonSummary[comparisonSummary$targetId %in% primaryTarCohortIds & comparisonSummary$comparatorId %in% primaryTarCohortIds, ]
# defaultAnalysisId <- 1
# 
# psTitles <- list()
# psFileNames <- list()
# 
# for (indicationGroup in indicationGroups) { # indicationGroup <- indicationGroups[[1]]
#   indication <- indicationGroup$indication
#   indicationCohortIds <- indicationGroup$cohortIds
# 
#   for (databaseId in databaseIds[databaseIds != "Meta-analysis"]) { # databaseId <- "CCAE"
#     fileName <- file.path(psFolder, sprintf("ps_plots_%s_%s.png", indication, databaseId))
#     fileName <- sub("THR/TKR", "THRTKR", fileName)
#     
#     if (!file.exists(fileName)) {
#       tcs <- tcds[tcds$databaseId == databaseId & tcds$targetId %in% indicationCohortIds & tcds$comparatorId %in% indicationCohortIds, ]
#       psPlots <- list()
#       
#       for (i in 1:nrow(tcs)) { # i=1
#         targetId <- tcs$targetId[i]
#         comparatorId <- tcs$comparatorId[i]
#         ps <- getPs(connection, targetId, comparatorId, defaultAnalysisId, databaseId)
#         if (!is.null(ps)) {
#           tcSizes <- getAttrition(connection, targetId, comparatorId, outcomeId, defaultAnalysisId, databaseId)
#           targetSize <- tcSizes$targetPersons[tcSizes$description == "Have at least 1 days at risk"] 
#           comparatorSize <- tcSizes$comparatorPersons[tcSizes$description == "Have at least 1 days at risk"]
#           psPlot <- plotReportPs(ps, targetSize = targetSize, comparatorSize = comparatorSize)
#         } else {
#           psPlot <- blankPlot
#         }
#         psPlots[[length(psPlots) + 1]] <- psPlot
#       }
#       row1 <- grid::textGrob("Rivaroxaban", rot = 90, gp = grid::gpar(fontsize = 16))
#       row2 <- grid::textGrob("Apixaban", rot = 90, gp = grid::gpar(fontsize = 16))
#       row3 <- grid::textGrob("Dabigatran", rot = 90, gp = grid::gpar(fontsize = 16))
#       row4 <- grid::textGrob("")
#       col1 <- grid::textGrob("Apixaban", gp = grid::gpar(fontsize = 16))
#       col2 <- grid::textGrob("Dabigatran", gp = grid::gpar(fontsize = 16))
#       col3 <- grid::textGrob("Warfarin", gp = grid::gpar(fontsize = 16))
#       plotGrid <- gridExtra::grid.arrange(row1, psPlots[[1]], psPlots[[2]], psPlots[[3]],
#                                           row2, blankPlot,    psPlots[[4]], psPlots[[5]],
#                                           row3, blankPlot,    blankPlot,    psPlots[[6]],
#                                           row4, col1, col2, col3,
#                                           nrow = 4,
#                                           heights = c(4, 4, 4, 0.25),
#                                           widths = c(0.25, 4, 4, 4))
#       ggplot2::ggsave(fileName, plotGrid, width = 12.5, height = 12.5, dpi = 400)
#     }
#     psFileNames[[length(psFileNames) + 1]] <- fileName
#     title <- sprintf("%s: %s", indication, databaseId)
#     title <- officer::fpar(officer::ftext(title, prop = titleFormat))  
#     psTitles[[length(psTitles) + 1]] <- title
#   }
# }
# psPlotPairs <- list(psTitles, psFileNames)
# doc <- officer::read_docx() %>% 
#   officer::body_add_fpar(section4, style = "heading 1")
# for(i in 1:length(psFileNames)) { #i=1
#   doc <- doc %>% 
#     officer::body_add_fpar(psPlotPairs[[1]][[i]], style = "heading 2") %>%
#     officer::body_add_img(psPlotPairs[[2]][[i]], width = 6, height = 6) %>%
#     officer::body_add_break()
# }
# print(doc, target = file.path(reportFolder, "psPlots.docx"))