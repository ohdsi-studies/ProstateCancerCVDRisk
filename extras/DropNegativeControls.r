# This program drops concepts from the set of negative controls (based on feedback from Laura Hester in PREP review)

# Original set of negative control concepts
    # Concept Set #9821 - [EPI_720] Negative Controls - Curated (contrast: abiraterone v. enzalutamide) - Deduped
# Set of concepts to drop from original negative control set
    # Concept Set #11020 - [EPI_720] Negative Controls - Curated (contrast: abiraterone v. enzalutamide) - Deduped (ROUND 2 VARS TO DROP 4/28/2020)
# New set of negative control concepts
    # Concept Set #11011 - [EPI_720] Negative Controls - Curated (contrast: abiraterone v. enzalutamide) - Deduped (ROUND 2 REDUCED 4/28/2020)

# negative controls csv
pathToCsv <- system.file("settings", "NegativeControls.csv", package = "EPI720")
negativeControls <- read.csv(pathToCsv)
write.csv(negativeControls, "inst/settings/NegativeControls_original.csv", row.names = FALSE)
drops <- negativeControls$outcomeId %in% c(321314,4305080,4009585,43531681,199866,440930,321318,312934,43020432,443344,4082311,198803,4127568,4169095,4216972,380094,4110192,432738,4201390,4079936,42872402,40481919,80816,435796,196523,433516,433440,30753,318566,194408,4214376,434610,436375,437390,4308509,76161,441417,4208390,4168222,4115991,197676,440577,4166126,317002,439926,442012,199075,45766207,321596,316429,4084229,436926,261236,75488,4068482,440360,4142875,4059290,199065,4275423,444070,441051,195590,200845)
negativeControls <- negativeControls[!drops, ]
write.csv(negativeControls, "inst/settings/NegativeControls.csv", row.names = FALSE)

# rebuild package

# all controls from each database
studyFolder <- "D:/StudyResults/EPI720_3"
timeAtRisks <- c("AsTreated", "IntentToTreat") 
outputFolders <- c(file.path(studyFolder, "MDCR"),
                   file.path(studyFolder, "Optum"))

for (outputFolder in outputFolders) { 
  for (timeAtRisk in timeAtRisks) { 
    pathToCsv <- file.path(outputFolder, timeAtRisk, "AllControls.csv")
    allControls <- read.csv(pathToCsv)
    write.csv(allControls, file.path(outputFolder, timeAtRisk, "AllControls_original.csv"), row.names = FALSE)
    drops <- allControls$oldOutcomeId %in% c(321314,4305080,4009585,43531681,199866,440930,321318,312934,43020432,443344,4082311,198803,4127568,4169095,4216972,380094,4110192,432738,4201390,4079936,42872402,40481919,80816,435796,196523,433516,433440,30753,318566,194408,4214376,434610,436375,437390,4308509,76161,441417,4208390,4168222,4115991,197676,440577,4166126,317002,439926,442012,199075,45766207,321596,316429,4084229,436926,261236,75488,4068482,440360,4142875,4059290,199065,4275423,444070,441051,195590,200845)
    allControls <- allControls[!drops, ]
    write.csv(allControls, pathToCsv, row.names = FALSE)
  }
}
# rerun export function for all databases for re-calibrate with 4 negative controls removed
# rerun prepare for shiny function