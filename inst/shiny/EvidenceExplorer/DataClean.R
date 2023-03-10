database$order <- match(database$databaseId, c("CCAE", "MDCR", "Optum"))
database <- database[order(database$order), ]
database$order <- NULL

#cohortMethodResult$i2 <- round(cohortMethodResult$i2, 2)

#cadPadExposureIds <- c(11408, 11409, 11410, 11411, 13207, 13208, 13209, 13236)

#drops <- tcos$targetId %in% cadPadExposureIds | tcos$comparatorId %in% cadPadExposureIds
#tcos <- tcos[!drops, ]

#drops <- cohortMethodResult$targetId %in% cadPadExposureIds | cohortMethodResult$comparatorId %in% cadPadExposureIds # drop CAD/PAD rows
#cohortMethodResult <- cohortMethodResult[!drops, ]

#drops <- attrition$targetId %in% cadPadExposureIds | attrition$comparatorId %in% cadPadExposureIds
#attrition <- attrition[!drops, ]

#exposureOfInterest$definition <- NULL
#exposureOfInterest <- exposureOfInterest[!duplicated(exposureOfInterest), ]
#drops <- exposureOfInterest$exposureId %in% cadPadExposureIds # drop CAD/PAD exposures
#exposureOfInterest <- exposureOfInterest[!drops, ]

#exposureOfInterest$exposureName <- gsub("\\[680\\] Female ", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(" new users with prior", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName[grep("non-valvular atrial fibrillation", exposureOfInterest$exposureName)] <- paste("NVAF:", exposureOfInterest$exposureName[grep("non-valvular atrial fibrillation", exposureOfInterest$exposureName)])
#exposureOfInterest$exposureName[grep("venous thromboembolism", exposureOfInterest$exposureName)] <- paste("VTE:", exposureOfInterest$exposureName[grep("venous thromboembolism", exposureOfInterest$exposureName)])
#exposureOfInterest$exposureName[grep("knee or hip replacement surgery", exposureOfInterest$exposureName)] <- paste("TKR/THR:", exposureOfInterest$exposureName[grep("knee or hip replacement surgery", exposureOfInterest$exposureName)])
#exposureOfInterest$exposureName <- gsub(" non-valvular atrial fibrillation", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(" venous thromboembolism", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(" coronary artery disease or peripheral artery disease", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(" knee or hip replacement surgery", "", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(", sens TAR1", " (30d gap)", exposureOfInterest$exposureName)
#exposureOfInterest$exposureName <- gsub(", sen TAR1", " (30d gap)", exposureOfInterest$exposureName)
#exposureOrder <- c("NVAF: rivaroxaban",
#                   "NVAF: apixaban",
#                   "NVAF: dabigatran",
#                   "NVAF: warfarin",
#                   "VTE: rivaroxaban",
#                   "VTE: apixaban",
#                   "VTE: dabigatran",
#                   "VTE: warfarin",
#                   "TKR/THR: rivaroxaban",
#                   "TKR/THR: apixaban",
#                   "TKR/THR: dabigatran",
#                   "TKR/THR: warfarin",
#                   "NVAF: rivaroxaban (30d gap)",
#                   "NVAF: apixaban (30d gap)",
#                   "NVAF: dabigatran (30d gap)",
#                   "NVAF: warfarin (30d gap)",
#                   "VTE: rivaroxaban (30d gap)",
#                   "VTE: apixaban (30d gap)",
#                   "VTE: dabigatran (30d gap)",
#                   "VTE: warfarin (30d gap)",
#                   "TKR/THR: rivaroxaban (30d gap)",
#                   "TKR/THR: apixaban (30d gap)",
#                   "TKR/THR: dabigatran (30d gap)",
#                   "TKR/THR: warfarin (30d gap)")
#exposureOfInterest$order <- match(exposureOfInterest$exposureName, exposureOrder)
#exposureOfInterest <- exposureOfInterest[order(exposureOfInterest$order), ]
#exposureOfInterest$order <- NULL

#outcomeOfInterest$definition <- NULL
#outcomeOfInterest <- outcomeOfInterest[!duplicated(outcomeOfInterest), ]
#outcomeOfInterest$outcomeName <- "Severe uterine bleed"

#cohortMethodAnalysis$definition <- NULL
#cohortMethodAnalysis <- cohortMethodAnalysis[!duplicated(cohortMethodAnalysis), ]
#cohortMethodAnalysis$description[cohortMethodAnalysis$description == "1:1 ratio PS matching on-treatment unconditional Cox PH"] <- "On-treatment, 1:1 PS match"
#cohortMethodAnalysis$description[cohortMethodAnalysis$description == "1:100 variable ratio PS matching on-treatment conditional Cox PH"] <- "On-treatment, 1:100 PS match"
#cohortMethodAnalysis$description[cohortMethodAnalysis$description == "1:1 ratio PS matching intent-to-treat unconditional Cox PH"] <- "Intent-to-treat, 1:1 PS match"
#cohortMethodAnalysis$description[cohortMethodAnalysis$description == "1:100 variable ratio PS matching intent-to-treat conditional Cox PH"] <- "Intent-to-treat, 1:100 PS match"

#primaryTarCohortIds <- c(11400, 11403, 11401, 11402, 11404, 11407, 11405, 11406, 11412, 11415, 11413, 11414)
#sensitivityTarCohortIds <- c(13033, 13234, 13199, 13200, 13203, 13235, 13204, 13205, 13211, 13237, 13212, 13213)

#primaryTar <- cohortMethodResult$targetId %in% primaryTarCohortIds & cohortMethodResult$comparatorId %in% primaryTarCohortIds & cohortMethodResult$analysisId %in% c(1, 3)
#sensitivityTar <- cohortMethodResult$targetId %in% sensitivityTarCohortIds & cohortMethodResult$comparatorId %in% sensitivityTarCohortIds & cohortMethodResult$analysisId %in% c(1, 3)
#ittTar <- cohortMethodResult$targetId %in% primaryTarCohortIds & cohortMethodResult$comparatorId %in% primaryTarCohortIds & cohortMethodResult$analysisId %in% c(2, 4)
#cohortMethodResult <- cohortMethodResult[primaryTar | sensitivityTar | ittTar, ] # this drops 1:1 SENS ITT amd 1:100 SENS ITT (since the same as 1:1 ITT and 1:100 ITT)



UnblindOutcomes <- cohortMethodResult$outcomeId[!(cohortMethodResult$outcomeId %in% c(13320,13323,13326,13327))]

#cohortMethodResult$outcomeId %in% ()
#outcomeOfInterest$outcomeId

#subset(D1, !(V1 %in% c('B','N','T')))

#13320 13321 13323 13324 13325 13326 13327




# This line only selects the AMI, ischemic stroke, and heart failure analyses for the un-stratified cohorts in the OptumDOD database (as-treated analysis using the KEEP-ALL approach and 1:100 variable ratio matching)
keeps <- (cohortMethodResult$databaseId %in% c("Optum") & cohortMethodResult$outcomeId %in% UnblindOutcomes & cohortMethodResult$targetId==13186
          & cohortMethodResult$comparatorId==13245 & cohortMethodResult$analysisId==5)

#test<- cohortMethodResult[keeps,]

#test[test$outcomeId==13324,]
#cohortMethodResult[cohortMethodResult$outcomeId==13324 & cohortMethodResult$analysisId==5,]

#cohortMethodResult[cohortMethodResult$outcomeId %in% c(13186,13245,13320,13323,13326,13326,13327) | cohortMethodResult$databaseId %in% c("CCAE","MDCR") | cohortMethodResult$analysisId != 5 | cohortMethodResult$targetId != 13248 | cohortMethodResult$comparatorId != 13247,]

  #(
  #    (cohortMethodResult$targetId %in% c(11400, 13033) & cohortMethodResult$comparatorId %in% c(11403, 13234) & cohortMethodResult$databaseId == "MDCD") |  # NVAF: rivaroxaban vs NVAF: warfarin
  #    (cohortMethodResult$targetId %in% c(11402, 13200) & cohortMethodResult$comparatorId %in% c(11403, 13234) & cohortMethodResult$databaseId == "MDCD") |  # NVAF: dabigatran  vs NVAF: warfarin
  #    (cohortMethodResult$targetId %in% c(11400, 13033) & cohortMethodResult$comparatorId %in% c(11402, 13200) & cohortMethodResult$databaseId == "MDCD") |  # NVAF: rivaroxaban vs NVAF: dabigatran
  #    (cohortMethodResult$targetId %in% c(11401, 13199) & cohortMethodResult$comparatorId %in% c(11402, 13200) & cohortMethodResult$databaseId == "MDCD") |  # NVAF: apixaban    vs NVAF: dabigatran
  #    
  #    (cohortMethodResult$targetId %in% c(11405, 13204) & cohortMethodResult$comparatorId %in% c(11407, 13235) & cohortMethodResult$databaseId == "MDCR") |                                # VTE: apixaban    vs VTE: warfarin
  #    (cohortMethodResult$targetId %in% c(11406, 13205) & cohortMethodResult$comparatorId %in% c(11407, 13235) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum")) |  # VTE: dabigatran  vs VTE: warfarin 
  #    (cohortMethodResult$targetId %in% c(11404, 13203) & cohortMethodResult$comparatorId %in% c(11406, 13205) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum")) |  # VTE: rivaroxaban vs VTE: dabigatran
  #    (cohortMethodResult$targetId %in% c(11405, 13204) & cohortMethodResult$comparatorId %in% c(11406, 13205) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum")) |  # VTE: apixaban    vs VTE: dabigatran
  #    
  #    (cohortMethodResult$targetId %in% c(11412, 13211) & cohortMethodResult$comparatorId %in% c(11415, 13212) & cohortMethodResult$databaseId == "MDCD") |                                # TKR/THR: rivaroxaban vs TKR/THR: warfarin
  #    (cohortMethodResult$targetId %in% c(11413, 13212) & cohortMethodResult$comparatorId %in% c(11415, 13237) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR")) |           # TKR/THR: apixaban    vs TKR/THR: warfarin
  #    (cohortMethodResult$targetId %in% c(11414, 13213) & cohortMethodResult$comparatorId %in% c(11415, 13237) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum")) |  # TKR/THR: dabigatran  vs TKR/THR: warfarin
  #    (cohortMethodResult$targetId %in% c(11412, 13211) & cohortMethodResult$comparatorId %in% c(11413, 13212) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR")) |           # TKR/THR: rivaroxaban vs TKR/THR: apixaban
  #    (cohortMethodResult$targetId %in% c(11412, 13211) & cohortMethodResult$comparatorId %in% c(11414, 13213) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum")) |  # TKR/THR: rivaroxaban vs TKR/THR: dabigatran
  #    (cohortMethodResult$targetId %in% c(11413, 13212) & cohortMethodResult$comparatorId %in% c(11414, 13213) & cohortMethodResult$databaseId %in% c("CCAE", "MDCD", "MDCR", "Optum"))    # TKR/THR: apixaban    vs TKR/THR: dabigatran
    
  #) & cohortMethodResult$outcomeId == outcomeOfInterest$outcomeId # don't blind negative control estimates

cohortMethodResult$rr[!keeps] <- NA
cohortMethodResult$ci95Ub[!keeps] <- NA
cohortMethodResult$ci95Lb[!keeps] <- NA
cohortMethodResult$logRr[!keeps] <- NA
cohortMethodResult$seLogRr[!keeps] <- NA
cohortMethodResult$p[!keeps] <- NA
cohortMethodResult$calibratedRr[!keeps] <- NA
cohortMethodResult$calibratedCi95Ub[!keeps] <- NA
cohortMethodResult$calibratedCi95Lb[!keeps] <- NA
cohortMethodResult$calibratedLogRr[!keeps] <- NA
cohortMethodResult$calibratedSeLogRr[!keeps] <- NA
cohortMethodResult$calibratedP[!keeps] <- NA
