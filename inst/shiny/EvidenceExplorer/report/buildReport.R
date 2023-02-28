start <- Sys.time()

# create report document

doc <- officer::read_docx()

# section 1: study population counts

doc <- doc %>% 
  officer::body_add_fpar(section1, style = "heading 1") %>%
  officer::body_add_fpar(exposureCountTitle, style = "heading 2") %>%
  flextable::body_add_flextable(countsByExposure) %>%
  officer::body_add_break() %>%
  officer::body_add_fpar(databaseCountTitle, style = "heading 2") %>%
  flextable::body_add_flextable(countsByDatabase) %>%
  officer::body_add_break()

# section 2: baseline characteristics

doc <- doc %>% 
  officer::body_add_fpar(section2, style = "heading 1") 
for(i in 1:length(table1s)) { #i=1
  if (i %in% c(1, 7, 13)) {
    doc <- doc %>% 
      officer::body_add_fpar(indicationTitles[[i]], style = "heading 2") %>%
      officer::body_add_fpar(table1Pairs[[1]][[i]], style = "heading 3") %>%
      flextable::body_add_flextable(table1Pairs[[2]][[i]]) %>%
      officer::body_add_break()
  } else {
    doc <- doc %>% 
      officer::body_add_fpar(table1Pairs[[1]][[i]], style = "heading 3") %>%
      flextable::body_add_flextable(table1Pairs[[2]][[i]]) %>%
      officer::body_add_break() 
  }
}

# section 3: crude IRs

doc <- doc %>% 
  officer::body_add_fpar(section3, style = "heading 1") %>%
  flextable::body_add_flextable(crudeIrTable) %>%
  officer::body_add_break() 

# section 4: after matching IRs

doc <- doc %>% 
  officer::body_add_fpar(section4, style = "heading 1")
for(i in 1:length(eventTables)) { #i=1
  doc <- doc %>% 
    officer::body_add_fpar(eventTablePairs[[1]][[i]], style = "heading 2") %>%
    flextable::body_add_flextable(eventTablePairs[[2]][[i]]) %>%
    officer::body_add_break()
}

# section 5: primary analysis diagnostics

doc <- doc %>% 
  officer::body_add_fpar(section5, style = "heading 1")
for(i in 1:length(primaryFileNames)) { #i=1
  if (i %in% c(1, 7, 13)) {
    doc <- doc %>%
      officer::body_add_fpar(indicationTitles[[i]], style = "heading 2") %>%
      officer::body_add_fpar(primaryPlotPairs[[1]][[i]], style = "heading 3") %>%
      officer::body_add_img(primaryPlotPairs[[2]][[i]], width = 6, height = 6) %>%
      officer::body_add_break()
  } else {
    doc <- doc %>%
      officer::body_add_fpar(primaryPlotPairs[[1]][[i]], style = "heading 3") %>%
      officer::body_add_img(primaryPlotPairs[[2]][[i]], width = 6, height = 6) %>%
      officer::body_add_break()
  }
}

# section 6: HR plots (need reordering/relabeling)

doc <- doc %>% 
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


# print full report

print(doc, target = file.path(reportFolder, "report.docx"))

delta <- Sys.time() - start
paste("Creating document took", signif(delta, 3), attr(delta, "units"))