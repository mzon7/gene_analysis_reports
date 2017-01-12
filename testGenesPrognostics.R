
testOvCancSurv = function(genesEntrez, geneSigCoef, cancerType, bestProbes = NULL, ovDataList = NULL, threeGroups = FALSE, infoString = NULL)
{
  #print(genesEntrez)
  #assign("geneTracker",c(), envir = .GlobalEnv)
  #geneTracker = c(geneTracker, genesEntrez)
  
  #sum data sets dont work with the subtyping functions, keep track of which data not to use
  cancerTypeOrig = cancerType
  forestPlotFrame = as.data.frame(NULL)
  idString = c()
  #infoString = c()
  sampleString = c()
  broken = c()
  #print(length(ovDataList))
  
  if(cancerType == "ovOs" || cancerType =="ovRec")
    broken = c()

  dataEvTime = list()
  dataEvStatus = list()
  #get valid ovarian esets sets from metaGXOvarian if they are not provided
  geneSig = genesEntrez
  if(length(ovDataList) == 0)
  {
    ovDataList = list()
    
    #setwd("D:\\PMH Data\\Breast Cancer Data\\Breast Cancer Data Haibe-Kains\\wetransfer-474706")
    #install.packages("MetaGxBreast_0.99.0.tar.gz", type ="source", repos = NULL)
    if(cancerType == "ovOs" || cancerType == "ovRec")
    {
      if(cancerType == "ovOs")
      {
      library(MetaGxOvarian)
      source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
      #is below line out okay?
      #source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
      #manually determined which data sets had necessary survival info, should automate this process in the
      #future as more data sets may be added
      #should adapt this to work like ovRec section, then just need a not manual way to load data sets into ovDataList and infoString
      ovDataList[[1]] = TCGA
      ovDataList[[2]] = E.MTAB.386
      ovDataList[[3]] = GSE13876
      ovDataList[[4]] = GSE14764
      ovDataList[[5]] = GSE17260
      ovDataList[[6]] = GSE18520
      #GSE19829 does not cooperate with subtyping so remove
      #ovDataList[[7]] = GSE19829
      ovDataList[[7]] = GSE26193
      ovDataList[[8]] = GSE26712
      #GSE30009 does not cooperate with subtyping so remove
      #ovDataList[[10]] = GSE30009
      ovDataList[[9]] = GSE30161
      ovDataList[[10]] = GSE32062
      ovDataList[[11]] = GSE32063
      ovDataList[[12]] = GSE49997
      ovDataList[[13]] = GSE51088
      ovDataList[[14]] = GSE8842
      ovDataList[[15]] = GSE9891
      ovDataList[[16]] = PMID17290060
      ovDataList[[17]] = PMID19318476
      infoString = c("TCGA", "E.MTAB.386" ,"GSE13876", "GSE14764", "GSE17260", "GSE18520", 
                     "GSE26193", "GSE26712", "GSE30161", "GSE32062", "GSE32063"
                     , "GSE49997", "GSE51088", "GSE8842", "GSE9891", "PMID17290060", "PMID19318476")
      
      }
      
      
    if(cancerType == "ovRec"){
        library(MetaGxOvarian)
        source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
        #is below line out okay?
        #source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
        #manually determined which data sets had necessary survival info, should automate this process in the
        #future as more data sets may be added
        ovDataList[[1]] = E.MTAB.386
        ovDataList[[2]] = GSE2109
        ovDataList[[3]] = GSE6008
        ovDataList[[4]] = GSE6822
        ovDataList[[5]] = GSE8842
        ovDataList[[6]] = GSE9891
        ovDataList[[7]] = GSE12418
        ovDataList[[8]] = GSE12470
        ovDataList[[9]] = GSE13876
        ovDataList[[10]] = GSE14764
        ovDataList[[11]] = GSE17260
        ovDataList[[12]] = GSE18520
        ovDataList[[13]] = GSE19829
        ovDataList[[14]] = GSE20565
        ovDataList[[15]] = GSE26193
        ovDataList[[16]] = GSE26712
        ovDataList[[17]] = GSE30009
        ovDataList[[18]] = GSE30161
        ovDataList[[19]] = GSE32062
        ovDataList[[20]] = GSE32063
        ovDataList[[21]] = GSE44104
        ovDataList[[22]] = GSE49997
        ovDataList[[23]] = GSE51088
        ovDataList[[24]] = PMID15897565
        ovDataList[[25]] = PMID17290060
        ovDataList[[26]] = PMID19318476
        ovDataList[[27]] = TCGA
        
        infoString = c("E.MTAB.386", "GSE2109" ,"GSE6008", "GSE6822", "GSE8842", "GSE9891", "GSE12418", 
                       "GSE12470", "GSE13876", "GSE14764", "GSE17260", "GSE18520", "GSE19829"
                       , "GSE20565", "GSE26193", "GSE26712", "GSE30009", "GSE30161", "GSE32062"
                       , "GSE32063", "GSE44104", "GSE49997", "GSE51088", "PMID15897565", "PMID17290060", "PMID19318476", "TCGA")
        
        #remove data without recurrence survival info
        
        remove = c()
        for(i in 1:length(ovDataList))
        {
          numMissing = sum(is.na(ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence))
          numPresent = length(ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence)
          if(numMissing == numPresent)
            remove = c(remove, i)
          numMissing = sum(is.na(ovDataList[[i]]@phenoData@data$recurrence_status))
          numPresent = length(ovDataList[[i]]@phenoData@data$recurrence_status)
          if(numMissing == numPresent)
            remove = c(remove, i)
          
        }
        
        remove = unique(remove)
        ovDataList[remove] = NULL
        infoString = infoString[-remove]
    }
      cancerType = "ovarian"
      }else if(cancerType == "brAll" || cancerType == "brMet"|| cancerType == "brMost" || cancerType == "brTcga"){
      if(mord == FALSE)
        setwd("C:\\Users\\micha\\Documents\\MetaGxBreast_0.99.0\\MetaGxBreast\\data")
      if(mord == TRUE)
        setwd("//mnt//work1//users//home2//mzon//sweaveBreastAll//data")
        
      #if(mord = TRUE)
      #  setwd("C:\\Users\\micha\\Documents\\MetaGxBreast_0.99.0\\MetaGxBreast\\data")
      #setwd("D:\\PMH Data\\Breast Cancer Data\\Breast Cancer Data Haibe-Kains\\wetransfer-474706\\MetaGxBreast_0.99.0\\MetaGxBreast\\data")
      fileInd = 1
      for(i in 1:length(list.files()))
      {
        #print(i)
        #i = 2 file is corrupted, think it prevents proper loading of metaGxBreast package
        if(i != 2)
        {
          load(list.files()[i])
          fileInd = fileInd + 1 
        }
      }
      if(cancerType == "brAll")
      {
        ovDataList[[1]] = CAL
        ovDataList[[2]] = DUKE
        ovDataList[[3]] = METABRIC
        ovDataList[[4]] = NKI
        ovDataList[[5]] = PNC
        ovDataList[[6]] = STNO2
        ovDataList[[7]] = TCGA
        ovDataList[[8]] = TRANSBIG
        ovDataList[[9]] = UCSF
        ovDataList[[10]] = UNC4 
        infoString = c("CAL", "DUKE", "METABRIC", "NKI", "PNC", "STNO2", "TCGA", "TRANSBIG", "UCSF", "UNC4")
      }
      if(cancerType == "brMet")
      {
        ovDataList[[1]] = METABRIC
        infoString = c("METABRIC")
      }
      if(cancerType == "brTcga")
      {
        ovDataList[[1]] = TCGA
        infoString = c("TCGA")
      }
      #brMost is recurrence free survival for all data sets minus tcga and metabric
      if(cancerType == "brMost")
      {
        ovDataList[[1]] = CAL
        ovDataList[[2]] = DFHCC
        ovDataList[[3]] = DFHCC2
        ovDataList[[4]] = DFHCC3
        ovDataList[[5]] = DUKE
        ovDataList[[6]] = DUKE2
        ovDataList[[7]] = EMC2
        ovDataList[[8]] = EORTC10994
        ovDataList[[9]] = EXPO
        ovDataList[[10]] = FNCLCC
        ovDataList[[11]] = GSE25066
        ovDataList[[12]] = GSE32646
        ovDataList[[13]] = GSE48091
        ovDataList[[14]] = GSE58644
        ovDataList[[15]] = HLP
        ovDataList[[16]] = IRB
        ovDataList[[17]] = KOO
        ovDataList[[18]] = LUND
        ovDataList[[19]] = LUND2
        ovDataList[[20]] = MAINZ
        ovDataList[[21]] = MAQC2
        ovDataList[[22]] = MCCC
        ovDataList[[23]] = MDA4
        ovDataList[[24]] = METABRIC
        ovDataList[[25]] = MSK
        ovDataList[[26]] = MUG
        ovDataList[[27]] = NCCS
        ovDataList[[28]] = NCI
        ovDataList[[29]] = NKI
        ovDataList[[30]] = PNC
        ovDataList[[31]] = STK
        ovDataList[[32]] = STNO2
        ovDataList[[33]] = TCGA
        ovDataList[[34]] = TRANSBIG
        ovDataList[[35]] = UCSF
        ovDataList[[36]] = UNC4
        ovDataList[[37]] = UNT
        ovDataList[[38]] = UPP
        ovDataList[[39]] = VDX
          
        infoString = c("CAL", "DFHCC", "DFHCC2", "DFHCC3", "DUKE", "DUKE2", "EMC2", "EORTC10994", "EXPO",
                       "FNCLCC", "GSE25066", "GSE32636","GSE48091","GSE58644","HLP","IRB", "KOO","LUND",
                       "LUND2","MAINZ","MAQC2","MCCC","MDA4","METABRIC","MSK","MUG","NCCS","NCI","NKI",
                       "PNC","STK","STNO2","TCGA","TRANSBIG","UCSF","UNC4","UNT","UPP","VDX")
      
        #since metabric and tcga dont have recurrence free survival, can just include all data and then remove
        remove = c()
        for(i in 1:length(ovDataList))
        {
          numMissing = sum(is.na(ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence))
          numPresent = length(ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence)
          if(numMissing == numPresent)
            remove = c(remove, i)
          numMissing = sum(is.na(ovDataList[[i]]@phenoData@data$recurrence_status))
          numPresent = length(ovDataList[[i]]@phenoData@data$recurrence_status)
          if(numMissing == numPresent)
            remove = c(remove, i)
          
        }
        
        remove = unique(remove)
        ovDataList[remove] = NULL
        infoString = infoString[-remove]
        
      }
      cancerType = "breast"
    }
    
  }

  if(cancerTypeOrig == "ovOs" || cancerTypeOrig == "brTcga" || cancerTypeOrig == "brMet")
  {
    dataStr = "os"
    for(i in 1:length(ovDataList))
    {
      #dataEvStatus[[i]] = ovDataList[[i]]@phenoData@data$vital_status;
      #dataEvTime[[i]] = ovDataList[[i]]@phenoData@data$days_to_death/365.25;
      
      events = ovDataList[[i]]@phenoData@data$vital_status;
      evTime = ovDataList[[i]]@phenoData@data$days_to_death/365.25;
      events = as.character(ovDataList[[i]]@phenoData@data$vital_status)
      events[events == "living"] = 0
      events[events == "deceased"] = 1
      events = as.numeric(events)
      newDat = censor.time(evTime, events, time.cens = 10)
      evTime = newDat[[1]]
      events = newDat[[2]]
      dataEvStatus[[i]] = events
      dataEvTime[[i]] = evTime
    }
  }
  if(cancerTypeOrig == "brMost" || cancerTypeOrig == "ovRec")
  {
    dataStr = "rec"
    for(i in 1:length(ovDataList))
    {
      #dataEvStatus[[i]] = ovDataList[[i]]@phenoData@data$recurrence_status;
      #dataEvTime[[i]] = ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence/365.25;
      
      events = ovDataList[[i]]@phenoData@data$recurrence_status;
      evTime = ovDataList[[i]]@phenoData@data$days_to_tumor_recurrence/365.25;
      events = as.character(ovDataList[[i]]@phenoData@data$recurrence_status)
      events[events == "norecurrence"] = 0
      events[events == "recurrence"] = 1
      events = as.numeric(events)
      newDat = censor.time(evTime, events, time.cens = 10)
      evTime = newDat[[1]]
      events = newDat[[2]]
      dataEvStatus[[i]] = events
      dataEvTime[[i]] = evTime
    }
  }
  
  if(mord == TRUE)
    setwd("//mnt//work1//users//home2//mzon//sweaveBreastAll")
  for(i in 1:length(ovDataList))
  {
    idString = c(idString, ovDataList[[i]]@experimentData@pubMedIds)
    sampleString = c(sampleString, dim(ovDataList[[i]]@assayData$exprs)[2])
  }
  forestPlotFrame = as.data.frame(cbind(infoString, idString, sampleString))
  
  #determine the best probes based for the esets if they are not provided
  ovDataListBestProbes = bestProbes
  if(length(bestProbes) == 0)
  {
    if(mord == FALSE)
      setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code");
    source("probeSelection.R");
    
    ovDataListBestProbes = list()
    for(i in 1:length(ovDataList))
    {
      #print(i)
      ovDataListBestProbes[[i]] = getBestProbes(ovDataList[[i]])
    }
  }
  #print(length(ovDataListFunc[[1]]$Verhaak.subtypes) == 0)
  #note the name is a relic from the past, should really just be called subtypes
  if(length(ovDataList[[1]]$Verhaak.subtypes) == 0)
  {
    #print("determining subtypes")
    if(cancerType == "ovarian")
    {
      if(mord == FALSE)
        setwd("C:\\Users\\micha\\Documents\\OVC Project\\Subtyping")
      #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Subtyping")
      source("verhaakSubtypeFunction.R")
      supplementary.data.sheet7 <- read.xls("JCI65833sd1.xls", sheet=7, skip=1,  perl = "C:\\Perl64\\bin\\perl.exe")
      supplementary.data.sheet1 <- read.xls("JCI65833sd1.xls", skip=1,  perl = "C:\\Perl64\\bin\\perl.exe")
      
      for(i in 1:length(ovDataList))
      {
        #print(i)
        if(sum(i == broken) == 0)
          ovDataList[[i]] = getVerhaakSubtypes(ovDataList[[i]], supplementary.data.sheet1, supplementary.data.sheet7)
      } 
    }
    else if(cancerType == "breast"){
      if(mord == FALSE)
        setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code")
      #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code")
      source("getBreastSubtypes.R")
     
      for(i in 1:length(ovDataList))
      {
        #print(i)
        if(sum(i == broken) == 0)
          ovDataList[[i]] = getBreastSubtypes(ovDataList[[i]])
      } 
    }
  }
  
  if(mord == FALSE)
    setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code");
  
  #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code");
  source("getSignatureSigscore.R");
  source("getSurvivalStatsAndPlot.R");
  
  #conduct a survival analysis on all of the patients, i.e no stratification by subtypye
  survInfoAll = list()
  survInfoAll[[1]] = as.data.frame(NULL)
  
  dInds = as.data.frame(NULL)
  dIndsLow = as.data.frame(NULL)
  dIndsHigh = as.data.frame(NULL)
  for(i in 1:length(ovDataList))
  {
    #print(i)
    data = ovDataList[[i]];
    dataVals = data@assayData$exprs;
    dataValsOrig = dataVals
    dataInfo = data@featureData@data;
    #print(as.character(rownames(dataInfo)[1000]))
    #print(as.character(dataInfo$probeset[1000]))
    eventTimes = dataEvTime[[i]];
    dataTimeToDeath = eventTimes
    #below is factor in TCGA and character vector in CAL data
    eventStatuses = dataEvStatus[[i]];
    dataVitalStat = eventStatuses
    
    #bestProbes = ovDataListBestProbes[[2]]
    #survInfo[[2]] has the d indexes
    survInfoAll = getSurvInfo(survInfoAll, ovDataListBestProbes[[i]], geneSig, geneSigCoef,  dataVals, dataInfo, eventTimes, eventStatuses)
  }
  #if NA gene entrez then cant name columns of empty frame
  if(!(length(survInfoAll[[1]]) == 0))
    colnames(survInfoAll[[1]]) = c("timeToDeath", "vitalStat", "groups", "scoreVals", "scoreValsOrig")
  
  dInds = cbind(survInfoAll[[2]])
  dIndsLow = cbind(survInfoAll[[2]]*exp(-qnorm(.05, lower.tail = FALSE)*survInfoAll[[3]]))
  dIndsHigh = cbind(survInfoAll[[2]]*exp(qnorm(.05, lower.tail = FALSE)*survInfoAll[[3]]))
  
  #rownames(forestPlotFrame) = survInfoAll[[5]]
  
  subtypes = unique(ovDataList[[1]]@phenoData@data$Verhaak.subtypes)
  numSubtypes = length(unique(ovDataList[[1]]@phenoData@data$Verhaak.subtypes))
  numExtras = 4
  infoPerSubtype = 3
  
  retList = list()
  for(k in 1:(numExtras + (numSubtypes + 1)*infoPerSubtype))
    retList[[k]] = "empty"
  retList[[1]] = genesEntrez
  retList[[2]] = ovDataListBestProbes
  retList[[3]] = ovDataList
  retList[[4]] = list(numExtras, infoPerSubtype, numSubtypes, subtypes)
  
  getSurvStatsGen = function(ovDataList, broken, subtypeName, retList, retListInd, dInds, dIndsLow, dIndsHigh, survInfoAll)
  {
    #retListLong = getSurvStatsGen(ovDataList, broken, subtypes[i], retList, retListInd, dInds, dIndsLow, dIndsHigh, survInfoAll)
    
    survInfoDat = list()
    survInfoDat[[1]] = as.data.frame(NULL)
    
    dindVec = c()
    dindSeVec = c()
    dindPvalVec = c()
    subIndVec = c()
    shift = 0
    #for(i in 1:length(ovDataList))
    #{
    #  #data set 10 and 7 do not work with the ovarian subtype code
    #  if(sum(i == broken) == 0)
    #  {
    #    #this was tested to ensure it still only grabs the subtype patients with valid data
    #    #remove first element called "start", was in place so empty lists still generate
    #    invalidPats = as.integer(survInfoAll[[6]][[i]][-1])
    #    subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes);
    #    if(length(invalidPats > 1))
    #      subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes[-invalidPats]);
    #    
    #    subInds = which(subtypes == subtypeName)
    #    data = ovDataList[[i]]
    #    dataVals = data@assayData$exprs[, subInds];
    #    dataInfo = data@featureData@data;
    #    dataTimeToDeath = data@phenoData@data$days_to_death[subInds]/365.25;
    #    dataVitalStat = data@phenoData@data$vital_status[subInds];
        
    #    dataSetInds = subInds + shift
    #    subIndVec = c(subIndVec, dataSetInds)
        
    #    shift = shift +length(subtypes)
        
    #    dindDat = survInfoAll[[1]][dataSetInds, ]
        #all the d indexes were the same yet, individual genes must have yielded NA values
        #should double check this later
    #    dind <- D.index(x=dindDat[-which(is.na(dindDat[,4])), 4], surv.time=dindDat[-which(is.na(dindDat[,1])), 1], surv.event=dindDat[-which(is.na(dindDat[,2])), 2])
    #    dindVec = c(dindVec, dind$d.index)
    #    dindSeVec = c(dindSeVec, dind$se)
    #    dindPvalVec = c(dindPvalVec, dind$p.value)
        
        #survInfoDat = getSurvInfo(survInfoDat, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, subInds], dataInfo, dataTimeToDeath[subInds], dataVitalStat[subInds])
        #print(i) 
    #  }
    #}
    #with shift variable this gives the same result as before/below commented out code
    #survInfoDat[[1]] = survInfoAll[[1]][subIndVec, ]
    #survInfoDat[[2]] = dindVec
    #survInfoDat[[3]] = dindSeVec
    #survInfoDat[[4]] = dindPvalVec
    #survInfoDat[[5]] = survInfoAll[[5]]
    
    for(i in 1:length(ovDataList))
    {
      #data set 10 and 7 do not work with the ovarian subtype code
      if(sum(i == broken) == 0)
      {
        subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes);
        
        data = ovDataList[[i]]
        dataVals = data@assayData$exprs;
        dataInfo = data@featureData@data;
        eventTimes = dataEvTime[[i]];
        eventStatuses = dataEvStatus[[i]];    
        
        subInds = which(subtypes == subtypeName)
        
        survInfoDat = getSurvInfo(survInfoDat, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, subInds], dataInfo, eventTimes[subInds], eventStatuses[subInds])
        #print(i) 
      }
    }

    dInds = cbind(dInds, survInfoDat[[2]])
    dIndsLow = cbind(dIndsLow, survInfoDat[[2]]*exp(-qnorm(.05, lower.tail = FALSE)*survInfoDat[[3]]))
    dIndsHigh = cbind(dIndsHigh, survInfoDat[[2]]*exp(qnorm(.05, lower.tail = FALSE)*survInfoDat[[3]]))
    
    if(!(length(survInfoDat[[1]]) == 0))
    {
      survStatsDat = getSurvivalStats(survInfoDat, FALSE, threeGroups, TRUE)
      pvalCombDAllStats = 2*pnorm(-abs(log(survStatsDat[[1]][1])/survStatsDat[[1]][2]))
      
      retList[[retListInd]] = survStatsDat[[1]]
      retList[[retListInd + 1]] = survStatsDat[[2]]
      #survStatsDat[[3]] is the same as survInfoDat, survInfoDat[[1]] and [[2]] are d.index and se for the data sets
      retList[[retListInd + 2]] = survStatsDat[[3]]
    }
    retList = list(retList, dInds, dIndsLow, dIndsHigh)
    return(retList)
  }
  
  if(!(length(survInfoAll[[1]]) == 0))
  {
    
    survStatsAll = getSurvivalStats(survInfoAll, FALSE, threeGroups, TRUE)
    pvalCombDAll = 2*pnorm(-abs(log(survStatsAll[[1]][1])/survStatsAll[[1]][2]))
    
    retList[[numExtras + 1]] = survStatsAll[[1]]
    retList[[numExtras + 2]] = survStatsAll[[2]]
    #survStatsAll[[3]] is the same as survInfoAll, survInfoAll[[1]] and [[2]] are d.index and se for the data sets
    retList[[numExtras + 3]] = survStatsAll[[3]]  
    
    for(i in 1:numSubtypes)
    {
      #print(i)
      #ovDataList[[19]] dInds are NA for Paul Signature
      retListInd = (numExtras + i*infoPerSubtype + 1)
      retListLong = getSurvStatsGen(ovDataList, broken, subtypes[i], retList, retListInd, dInds, dIndsLow, dIndsHigh, survInfoAll)
      retList = retListLong[[1]]
      
      dInds = retListLong[[2]]
      dIndsLow = retListLong[[3]]
      dIndsHigh = retListLong[[4]]
    }
  }
  numGenesPresent = survInfoAll[[5]]
  forestPlotFrame = cbind(forestPlotFrame, dInds, dIndsLow, dIndsHigh, numGenesPresent)
  
  #if(!(length(survInfoAll[[1]]) == 0))
  #{
  #  survInfoImr = list()
  #  survInfoImr[[1]] = as.data.frame(NULL)
  #  survInfoPro = list()
  #  survInfoPro[[1]] = as.data.frame(NULL)
  #  survInfoDif = list()
  #  survInfoDif[[1]] = as.data.frame(NULL)
  #  survInfoMes = list()
  #  survInfoMes[[1]] = as.data.frame(NULL)
    
  #  for(i in 1:length(ovDataList))
  #  {
      #data set 10 and 7 do not work with the ovarian subtype code
  #    if(sum(i == broken) == 0)
  #    {
  #      subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes);
  #      
  #      data = ovDataList[[i]]
  #      dataVals = data@assayData$exprs;
  #      dataInfo = data@featureData@data;
  #      dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
  #      dataVitalStat = data@phenoData@data$vital_status;
  #      
  #      imrInds = which(subtypes == "IMR")
  #      proInds = which(subtypes == "PRO")
  #      difInds = which(subtypes == "DIF")
  #      mesInds = which(subtypes == "MES")
  #      
  #      survInfoImr = getSurvInfo(survInfoImr, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, imrInds], dataInfo, dataTimeToDeath[imrInds], dataVitalStat[imrInds])
  #      survInfoPro = getSurvInfo(survInfoPro, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, proInds], dataInfo, dataTimeToDeath[proInds], dataVitalStat[proInds])
  #      survInfoDif = getSurvInfo(survInfoDif, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, difInds], dataInfo, dataTimeToDeath[difInds], dataVitalStat[difInds])
  #      survInfoMes = getSurvInfo(survInfoMes, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, mesInds], dataInfo, dataTimeToDeath[mesInds], dataVitalStat[mesInds])
        
        #print(i) 
  #    }
  #  }  
  #  
  #  if(!(length(survInfoImr[[1]]) == 0))
  #  {
  #    survStatsImr = getSurvivalStats(survInfoImr, FALSE, threeGroups, TRUE)
  #    pvalCombDAllImr = 2*pnorm(-abs(log(survStatsImr[[1]][1])/survStatsImr[[1]][2]))
  #    
  #    retList[[7]] = survStatsImr[[1]]
  #    retList[[8]] = survStatsImr[[2]]
  #    retList[[9]] = survStatsImr[[3]]
  #  }
  #  if(!(length(survInfoPro[[1]]) == 0))
  #  {
  #    survStatsPro = getSurvivalStats(survInfoPro, FALSE, threeGroups, TRUE)
  #    pvalCombDAllPro = 2*pnorm(-abs(log(survStatsPro[[1]][1])/survStatsPro[[1]][2]))
  #    
  #    retList[[10]] = survStatsPro[[1]]
  #    retList[[11]] = survStatsPro[[2]]
  #    retList[[12]] = survStatsPro[[3]]
  #  }
  #  if(!(length(survInfoDif[[1]]) == 0))
  #  {
  #    survStatsDif = getSurvivalStats(survInfoDif, FALSE, threeGroups, TRUE)
  #    pvalCombDAllDif = 2*pnorm(-abs(log(survStatsDif[[1]][1])/survStatsDif[[1]][2]))
  #    
  #    retList[[13]] = survStatsDif[[1]]
  #    retList[[14]] = survStatsDif[[2]]
  #    retList[[15]] = survStatsDif[[3]]
  #  }
  #  if(!(length(survInfoMes[[1]]) == 0))
  #  {
  #    survStatsMes = getSurvivalStats(survInfoMes, FALSE, threeGroups, TRUE)
  #    pvalCombDAllMes = 2*pnorm(-abs(log(survStatsMes[[1]][1])/survStatsMes[[1]][2]))
  #    
  #    retList[[16]] = survStatsMes[[1]]
  #    retList[[17]] = survStatsMes[[2]]
  #    retList[[18]] = survStatsMes[[3]]
  #  }
  #  
  #}
  dataTabInfo = getDataInfo(ovDataList, infoString, dataStr, ovDataListBestProbes, dataEvStatus, dataEvTime)
  returnList = list(retList, forestPlotFrame, dataTabInfo)
  return(returnList)
  
}

getDataInfo = function(dataList, infoString, dataStr, bestProbesList, evStatList, evTimeList)
{
  #dont censor data for table (Haibe-Kains)
  dataInfoFrame = as.data.frame(NULL)
  
  rowInfo = c()
  for(i in 1:length(dataList))
  {
    if(dataStr == "os")
    {
      numEvStr = "# Deceased"
      evTimeStr = "Median Survival Time (Years)"
      
      events = as.character(dataList[[i]]@phenoData@data$vital_status)
      events[events == "living"] = 0
      events[events == "deceased"] = 1
      events = as.numeric(events)
      #events = evStatList[[i]]
      #if(length(which(is.na(events))) > 0)
      #  numEvents = sum(events[-which(is.na(events))])
      #if(length(which(is.na(events))) == 0)
      #  numEvents = sum(events)
      
      evTime = dataList[[i]]@phenoData@data$days_to_death/365.25
      #if(length(which(is.na(evTime))) > 0)
      #  medEvTime = median(evTime[-which(is.na(evTime))])
      #if(length(which(is.na(evTime))) == 0)
      #  medEvTime = median(evTime)
      #evTime = evTimeList[[i]]
      #evTime = censor.time(evTime, events, time.cens = 10)[1]
      #events = censor.time(evTime, events, time.cens = 10)[2]
      
    }
    if(dataStr == "rec")
    {
      numEvStr = "# Tumor Recurrences"
      evTimeStr = "Median Recurrence Time (Years)"
      
      events = as.character(dataList[[i]]@phenoData@data$recurrence_status)
      events[events == "norecurrence"] = 0
      events[events == "recurrence"] = 1
      events = as.numeric(events)
      
      #if(length(which(is.na(events))) > 0)
      #  numEvents = sum(events[-which(is.na(events))])
      #if(length(which(is.na(events))) == 0)
      #  numEvents = sum(events)
      
      evTime = dataList[[i]]@phenoData@data$days_to_tumor_recurrence/365.25
      #if(length(which(is.na(evTime))) > 0)
      #  medEvTime = median(evTime[-which(is.na(evTime))])
      #if(length(which(is.na(evTime))) == 0)
      #  medEvTime = median(evTime)
      #evTime = censor.time(evTime, events, time.cens = 10)[1]
      #events = censor.time(evTime, events, time.cens = 10)[2]
      
    }
    survInf = survfit(formula = Surv(evTime, events) ~ 1)
    medEvTime = as.numeric(summary(survInf)$table["median"])
    numPatients = as.integer(summary(survInf)$table["records"])
    numEvents = as.integer(summary(survInf)$table["events"])
    #rowInfo = c(data name, number patients, number of genes, number of events, median event time, platform)
    platformStr = as.character(dataList[[i]]@experimentData@other$platform_shorttitle)
    if(length(platformStr) == 0)
      platformStr = "Unknown"
    rowInfo = c(infoString[i], numPatients, length(bestProbesList[[i]]), numEvents, sprintf("%.3g",medEvTime), platformStr)
    dataInfoFrame = rbind(dataInfoFrame, rowInfo, stringsAsFactors = FALSE)
  }
  colnames(dataInfoFrame) = c("Data Name", "# Patients", "# Genes", numEvStr, evTimeStr, "platform")
  #colNameVec = c("Data Name", "# Patients", "# Genes", numEvStr, evTimeStr, "platform")
  
  return(list(dataInfoFrame, infoString))
}