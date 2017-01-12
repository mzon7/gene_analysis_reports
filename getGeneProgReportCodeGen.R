

setupReport = function(neel)
{
  
  #Obtain/install the MetaGXOvarian package prior to runnin the analysis
  #install.packages("C:\\Users\\micha\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\MetaGxOvarian_0.99.0.tar.gz", repos = NULL, type="source")
  
  #may want to check this infrequent warning, likely from coxph
  #Error in fitter(X, Y, strats, offset, init, control, weights = weights,  : 
  #  (converted from warning) Ran out of iterations and did not converge
  #In addition: Warning message:
  #In SuppDists::normOrder(N = length(sx)) : 
  #Values may be inaccurate because of the size of N
  
  #install.packages('knitr', dependencies = TRUE)
  Sys.setenv(JAVA_HOME='C:\\Users\\micha\\Documents\\jre1.8.0_111')
  library(XLConnect)
  library(gdata)
  library(xlsx)
  
  library(limma);
  library(genefu)
  library(illuminaHumanv4.db)
  #library(TCGAbiolinks)
  library(matrixStats)
  library(ROCR)
  library(org.Hs.eg.db)
  library(parallel)
  #library(jetset)
  library(Biobase)
  library(gdata)
  library(GSVA)
  library(org.Hs.eg.db)
  library(forestplot)
  
  #To Fix
  # - map gene names and try to fill in ones not provided by accepting file gene names
  # - in report if no gene name use symbol and then entrez id so no NA heading
  # - if no name/data section shouldnt be gene not present, should put entrez or gene name/symbol/entrez
  #   and than have a note that the gene was not present in the data sets

  #Q1: why is gene name list 4 elements not 5
  #Note: gene name list elements become section names so replace NA wwith whatever one wants
  #need to automate finding the number of sheets
  
  #need to add option to use entrez ids from file first as 2/5 didnt map and just better to let user specify
  geneSigFrame = as.data.frame(NULL)
  
  
  if(neel == "paul")
  {
    numSheets = 1
    setwd("C:\\Users\\micha\\Documents\\PMH Research\\Michael Paul Data")
    fileId = "Paul Gene Signature.xls" 
  }
  
  #note: skipping last 2 sheets, last sheet is sum of 1:16 sheets, raw data just a bunch of genes, 5 times as many as the other sheets combines , about 2500
  if(neel == TRUE)
  {
    numSheets = 1
    if(mord == FALSE)
    {
      setwd("C:\\Users\\micha\\Documents\\PMH Research\\BenNeelCmap\\Data")
      fileId = "Neel Gene Signature adjP less 05 and abs(logFC) big 1.xls" 
    }
  }
  if(neel == FALSE)
  {
    numSheets = 16
    if(mord == FALSE)
    {
      setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis")
      fileId = "Normal Breast Tissue (29) VS Basal Type Brest Tissue  (59).xls"
    }
  }
    
  geneSymbolList = list()
  geneNameList = list()
  geneDescripList = list()
  geneEntrezList = list()
  for(i in 1:numSheets)
  {
    geneSymbolList[[i]] = c()
    geneEntrezList[[i]] = c()
    geneNameList[[i]] = c()
    geneDescripList[[i]] = c()
    #cant have sheet = i argument if only 1 sheet is present
    if(mord == FALSE)
    {
      if(numSheets > 1)
        genesOfInt <- read.xls(fileId, sheet=i,  perl = "C:\\Perl64\\bin\\perl.exe")
      if(numSheets == 1)
        genesOfInt <- read.xls(fileId, perl = "C:\\Perl64\\bin\\perl.exe") 
    }
    if(mord == TRUE)
    {
      if(numSheets > 1)
        genesOfInt <- read.xls(fileId, sheet=i)
      if(numSheets == 1)
        genesOfInt <- read.xls(fileId)
    }
    
    #if(numSheets > 1)
    #  genesOfInt <- read.xls(fileId, sheet=i)
    #if(numSheets == 1)
    #  genesOfInt <- read.xls(fileId)
    
    geneSymbolList[[i]] = unique(as.character(genesOfInt$Gene.Symbol))
    if(!is.null(genesOfInt$entrezid))
    {
      geneEntrezList[[i]] = unique(as.character(genesOfInt$entrezid))
      #if(is.null(genesOfInt$Gene.Symbol))
      #{
        mapping = select(org.Hs.eg.db, keytype="ENTREZID", keys=as.character(genesOfInt$entrezid), columns=c("SYMBOL"))
        mapping = mapping[!duplicated(mapping$ENTREZID),]
        if(length(geneSymbolList[[i]]) == 0)
          geneSymbolList[[i]] = mapping[, 2]
        if(length(geneSymbolList[[i]]) > 0)
          missingSymbs = which(is.na(geneSymbolList[[i]]))
        if(length(missingSymbs) > 0)
          geneSymbolList[[i]][missingSymbs] = mapping[missingSymbs, 2]
      #}
    }else if(!is.null(genesOfInt$ensembleid)){
      ensIds = as.character(genesOfInt$ensembleid)
      mapping = select(org.Hs.eg.db, keytype="ENSEMBL", keys=ensIds, columns=c("ENTREZID"))
      mapping = mapping[!duplicated(mapping$ENSEMBL),]
      geneEntrezList[[i]] = mapping[, 2]
      #if(is.null(genesOfInt$Gene.Symbol))
      #{
        mapping = select(org.Hs.eg.db, keytype="ENSEMBL", keys=ensIds, columns=c("SYMBOL"))
        mapping = mapping[!duplicated(mapping$ENSEMBL),]
        if(length(geneSymbolList[[i]]) == 0)
          geneSymbolList[[i]] = mapping[, 2]
        if(length(geneSymbolList[[i]]) > 0)
          missingSymbs = which(is.na(geneSymbolList[[i]]))
        if(length(missingSymbs) > 0)
          geneSymbolList[[i]][missingSymbs] = mapping[missingSymbs, 2]
      #}
    }else if(!is.null(genesOfInt$Gene.Symbol)){
      geneSymbolList[[i]] = unique(as.character(genesOfInt$Gene.Symbol))
      mapping = select(org.Hs.eg.db, keytype="SYMBOL", keys=geneSymbolList[[i]], columns=c("ENTREZID"))
      mapping = mapping[!duplicated(mapping$SYMBOL),]
      geneEntrezList[[i]] = mapping[, 2]
      #if(is.null(genesOfInt$ensembleid))
      #{
        mapping = select(org.Hs.eg.db, keytype="SYMBOL", keys=geneSymbolList[[i]], columns=c("ENSEMBL"))
        mapping = mapping[!duplicated(mapping$SYMBOL),]
        #store ensemble ids if need be, above gets them though
        if(length(geneSymbolList[[i]]) == 0)
          geneSymbolList[[i]] = mapping[, 2]
        if(length(geneSymbolList[[i]]) > 0)
          missingSymbs = which(is.na(geneSymbolList[[i]]))
        if(length(missingSymbs) > 0)
          geneSymbolList[[i]][missingSymbs] = mapping[missingSymbs, 2]
      #}
    }
    
    #adjust not to take user inputs as above code removes non unique elements and one must do the same
    #for geneNameList to have the report put the right names on the right curves
    #geneNameList[[i]] = as.character(genesOfInt$GENE.NAME)
    #if(is.null(genesOfInt$GENE.NAME))
    #{
      mapping = select(org.Hs.eg.db, keytype="SYMBOL", keys=geneSymbolList[[i]], columns=c("GENENAME"))
      mapping = mapping[!duplicated(mapping$SYMBOL),]
      geneNameList[[i]] = mapping[, 2]
    #}
    #geneDescripList[[i]] = as.character(genesOfInt$GENE.NAME)
    #if(is.null(genesOfInt$Description))
    #{
      mapping = select(org.Hs.eg.db, keytype="SYMBOL", keys=geneSymbolList[[i]], columns=c("GENENAME"))

      mapping = mapping[!duplicated(mapping$SYMBOL),]
      geneDescripList[[i]] = mapping[, 2]
    #}
      missName = which(is.na(geneNameList[[i]]))
      if(length(missName) > 0)
      {
        for(k in 1:length(missName))
        {
          missInd = missName[k]
          geneNameList[[i]][missInd] = geneSymbolList[[i]][missInd]
          geneDescripList[[i]][missInd] = geneSymbolList[[i]][missInd]
          if(is.na(geneSymbolList[[i]][missInd]))
          {
            geneNameList[[i]][missInd] = paste("Entrez ID:", geneEntrezList[[i]][missInd])
            geneDescripList[[i]][missInd] = paste("Entrez ID: ",geneEntrezList[[i]][missInd])
          }
          
        }
      }
    
    if(length(genesOfInt$geneInSig) > 0)
    {
      geneSigInds = which(genesOfInt$geneInSig == TRUE)
      geneSigFrame = cbind(geneEntrezList[[i]][geneSigInds], genesOfInt$t[geneSigInds], geneSymbolList[[i]][geneSigInds], geneNameList[[i]][geneSigInds], geneDescripList[[i]][geneSigInds])
      colnames(geneSigFrame) = c("entrez", "tstat", "geneSymbol", "geneName", "geneDescription")
    }
    
  }
  
  #x <- org.Hs.egSYMBOL2EG
  # Get the entrez gene identifiers that are mapped to a gene symbol
  #mapped_genes <- mappedkeys(x)
  # Convert to a list
  #xx <- as.list(x[mapped_genes])
  
  #geneEntrezList = list()
  #for (i in 1:length(geneSymbolList))
  #{
  #  geneEntrezList[[i]] = c()
  #  geneEntrezs = c()
  #  for(j in 1:length(geneSymbolList[[i]]))
  #  {
  #    index = which(names(xx) == geneSymbolList[[i]][j])
  #    geneEntrezs = c(geneEntrezs, as.numeric(xx[[index]]))
  #  }
  #  geneEntrezs = unique(geneEntrezs)
  #  geneEntrezList[[i]] = geneEntrezs
  #}
  
  
  library(knitr)
  #library(MetaGxOvarian)
  #source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
  #source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
  
  if(mord == FALSE)
    setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code")
  #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project\\R Code")
  source("testGenesPrognostics.R")
  source("getSurvPlots.R")
  
  return(list(geneEntrezList, geneNameList, geneDescripList, geneSigFrame, geneSymbolList))
}

getSigResults = function(cancerType, neel)
{
  geneFileInfo = setupReport(neel)
  geneEntrezList = geneFileInfo[[1]]
  geneNameList = geneFileInfo[[2]]
  geneDescripList = geneFileInfo[[3]]
  geneSigFrame = geneFileInfo[[4]]
  geneSymbList = geneFileInfo[[5]]
  
  pathwayNames = c("PPAR Signaling","Complement and Coag Cascades","Regulation of lipolysis in adip", "Leukocy transendothelial migra", "
                   P4-mediated oocyte mat","Adrenergic signaling in cardiom","ECM-receptor interaction", "cGMP-PKG signaling ",
                   " Calcium signaling ", "Cell adhesion molecules (CAMs)", " Oocyte meiosis", "cAMP signaling ", "Proteoglycans in cancer","Focal adhesion", "Cell cycle", " PI3K-Akt")
  
  x <- org.Hs.egSYMBOL2EG
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  
  geneInfo = list()
  dataTabInfo = list()
  geneInfo[[1]] = as.data.frame(NULL)
  geneInfo[[2]] = as.data.frame(NULL)
  geneInfo[[3]] = as.data.frame(NULL)
  infString = NULL
  
  geneSigForestPlotFrame = as.data.frame(NULL)
  geneSigInfo = list()
  geneSigInfoList = list()
  if(dim(geneSigFrame)[1] > 0)
  {
    geneSigInfoMult = analyzeSignature(cancerType, geneSigFrame)
    geneSigInfo = geneSigInfoMult[[1]]
    geneSigForestPlotFrame = geneSigInfoMult[[2]]
    dataTabInfo = geneSigInfoMult[[3]]
    infString = dataTabInfo[[2]]
    geneInfo[[2]] = geneSigInfo[[2]]
    geneInfo[[3]] = geneSigInfo[[3]]
    
    
    for(k in 1:length(geneSigInfo))
    {
      #curSize = length(geneInfoList)
      geneSigInfoList[[length(geneSigInfoList) + 1]] = geneSigInfo[[k]]
      #geneInfoList[[curSize + 1]] = list(geneInfoList[[curSize + 1]], geneInfo[[k]])
    }
  }
  
  #remove these lines in the future as data sets may be added to metGx
  #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project")
  #load("bestProbesFunc.RData")
  #load("ovDataListFunc.RData")
  #geneInfo[[2]] = bestProbesFunc 
  #geneInfo[[3]] = ovDataListFunc 
  
  #cat("\\tableofcontents")
  
  iEnd = length(geneEntrezList)
  #iEnd = 1
  
  geneInfoList = list()
  geneSymbVec = c()
  pathwayNameVec = c()
  statTableList = list()
  #summaryTabList = list()

  for(i in 1:iEnd)
  {
    pathwayName = paste(pathwayNames[i], " Pathway")
    pathwayNameVec = c(pathwayNameVec, pathwayName)
    #cat(paste("\\section{",pathwayName,"}", sep="")) 
    
    #jEnd = 3
    #for(j in 1:length(geneEntrezList[[i]]))
    jEnd = length(geneEntrezList[[i]])
    #jEnd = 2
    for(j in 1:jEnd)
    {
      #print(j)
      #knit_print("test1")
      #kable(statTable, digits = 5, caption = paste("survival statistics summary for gene ", geneSymb))
      
      statTabRow = c()
      statTable = as.data.frame(NULL)
      #print(geneEntrezList[[i]][j])
      #geneInfo[[7]] has the info needed for the forest plot from the analysis on ALL the patients
      #print("hi")
      #if running this manually change cancerType back to normal
      geneData = testOvCancSurv(geneEntrezList[[i]][j], 1, cancerType, bestProbes = geneInfo[[2]], ovDataList = geneInfo[[3]], infoString = infString)
      #print("bye")
      #get list of dataset names as only generate it during first iteration
      dataTabInfo = geneData[[3]]
      infString = dataTabInfo[[2]]
      
      geneInfo = geneData[[1]]
      forestPlotFrame = geneData[[2]]
      dataTabInfo = geneData[[3]]
      #**this number 4 below is the only variables that has to be changed upon making changes to testGenePrognostics.R
      indexInfo = geneInfo[[4]]
      numExtras = indexInfo[[1]] 
      infoPerSubtype = indexInfo[[2]] 
      numSubtypes = indexInfo[[3]] 
      subtypes = indexInfo[[4]] 
      
      if(i == 1 & j == 1)
      {
        summaryTabList = list()
        geneSymbRowsList = list()
        
        for(k in 1:(numSubtypes + 1))
        {
          summaryTabList[[k]] = as.data.frame(NULL)
          geneSymbRowsList[[k]] = c()
        }
        geneSymbRowsList[[numSubtypes+2]] = "this sets previous elements = null and is needed for appending to prior elements"
        
      }
      
      for(k in 1:length(geneInfo))
      {
        #curSize = length(geneInfoList)
        geneInfoList[[length(geneInfoList) + 1]] = geneInfo[[k]]
        #geneInfoList[[curSize + 1]] = list(geneInfoList[[curSize + 1]], geneInfo[[k]])
      }
      geneSymb = names(which(xx == as.character(geneEntrezList[[i]][j])))
      geneSymbVec = c(geneSymbVec, geneSymb)
      
      #geneSymbRows = c(geneSymbRows, geneSymb)
      #print(geneSymb)
      pathwayGene = geneSymb
      #cat(paste("\\subsection{",geneNameList[[i]][j],"}", sep=""))
      #cat("\\subsection{",geneNameList[[i]][j],"}")
      #cat("\\vbox{Gene Description:",geneDescripList[[i]][j],"}") 
      
      if(identical(geneInfo[[numExtras + 1]], "empty"))
        cat("\\vbox{Gene not present in the expression data of the external data sets}")
      
      if(!identical(geneInfo[[numExtras + 1]], "empty"))
      {
        for(k in 1:(numSubtypes + 1))
        {
          geneInfInd = numExtras + 1 + infoPerSubtype*(k-1)
          if(!identical(geneInfo[[geneInfInd]], "empty"))
          {
            if(k == 1)
              statTabRow = c("All Data")
            if(k > 1)
              statTabRow = c(statTabRow, as.character(subtypes[k - 1])) #a lot of this +/- 1 business as the all data is done and then subtypes are done
            #d.index*exp(-alpha/2*standardError)
            seLow = sprintf("%.2g", geneInfo[[geneInfInd]][1]*exp(-qnorm(.05, lower.tail = FALSE)*geneInfo[[geneInfInd]][2]))
            seHigh = sprintf("%.2g",geneInfo[[geneInfInd]][1]*exp(qnorm(.05, lower.tail = FALSE)*geneInfo[[geneInfInd]][2]))
            statTabVec = data.frame(sprintf("%.2g", geneInfo[[geneInfInd]][1]), paste("(", seLow, ", ", seHigh, ")", sep=""), format(geneInfo[[geneInfInd]][4], scientific =  TRUE, digits = 2), format(geneInfo[[geneInfInd]][3], scientific =  TRUE, digits = 2))
            #statTable = rbind(statTable, c(geneInfo[[geneInfInd]][4], geneInfo[[geneInfInd]][1], geneInfo[[geneInfInd]][2], geneInfo[[geneInfInd]][3]))
            #statTabVec = c(sprintf("%.4g", geneInfo[[geneInfInd]][1]), sprintf("%.4g", geneInfo[[geneInfInd]][4]), sprintf("%.4g", geneInfo[[geneInfInd]][3]))
            colnames(statTabVec) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P")
            statTable = rbind(statTable, statTabVec)
            #colnames(statTable) = c("D Index P Val", " D Index", "D Index Standard Error")
            colnames(statTable) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P")
            summaryTabList[[k]] = rbind(summaryTabList[[k]], statTabVec)
            #print(statTabVec)
            #summaryTabList[[k]] = rbind(summaryTabList[[k]], c(geneInfo[[geneInfInd]][4], geneInfo[[geneInfInd]][1], geneInfo[[geneInfInd]][2], geneInfo[[geneInfInd]][3]))
            geneSymbRowsList[[k]] = c(geneSymbRowsList[[k]], geneSymb)
            
          }
        }
      }
      
      if(length(statTable) > 0)
      {
        #below line broken, NAs in statTabRow :s? just needed as.character around subtypes[k] above 
        rownames(statTable) = statTabRow
        #print(xtable::xtable(statTable, digits = 5, caption = paste("survival statistics summary for gene ", geneSymb)))
      }
      statTableList[[length(statTableList)+1]] = statTable
    }
    #cat("\\clearpage") 
  }
  
  setupTable = function(table, geneRowNames)
  {
    table = cbind(table, geneRowNames)
    table = unique(table)
    table = table[,1:4]
    
    table = cbind(table, format(p.adjust(as.vector(table[,3]), method = "fdr"), scientific =  TRUE, digits = 2))
    table = cbind(table, format(p.adjust(as.vector(table[,4]), method = "fdr"), scientific =  TRUE, digits = 2))
    
    colnames(table) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P", "D Index FDR", "Log Rank FDR")
    #colnames(table) = c("D Index P Val", " D Index", "D Index Standard Error", "Log Rank Test P Val")
    rownames(table) = unique(geneRowNames)
    sortInd = order(as.vector(table["Log Rank Test P"]))
    table = table[sortInd, ]
    return(table)
  }
  
  for(k in 1:(numSubtypes + 1))
  {
    summaryTabList[[k]] = setupTable(summaryTabList[[k]], geneSymbRowsList[[k]])
  }
  
  #summaryTabList = list(statTabAll, statTabMes, statTabDif, statTabPro, statTabImr)
  if(dim(geneSigFrame)[1] > 0)
    forestPlotFrame = geneSigForestPlotFrame
  
  geneInfoLists = list(geneInfoList, list(geneSymbVec, pathwayNameVec, geneNameList, geneDescripList, geneSymbList), geneFileInfo, statTableList, summaryTabList, indexInfo, geneSigInfoList, geneSigFrame, forestPlotFrame, dataTabInfo)
  return(geneInfoLists)
}

analyzeSignature = function(cancerType, geneSigFrame)
{
  geneSigInfo = list()
  geneSigInfo[[1]] = as.data.frame(NULL)
  geneSigInfo[[2]] = as.data.frame(NULL)
  geneSigInfo[[3]] = as.data.frame(NULL)
  
  #whats going on with geneSigInfo[[2]] for paul signature ovOs
  geneSigInfo = testOvCancSurv(geneSigFrame[, 1], geneSigFrame[, 2], cancerType, bestProbes = geneSigInfo[[2]], ovDataList = geneSigInfo[[3]], infoString = geneSigInfo[[4]])
  return(geneSigInfo)
}

createForestPlot = function(dInds, dIndsLow, dIndsHigh, tabText)
{
  # plotInfo = createForestPlot(dInds, dIndsLow, dIndsHigh, tabText)
  #need standard errors for final d index, solve for them below
  dIndsSe = log(dIndsHigh/dInds)/qnorm(.05, lower.tail = FALSE)
  dindComb = combine.est(dInds, dIndsSe, hetero = TRUE)
  dindCombLow = combine.est(dIndsLow, dIndsSe, hetero = TRUE)
  dindCombHigh = combine.est(dIndsHigh, dIndsSe, hetero = TRUE)
  # Cochrane data from the 'rmeta'-package

  plotObj <- 
    structure(list(
      mean  = c(NA, dInds, dindComb$estimate), 
      lower = c(NA, dIndsLow, dindCombLow$estimate),
      upper = c(NA, dIndsHigh, dindCombHigh$estimate)),
      .Names = c("mean", "lower", "upper"), 
      row.names = c(NA, -(length(dInds) + 2)), 
      class = "data.frame")
  
  tabletext<-cbind(
    c("Study", as.character(tabText[,1]), "Summary"),
    c("Patients", as.character(tabText[,3]), NA),
    c("Genes", as.character(tabText[,4]), NA),
    c("D Index", as.character(format(dInds, digits = 2)), as.character(format(dindComb$estimate, digits = 2))))
  
  #forestplot(tabletext, 
  #           plotObj,new_page = TRUE,
  #           is.summary=c(TRUE,rep(FALSE,10),TRUE),
  #           clip=c(0.1,2.5), 
  #           xlog=TRUE, 
  #           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
  
  return(list(tabletext, plotObj))
  
}




