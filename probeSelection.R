getBestProbes = function(data)
{
  dataValsOrig = data@assayData$exprs;
  #Result is much worse if the below line is included
  #dataVals <- normalizeBetweenArrays(dataVals,method="quantile")
  dataInfo = data@featureData@data;
  dataInfoEntrezGene.ID = dataInfo$EntrezGene.ID
  

  dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
  #some missing ids that are "///"
  slashStarts = which(substr(dataInfoEntrezGene.ID, 1, 1) == "/")
  slashStLeng = nchar(dataInfoEntrezGene.ID[slashStarts])
  slashRem = substr(dataInfoEntrezGene.ID[slashStarts], 4, slashStLeng)
  dataInfoEntrezGene.ID[slashStarts] = slashRem
  missingEnts = which(dataInfoEntrezGene.ID == "")
  #set the missing ents to -100 so the probe selection works
  #ont end p using any of them anyways
  dataInfoEntrezGene.ID[missingEnts]  = -100
  
  #multiple entrez ids have ent1///ent2///....///entN
  multEnts = which(grepl("/", dataInfoEntrezGene.ID))
  multStrings = dataInfoEntrezGene.ID[multEnts]
  firstSlashes = rapply(gregexpr(pattern = "/", multStrings), function(x) head(x, 1))
  firstEnts = substr(multStrings, 1, firstSlashes - 1)
  dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  multEnts = which(grepl(",", dataInfoEntrezGene.ID))
  multStrings = dataInfoEntrezGene.ID[multEnts]
  firstSlashes = rapply(gregexpr(pattern = ",", multStrings), function(x) head(x, 1))
  firstEnts = substr(multStrings, 1, firstSlashes - 1)
  dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  dataInfoEntrezGene.ID = as.numeric(dataInfoEntrezGene.ID)
  
  dataInfoProbeset = as.vector(dataInfo$probeset)
  dataInfoGene = as.vector(dataInfo$gene)
  
  #repeatEnt = duplicated(dataInfoEntrezGene.ID)
  #repeatEnt = which(repeatEnt == TRUE)
  
  #rmList = c()
  #for(i in 1:length(repeatEnt))
  #{
  #  repInd = repeatEnt[i]
  #  entId = dataInfoEntrezGene.ID[repInd]
  #  probesWithId = which(dataInfoEntrezGene.ID == entId)
  #  probes = dataValsOrig[probesWithId,]
  #  iqrs = rowIQRs(probes);
  #  winnerInd = which(iqrs == max(iqrs))
  #  keeper = probesWithId[winnerInd]
  #  probesRm = probesWithId[-winnerInd]
  #  rmList = c(rmList, probesRm)
  #  bestProbes = c(bestProbes, keeper)
  #}
  #above code checks more then once per repeat id, remove non unique
  #rmList = unique(rmList)
  #bestProbes is now a just a list of all the rows to keep in the data set
  #bestProbes = c(1:length(dataInfoEntrezGene.ID))
  #bestProbes = bestProbes[-rmList]
  
  #return(bestProbes)
  
  #old method below that is faster, should probably switch back to this and add a reordering 
  
  sortInd = sort(dataInfoEntrezGene.ID, decreasing = FALSE, index.return=TRUE)$ix;
  #definitly sorts least to greatest
  dataInfoEntrezGene.IDPreSort = dataInfoEntrezGene.ID
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[sortInd];
  dataVals = dataValsOrig[sortInd, ];
  dataInfoProbeset = dataInfoProbeset[sortInd];
  dataInfoGene = dataInfoGene[sortInd];

  bestProbes = c()
  i = 1;
  while(i < dim(dataVals)[1])
  {
    entrezId = dataInfoEntrezGene.ID[i];
    origInd = i;
    probesWithId = c();
    while(dataInfoEntrezGene.ID[i] == entrezId)
    {
      probesWithId = c(probesWithId, i);
      i = i + 1;
      if(i == dim(dataVals)[1])
      {
        if(dataInfoEntrezGene.ID[i] == entrezId)
          probesWithId = c(probesWithId, i)
        if(dataInfoEntrezGene.ID[i] != entrezId)
          bestProbes = c(bestProbes, i)
        
        break
      }
    }
    if(length(probesWithId) > 1)
    {
      probes = dataVals[probesWithId,]
      iqrs = rowIQRs(probes);
      keepProbe = which(iqrs == max(iqrs)) + (origInd - 1);
      bestProbes = c(bestProbes, keepProbe)
    }
    if(length(probesWithId) == 1)
    {
      bestProbes = c(bestProbes, origInd)
    }
  }
  #get the rownames of the probes to keep and then match them with the rows from the unsorted data
  keepRows = rownames(dataVals)[bestProbes]
  bestProbes = match(keepRows, rownames(dataValsOrig))
  return(bestProbes)  
}
  