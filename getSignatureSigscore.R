getSurvInfo <- function(survObjData, bestProbes, geneSignature, geneSigCoef, dataValsOrig, dataInfo, dataTimeToDeath, dataVitalStat)
{

  #new start
  dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
  #some missing ids that are "///"
  slashStarts = which(substr(dataInfoEntrezGene.ID, 1, 1) == "/")
  slashStLeng = nchar(dataInfoEntrezGene.ID[slashStarts])
  slashRem = substr(dataInfoEntrezGene.ID[slashStarts], 4, slashStLeng)
  dataInfoEntrezGene.ID[slashStarts] = slashRem
  missingEnts = which(dataInfoEntrezGene.ID == "")
  #set the missing ents to -100 so the probe selection works
  #dont end up using any of them anyways
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
  #new end
  
  #old line
  #dataInfoEntrezGene.ID = as.numeric(dataInfo$EntrezGene.ID)
  
  dataVals = dataValsOrig[bestProbes, ]
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[bestProbes]
  #print("hello")
  dataInfoProbeset = as.vector(dataInfo$probeset[bestProbes])
  dataInfoGene = as.vector(dataInfo$gene[bestProbes])
  
  geneSigInds = c()
  geneSigEntrez = c()
  geneSigNames = c()
  dataGeneSigInds = c()
  genesFoundInd = c()
  for(i in 1:length(geneSignature))
  {
    geneSigEntrez = geneSignature[i]
    #geneSigNames = c(geneSigNames, xOrig$geneNames[geneSigInd])
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
    if(length(which(dataInfoEntrezGene.ID == geneSigEntrez)) > 0)
      genesFoundInd = c(genesFoundInd, i)
    dataGeneSigInds = c(dataGeneSigInds, which(dataInfoEntrezGene.ID == geneSigEntrez));
    #print(i)
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
  }
  #geneSigInds = which(xOrig$probeIds %in% geneSig);
  #geneSigEntrez = xOrig$entrezIds[geneSigInds]
  
  #tcgaGeneSigInds = which(tcgaInfo$EntrezGene.ID %in% geneSigEntrez);
  
  scores = NULL
  topGenes = dataInfoProbeset[dataGeneSigInds]
  
  #if none of the randomly selected genes are in the data set then move on
  if(length(topGenes) < 1)
  {
    if(length(survObjData) == 1)
    {
      survObjData[[2]] = c(NA)
      survObjData[[3]] = c(NA)
      survObjData[[4]] = c(NA)
      survObjData[[5]] = c(length(topGenes))
      survObjData[[6]] = list()
      survObjData[[6]][[length(survObjData[[6]]) + 1]] = FALSE
    }else if(length(survObjData) > 1){
      survObjData[[2]] = c(survObjData[[2]], NA)
      survObjData[[3]] = c(survObjData[[3]], NA)
      survObjData[[4]] = c(survObjData[[4]] ,NA)
      survObjData[[5]] = c(survObjData[[5]] ,length(topGenes))
      survObjData[[6]][[length(survObjData[[6]]) + 1]] = FALSE
    }
    return(survObjData) 
  }
  
  #print(paste("number of genes in data set from signature found = ", length(topGenes)))
  
  #if even one gene expression value is NA than the score is NA, for now just remove NAs at end
  # may be better to check how many of the genes are NA and than remove NA genes while getting
  # the score with the remaining genes if enough genes remain 
  #print("hello")
  #
  #top genes are the ensemble/probe id for the genes in the signature, found using there entrez gene id
  #have verified by viewing dataInfo frame that the returned topGenes correspond to the genes of interest
  #as one can see the gene names in the frame
  
  #addition of below line fixes issue, cant assume rownames of dataVals are probesets. Usually are,
  #but when they are not sig.score doesnt find correct columns
  rownames(dataVals) = dataInfoProbeset
  for(i in 1:dim(dataVals)[2])
  {
    testSamp = t(dataVals[, i, drop = FALSE]);
    mysig <- cbind("probe"=topGenes, "EntrezGene.ID"=NA, "coefficient"= (as.numeric(geneSigCoef[genesFoundInd])/sum(abs(as.numeric(geneSigCoef[genesFoundInd])))))
    rownames(mysig) = topGenes;
    scores <- rbind(scores, cbind("score"=genefu::sig.score(x=mysig, testSamp, annot=NULL, do.mapping=FALSE, signed=TRUE)$score, "fold"=i))
  }
  #print("hi")
  scoreVals = as.numeric(scores[,1])
  
  #probably should adjust to not assume NA in time to death is the same as NA in vital stat
  #dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
  #dataVitalStat = data@phenoData@data$vital_status;
  invalidPatVec = c("start")
  
  unknown = which(is.na(dataTimeToDeath) == TRUE)
  if(length(unknown) > 0)
  {
    dataTimeToDeath = dataTimeToDeath[-unknown];
    dataVitalStat = dataVitalStat[-unknown]
    scoreVals = scoreVals[-unknown]
    invalidPatVec = c(invalidPatVec, unknown)
  }
  
  #Is this dirupting something else, why didnt I add this earlier if it was so easy --> double check
  unknown = which(is.na(dataVitalStat) == TRUE)
  if(length(unknown) > 0)
  {
    dataTimeToDeath = dataTimeToDeath[-unknown];
    dataVitalStat = dataVitalStat[-unknown]
    scoreVals = scoreVals[-unknown]
    invalidPatVec = c(invalidPatVec, unknown)
  }
  
  missingGenes = which(is.na(scoreVals) == TRUE)
  if(length(missingGenes > 0))
  {
    dataTimeToDeath = dataTimeToDeath[-missingGenes];
    dataVitalStat = dataVitalStat[-missingGenes]
    scoreVals = scoreVals[-missingGenes] 
    invalidPatVec = c(invalidPatVec, missingGenes)
  }
  
  threeGroups = FALSE
  if(threeGroups == TRUE)
  {
    sortInd = sort(scoreVals, decreasing = TRUE, index.return=TRUE)$ix;
    thirdData = round(length(scoreVals)/3)
    dataGroups = matrix(0, length(scoreVals))
    dataGroups[sortInd[1:thirdData]] = 0;
    dataGroups[sortInd[thirdData:(2*thirdData)]] = 1;
    dataGroups[sortInd[(2*thirdData + 1):length(sortInd)]] = 2;
  }
  if(threeGroups == FALSE)
  {
    scoreMed = median(scoreVals)
    dataGroups = as.integer((scoreVals <= scoreMed));
  }
  #print(sum(is.na(dataGroups)))
  
  #dataVitalStat = as.character(dataVitalStat)
  #dataVitalStat[dataVitalStat == "living"] = 0
  #dataVitalStat[dataVitalStat == "deceased"] = 1
  #dataVitalStat[dataVitalStat == "norecurrence"] = 0
  #dataVitalStat[dataVitalStat == "recurrence"] = 1
  #dataVitalStat = as.numeric(dataVitalStat)
  
  scoreValsOrig = scoreVals
  quantVals = quantile(scoreVals, c(.025, .975))
  oldLow = quantVals[1]
  oldHigh = quantVals[2]
  newLow = -1
  newHigh = 1
  scoreVals = newLow*(1 - (scoreVals - oldLow)/(oldHigh - oldLow)) + newHigh*((scoreVals - oldLow)/(oldHigh - oldLow))
  
  dataMatrix = matrix(c(dataTimeToDeath, dataVitalStat, dataGroups, scoreVals, scoreValsOrig), nrow = length(dataVitalStat), ncol = 5)
  survObjData[[1]] = rbind(survObjData[[1]], dataMatrix)
  
  dind <- D.index(x=scoreVals, surv.time=dataTimeToDeath, surv.event=dataVitalStat)
  
  if(length(survObjData) == 1)
  {
    survObjData[[2]] = c(dind$d.index)
    survObjData[[3]] = c(dind$se)
    survObjData[[4]] = c(dind$p.value)
    survObjData[[5]] = c(length(topGenes))
    survObjData[[6]] = list()
    survObjData[[6]][[length(survObjData[[6]]) + 1]] = invalidPatVec
  }else if(length(survObjData) > 1){
    survObjData[[2]] = c(survObjData[[2]], dind$d.index)
    survObjData[[3]] = c(survObjData[[3]], dind$se)
    survObjData[[4]] = c(survObjData[[4]] ,dind$p.value)
    survObjData[[5]] = c(survObjData[[5]] ,length(topGenes))
    survObjData[[6]][[length(survObjData[[6]]) + 1]] = invalidPatVec
  }
  
  return(survObjData)
}
