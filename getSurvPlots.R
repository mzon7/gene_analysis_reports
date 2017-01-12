makeSurvPlot = function(survStats, survObj, survInfo, threeGroups, dindStr = "", titleStr = "Kaplan-Meier Survival Curve", xLoc = 2.)
{
  temp = NULL
  if(threeGroups == FALSE)
  {
    timesToDeath = survInfo[,1]
    vitalStats = survInfo[,2]
    groups = survInfo[,3]
    scores = survInfo[,4]
    
    naList = c(which(is.na(scores)), which(is.na(vitalStats)), which(is.na(timesToDeath)))
    naList = unique(naList)
    
    if(length(naList > 0))
      dind <- D.index(x=scores[-naList], surv.time=timesToDeath[-naList], surv.event=vitalStats[-naList])
    if(length(naList) == 0)
      dind <- D.index(x=scores, surv.time=timesToDeath, surv.event=vitalStats)
    
    bb <- survdiff(Surv(timesToDeath, vitalStats) ~ groups,rho=0)
    logRankP = (1 - pchisq(bb$chisq, 1))
    
    datFrame = data.frame("surv.time" = timesToDeath, "surv.event"=vitalStats, "strat"=groups)
    ddweights <- array(1, dim=nrow(datFrame))
    #km.coxph.plot(formula.s = Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
    km.coxph.plot(formula.s= Surv(timesToDeath, vitalStats) ~ groups, data.s=datFrame,
                  weight.s=ddweights, x.label="Time (years)", y.label="Probability of survival",
                  main.title=titleStr, leg.text=paste(c("High Score Group", "Low Score Group"), "   ", sep=""),
                  leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkred"),
                  .lty=c(1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE, leg.bty = "n", o.text = "")
    
    if(logRankP == 0)
      logRankP = 1.00*10^-16
    
    text(xLoc, 0.115, paste("D index = ", dindStr))
    text(xLoc + 0.5, 0.03, location = "bottomright", paste("Log Rank P = ", format(logRankP, scientific =  TRUE, digits = 2)))
    
    #text(temp, paste("D index = ", sprintf("%.2g", dind$d.index)))
    #plot(survObj, lty = 2:3,col=3:4, main=titleStr, xlab="time (years)", ylab="survival function")
    #high sigscore class = 0, low sigscore class = 1
    #legend("top",legend=c(paste("Group 1 - High Sigscore, n =", sum(survInfo$groups == 0)), paste("Group 2 - Low Sigscore, n =", sum(survInfo$groups == 1))), col=3:4, lty=2:3, horiz=FALSE, bty='n')
  }
  if(threeGroups == TRUE)
  {
    timesToDeath = survInfo[,1]
    vitalStats = survInfo[,2]
    groups = survInfo[,3]
    datFrame = data.frame("surv.time" = timesToDeath, "surv.event"=vitalStats, "strat"=groups)
    ddweights <- array(1, dim=nrow(datFrame))
    #km.coxph.plot(formula.s = Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
    km.coxph.plot(formula.s= Surv(timesToDeath, vitalStats) ~ groups, data.s=datFrame,
                  weight.s=ddweights, x.label="Time (years)", y.label="Probability of survival",
                  main.title=titleStr, leg.text=paste(c("High Score Group", "Intermediate Score Group", "Low Score Group"), "   ", sep=""),
                  leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen", "darkred"),
                  .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE, leg.bty = "n")
    }

  #temp$x = 4000.5
  #temp$y = 0.5
  #text(temp, cat("D.index p value", dind$p.value))
  #Corner_text <- function(text, location="topright"){
  #  legend(location,legend=text, bty ="n", pch=NA) 
  ##textVal = paste("Combined D.index ", round(survStats[1], digits = 6));
  #Corner_text(text= textVal, location = "bottomleft")
    
  #logRankP = survStats[3]
  #temp$y = 0.6
  #textVal = paste("log rank test p value", round(logRankP, digits = 5));
  #Corner_text(text= textVal, location = "bottomright")
    
  #get gene names for significant genes
  #sigGeneNames = xOrig$geneNames[geneSigInds]
  #xOrig$geneNames[which(xOrig$probeIds %in% geneListSecond500[[472]])]
    
  #for survivak curve, strat is group number, survival time obvious,x is my sigscores,
  #log rank test is enough for now, D.index later maybe
  
}
