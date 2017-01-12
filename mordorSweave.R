
args = commandArgs(trailingOnly = TRUE)
print(length(args))
args

## n shows the number of up and down regulated genes. can be 30, 50 or 70
mord = as.numeric(args[1])
if(is.na(mord))
{
  fileString = "Normal Breast Tissue (29) VS Basal Type Brest Tissue  (59).xls"
  #assign("mord", FALSE, envir = .GlobalEnv)
}
if(mord == TRUE)
{
  fileString = as.character(args[2])
}
mord = TRUE

install.packages("forestplot_1.6.tar.gz", type ="source", repos = NULL)
library(forestplot)


#if(mord == TRUE)
#{
  assign("mord",TRUE, envir = .GlobalEnv)
  assign("fileId",fileString, envir = .GlobalEnv)
  .libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
#}
library(knitr)

knit2pdf("getOvarianGenePrognosisReport.Rnw")

