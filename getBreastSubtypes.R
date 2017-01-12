
getBreastSubtypes <- function(eset, supplementary.data.sheet1, supplementary.data.sheet7)
{
  
  subtypePredictions = molecular.subtyping(sbt.model = "scmod2",data = t(eset@assayData$exprs),
                                           annot = eset@featureData@data, do.mapping = TRUE)
  
  #confirmed that order of data columns is the same as the order of the subtype results, i.e subtye result 1 is for patient in column 1
  eset$Verhaak.subtypes = subtypePredictions$subtype
  
  return(eset)
  
}
