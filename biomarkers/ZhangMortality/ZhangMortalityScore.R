############################
#######Zhang Mortality Score Zhang et al 2017: PMID: 28303888
############################
ZhangMortalityScore = function(dat0, path=NULL){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #path: directory where files of Zhang et al are store (not needed, weights are below)
  
  #weights to calculate Zhang Mortality Score
  Weights = data.frame(CpG = c("cg01612140", "cg05575921","cg06126421","cg08362785","cg10321156","cg14975410","cg19572487","cg23665802","cg24704287","cg25983901"),
                       Coefficient =c(-0.38253, -0.92224, -1.70129, 2.71749, -0.02073, -0.04156, -0.28069, -0.89440, -2.98637, -1.80325))
  
  selectCpGsClock = dat0$ProbeID %in% Weights$CpG
  
  #print how many probes are missing
  print("calculating score for.. Zhang Mortality Score")
  print(paste0("number of probes missing.. ",(length(Weights$CpG)) - sum(selectCpGsClock)))
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with datClock
  Weights2 = Weights[Weights$CpG %in% c(colnames(datMethClock0)),]
  datMethClock= datMethClock0[,match(Weights2$CpG,colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  score=as.numeric((as.matrix(datMethClock) %*% as.numeric(Weights2$Coefficient)))
  
  return(score)
}
