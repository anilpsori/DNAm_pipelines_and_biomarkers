############################
#######Weidner et al 2014: PMID: 24490752
############################
WeidnerAge = function(dat0, path=NULL){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #path: directory where files of Weidner et al are store (not needed, weights are below)
  
  #weights to calculate Weidner 3 CpG DNAm Age
  Weights = data.frame(CpG = c("Intercept", "cg02228185", "cg25809905", "cg17861230"),
                       Coefficient =c(38.0, -26.4, -23.7, 164.7))
  
  selectCpGsClock = dat0$ProbeID %in% Weights$CpG[-1]
  
  #print how many probes are missing
  print("calculating score for.. Weidner DNAm Age")
  print(paste0("number of probes missing.. ",(length(Weights$CpG[-1])) - sum(selectCpGsClock)))
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with datClock
  Weights2 = Weights[Weights$CpG %in% c("Intercept", colnames(datMethClock0)),]
  datMethClock= datMethClock0[,match(Weights2$CpG[-1],colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  WeidnerAge=as.numeric(Weights2$Coefficient[1] + (as.matrix(datMethClock) %*% as.numeric(Weights2$Coefficient[-1])))
  
  return(WeidnerAge)
}