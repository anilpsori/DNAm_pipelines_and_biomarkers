############################
#######DNAm Age Lin et al 2015: PMID: 26928272
############################
LinAge = function(dat0, path=NULL){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #directory is path where Lin_Weights.csv file is stored
  Weights=read.csv(paste0(path,"/Lin_Weights.csv"))
  Weights$probeID = as.character(Weights$probeID)
  Weights$Coefficient = as.numeric(as.character(Weights$Coefficient))
  
  selectCpGsClock = dat0$ProbeID %in% Weights$probeID #logical for x %in% y
  
  #print how many probes are missing
  print("calculating score for.. Lin DNAm Age")
  print(paste0("number of probes missing.. ",print((length(Weights$probeID)-1) - sum(selectCpGsClock)))) #+1 to account for intercept term
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with Weights
  Weights2 = Weights[Weights$probeID %in% c("Intercept",colnames(datMethClock0)),]
  datMethClock= datMethClock0[,match(Weights2$probeID[-1],colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  LinAge=as.numeric(Weights2$Coefficient[1] + (as.matrix(datMethClock) %*% as.numeric(Weights2$Coefficient[-1])))
  
  return(LinAge)
}