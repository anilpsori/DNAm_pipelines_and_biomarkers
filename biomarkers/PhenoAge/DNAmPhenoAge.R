############################
#######Levine PhenoAge
############################
PhenoAge = function(dat0,path){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #directory is path where Horvath_phenoAge_weights.csv file is stored
  datClock=read.csv(paste0(path,"/Horvath_phenoAge_weights.csv"))
  
  selectCpGsClock = dat0$ProbeID %in% datClock$CpG[-1] #logical for x %in% y
  
  #print how many probes are missing
  print("calculating score for.. Levine Pheno Age")
  print(paste0("number of probes missing..",print((length(datClock$CpG)-1) - sum(selectCpGsClock)))) #-1 to account for intercept term
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with datClock
  datClock2 = datClock[datClock$CpG %in% c("Intercept",colnames(datMethClock0)),]
  datMethClock= datMethClock0[,match(datClock2$CpG[-1],colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  PhenoAge=as.numeric(datClock2$Coefficient[1] + (as.matrix(datMethClock) %*% as.numeric(datClock2$Coefficient[-1])))
  
  return(PhenoAge)
}