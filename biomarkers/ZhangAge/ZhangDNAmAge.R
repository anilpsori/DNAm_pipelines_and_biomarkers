############################
#######ZHANG DNAmAGE: BiorXiv:327890
############################
ZhangDNAmAge = function(dat0,path){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #directory is path where Zhang_elastic_weights.csv file is stored
  datClock=read.csv(paste0(path,"/Zhang_elastic_weights.csv"))
  
  selectCpGsClock = dat0$ProbeID %in% datClock$probe.ID #logical for x %in% y
  
  #print how many probes are missing
  print("calculating score for.. Zhang DNAm age")
  print(paste0("number of probes missing..",print((length(datClock$probe.ID)-1) - sum(selectCpGsClock)))) #+1 to account for intercept term
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with datClock
  datClock2 = datClock[datClock$probe.ID %in% c("Intercept",colnames(datMethClock0)),]
  datMethClock= datMethClock0[,match(datClock2$probe.ID[-1],colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  ZhangDNAmAge=as.numeric(datClock2$Coefficient[1] + (as.matrix(datMethClock) %*% as.numeric(datClock2$Coefficient[-1])))
  
  return(ZhangDNAmAge)
}