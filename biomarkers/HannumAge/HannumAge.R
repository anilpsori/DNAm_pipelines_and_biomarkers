############################
#######Hannum DNAmAge
############################
##The resulting age estimate is based on the 71 CpGs and coefficient values from the third supplementary table
##The authors developed this age prediction method by using an elastic net regression model for predicting chronological age based on DNA methylation levels from whole blood (450K).

HannumAge = function(dat0,path){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #path is path where Hannum_TableS2_71probes.txt file is stored
  datClock = read.table(paste0(path,"/Hannum_TableS2_71probes.txt"),header=T,sep="\t")
  
  selectCpGsClock = dat0$ProbeID %in% datClock$Marker #logical for x %in% y
  
  #print how many probes are missing
  print("calculating score for.. Hannum DNAm age")
  print(paste0("number of probes missing..",print((length(datClock$Marker)-1) - sum(selectCpGsClock)))) #+1 to account for intercept term
  
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
  
  #match order of CpGs with datClock
  datClock2 = datClock[datClock$Marker %in% colnames(datMethClock0),]
  datMethClock= datMethClock0[,match(datClock2$Marker,colnames(datMethClock0))]
  
  #Output DNAm age estimator for the skin & blood clock
  DNAmAgeHannum=as.numeric(as.matrix(datMethClock) %*% as.numeric(datClock2$Coefficient))
  
  return(DNAmAgeHannum)
}