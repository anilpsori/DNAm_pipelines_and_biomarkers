############################
#######SKIN AND BLOOD CLOCK
############################
DNAmAgeSkinClock = function(dat0,path){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #directory is path where datSkinClock.csv file is stored
  
  #R functions for transforming age
  adult.age1=20
  trafo= function(x,adult.age=adult.age1){x=(x+1)/(1+adult.age); y=ifelse(x<=1,log( x),x-1);y}
  anti.trafo= function(x,adult.age=adult.age1) {ifelse(x<0, (1+adult.age)*exp(x)-1,(1+adult.age)*x+adult.age) }
  
  datClock=read.csv(paste0(path,"/datSkinClock.csv"))
  
  selectCpGsClock=is.element(dat0[,1], as.character(datClock[-1,1])) #logical for x %in% y
  datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
  
  colnames(datMethClock0)=as.character(dat0[selectCpGsClock ,1])
  
  # Reality check: the following output should only contain numeric values.
  # Further, the column names should be CpG identifiers (cg numbers).
  #datMethClock0[1:5,1:5]
  datMethClock= data.frame(datMethClock0[as.character(datClock[-1,1])])
  
  # The number of rows should equal the number of samples (Illumina arrays)
  #dim(datMethClock)
  
  #Output DNAm age estimator for the skin & blood clock
  DNAmAgeSkin=as.numeric(anti.trafo(datClock$Coef[1]+as.matrix(datMethClock) %*% as.numeric(datClock$Coef[-1])))
  
  return(DNAmAgeSkin)
}

