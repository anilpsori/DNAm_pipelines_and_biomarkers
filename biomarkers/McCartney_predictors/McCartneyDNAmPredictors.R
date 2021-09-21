############################
#######ZHANG DNAmAGE: BiorXiv:327890
############################
McCartneyDNAmPredictors = function(dat0,path){
  
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #Path is where files with weights are stored
  alcohol = read.csv(paste0(path,"/McCartney_Alcohol.csv"))
  bmi = read.csv(paste0(path,"/McCartney_BMI.csv"))
  bodyfat = read.csv(paste0(path,"/McCartney_BodyFat.csv"))
  cholesterol = read.csv(paste0(path,"/McCartney_cholesterol.csv"))
  edu = read.csv(paste0(path,"/McCartney_EDU.csv"))
  hdl = read.csv(paste0(path,"/McCartney_HDL.csv"))
  hdlratio = read.csv(paste0(path,"/McCartney_HDLratio.csv"))
  ldl = read.csv(paste0(path,"/McCartney_LDL.csv"))
  smoking = read.csv(paste0(path,"/McCartney_Smoking.csv"))
  whr = read.csv(paste0(path,"/McCartney_WHR.csv"))
  
  predictors = list(alcohol,bmi,bodyfat,cholesterol,edu,hdl,hdlratio,ldl,smoking,whr)
  names(predictors) = c("mccartney.alcohol","mccartney.bmi","mccartney.bodyfat","mccartney.cholesterol","mccartney.edu",
                        "mccartney.hdl","mccartney.hdlratio","mccartney.ldl","mccartney.smoking","mccartney.whr")
  
  DNAmPheno = data.frame(SampleID = colnames(dat0)[-1])
  
  for(n in 1:length(predictors)){
    #get CpGs needed to compute biomarker
    datClock = predictors[[n]]
    selectCpGsClock = dat0$ProbeID %in% datClock$CpG #logical for x %in% y
    
    #print how many probes are missing
    print(paste0("calculating score for..",names(predictors)[n]))
    print(paste0("number of probes missing..",print(length(datClock$CpG) - sum(selectCpGsClock))))
    
    datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1]))
    colnames(datMethClock0) = as.character(dat0$ProbeID[selectCpGsClock])
    
    #match order of CpGs with datClock
    datClock2 = datClock[datClock$CpG %in% colnames(datMethClock0),]
    datMethClock= datMethClock0[,match(datClock2$CpG,colnames(datMethClock0))]
    
    #calculate estimate
    DNAmPheno$Temp = as.numeric(as.matrix(datMethClock) %*% as.numeric(datClock2$Beta))
    names(DNAmPheno)[1+n]=names(predictors)[n]
  }
  
  return(DNAmPheno)
}
