############################
#######MiAge Youn et al 2018: PMID: 29160179
############################
MiAge = function(dat0,path){
  #dat0: is a dataframe of beta values with colnames; ProbeID, Sample1, Sample2, etc
  #directory is path where files for MiAge calculations are stored
  
  #load functionals and datasets
  source(paste0(path,"/function_library.r")) ### library of all functions used in the calculation of mitotic ages    
  load(paste0(path,"/site_specific_parameters.Rdata")) #### site-specific parameters
  clocksitesID=as.matrix(read.csv(paste0(path,"/Additional_File1.csv"),header=T))[,1]   ### epigenetic clock CpG sites: 268 (of which 241 are on EPIC)
  
  #print how many probes are missing
  selectCpGsClock = clocksitesID %in% dat0$ProbeID #logical for x %in% y
  print("calculating score for.. DNAm MiAge")
  print(paste0("number of probes missing is.. ",(length(selectCpGsClock) - sum(selectCpGsClock))))
  
  #transform dataframe of beta values to matrix with probe IDs as rownames and sample ID as colnames
  rownames(dat0) = dat0$ProbeID
  dat1 = as.matrix(dat0[,-1])
  
  #get beta values
  beta=  dat1[ match(clocksitesID[selectCpGsClock],rownames(dat1)),]##select clock CpG sites 
  
  #estimate MiAge
  b=methyl.age[[1]][selectCpGsClock];c=methyl.age[[2]][selectCpGsClock];d=methyl.age[[3]][selectCpGsClock]
  n=mitotic.age(beta,b,c,d) ### estimated mitotic age
  
  return(n)
}