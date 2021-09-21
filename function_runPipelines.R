runPipelines = function(RGset, method=c(1:101),selectCpG,returnAll=FALSE){
  ##RGset: RGChannelSetExtended object of DNAm data
  ##this function can call different DNAm data processing pipelines from the wateRmelon, minfi, en ENmix package
  ##method variable: a numeric single variable or vector that says which pipeline to run
  ##selectCpG: CpGs you would like to subset and write to .csv as output (e.g. those in datMiniAnnotation3.csv)
  ##returnALL: if processed values for ALL probes should be returned and no subsetting be done
  
  ##define methylation data
  M=RGset #this should be an RGChannelSet(Extended) object
  
  #########################
  #### ENmix methods
  #########################
  #apply normalization and background correction
  #bgParaEst: oob, est, neg
  #dyeCorr: mean, RELIC, none
  #normalization: q1,q2,q3, none
  #RCP (dyebias): TRUE or FALSE
  
  if(method == 1){
    pipeline = "enmix_oob_mean_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 2){
    pipeline = "enmix_est_mean_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 3){
    pipeline = "enmix_neg_mean_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 4){
    pipeline = "enmix_oob_relic_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 5){
    pipeline = "enmix_est_relic_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 6){
    pipeline = "enmix_neg_relic_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 7){
    pipeline = "enmix_oob_nodye_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 8){
    pipeline = "enmix_est_nodye_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 9){
    pipeline = "enmix_neg_nodye_nonorm_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 10){
    pipeline = "enmix_oob_mean_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 11){
    pipeline = "enmix_est_mean_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 12){
    pipeline = "enmix_neg_mean_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 13){
    pipeline = "enmix_oob_relic_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 14){
    pipeline = "enmix_est_relic_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 15){
    pipeline = "enmix_neg_relic_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 16){
    pipeline = "enmix_oob_nodye_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 17){
    pipeline = "enmix_est_nodye_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 18){
    pipeline = "enmix_neg_nodye_q1_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 19){
    pipeline = "enmix_oob_mean_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 20){
    pipeline = "enmix_est_mean_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 21){
    pipeline = "enmix_neg_mean_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 22){
    pipeline = "enmix_oob_relic_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 23){
    pipeline = "enmix_est_relic_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 24){
    pipeline = "enmix_neg_relic_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 25){
    pipeline = "enmix_oob_nodye_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 26){
    pipeline = "enmix_est_nodye_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 27){
    pipeline = "enmix_neg_nodye_q2_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 28){
    pipeline = "enmix_oob_mean_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 29){
    pipeline = "enmix_est_mean_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 30){
    pipeline = "enmix_neg_mean_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 31){
    pipeline = "enmix_oob_relic_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 32){
    pipeline = "enmix_est_relic_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 33){
    pipeline = "enmix_neg_relic_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 34){
    pipeline = "enmix_oob_nodye_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 35){
    pipeline = "enmix_est_nodye_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 36){
    pipeline = "enmix_neg_nodye_q3_norcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = getBeta(normalized)
  }
  
  if(method == 37){
    pipeline = "enmix_oob_mean_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 38){
    pipeline = "enmix_est_mean_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 39){
    pipeline = "enmix_neg_mean_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 40){
    pipeline = "enmix_oob_relic_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    #no normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 41){
    pipeline = "enmix_est_relic_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 42){
    pipeline = "enmix_neg_relic_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 43){
    pipeline = "enmix_oob_nodye_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 44){
    pipeline = "enmix_est_nodye_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 45){
    pipeline = "enmix_neg_nodye_nonorm_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    #normalization
    normalized = background
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 46){
    pipeline = "enmix_oob_mean_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 47){
    pipeline = "enmix_est_mean_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 48){
    pipeline = "enmix_neg_mean_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 49){
    pipeline = "enmix_oob_relic_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 50){
    pipeline = "enmix_est_relic_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 51){
    pipeline = "enmix_neg_relic_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 52){
    pipeline = "enmix_oob_nodye_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 53){
    pipeline = "enmix_est_nodye_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 54){
    pipeline = "enmix_neg_nodye_q1_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile1: will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile1")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 55){
    pipeline = "enmix_oob_mean_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 56){
    pipeline = "enmix_est_mean_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 57){
    pipeline = "enmix_neg_mean_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 58){
    pipeline = "enmix_oob_relic_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 59){
    pipeline = "enmix_est_relic_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 60){
    pipeline = "enmix_neg_relic_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 61){
    pipeline = "enmix_oob_nodye_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 62){
    pipeline = "enmix_est_nodye_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 63){
    pipeline = "enmix_neg_nodye_q2_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile2: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes
    normalized = norm.quantile(background, method="quantile2")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 64){
    pipeline = "enmix_oob_mean_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 65){
    pipeline = "enmix_est_mean_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 66){
    pipeline = "enmix_neg_mean_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #mean: dye bias correction based on averaged red/green ratio
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="mean",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 67){
    pipeline = "enmix_oob_relic_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 68){
    pipeline = "enmix_est_relic_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 69){
    pipeline = "enmix_neg_relic_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #RELIC: REgression on Logarithm of Internal Control probes 
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 70){
    pipeline = "enmix_oob_nodye_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #oob: uses out-of-band Infinium I intensities to estimate normal distribution parameters to model background noise
    #no dye correction
    background = preprocessENmix(M, bgParaEst="oob", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 71){
    pipeline = "enmix_est_nodye_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #est (needs MethylSet?): use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type
    #no dye correction
    background = preprocessENmix(M, bgParaEst="est", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  if(method == 72){
    pipeline = "enmix_neg_nodye_q3_rcp"
    print("Now running background:  ");print(pipeline)
    
    #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
    #no dye correction
    background = preprocessENmix(M, bgParaEst="neg", dyeCorr="none",nCores=1)
    
    ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
    normalized = norm.quantile(background, method="quantile3")
    
    #Probe design type bias correction using Regression on Correlated Probes (RCP) method
    probecorrected = rcp(normalized)
  }
  
  #########################
  #### WateRmelon Methods
  #########################
  if(method == 73){
    pipeline = "watermelon_dasen"
    print("Now running method:  ");print(pipeline)
    
    #dasen: same as nasen but type I and type II backgrounds are normalized first. This is our recommended   method
    probecorrected = dasen(M)
  }
  
  if(method == 74){
    pipeline = "watermelon_naten"
    print("Now running method:  ");print(pipeline)
    
    #naten: quantile normalizes methylated and unmethylated intensities separately, then calculates betas
    probecorrected = naten(M)
  }
  
  if(method == 75){
    pipeline = "watermelon_nanet"
    print("Now running method:  ");print(pipeline)
    
    #nanet: quantile normalizes methylated and unmethylated intensities together, then calculates betas. This should equalize dye bias.
    probecorrected = nanet(M)
  }
  
  if(method == 76){
    pipeline = "watermelon_nanes"
    print("Now running method:  ");print(pipeline)
    
    #nanes: quantile normalizes methylated and unmethylated intensities separately, except for type II probes where methylated and unmethylated are normalized together
    #This should equalize dye bias without affecting type I probes which are not susceptible
    probecorrected = nanes(M)
  }
  
  if(method == 77){
    pipeline = "watermelon_danes"
    print("Now running method:  ");print(pipeline)
    
    #danes: same as nanes, except typeI and type II background are equalised first
    probecorrected = danes(M)
  }
  
  if(method == 78){
    pipeline = "watermelon_danet"
    print("Now running method:  ");print(pipeline)
    
    #danet: same as nanet, except typeI and type II background are equalised first
    probecorrected = danet(M)
  }
  
  if(method == 79){
    pipeline = "watermelon_danen"
    print("Now running method:  ");print(pipeline)
    
    #danen: background equalisation only, no normalization
    probecorrected = danen(M)
  }
  
  if(method == 80){
    pipeline = "watermelon_daten1"
    print("Now running method:  ");print(pipeline)
    
    #daten1: same as naten, except typeI and type II background are equalised first (smoothed only for methylated)
    probecorrected = daten1(M)
  }
  
  if(method == 81){
    pipeline = "watermelon_daten2"
    print("Now running method:  ");print(pipeline)
    
    #daten2: same as naten, except typeI and type II background are equalised first (smoothed for methylated an unmethylated)
    probecorrected = daten2(M)
  }
  
  if(method == 82){
    pipeline = "watermelon_nasen"
    print("Now running method:  ");print(pipeline)
    
    #nasen: same as naten but typeI and typeII intensities quantile normalized separately
    temp = nasen(M)
    probecorrected = temp$beta
  }
  
  if(method == 83){
    pipeline = "watermelon_raw_BMIQ"
    print("Now running method:  ");print(pipeline)
    
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    raw = preprocessRaw(M)
    probecorrected = BMIQ(raw)
  }
  
  if(method == 84){
    pipeline = "watermelon_raw_PBC"
    print("Now running method:  ");print(pipeline)
    
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    raw = preprocessRaw(M)
    probecorrected = fuks(raw)
  }
  
  #########################
  #### Minfi methods
  #########################
  if(method == 85){
    pipeline = "minfi_raw"
    print("Now running method:  ");print(pipeline)
    
    #preprocessRaw: Converts the Red/Green channel for an Illumina methylation array into methylation signal, without using any normalization.
    temp = preprocessRaw(M)
    probecorrected = getBeta(temp)
  }
  
  if(method == 86){
    pipeline = "minfi_illumina_bg_nonorm"
    print("Now running method:  ");print(pipeline)
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction
    temp = preprocessIllumina(M, bg.correct=T,normalize="no")
    probecorrected = getBeta(temp)
  }
  
  if(method == 87){
    pipeline = "minfi_illumina_bg_normcontrol"
    print("Now running method:  ");print(pipeline)
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction and normalization
    temp = preprocessIllumina(M, bg.correct=T,normalize="controls",reference=1)
    probecorrected = getBeta(temp)
  }
  
  if(method == 88){
    pipeline = "minfi_illumina_nobg_normcontrol"
    print("Now running method:  ");print(pipeline)
    
    #bgcorrect.illumina: implements Illumina GenomeStudio background correction and normalization
    temp = preprocessIllumina(M, bg.correct=F,normalize="controls",reference=1)
    probecorrected = getBeta(temp)
  }
  
  if(method == 89){
    pipeline = "minfi_noob_dyecorr"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    temp = preprocessNoob(M,dyeCorr=T,dyeMethod="single")
    probecorrected = getBeta(temp)
  }
  
  if(method == 90){
    pipeline = "minfi_noob_nodyecorr"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    temp = preprocessNoob(M,dyeCorr=F)
    probecorrected = getBeta(temp)
  }
  
  if(method == 91){
    pipeline = "minfi_funnorm_nobg_nodyecorr"
    print("Now running method:  ");print(pipeline)
    
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform
    ##It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=F,dyeCorr=F)
    probecorrected = getBeta(temp)
  }
  
  if(method == 92){
    pipeline = "minfi_funnorm_bg_dyecorr"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform. It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=T,dyeCorr=T)
    probecorrected = getBeta(temp)
  }
  
  if(method == 93){
    pipeline = "minfi_funnorm_bg_nodyecorr"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    ##preprocessFunnorm: a between-array normalization method for the Illumina Infinium HumanMethylation450 platform. It removes unwanted variation by regressing out variability explained by the control probes present on the array.
    temp = preprocessFunnorm(M, bgCorr=T,dyeCorr=F)
    probecorrected = getBeta(temp)
  }
  
  if(method == 94){
    pipeline = "minfi_raw_quantile_strat"
    print("Now running method:  ");print(pipeline)
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays.
    ##Probes are stratified by region (CpG island, shore, etc.)
    temp = preprocessRaw(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = T)
    probecorrected = getBeta(temp2)
  }
  
  if(method == 95){
    pipeline = "minfi_raw_quantile_nostrat"
    print("Now running method:  ");print(pipeline)
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays.
    ##Probes are stratified by region (CpG island, shore, etc.)
    temp = preprocessRaw(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = F)
    probecorrected = getBeta(temp2)
  }
  
  if(method == 96){
    pipeline = "minfi_illumina_bg_quantile_strat"
    print("Now running method:  ");print(pipeline)
    
    ##PreprocessQuantile:This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays. Probes are stratified by region (CpG island, shore, etc.)
    temp = bgcorrect.illumina(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = T)
    probecorrected = getBeta(temp2)
  }
  
  if(method == 97){
    pipeline = "minfi_illumina_bg_quantile_nostrat"
    print("Now running method:  ");print(pipeline)
    
    ##preprocessQuantile: This function implements stratified quantile normalization preprocessing for Illumina methylation microarrays. Probes are stratified by region (CpG island, shore, etc.)
    temp = bgcorrect.illumina(M)
    temp2 = preprocessQuantile(temp, quantileNormalize = T,stratified = F)
    probecorrected = getBeta(temp2)
  }
  
  if(method == 98){
    pipeline = "minfi_raw_SWAN"
    print("Now running method:  ");print(pipeline)
    
    ##preprocessSWAN: Subset-quantile Within Array Normalisation (SWAN) is a within array normalisation method 
    temp = preprocessSWAN(M)
    probecorrected = getBeta(temp)
  }
  
  if(method == 99){
    pipeline = "minfi_illumina_bg_SWAN"
    print("Now running method:  ");print(pipeline)
    
    ##preprocessSWAN: Subset-quantile Within Array Normalisation (SWAN) is a within array normalisation method 
    temp = bgcorrect.illumina(M)
    temp2 = preprocessSWAN(temp)
    probecorrected = getBeta(temp2)
  }
  
  if(method == 100){
    pipeline = "cross_noob_dyecorr_BMIQ"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    temp = preprocessNoob(M,dyeCorr=T,dyeMethod="single")
    probecorrected = BMIQ(temp)
  }
  
  if(method == 101){
    pipeline = "cross_noob_nodyecorr_BMIQ"
    print("Now running method:  ");print(pipeline)
    
    #preprocessNoob: Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays
    #BMIQ: is an intra-sample normalisation procedure, correcting the bias of type-2 probe values
    temp = preprocessNoob(M,dyeCorr=F,dyeMethod="single")
    probecorrected = BMIQ(temp)
  }
  
  ##get betas and write 30K probes from datProbes_30491.csv
  dat0 = data.frame(probecorrected)
  probes = data.frame(ProbeID=rownames(dat0))
  dat1=cbind(probes,dat0)
  
  if(returnAll){
    output = list(pipeline=pipeline,probes = dat1)
    return(output)
  }
  else{
    dat2=dat1[dat1$ProbeID %in% selectCpG,]
    output = list(pipeline=pipeline,probes = dat2)
    return(output)
  }
}
