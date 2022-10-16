rm(list = ls())

library(TwoSampleMR)

IQ <-extract_instruments(outcomes='ebi-a-GCST006250',access_token = NULL) 
dim(IQ) 

pd <- extract_outcome_data(
  snps=IQ$SNP,
  outcomes='ieu-b-7',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL)

dim(pd)
mydata <- harmonise_data(exposure_dat=IQ,outcome_dat=pd,action= 2)
res <- mr(mydata)
res
res2=generate_odds_ratios(res)
res2
het <- mr_heterogeneity(mydata)
het
mr(mydata,method_list=c('mr_ivw_mre'))
pleio <- mr_pleiotropy_test(mydata)
pleio
single <- mr_leaveoneout(mydata)
single<-single[!duplicated(single$SNP),]
mr_leaveoneout_plot(single)

mr_scatter_plot(res,mydata)

res_single <- mr_singlesnp(mydata)
res_single<-res_single[!duplicated(res_single$SNP),]
mr_forest_plot(res_single)

mr_funnel_plot(res_single)


library(MRPRESSO)

data(mydata)
head(mydata)

mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata, NbDistribution = 1000,   
                SignifThreshold = 0.05)
