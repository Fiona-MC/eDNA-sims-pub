library(stats)

schliepXdata<-readRDS("/home/fiona_callahan/2sp_2cov_noInteractions/SchliepRes/trial1X/dataSchliep1Burn500_t100.Rdata")

logisticDatafromSchliepX<-function(schliepXdata){
  sitetab<-data.frame()
  for(t in 1:schliepXdata$t){
    for(loc_index in 1:schliepXdata$n){
      for(species in 1:schliepXdata$S){
        cov1<-schliepXdata$X[[t]][loc_index,2]
        cov2<-schliepXdata$X[[t]][loc_index,3]
        presence<-schliepXdata$Y[[t]][[loc_index]][1,species]
        sitetab<-rbind(sitetab, c(loc_index, t, cov1, cov2, species, presence))
      }
    }
  }
  names(sitetab)=c("location","time_index","Cov1", "Cov2", "Species", "Presence")
  return(sitetab)
}

sitetab_schliepXData<-logisticDatafromSchliepX(schliepXdata)

model <- glm(Presence ~ Cov1+Cov2, family=binomial(link='logit'), data=sitetab_schliepXData[sitetab_schliepXData$Species == 1,])
summary(model)
model <- glm(Presence ~ Cov1+Cov2, family=binomial(link='logit'), data=sitetab_schliepXData[sitetab_schliepXData$Species == 2,])
summary(model)