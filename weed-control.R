DK<-read.csv("DKWA.csv",header=T)
attach(DK)

library(lme4)
library(MuMIn)

Global.mod<-glmer(ALL~JulDaySurv+TimeElapsed+VegRemoved+Type+Visibility+Weather+Temp+(1|Site/SurvCat/DepthCat),family=poisson,na.action="na.fail")  #Insert response
summary(Global.mod)


mod.sel<-dredge(Global.mod,rank="AICc",REML = FALSE,fixed=c("Type","(1|Site/SurvCat/DepthCat)"),subset=!("Type" && "VegRemoved"))
(top.mods<-subset(mod.sel,subset = delta<4))
(avgm<-model.avg(top.mods))
summary(avgm)
