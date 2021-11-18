DK<-read.csv("DKWA.csv",header=T)
attach(DK)

library(lme4)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(lattice)
library(visreg)

View(DK)

##scale n center data due to error in reduced model
DK.s<-DK
DK.s[,names(DK) %in% c("JulDaySurv","VegRemoved","Visibility")]<-scale(DK[,names(DK) %in% c("JulDaySurv","VegRemoved","Visibility")])

## summarize
library(tidyverse)
DK2<-DK %>%
  group_by(Site,DepthCat,SurvCat,Type) %>%
  summarize(mean=mean(ALL), SD=sd(ALL), n=n())

library(data.table)

setDT(DK2)
DK.sum<-dcast(DK2, Site+DepthCat+SurvCat ~ Type, value.var = "mean")
(DK.sum$TCRatio<-DK.sum$Treat/DK.sum$Control)
mean(DK.sum$TCRatio)
sd(DK.sum$TCRatio)

shapiro.test(DK.sum$TCRatio)
library("ggpubr")
ggqqplot(DK.sum$TCRatio, ylab = "Treat:Control Ratio",
         ggtheme = theme_minimal())
#both suggest at least one outlier, that is removed. analyses repeated [-3]
shapiro.test(DK.sum$TCRatio[-3])

ggqqplot(DK.sum$TCRatio[-3], ylab = "Treat:Control Ratio",
         ggtheme = theme_minimal())
#mean and SD of T:C Ratio with outlier removed
mean(DK.sum$TCRatio[-3])
sd(DK.sum$TCRatio[-3])

#one sample t-test to show if H0=1
t.test(DK.sum$TCRatio[-3], mu = 1)

#t-test to compare ratios between 11 and 22 ft depths. outlier removed
var.test(TCRatio ~ DepthCat, data = DK.sum[-3,])
t.test(DK.sum[DepthCat=="Eleven"]$TCRatio, DK.sum[DepthCat=="TwentyTwo"][-1]$TCRatio, var.equal = FALSE)
mean(DK.sum[DepthCat=="Eleven"]$TCRatio)
sd(DK.sum[DepthCat=="Eleven"]$TCRatio)
mean(DK.sum[DepthCat=="TwentyTwo"]$TCRatio[-1])
sd(DK.sum[DepthCat=="TwentyTwo"]$TCRatio[-1])

t.test(DK.sum[DepthCat=="Eleven"]$TCRatio, mu = 1)

t.test(DK.sum[DepthCat=="TwentyTwo"]$TCRatio[-1], mu = 1)

## Variable Correlations
corr<-cor(DK[sapply(DK,is.numeric)])
(corr.tbl<-as.data.frame(ifelse(abs(corr[c(1:8),c(1:8)])>0.6,round(corr[c(1:8),c(1:8)],3),"*")))

##Exploration of 11 ft and survey 1 only
DK.11.1<-subset(DK.s, DepthCat=="Eleven" & SurvNum==1)

m1<-lm(log(ALL+1)~Type,data=DK.11.1)
m1<-glmer(ALL~Type+(1|Site), data=DK.11.1, family=poisson,na.action="na.fail")
summary(m1)
visreg(m1,"Type",type='conditional', band=TRUE)
visreg(m1,"Type",type='contrast')
(m1.rand.eff<-dotplot(ranef(m1))$Site)

## Site Random Effect Necessary?
gm1<-glm(ALL~JulDaySurv+Type*DepthCat+Visibility, data=DK.s, family=poisson, na.action="na.fail")
gm1.site.reff<-glmer(ALL~JulDaySurv+Type*DepthCat+Visibility+(1|Site),
                     data=DK.s, family=poisson, na.action="na.fail")
gm2.site.surv.reff<-glmer(ALL~JulDaySurv+Type*DepthCat+Visibility+(1|Site)+(1|SurvCat),
                     data=DK.s, family=poisson, na.action="na.fail")
#gm3.site.surv_reff<-glmer(ALL~JulDaySurv+Type*DepthCat+Visibility+(1|Site/SurvCat),
                     #data=DK.s, family=poisson, na.action="na.fail")  #Is singular
AIC(gm1, gm1.site.reff, gm2.site.surv.reff,gm3.site.surv_reff)#AIC,  suggests random effect needed


## Global model
Global.mod<-glmer(ALL~JulDaySurv+Type*DepthCat+Visibility+(1|Site)+(1|SurvCat),
                data=DK.s,family=poisson, na.action="na.fail")  #Insert response
summary(Global.mod)
plot(Global.mod)

mod.sel<-dredge(Global.mod,rank="AICc", REML = FALSE, fixed=c("Type"))#, subset=!("Type" && "VegRemoved"))
(top.mods<-subset(mod.sel,subset = delta<4))
(avgm<-model.avg(top.mods))
summary(avgm)
#drop Julian Day of Survey because Survey Category captures that effect, and parameter estimate is imprecise
#drop Visibility because it is non-sensical (less fish with more visibility)

##Reduced model
r.mod<-glmer(ALL~Type*DepthCat+(1|Site)+(1|SurvCat), 
             data=DK.s, family=poisson,na.action="na.fail")  
summary(r.mod)
plot(r.mod)
## Plots of species/group vs. treatment

visreg(r.mod, "Type", by="DepthCat",
       type='conditional',ylab="log(Number of Fish)")

randoms <- ranef(r.mod)
(srvy.rand.eff<-dotplot(randoms)$SurvCat)
(site.rand.eff<-dotplot(randoms)$Site)
#rand.eff.fig<-ggarrange(srvy_site,site, labels = c("A)", "B)"),
                  #ncol = 1, nrow = 2, heights=c(1,1), common.legend = TRUE)
#rand.eff.fig
#ggsave("Random_effects.tiff", srvy_site.rand.eff, width = 10, height = 5.5)

###############SAME BUT WITH 22Ft treatment removed

DK.sub<-subset(DK,DepthCat=="Eleven")
View(DK.sub)

##Global model
Global.mod<-glmer(ALL~JulDaySurv+TimeElapsed+VegRemoved+Type+Visibility+Weather+Temp+(1|Site/SurvCat),
                  data=DK.sub,family=poisson,na.action="na.fail")  #Insert response
summary(Global.mod)

mod.sel<-dredge(Global.mod,rank="AICc",REML = FALSE,fixed=c("Type","(1|Site/SurvCat)"),subset=!("Type" && "VegRemoved"))
top.mods<-subset(mod.sel,subset = delta<4)
avgm<-model.avg(top.mods)
summary(avgm)

##Reduced model
r.mod<-glmer(ALL~Type+(1|Site/SurvCat),data=DK.sub,family=poisson,na.action="na.fail")  
summary(r.mod)

## Plots of species/group vs. treatment
par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(3,3,1,1))

visreg(r.mod, "Type",type='conditional')

randoms <- ranef(r.mod)
tiff("randeffs.tiff",width=7,height=10,units='in',res=300)
par(mfrow=c(1,3),mar=c(2,2,2,2),oma=c(3,3,1,1))
(depthsurvsite.rand.eff<-dotplot(randoms)$`DepthCat:(SurvCat:Site)`)
(srvy_site.rand.eff<-dotplot(randoms)$`SurvCat:Site`)
(site.rand.eff<-dotplot(randoms)$Site)
dev.off()

