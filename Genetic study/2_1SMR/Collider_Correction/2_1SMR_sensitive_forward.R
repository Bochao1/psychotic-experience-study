# Overview
# This file contains the R script of the 1SMR sensitivity analyses (1SMRsensitivity1 (IVW), 1SMRsensitivity2 (MR-Egger), 1SMRsensitivity3 (LADreg))
# Here we take sleep duration and insomnia vs PE as examples  (R codes of other sleep traits are similar)
# To be noticed, these estimates are in log unit which haven't be transformed into SD unit 
# For PE, you need to divide the estimates by 0.15 log mmol/mol to obtain the SD estimates; for glucose, you need to divide the estimates by 0.17 log mmol/l to obtain the SD unit

library(foreign)
library(simex)
library(L1pack)
library(mr.raps)
library(coda)
library(xtable)
library(data.table)
rm(list = ls())

setwd('/Users/bochaolin/Documents/UM/PE/MR/res/onesample_NEW')
data<-read.csv('./pheno/One_Sample_pheno_130363.csv')[,2:83]
head(data)
vars_names<-colnames(data)[3:38]
vars_names_s<-gsub("\\_(_*?)(?=f\\d)"," ",vars_names,perl=T)
fid<-stringr::str_split_fixed(vars_names_s,' f',2)[,2]
FID<-c('20529','20531','2040')

vars_names<- vars_names[c(16,21,25,35)]
res<-data.frame(outcome=vars_names,Fstatistics=NA,Qexact=NA,Qp=NA, I2=NA, OR.WM=NA,lower.CI.WM=NA,upper.CI.WM=NA,P.WM=NA,
                OR.Egger=NA,lower.CI.Egger=NA,upper.CI.Egger=NA,P.Egger=NA,OR.LAD=NA,lower.CI.LAD=NA,upper.LAD=NA,P.LAD=NA)
#for (i in 1:length(vars_names)){


df <- read.csv(paste("./sums/", FID[3],"_PE.csv",sep = ''))
set.seed(888)
dim(df)
source("ExactQ.R")  #inport "ExactQ.R" function file

# function of generating the IGx2 
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

BetaXG<-df$BetaXG 
seBetaXG<-df$seBetaXG
BetaYG<-df$BetaYG 
seBetaYG<-df$seBetaYG
alphahatstar<-df$alphahatstar 
se.alphahatstar<-df$se.alphahatstar  
betastar<-df$betastar 
se.betastar<-df$se.betastar

##############################
F = BetaXG^2/seBetaXG^2
res$Fstatistics[var]<-mean(F)

# 1SMRsensitivity1
# Implement IVW using Modified weights 
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00000001)
IVW_mw  = Results$RESULTS[5,1]
seIVW= Results$RESULTS[5,2]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]
betaIVW = betastar[1] + IVW_mw
names(Results)

IsqGX   
#1SMRsensitivity2
# MR-Egger   
IsqGX         = Isq(BetaXG,seBetaXG)
res$I2[var]<-IsqGX

betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef    
MREgger       = betastar[1] + betahat[2,1]  
# Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                              
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    
Egger_intercept = bSIMEX[1,4]


# 1SMRsensitivity3
# LAD regression
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat  

# Collider correction + SiMEX 
Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,       
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX
# obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = df[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  # LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                              
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}

seLAD = sd(Ests) 

# sort the estimates
BetaXG          = df$BetaXG 
seBetaXG        = df$seBetaXG
alphahatstar    = df$alphahatstar
se.alphahatstar = df$se.alphahatstar

Fbar=mean(F)

#betaIVW<-seIVW<-Qexact<-Qp<-NA

Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)

ColliderCorrections = Estimates-betastar[1]
SEs                 = sqrt(c(seIVW,EggerSE,seLAD)^2 + (se.betastar[1])^2)
LCIs                =exp( Estimates - 1.96 * SEs)
UCIs                = exp(Estimates + 1.96 * SEs)
pval                = 2*(1-pnorm(abs(Estimates/SEs)))

ORs<-exp(Estimates)

Final = data.frame(ORs,LCIs ,UCIs  ,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))
Final 
vars_names[var]
res[var,6:9]<-Final[var,]; res[var,10:13]<-Final[2,]; res[var,14:17]<-Final[3,]
res[var,3]<- Qexact; res[var,4]<-Qp



WriteXLS::WriteXLS(res,'./res/Stable_8B_forward_CCed.xls')




