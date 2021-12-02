# Overview
# This file contains the R script of generating the summary statistics for the 1SMR sensitivity analyses (i.e., collider-corrected estimates)
library(dplyr)
library(data.table)
library(stringr)

# Continuous exposure (e.g., sleep duration)
# Summary statistics of genetic associations of with sleep duration and HbA1c in the UKB 
# (adjust for sex, age, assessment centre, chip, 40 genetic principal components)
rm(list = ls())

df <- read.csv('./pheno/One_Sample_pheno_130363.csv')[,2:83]

vars_names<-colnames(df)[3:38]
vars_names_Q<-vars_names[c(18,21,27,31)]
vars_names_B<-vars_names[-c(18,21,27,31)]
covs_names<-colnames(df)[39:82]

df$PE<-ifelse(df$PE=='Yes',1,0)

head(df)
vars<-'PE'
geno<-read.table(paste('./raw/',vars,'.raw',sep=''),header = T)
geno$eid<-geno$FID;geno<-geno[,-(1:6)]
genonames<-colnames(geno)
geno<-plyr::join(df, geno)
geno<- geno[,genonames]
geno<-geno[,1:(ncol(geno)-1)]; size<-ncol(geno)

n <- nrow(df)
colnames(df)[41]
G1.2 <- as.matrix(geno)
head(G1.2)
for (i in 1:length(vars_names_B)){
fid<-vars_names_B[i]
grep(fid,names(df))
names(df)[grep(fid,names(df))]
Y<- as.matrix(df[,grep(fid,names(df))])
print(head(Y))
Y<-ifelse(Y=='Yes',1,0)

X.2 <- as.matrix(df$PE) 

# covariates: sex, age, assessment centre, chip, 40 genetic principal components
sex<-as.matrix(df[,39])
age<-as.matrix(df[,40])
centre<-as.matrix(df$uk_biobank_assessment_centre_f54_0_0)

chip<-as.matrix(df$chip)

pc1 <- as.matrix(df[,41])     
pc2 <- as.matrix(df[,42]) 
pc3 <- as.matrix(df[,43]) 
pc4 <- as.matrix(df[,44]) 
pc5 <- as.matrix(df[,45]) 
pc6 <- as.matrix(df[,46]) 
pc7 <- as.matrix(df[,47]) 
pc8 <- as.matrix(df[,48]) 
pc9 <- as.matrix(df[,49]) 
pc10 <- as.matrix(df[,50]) 

pc11 <- as.matrix(df[,51])     
pc12 <- as.matrix(df[,52]) 
pc13 <- as.matrix(df[,53]) 
pc14 <- as.matrix(df[,54]) 
pc15 <- as.matrix(df[,55]) 
pc16 <- as.matrix(df[,56]) 
pc17 <- as.matrix(df[,57]) 
pc18 <- as.matrix(df[,58]) 
pc19 <- as.matrix(df[,59]) 
pc20 <- as.matrix(df[,60]) 

pc21 <- as.matrix(df[,61])     
pc22 <- as.matrix(df[,62]) 
pc23 <- as.matrix(df[,63]) 
pc24 <- as.matrix(df[,64]) 
pc25 <- as.matrix(df[,65]) 
pc26 <- as.matrix(df[,66]) 
pc27 <- as.matrix(df[,67]) 
pc28 <- as.matrix(df[,68]) 
pc29 <- as.matrix(df[,69]) 
pc30 <- as.matrix(df[,70]) 

pc31 <- as.matrix(df[,71])     
pc32 <- as.matrix(df[,72]) 
pc33 <- as.matrix(df[,73]) 
pc34 <- as.matrix(df[,74]) 
pc35 <- as.matrix(df[,75]) 
pc36 <- as.matrix(df[,76]) 
pc37 <- as.matrix(df[,77]) 
pc38 <- as.matrix(df[,78]) 
pc39 <- as.matrix(df[,79]) 
pc40 <- as.matrix(df[,80]) 

# X~G
XGdata2        = data.frame(X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   
head(XGdata2)
FIT2           = summary(glm(XGdata2, family = 'binomial') )   
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= size)                         
seBetaXG       = FIT2$coef[-1,2]                               
seBetaXG         = head(seBetaXG, n= size)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                


# Y~G 
YGdata        = data.frame(Y, G1.2,
                           sex, age, as.factor(centre), chip, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                           pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   


FIT3           = summary(glm(YGdata,family = 'binomial'))                          
BetaYG         = FIT3$coef[-1,1]                            
BetaYG         = head(BetaYG, n= size)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG         = head(seBetaYG, n= size)                          


# Y~G+X 
YXGdata        = data.frame(Y, X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40) 


FIT4             = summary(glm(YXGdata,family = 'binomial'))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar         = head(alphahatstar, n= size)                           
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar         = head(se.alphahatstar, n= size)                           
betastar         = FIT4$coef[2,1] 
se.betastar  =  FIT4$coef[2,2] 


#summary statistics
# X- G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")
write.csv(data, paste('sums/PE_',fid,".csv",sep=''))
}


for (i in 1:length(vars_names_Q)){
  fid<-vars_names_Q[i]
  grep(fid,names(df))
  names(df)[grep(fid,names(df))]
  Y<- as.matrix(df[,grep(fid,names(df))])
  X.2 <- as.matrix(df$PE) 
  
  # covariates: sex, age, assessment centre, chip, 40 genetic principal components
  sex<-as.matrix(df[,39])
  age<-as.matrix(df[,40])
  centre<-as.matrix(df$uk_biobank_assessment_centre_f54_0_0)
  
  chip<-as.matrix(df$chip)
  
  pc1 <- as.matrix(df[,41])     
  pc2 <- as.matrix(df[,42]) 
  pc3 <- as.matrix(df[,43]) 
  pc4 <- as.matrix(df[,44]) 
  pc5 <- as.matrix(df[,45]) 
  pc6 <- as.matrix(df[,46]) 
  pc7 <- as.matrix(df[,47]) 
  pc8 <- as.matrix(df[,48]) 
  pc9 <- as.matrix(df[,49]) 
  pc10 <- as.matrix(df[,50]) 
  
  pc11 <- as.matrix(df[,51])     
  pc12 <- as.matrix(df[,52]) 
  pc13 <- as.matrix(df[,53]) 
  pc14 <- as.matrix(df[,54]) 
  pc15 <- as.matrix(df[,55]) 
  pc16 <- as.matrix(df[,56]) 
  pc17 <- as.matrix(df[,57]) 
  pc18 <- as.matrix(df[,58]) 
  pc19 <- as.matrix(df[,59]) 
  pc20 <- as.matrix(df[,60]) 
  
  pc21 <- as.matrix(df[,61])     
  pc22 <- as.matrix(df[,62]) 
  pc23 <- as.matrix(df[,63]) 
  pc24 <- as.matrix(df[,64]) 
  pc25 <- as.matrix(df[,65]) 
  pc26 <- as.matrix(df[,66]) 
  pc27 <- as.matrix(df[,67]) 
  pc28 <- as.matrix(df[,68]) 
  pc29 <- as.matrix(df[,69]) 
  pc30 <- as.matrix(df[,70]) 
  
  pc31 <- as.matrix(df[,71])     
  pc32 <- as.matrix(df[,72]) 
  pc33 <- as.matrix(df[,73]) 
  pc34 <- as.matrix(df[,74]) 
  pc35 <- as.matrix(df[,75]) 
  pc36 <- as.matrix(df[,76]) 
  pc37 <- as.matrix(df[,77]) 
  pc38 <- as.matrix(df[,78]) 
  pc39 <- as.matrix(df[,79]) 
  pc40 <- as.matrix(df[,80]) 
  
  # X~G
  XGdata2        = data.frame(X.2, G1.2,
                              sex, age, as.factor(centre), chip, 
                              pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                              pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                              pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                              pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   
  head(XGdata2)
  FIT2           = summary(glm(XGdata2, family = 'binomial') )   
  BetaXG         = FIT2$coef[-1,1]                              
  BetaXG         = head(BetaXG, n= size)                         
  seBetaXG       = FIT2$coef[-1,2]                               
  seBetaXG         = head(seBetaXG, n= size)                           
  Fbar           = mean((BetaXG^2)/(seBetaXG^2))                
  
  
  # Y~G 
  YGdata        = data.frame(Y, G1.2,
                             sex, age, as.factor(centre), chip, 
                             pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                             pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                             pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                             pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   
  
  
  FIT3           = summary(lm(YGdata))                          
  BetaYG         = FIT3$coef[-1,1]                            
  BetaYG         = head(BetaYG, n= size)                           
  seBetaYG       = FIT3$coef[-1,2]                              
  seBetaYG         = head(seBetaYG, n= size)                          
  
  
  # Y~G+X 
  YXGdata        = data.frame(Y, X.2, G1.2,
                              sex, age, as.factor(centre), chip, 
                              pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                              pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                              pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                              pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40) 
  
  
  FIT4             = summary(lm(YXGdata))                       
  alphahatstar     = FIT4$coef[-c(1,2),1]                       
  alphahatstar         = head(alphahatstar, n= size)                           
  se.alphahatstar  = FIT4$coef[-c(1,2),2]
  se.alphahatstar         = head(se.alphahatstar, n= size)                           
  betastar         = FIT4$coef[2,1] 
  se.betastar  =  FIT4$coef[2,2] 
  
  
  #summary statistics
  # X- G
  BetaXG_data<-as.matrix(BetaXG)
  colnames(BetaXG_data) <- c("BetaXG")
  seBetaXG_data<-as.matrix(seBetaXG)
  colnames(seBetaXG_data) <- c("seBetaXG")
  
  #Y - G
  BetaYG_data<-as.matrix(BetaYG)
  colnames(BetaYG_data) <- c("BetaYG")
  seBetaYG_data<-as.matrix(seBetaYG)
  colnames(seBetaYG_data) <- c("seBetaYG")
  
  #Y - X - G
  alphahatstar_data<-as.matrix(alphahatstar)
  colnames(alphahatstar_data) <- c("alphahatstar")
  se.alphahatstar_data<-as.matrix(se.alphahatstar)
  colnames(se.alphahatstar_data) <- c("se.alphahatstar")
  
  betastar_data<-as.matrix(betastar)
  colnames(betastar_data) <- c("betastar")
  se.betastar_data<-as.matrix(se.betastar)
  colnames(se.betastar_data) <- c("se.betastar")
  
  data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
  
  SNP<-as.matrix(row.names(data))
  colnames(SNP) <- c("SNP_ukb")
  write.csv(data, paste('sums/PE_',fid,".csv",sep=''))
}

