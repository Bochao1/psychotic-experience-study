setwd('/Users/bochaolin/Documents/PE/MR/onesample/')
rm(list = ls())
PE<-read.csv('./pheno/PE_and_covariates.csv')
system('ls ./pheno/*')
head(pheno)

#---------------------------reverse One_Sample MR legge 
data<-read.csv('One_Sample_pheno_130363.csv')[,2:83]
head(data)

vars_names<-colnames(data)[3:38]
vars_names_Q<-vars_names[c(18,21,27,31)]
vars_names_B<-vars_names[-c(18,21,27,31)]
covs_names<-colnames(data)[39:82]

#factors:

data$uk_biobank_assessment_centre_f54_0_0<-as.factor(data$uk_biobank_assessment_centre_f54_0_0)
data$chip<-as.factor(data$chip)
data$sex_f31_0_0<-as.factor(data$sex_f31_0_0)

system('ls ./score/')
score<-read.table('./score/PE_legge.profile',header = T)[,c(1,6)]
head(score)
data$PE<-ifelse(data$PE=='Yes',0,1)
data<-merge(data, score, by.x='eid',by.y='FID',all=F)

### First stage of regression 
data_reg1<-data[,c('PE', covs_names, 'SCORESUM')]
olsreg<- glm( PE~.,data=data_reg1,family = 'binomial')
R2<- rcompanion::nagelkerke(olsreg)
data$Y2hat<- scale(fitted(olsreg))
summary(data$Y2hat)
####second stage of regression 

reverse_oneMR_B <- data.frame(Predictor=vars_names_B,beta=NA, SE=NA, Z=NA,P=NA)
for(i in 1:length(vars_names_B)){
  data[,vars_names_B[i]]<-ifelse(data[,vars_names_B[i]]=='Yes',1,0)
  data_reg2<- data[,c(vars_names_B[i], covs_names, 'Y2hat')]
  colnames(data_reg2)[1]<-'outcome'
  olsreg2 <- glm(outcome~., data=data_reg2, family = 'binomial')
  results<-summary(olsreg2)
  reverse_oneMR_B[i,2:5]<-results$coefficients['Y2hat',]
}

reverse_oneMR_Q <- data.frame(Predictor=vars_names_Q,beta=NA, SE=NA, Z=NA,P=NA)
for(i in 1:length(vars_names_Q)){
  data_reg2<- data[,c(vars_names_B[i], covs_names, 'Y2hat')]
  colnames(data_reg2)[1]<-'outcome'
  olsreg2 <- lm(outcome~., data=data_reg2)
  results<-summary(olsreg2)
  reverse_oneMR_Q[i,2:5]<-results$coefficients['Y2hat',]
}

reverse_oneMR<-rbind(reverse_oneMR_B,reverse_oneMR_Q)
table(reverse_oneMR$P<0.05/36)

reverse_oneMR$OR<- exp(reverse_oneMR$beta)
reverse_oneMR$lower.CI<-  exp(reverse_oneMR$beta-1.96*reverse_oneMR$SE)
reverse_oneMR$upper.CI<-  exp(reverse_oneMR$beta+1.96*reverse_oneMR$SE)

reverse_oneMR<-subset(reverse_oneMR,select = -c(beta,SE,Z))
WriteXLS::WriteXLS(reverse_oneMR,'Stable_B_reverse_1SMR_legge.xls')


#---------------------------reverse One_Sample MR legge 
data<-read.csv('One_Sample_pheno_130363.csv')[,2:83]
head(data)

vars_names<-colnames(data)[3:38]
vars_names_Q<-vars_names[c(18,21,27,31)]
vars_names_B<-vars_names[-c(18,21,27,31)]
covs_names<-colnames(data)[39:82]

system('ls ./score/')
score<-read.table('./score/PE_ukb.raw',header = T)[,c(1,7)]
head(score)
colnames(score)[2]<-'SCORESUM'
data$PE<-ifelse(data$PE=='Yes',0,1)
data<-merge(data, score, by.x='eid',by.y='FID',all=F)

### First stage of regression 
data_reg1<-data[,c('PE', covs_names, 'SCORESUM')]
data_reg1<-tidyr::complete(data_reg1)
olsreg<- glm( PE~.,data=data_reg1,family = 'binomial',na.action=na.exclude)
R2<- rcompanion::nagelkerke(olsreg)
head(data_reg1)
data$Y2hat<- scale(fitted(olsreg))
summary(data$Y2hat)
####second stage of regression 

reverse_oneMR_B <- data.frame(Predictor=vars_names_B,beta=NA, SE=NA, Z=NA,P=NA)
for(i in 1:length(vars_names_B)){
  data[,vars_names_B[i]]<-ifelse(data[,vars_names_B[i]]=='Yes',1,0)
  data_reg2<- data[,c(vars_names_B[i], covs_names, 'Y2hat')]
  colnames(data_reg2)[1]<-'outcome'
  olsreg2 <- glm(outcome~., data=data_reg2, family = 'binomial')
  results<-summary(olsreg2)
  reverse_oneMR_B[i,2:5]<-results$coefficients['Y2hat',]
}

reverse_oneMR_Q <- data.frame(Predictor=vars_names_Q,beta=NA, SE=NA, Z=NA,P=NA)
for(i in 1:length(vars_names_Q)){
  data_reg2<- data[,c(vars_names_B[i], covs_names, 'Y2hat')]
  colnames(data_reg2)[1]<-'outcome'
  olsreg2 <- lm(outcome~., data=data_reg2)
  results<-summary(olsreg2)
  reverse_oneMR_Q[i,2:5]<-results$coefficients['Y2hat',]
}

reverse_oneMR<-rbind(reverse_oneMR_B,reverse_oneMR_Q)
table(reverse_oneMR$P<0.05/36)

reverse_oneMR$OR<- exp(reverse_oneMR$beta)
reverse_oneMR$lower.CI<-  exp(reverse_oneMR$beta-1.96*reverse_oneMR$SE)
reverse_oneMR$upper.CI<-  exp(reverse_oneMR$beta+1.96*reverse_oneMR$SE)

reverse_oneMR<-subset(reverse_oneMR,select = -c(beta,SE,Z))
WriteXLS::WriteXLS(reverse_oneMR,'Stable_C_reverse_1SMR_PEUKB.xls')



