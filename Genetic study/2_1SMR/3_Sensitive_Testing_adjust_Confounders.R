
setwd('/Users/bochaolin/Documents/PE/onesample/')
rm(list = ls())

data<-read.csv('./pheno/One_Sample_pheno_130363.csv')[,2:83]

###extracting signfiicant variables 
FID<-c('20529','2040','20531')

vars_names<-colnames(data)[3:38]
vars_names_s<-gsub("\\_(_*?)(?=f\\d)"," ",vars_names,perl=T)
fid<-stringr::str_split_fixed(vars_names_s,' f',2)[,2]
vars_names_Q<-vars_names[c(18,21,27,31)]; fid_Q<- fid[c(18,21,27,31)]
vars_names_B<-vars_names[-c(18,21,27,31)];fid_B<- fid[-c(18,21,27,31)]
covs_names<-colnames(data)[39:82]

#factors:
table(data$victim_of_sexual_assault_f20531, data$PE)
table(data$risk_taking_f2040,data$PE)
data[data=='Yes']<-1
data[data=='No']<-0
for( i in 1:ncol(data)){
  data[,i]<-as.numeric(data[,i])
}

data$uk_biobank_assessment_centre_f54_0_0<-as.factor(data$uk_biobank_assessment_centre_f54_0_0)
data$chip<-as.factor(data$chip)
data$sex_f31_0_0<-as.factor(data$sex_f31_0_0)
data$PE<-ifelse(data$PE=='Yes',0,1)
data_org<-data

system('ls ./score/')
forward_oneMR_B <- data.frame(Predictor=vars_names_B,R2=NA, Fstats=NA,beta=NA, SE=NA, Z=NA,P=NA)
  score<-read.table(paste('./score/',fid_B[i],'.profile',sep=''),header = T)[,c(1,6)]
 
  data<-data_org

  colnames(score)<-c('eid','SCORESUM')
  Confounder<-read.table(paste('./pheno/Conf_',fid_B[i],'.txt'),header = T)
  head(Confounder)
  Confounder<-Confounder[,-10]
  grep(fid_B[i],colnames(Confounder))
  
  score<-merge(score,Confounder,by='eid')
  data<-data[,c('eid','PE',vars_names_B[i],  covs_names)]
  data<-merge(data, score, by='eid',all=F)
  colnames(data)
  score_name<-colnames(score)[2:ncol(score)]
  ### First stage of regression 
  data_reg1<-data[,c(vars_names_B[i],  covs_names,  score_name)]
  colnames(data_reg1)[1]<-'exposure'
  data_reg1$exposure<-ifelse(data_reg1$exposure=='Yes',0,1)
  olsreg<- glm(exposure ~ .,data=data_reg1,family = 'binomial',na.action=na.exclude)
  R2<- rcompanion::nagelkerke(olsreg)
  data_null<-subset(data_reg1,select = -c(SCORESUM))
  olsnull<-glm(exposure ~ .,data=data_null,family = 'binomial',na.action=na.exclude)
  R2.null<- rcompanion::nagelkerke(olsnull)
  data$Y2hat<- scale(fitted(olsreg))
  olsreg.sum<-summary(olsreg)
  forward_oneMR_B$R2[i]<- R2$Pseudo.R.squared.for.model.vs.null[3]*100-R2.null$Pseudo.R.squared.for.model.vs.null[3]*100
  forward_oneMR_B$Fstats[i]<- (olsreg.sum$coefficients['SCORESUM',1]^2)/(olsreg.sum$coefficients['SCORESUM',2]^2)
  ####second stage of regression 
  data_reg2<- data[,c('PE', covs_names,  score_name[2:length(score_name)],'Y2hat')]
  olsreg2 <- glm(PE~., data=data_reg2, family = 'binomial')
  results<-summary(olsreg2)
  
  forward_oneMR_B[i,4:7]<-results$coefficients['Y2hat',]

  
forward_oneMR_Q <- data.frame(Predictor=vars_names_Q,R2=NA,Fstats=NA,beta=NA, SE=NA, Z=NA,P=NA)
i=grep(FID[2], fid_Q) 
  score<-read.table(paste('./score/',fid_Q[i],'.profile',sep=''),header = T)[,c(1,6)]
  data<-data_org
  data<-data[,c('eid','PE',vars_names_Q[i],  covs_names)]
  colnames(score)<-c('eid','SCORESUM')
  Confounder<-read.table(paste('./pheno/Conf_',fid_Q[i],'.txt'),header = T)
  head(Confounder)
 head(score)
  score<-merge(score,Confounder,by='eid')
  #score_reg<-subset(score,select = -c(eid))
  #cof_reg<- lm(SCORESUM ~., data=score_reg,na.action=na.exclude)
  #score$SCORESUM<-scale(fitted( cof_reg))
  #score<-subset(score,select = c(eid,SCORESUM ))
  data<-data[,c('eid','PE',vars_names_Q[i],  covs_names)]
  data<-merge(data, score, by='eid',all=F)
  
  score_name<-colnames(score)[2:ncol(score)]
  
  ### First stage of regression 
  data_reg1<-data[,c(vars_names_Q[i], covs_names, score_name)]
  colnames(data_reg1)[1]<-'exposure'
  olsreg<- lm(exposure ~ .,data=data_reg1,na.action=na.exclude)
  olsreg.sum<-summary(olsreg)
  forward_oneMR_Q$Fstats[i]<- (olsreg.sum$coefficients['SCORESUM',1]^2)/(olsreg.sum$coefficients['SCORESUM',2]^2)
  olsreg<- lm(exposure ~ .,data=data_reg1,na.action=na.exclude)
  data_null<-subset(data_reg1,select = -c(SCORESUM))
  olsnull.sum<- summary(lm(exposure ~ .,data=data_null,na.action=na.exclude))
  forward_oneMR_Q$R2[i]<- (olsreg.sum$r.squared- olsnull.sum$r.squared)*100
  data$Y2hat<- scale(fitted(olsreg))
  data_reg2<- data[,c("PE", covs_names, 'Y2hat')]
  olsreg2 <- glm(PE~., data=data_reg2, family = 'binomial')
  results<-summary(olsreg2)
  forward_oneMR_Q[i,4:7]<-results$coefficients['Y2hat',]
  
  forward_oneMR<-rbind(forward_oneMR_B,forward_oneMR_Q)
  forward_oneMR<-subset(forward_oneMR,forward_oneMR$P<1)
  table(forward_oneMR$P<0.05/36)
  
  
  forward_oneMR$OR<- exp(forward_oneMR$beta)
  forward_oneMR$lower.CI<-  exp(forward_oneMR$beta-1.96*forward_oneMR$SE)
  forward_oneMR$upper.CI<-  exp(forward_oneMR$beta+1.96*forward_oneMR$SE)
  
  head(forward_oneMR)
  forward_oneMR<-subset(forward_oneMR,select = c(Predictor,R2,Fstats,OR,lower.CI,upper.CI,P ))
  WriteXLS::WriteXLS(forward_oneMR,'Stable_B_forward_1SMR_ajustedconfounder.xls')

