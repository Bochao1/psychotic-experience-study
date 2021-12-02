###Author: HSSarac 
###Date: 22.11.2021
###Script for detecting the non-genetic correlates of Psychosis. See UK Biobank open-public data showcase for the non-genetic variables
#df is prefered as a generic name for dataframe

#Split data and conduct XWAS script

#Randomly split the dataset (N = 155247, 247 variables, age and sex) into train (N=77624) and test (N = 77623) datasets
##Split it proportionately to the outcome measure(PE) 
set.seed(123)
discovery.index <- 
  createDataPartition(df$PE,p = 0.5, list = FALSE)
df_discovery <- df[discovery.index, ]
df_replication <- df[-discovery.index, ]

#XWAS using logistic regression on each 247 variable seperately.
#a) In the training dataset, analyze each predictor's association with PE with the covariates: age and sex.

results_discovery <- data.frame(Predictor=rownames(cols_predictors),Discovery_Pval=NA,Discovery_Coefficient=NA,Discovery_SD=NA,Discovery_Oddsratios=NA,Discovery_R2=NA)

for(i in 1:nrow(cols_predictors)){
  modelA <- summary(glm(PE ~ df_discovery[,i]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_discovery, family = binomial))
  results_discovery[i,2] <- modelA$coefficients[2,4] ##2.line 4.column
  results_discovery[i,3] <- modelA$coefficients[2,1]
  results_discovery[i,4] <- modelA$coefficients[2,2]
}

#To get the pseudo R2 for each model
R2_glms_discovery<-lapply(1:nrow(cols_predictors),function(x)nagelkerke(glm(PE ~ df_discovery[,x]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_discovery, family = binomial),null=glm(PE ~ 1, binomial, df_discovery[!is.na(df_discovery[,x]),])))
Nagelkerke_Cragg_Uhler_discovery<- lapply(R2_glms_discovery,function(x) x$Pseudo.R.squared.for.model.vs.null[3,1])
Nagelkerke_Cragg_Uhler_discovery<-unlist(Nagelkerke_Cragg_Uhler_discovery)
results_discovery$Discovery_R2<-Nagelkerke_Cragg_Uhler_discovery

results_discovery$Discovery_Oddsratios<- exp(results_discovery$Discovery_Coefficient)
results_discovery$Discovery_significance<-ifelse(results_discovery$Discovery_Pval<=0.05/nrow(cols_predictors),1,0)
results_discovery$Discovery_ORlci <- exp(results_discovery$Discovery_Coefficient + qnorm(c(0.025)) * results_discovery$Discovery_SD)
results_discovery$Discovery_ORhci <- exp(results_discovery$Discovery_Coefficient + qnorm(c(0.975)) * results_discovery$Discovery_SD)

results_discovery<-as.data.frame(results_discovery)
write.csv(results_discovery,file="results_discovery.csv",row.names = FALSE)

#The same procedure was applied for the replication dataset.