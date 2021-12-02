###Author: HSSarac 
###Date: 22.11.2021
###Script for detecting the non-genetic correlates of Psychosis. See UK Biobank open-public data showcase for the non-genetic variables
#df is prefered as a generic name for dataframe

# Multivariate model R script

# Sample: 148 replicated variables + age + sex, N = 57702 participants with full data on all 148 variables

# final_multivariatemodel<-summary(glm(PE ~., data = df, family = binomial))
coeffinal_mul<-final_multivariatemodel$coefficients

rownames(coeffinal_mul)<-lapply(rownames(coeffinal_mul),function(x)gsub("[_]"," ",x))

final_multivariateresults <- data.frame(Predictor=rownames(coeffinal_mul),Final_Pval=NA,Final_Coefficient=NA,Final_SD=NA,Final_Oddsratios=NA)

#Evaluating the nagelkerke R2 of the final model. 
Final_R2<-nagelkerke(glm(PE ~., data = df, family = binomial))
R2_final=Final_R2$Pseudo.R.squared.for.model.vs.null[3,1]
Final_Nagelkerke_R2<-as.data.frame(R2_final)
write.csv(Final_Nagelkerke_R2,file="final_model_nagelkerke_pseudo_R2.csv",row.names = FALSE)

for(i in 1:nrow(final_multivariateresults)){
  final_multivariateresults[i,2] <- final_multivariatemodel$coefficients[i,4] ##line i,4.column
  final_multivariateresults[i,3] <- final_multivariatemodel$coefficients[i,1]
  final_multivariateresults[i,4] <- final_multivariatemodel$coefficients[i,2]
}
final_multivariateresults$Final_Oddsratios<- exp(final_multivariateresults$Final_Coefficient)
final_multivariateresults$Final_log10Pval<--log10(final_multivariateresults$Final_Pval)
final_multivariateresults$Final_significance<-ifelse(final_multivariateresults$Final_Pval<=(0.05),1,0)
final_multivariateresults$Final_ORlci <- exp(final_multivariateresults$Final_Coefficient + qnorm(c(0.025)) * final_multivariateresults$Final_SD)
final_multivariateresults$Final_ORhci <- exp(final_multivariateresults$Final_Coefficient + qnorm(c(0.975)) * final_multivariateresults$Final_SD)

write.csv(final_multivariateresults,file="final_multivariateresults.csv",row.names=FALSE)
