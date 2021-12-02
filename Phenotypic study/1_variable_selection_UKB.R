###Author: HSSarac 
###Date: 22.11.2021
###Script for detecting the non-genetic correlates of Psychosis. See UK Biobank open-public data showcase for the non-genetic variables
#df is prefered as a generic name for dataframe

#QC Steps Script

#1 Removing variables for Stable 1: Out of 23,091 variables, 16,171 have been removed from all participants, N ~ 500,000

df <- fread("Data_Dictionary_Showcase.csv", stringsAsFactors = F,data.table = F,header=T) 

#Keep a dataframe with all the variables for later to compare with the removed dataset.
all_columns<-rownames(df)

Strata_to_remove<- c("auxiliary","supporting") 
df<- df[-which(df$Strata %in% c(Strata_to_remove)),]
Itemtype_to_remove<-c("bulk","samples","records")
df<- df[-which(df$ItemType %in% c(Itemtype_to_remove)),]
Valuetype_to_remove<-c("compound","date","time","binary_object","records")
df<- df[-which(df$ValueType %in% c(Valuetype_to_remove)),]
Sexed_to_remove<-c("female","male")
df<- df[-which(df$Sexed %in% c(Sexed_to_remove)),]

#filter follow-up questions
df_followup<- df %>% 
  filter(str_detect(Notes, "except_those|asked_only|were_asked_if_the_participant_indicated|only_asked|as_defined_by_their_answer|was_only_asked_to_participants_who|question_was_asked_when"))
#remove these questions by filtering out the rows 
df<-df[!(rownames(df) %in% rownames(df_followup)), ]

#the list of variables to be removed.
remove<- all_columns [!(all_columns%in% rownames(df))] 

#remove the auxiliary, bulk, compound, follow-up, time, date and raw data variables.
df <- data.table::fread("ukb_data.csv",data.table = F,stringsAsFactors = F)

for(i in 1:nrow(remove)){
  removea<-remove[i,]
  df<-df[,!grepl(paste(removea,collapse="|"),names(df))]
}

#2 Calculating missingness and remove variables based on missingness. N = 155,247, only the participants who provided an valid answer for the Psychosis Expression (PE), out of 4678 variables, 4375 have been removed.

#missing rate calculation
missing_rate1_list_num <- sapply(df,function(x) sum(is.na(x)/nrow(df)))
dfmissing<-as.data.frame(missing_rate1_list_num)
dfmissing$missing_rate1_list_num<-dfmissing$missing_rate1_list_num*100

write.csv(dfmissing,file="missingratetable.csv")  

#Missing Rate Cut
df<-df %>%  purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=10)
remove(dfmissing,missing_rate1_list_num)
write.csv(df,file="variables_after_missing_cut.csv",row.names=FALSE)

#3 Calculate correlation and remove variables based on high correlation, N = 155247, out of 303 variables, 55 of them have been removed.

#Correlation analysis for detecting highly collinear variables
df_with_eid<-df

#We remove "eid" from the dataset for the collinearity analysis
df<-subset(df,select=-c(eid))

#We computed a heterogenous correlation matrix, consisting of Pearson product-moment correlations between numeric variables, polyserial correlations between numeric and ordinal variables, and polychoric correlations between ordinal variables.
correlation_model<-hetcor(df,std.err=F,use="pairwise.complete.obs",pd=F) 
correlation_matric<-correlation_model$correlations
df_correlation<-as.data.frame(correlation_matric)

write.csv(df_correlation,file="Correlation-matrice.csv")

#Find correlations higher than or equal to .90
correlation_high<-ifelse(df_correlation>=.90,1,0)

df_correlation_high<-as.data.frame(correlation_high)

high_cors<-sapply(df_correlation_high,function(x)length(which(x==1)))
high_cors_variables<-high_cors[high_cors>1]

#Find correlations lower than or equal to -.90
correlation_low<-ifelse(df_correlation<=-.90,1,0)
df_correlation_low<-as.data.frame(correlation_low)
low_cors<-sapply(df_correlation_low,function(x)length(which(x==1)))
high_neg_cors_variables<-low_cors[low_cors>0]

#Make a list of variables that are negatively and positively highly(<=-.90,>=.90) correlated
correlated_variables<-c(names(high_cors_variables),names(high_neg_cors_variables))
df_correlated_variables<-as.data.frame(correlated_variables)
corvars<-df_correlated_variables$correlated_variables

#Make a dataset of these highly correlated variables
df_high_cor=df[,grep(paste(corvars, collapse = "|"), names(df))]

remove(df,correlated_variables)

#Correlation matrix of highly correlated variables is saved. 
correlation_matrix_high_cors<-hetcor(df_high_cor,std.err=F,use="pairwise.complete.obs",pd=F)
highcors<-correlation_matrix_high_cors$correlations
highcors<-as.data.frame(highcors)
remove(df_high_cor,correlation_matrix_high_cors)
write_xlsx(highcors,"Supplementary_table_3_highcorvars.xlsx")
remove(highcors)

#Removal of a variable from each highly collinear pair through findCorrelation command that removes variables based on the idea to keep it collinearity of the whole dataset as low as possible. 
highly_collinear_variables_to_remove<-findCorrelation(correlation_matric,cutoff=0.9,names=TRUE)
df_with_eid=df_with_eid[,!grepl(paste(highly_collinear_variables_to_remove, collapse = "|"), names(df_with_eid))]

write.csv(df_with_eid,file="final_Prepared_Dataset.csv",row.names=FALSE)

#but give a clear of descripition of your input data (XX samples XX variables, have been recoded).