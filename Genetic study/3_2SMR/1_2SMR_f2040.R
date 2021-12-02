rm(list=ls(all.names = T))
setwd('/Users/bochaolin/Documents/UM/PE/MR/res/twosample/Sep/')
source('/Users/bochaolin/Documents/script/gsmr_plot.r')
library('TwoSampleMR')

## read sums stats for 
exp<-'f2040'
exposure_data<-read_exposure_data(
  paste('./sums',exp,'.ma.gz',sep=''),
  clump = FALSE,
  sep = " ",
  phenotype_col = "Pheno",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "N",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "POS"
)
out<-'PE_aD' ### for PE
#out<- 'SCZ'  ### for schizophrenia

outcome_data<-read_outcome_data(
  paste('/sums',out,'.ma.gz',sep=''),
  snps = NULL,
  sep = " ",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "N",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "BP"
)

exposure_data$exposure<-exp
outcome_data$outcome<-out

exposure_data_small<-clump_data(
  exposure_data_small,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-6,
  clump_p2 = 1,
  pop = "EUR"
)

###harmonise data

dat<-harmonise_data(exposure_data_small, outcome_data, action = 2)
dim(dat)

## get F statistic 

Fstats<- mean((dat$beta.exposure)^2/(dat$se.exposure)^2)

## conduct IVW, MW and Egger models:
mr_results <- mr(dat,method_list = c('mr_ivw_fe','mr_weighted_median','mr_egger_regression'))
mr_results

## get intercept of Egger model
egger<-mr_pleiotropy_test(dat)
egger

## get I2 for Egger model 
library(MendelianRandomization)
object<- subset(dat,select = c( beta.exposure,se.exposure, beta.outcome ,se.outcome ))
colnames(object)<-c('bx','bxse','by','byse')
object_in<-mr_input( bx=dat$beta.exposure, bxse=dat$se.exposure,by=dat$beta.outcome, byse= dat$se.outcome)
res<-mr_egger(object_in)
res

## get Q statitsics 
mr_heterogeneity(dat)

## conduct PRESSO model 
presso<-run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
presso

## leave one out analysis 
pdf(paste('Figure_LOO_',exp,'_',out,'.pdf',sep=''),width = 5,height = 15)
res_loo<-mr_leaveoneout(dat, parameters = default_parameters(), )
res_loo
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
dev.off()

mr_results$OR<-exp(mr_results$b)
mr_results$se<-as.numeric(mr_results$se)
mr_results$lower.ci<-exp(mr_results$b-1.96*mr_results$se )
mr_results$upper.ci<-exp(mr_results$b+1.96*mr_results$se)

mr<-subset(mr_results, select = c(nsnp,OR,lower.ci, upper.ci,pval))
round(mr_results[,2:4],3)

WriteXLS::WriteXLS(mr,paste(exp,out,'_forward.xls',sep=''))
