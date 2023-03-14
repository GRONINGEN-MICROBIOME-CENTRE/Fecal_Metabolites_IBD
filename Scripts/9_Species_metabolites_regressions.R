Phenotype association to metabolites
===

#load("~/Desktop/IBD_metabolomics_2022/3.Workspace/General_and_regressions.RData")
library(emmeans)
library(tidyverse)


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef",colnames(run_day))
my_phenos=controls_phenos[,c(regressors)]
my_phenos$clinical_ResectionAny=1
my_phenos=my_phenos[,c(1:7,28,8:27)]
controls_species2=merge(my_phenos,taxa3, by="row.names")
row.names(controls_species2)=controls_species2$Row.names
controls_species2$Row.names=NULL
q_mtbx2=q_mtbx
colnames(q_mtbx2)=make.names(colnames(q_mtbx2))
controls_species2=merge(controls_species2,q_mtbx2, by="row.names")
row.names(controls_species2)=controls_species2$Row.names
controls_species2$Row.names=NULL
controls_species2$cohort="aControl"
cd_specie2=cd_specie
cd_specie2$cohort="CD"
uc_specie2=uc_specie
uc_specie2$cohort="UC"

models_for_interac=bind_rows(controls_species2, uc_specie2, cd_specie2)
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day), "cohort3")

models_for_interac$Amount_sample_gram=scale(models_for_interac$Amount_sample_gram)[,1]
models_for_interac$metabolon_Month_in_freezer=scale(models_for_interac$metabolon_Month_in_freezer)[,1]
models_for_interac$host_Age=scale(models_for_interac$host_Age)[,1]
models_for_interac$host_BMI=scale(models_for_interac$host_BMI)[,1]
models_for_interac$clinical_BowelMovementADayDef=scale(models_for_interac$clinical_BowelMovementADayDef)[,1]


test_metabolites=colnames(models_for_interac)[138:999]
test_bugs=colnames(models_for_interac)[29:137]

####################

# Test relation between metabolite and bacteria presence/absence 
# Model 0: y ~ covariates + bacteria + dysbiosis
# Model 1: y ~ covariates + bacteria * dysbiosis

###################

dysbiosis_score3=dysbiosis_score_pheno[,c("dysbiotic"), drop=F]
models_for_interac=merge(dysbiosis_score3,models_for_interac,by="row.names")
row.names(models_for_interac)=models_for_interac$Row.names
models_for_interac$Row.names=NULL

flag=1
for ( i in test_bugs) {
  my_pheno=i
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in test_metabolites){
      my_trait=a
      my_uni_test=models_for_interac[,c(regressors,my_pheno,my_trait, "dysbiotic")]
      my_preds=c(regressors,my_pheno, "dysbiotic")
      interaction_factor=paste(my_pheno,"dysbiotic",sep="*")
      my_f0=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_lm0=summary(lm(my_f0,data = my_uni_test ))
      my_lm0_coef=as.data.frame(my_lm0$coefficients)
      my_f=as.formula(paste(my_trait, paste(paste(regressors, collapse = " + "), interaction_factor, sep="+"), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      my_lm=lm(my_f,data = my_uni_test )
      my_lm_sum=summary(my_lm)
      rg.nuis = ref_grid(my_lm, non.nuisance = c("dysbiotic", my_pheno))
      my_f2=as.formula(paste("pairwise", paste("dysbiotic",my_pheno, sep = "*"), sep = "~"))
      my_em=emmeans(rg.nuis, specs = my_f2, adjust="none")
      my_contrasts=as.data.frame(my_em$contrasts)
      my_lm_coef=as.data.frame(my_lm_sum$coefficients)
      my_first_result=data.frame(Metabolite=my_trait, 
                                 My_Bug=my_pheno, 
                                 combined_beta=my_lm0_coef[nrow(my_lm0_coef)-1,1],
                                 combined_pval=my_lm0_coef[nrow(my_lm0_coef)-1,4],
                                 beta_interaction=my_lm_coef[nrow(my_lm_coef),1],
                                 p_val_interaction=my_lm_coef[nrow(my_lm_coef),4],
                                 direction_eubiosis=my_contrasts[2,2]*-1,
                                 p_val_eubiosis=my_contrasts[2,6],
                                 direction_dysbiosis=my_contrasts[5,2]*-1,
                                 p_val_dysbiosis=my_contrasts[5,6])
      if (flag!=1){
        my_final_result=rbind(my_final_result,my_first_result)
      }else{
        my_final_result=my_first_result
        flag=5
      }
    }
  }
}

# Multiple testing correction and save
my_final_result$same_direction="no"
my_final_result$same_direction[my_final_result$Beta_bug>0 & my_final_result$direction_eubiosis>0 & my_final_result$direction_dysbiosis>0]="yes"
my_final_result$same_direction[my_final_result$Beta_bug<0 & my_final_result$direction_eubiosis<0 & my_final_result$direction_dysbiosis<0 ]="yes"
my_final_result$FDR_combined=p.adjust(my_final_result$p_val_bug, method = "BH")
my_final_result$FDR_eubiosis=p.adjust(my_final_result$p_val_eubiosis, method = "BH")
my_final_result$FDR_dysbiosis=p.adjust(my_final_result$p_val_dysbiosis, method = "BH")
my_final_result$FDR_interaction=p.adjust(my_final_result$p_val_interaction, method = "BH")
write.table(my_final_result, "~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/Bug-metabolite_interaction_bi.txt", sep = "\t", quote = F)



####################

# Test relation between metabolite and bacteria abundances (zeroes are coded as NA's)
# Model 0: y ~ covariates + bacteria + dysbiosis
# Model 1: y ~ covariates + bacteria * dysbiosis

###################


controls_species_lin2=controls_species_lin
controls_species_lin2$clinical_ResectionAny=1
controls_species_lin2=controls_species_lin2[,c(1:7,979,8:978)]
controls_species_lin2$cohort="aControl"
cd_specie_lin2=cd_specie_lin
cd_specie_lin2$cohort="CD"
uc_specie_lin2=uc_specie_lin
uc_specie_lin2$cohort="UC"
models_for_interac_lin=bind_rows(controls_species_lin2, uc_specie_lin2, cd_specie_lin2)
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "cohort")

models_for_interac_lin$Amount_sample_gram=scale(models_for_interac_lin$Amount_sample_gram)[,1]
models_for_interac_lin$metabolon_Month_in_freezer=scale(models_for_interac_lin$metabolon_Month_in_freezer)[,1]
models_for_interac_lin$host_Age=scale(models_for_interac_lin$host_Age)[,1]
models_for_interac_lin$host_BMI=scale(models_for_interac_lin$host_BMI)[,1]
models_for_interac_lin$clinical_BowelMovementADayDef=scale(models_for_interac_lin$clinical_BowelMovementADayDef)[,1]
models_for_interac_lin=merge(dysbiosis_score3,models_for_interac_lin,by="row.names")
row.names(models_for_interac_lin)=models_for_interac_lin$Row.names
models_for_interac_lin$Row.names=NULL


flag=1
for ( i in test_bugs) {
  my_pheno=i
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in test_metabolites){
      my_trait=a
      my_uni_test=models_for_interac_lin[,c(regressors,my_pheno,my_trait, "dysbiotic")]
      my_preds=c(regressors,my_pheno, "dysbiotic")
      interaction_factor=paste(my_pheno,"dysbiotic",sep="*")
      my_f0=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
      my_lm0=summary(lm(my_f0,data = my_uni_test ))
      my_lm0_coef=as.data.frame(my_lm0$coefficients)
      my_f=as.formula(paste(my_trait, paste(paste(regressors, collapse = " + "), interaction_factor, sep="+"), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      my_lm=lm(my_f,data = my_uni_test )
      my_lm_sum=summary(my_lm)
      rg.nuis = ref_grid(my_lm, non.nuisance = c("dysbiotic", my_pheno))
      my_em=test(emtrends(my_lm, pairwise ~ dysbiotic, var = my_pheno, adjust="none"))
      my_trends=as.data.frame(my_em$emtrends)
      my_lm_coef=as.data.frame(my_lm_sum$coefficients)
      #Warning, check if betas and pvalues are extracted corrected before running
      my_first_result=data.frame(Metabolite=my_trait, 
                                 My_Bug=my_pheno, 
                                 combined_beta=my_lm0_coef[nrow(my_lm0_coef)-1,1],
                                 combined_pval=my_lm0_coef[nrow(my_lm0_coef)-1,4],
                                 beta_interaction=my_lm_coef[nrow(my_lm_coef),1],
                                 p_val_interaction=my_lm_coef[nrow(my_lm_coef),4],
                                 direction_eubiosis=my_trends[1,2],
                                 p_val_eubiosis=my_trends[1,6],
                                 direction_dysbiosis=my_trends[2,2],
                                 p_val_dysbiosis=my_trends[2,6])
      if (flag!=1){
        my_final_result_lin=rbind(my_final_result_lin,my_first_result)
      }else{
        my_final_result_lin=my_first_result
        flag=5
      }
    }
  }
}

my_final_result_lin$same_direction="no"
my_final_result_lin$same_direction[my_final_result_lin$Beta_bug>0 & my_final_result_lin$direction_eubiosis>0 & my_final_result_lin$direction_dysbiosis>0]="yes"
my_final_result_lin$same_direction[my_final_result_lin$Beta_bug<0 & my_final_result_lin$direction_eubiosis<0 & my_final_result_lin$direction_dysbiosis<0 ]="yes"
my_final_result_lin$FDR_combined=p.adjust(my_final_result_lin$p_val_bug, method = "BH")
my_final_result_lin$FDR_eubiosis=p.adjust(my_final_result_lin$p_val_eubiosis, method = "BH")
my_final_result_lin$FDR_dysbiosis=p.adjust(my_final_result_lin$p_val_dysbiosis, method = "BH")
my_final_result_lin$FDR_interaction=p.adjust(my_final_result_lin$p_val_interaction, method = "BH")
write.table(my_final_result_lin, "~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/Bug-metabolite_interaction_li.txt", sep = "\t", quote = F)



####################

# Test relation between metabolite and metabolic gene clusters (zeroes are coded as NA's)
# Model 0: y ~ covariates + MGS + dysbiosis
# Model 1: y ~ covariates + MGS * dysbiosis

###################

#Remove bugs from the input file in the bug-metabolite association
models_for_interac_mgc=models_for_interac[,c(1:29,139:1001)]
mgc_fil_trans2=bgc_rpkm4
mgc_fil_trans2[mgc_fil_trans2==0]=NA
mgc_fil_trans3=transform_and_filter_mtb(mgc_fil_trans2, samples_row = T, method = "clr",missing_filter = 20)
models_for_interac_mgc=merge(models_for_interac_mgc, mgc_fil_trans3, by="row.names")
row.names(models_for_interac_mgc)=models_for_interac_mgc$Row.names
models_for_interac_mgc$Row.names=NULL
run_day_batch=cc_pheno6[,"run_day_cat", drop=F]
models_for_interac_mgc=merge(run_day_batch,models_for_interac_mgc, by="row.names")
row.names(models_for_interac_mgc)=models_for_interac_mgc$Row.names
models_for_interac_mgc$Row.names=NULL
test_metabolites=colnames(q_mtbx2)
test_bugs=colnames(mgc_fil_trans3)
total=ncol(mgc_fil_trans3)
#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day), "cohort")
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny","cohort" ,"run_day_cat")

flag=1
count=0
for ( i in test_bugs) {
  count=count+1
  my_pheno=i
  print (paste(paste(paste("Testing:",count), "/", total)))
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in test_metabolites){
      my_trait=a
      my_uni_test=models_for_interac_mgc[,c(regressors,my_pheno,my_trait, "dysbiotic")]
      my_preds=c(regressors,"dysbiotic",my_pheno)
      my_f0=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
      my_lm0=summary(lm(my_f0,data = my_uni_test ))
      my_lm0_coef=as.data.frame(my_lm0$coefficients)
      interaction_factor=paste(my_pheno,"dysbiotic",sep="*")
      my_preds=c(regressors)
      my_f=as.formula(paste(my_trait, paste(paste(my_preds, collapse = " + "), interaction_factor, sep="+"), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      my_lm=lm(my_f,data = my_uni_test )
      my_lm_sum=summary(my_lm)
      my_em=test(emtrends(my_lm, pairwise ~ dysbiotic, var = my_pheno, adjust="none"))
      my_trends=as.data.frame(my_em$emtrends)
      my_lm_coef=as.data.frame(my_lm_sum$coefficients)
       #Warning, check if betas and pvalues are extracted corrected before running
      my_first_result=data.frame(Metabolite=my_trait, 
                                 My_Bug=my_pheno, 
                                 combined_beta=my_lm0_coef[nrow(my_lm0_coef)-1,1],
                                 combined_pval=my_lm0_coef[nrow(my_lm0_coef)-1,4],
                                 beta_interaction=my_lm_coef[nrow(my_lm_coef),1],
                                 p_val_interaction=my_lm_coef[nrow(my_lm_coef),4],
                                 direction_eubiosis=my_trends[1,2],
                                 p_val_eubiosis=my_trends[1,6],
                                 direction_dysbiosis=my_trends[2,2],
                                 p_val_dysbiosis=my_trends[2,6])

      if (flag!=1){
        my_final_result_mgc=rbind(my_final_result_mgc,my_first_result)
      }else{
        my_final_result_mgc=my_first_result
        flag=5
      }
    }
  }
}


my_final_result_mgc$same_direction="no"
my_final_result_mgc$same_direction[my_final_result_mgc$Beta_bug>0 & my_final_result_mgc$direction_eubiosis>0 & my_final_result_mgc$direction_dysbiosis>0]="yes"
my_final_result_mgc$same_direction[my_final_result_mgc$Beta_bug<0 & my_final_result_mgc$direction_eubiosis<0 & my_final_result_mgc$direction_dysbiosis<0 ]="yes"
my_final_result_mgc$FDR_combined=p.adjust(my_final_result_mgc$p_val_bu, method = "BH")
my_final_result_mgc$FDR_eubiosis=p.adjust(my_final_result_mgc$p_val_eubiosis, method = "BH")
my_final_result_mgc$FDR_dysbiosis=p.adjust(my_final_result_mgc$p_val_dysbiosis, method = "BH")
my_final_result_mgc$FDR_interaction=p.adjust(my_final_result_mgc$p_val_interaction, method = "BH")


write.table(my_final_result_mgc, "~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/Bug-metabolite_dysbiosis_interaction_mgc_final2.txt", sep = "\t", quote = F)


####################

# Test relation between metabolite and METACYC pathways (zeroes are coded as NA's)
# Model 0: y ~ covariates + PATH + dysbiosis
# Model 1: y ~ covariates + PATH * dysbiosis

###################

pathways_raw <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/Metacyc/Merged_metacyc_uniref90_filtered_renames.txt", row.names = 1)
pathways_samples=subset(pathways_raw,rownames(pathways_raw) %in% rownames(cc_pheno2))
pathways_samples[pathways_samples==0]=NA
pathways_filt=transform_and_filter_mtb(pathways_samples,samples_row = T,method = "clr",missing_filter = 20)
pathways_filt$UNMAPPED=NULL
pathways_filt$UNINTEGRATED=NULL
models_for_interac_pwy=models_for_interac[,c(1:29,139:1001)]
models_for_interac_pwy=merge(models_for_interac_pwy, pathways_filt, by="row.names")
row.names(models_for_interac_pwy)=models_for_interac_pwy$Row.names
models_for_interac_pwy$Row.names=NULL
run_day_batch=cc_pheno6[,"run_day_cat", drop=F]
models_for_interac_pwy=merge(run_day_batch,models_for_interac_pwy, by="row.names")
row.names(models_for_interac_pwy)=models_for_interac_pwy$Row.names
models_for_interac_pwy$Row.names=NULL

test_metabolites=colnames(q_mtbx2)
test_bugs=colnames(pathways_filt)
total=ncol(pathways_filt)
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "cohort","run_day_cat")



flag=1
count=0
for ( i in test_bugs) {
  count=count+1
  my_pheno=i
  print (paste(paste(paste("Testing:",count), "/", total)))
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in test_metabolites){
      my_trait=a
      my_uni_test=models_for_interac_pwy[,c(regressors,my_pheno,my_trait, "dysbiotic")]
      my_preds=c(regressors,"dysbiotic",my_pheno)
      my_f0=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
      my_lm0=summary(lm(my_f0,data = my_uni_test ))
      my_lm0_coef=as.data.frame(my_lm0$coefficients)
      interaction_factor=paste(my_pheno,"dysbiotic",sep="*")
      my_preds=c(regressors)
      my_f=as.formula(paste(my_trait, paste(paste(my_preds, collapse = " + "), interaction_factor, sep="+"), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      my_lm=lm(my_f,data = my_uni_test )
      my_lm_sum=summary(my_lm)
      my_em=test(emtrends(my_lm, pairwise ~ dysbiotic, var = my_pheno, adjust="none"))
      my_trends=as.data.frame(my_em$emtrends)
      my_lm_coef=as.data.frame(my_lm_sum$coefficients)
      #Warning, check if betas and pvalues are extracted corrected before running
      my_first_result=data.frame(Metabolite=my_trait, 
                                 My_Bug=my_pheno, 
                                 combined_beta=my_lm0_coef[nrow(my_lm0_coef)-1,1],
                                 combined_pval=my_lm0_coef[nrow(my_lm0_coef)-1,4],
                                 beta_interaction=my_lm_coef[nrow(my_lm_coef),1],
                                 p_val_interaction=my_lm_coef[nrow(my_lm_coef),4],
                                 direction_eubiosis=my_trends[1,2],
                                 p_val_eubiosis=my_trends[1,6],
                                 direction_dysbiosis=my_trends[2,2],
                                 p_val_dysbiosis=my_trends[2,6])
      if (flag!=1){
        my_final_result_pwy=rbind(my_final_result_pwy,my_first_result)
      }else{
        my_final_result_pwy=my_first_result
        flag=5
      }
    }
  }
}


my_final_result_pwy2=my_final_result_pwy
my_final_result_pwy2$same_direction="no"
my_final_result_pwy2$same_direction[my_final_result_pwy2$Beta_bug>0 & my_final_result_pwy2$direction_eubiosis>0 & my_final_result_pwy2$direction_dysbiosis>0]="yes"
my_final_result_pwy2$same_direction[my_final_result_pwy2$Beta_bug<0 & my_final_result_pwy2$direction_eubiosis<0 & my_final_result_pwy2$direction_dysbiosis<0 ]="yes"
my_final_result_pwy2$FDR_combined=p.adjust(my_final_result_pwy2$p_val_bug, method = "BH")
my_final_result_pwy2$FDR_eubiosis=p.adjust(my_final_result_pwy2$p_val_eubiosis, method = "BH")
my_final_result_pwy2$FDR_dysbiosis=p.adjust(my_final_result_pwy2$p_val_dysbiosis, method = "BH")
my_final_result_pwy2$FDR_interaction=p.adjust(my_final_result_pwy2$p_val_interaction, method = "BH")
write.table(my_final_result_pwy2, "~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/Bug-metabolite_dysbiosis_interaction_pwy.txt", sep = "\t", quote = F)
