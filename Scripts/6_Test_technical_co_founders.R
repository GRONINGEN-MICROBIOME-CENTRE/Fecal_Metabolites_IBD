Evaluating the effect of technical factors and host age, sex, BMI, Bowel movements per day
===

rownames(cc_pheno5)=cc_pheno5$Row.names
cc_pheno5$Row.names=NULL
cc_pheno6=cc_pheno5
table(cc_pheno6$run_day_cat, cc_pheno6$COMMENT)
cc_pheno6=subset(cc_pheno6,select = -c(Preparation_batch,Preparation_day,Preparation_order,CLIENT.IDENTIFIER,Box,Preparation_day_cat, RUN.DAY,SPL_AMNT_GRAM))

#Correction factors: day measurment was done, LC column, ammount sample per gram, month sample in a freezer, sex and age. 
my_correction_factors=cc_pheno6[,c("run_day_cat","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling", "host_BMI", "clinical_BowelMovementADayDef", "ibd_Diagnosis")]
my_correction_factors_test=merge(my_correction_factors,q_mtb,by="row.names")
rownames(my_correction_factors_test)=my_correction_factors_test$Row.names
my_correction_factors_test$Row.names=NULL
my_correction_factors_test=merge(my_correction_factors_test,SCFA3,by="row.names")
rownames(my_correction_factors_test)=my_correction_factors_test$Row.names
my_correction_factors_test$Row.names=NULL
colnames(my_correction_factors_test)=make.names(colnames(my_correction_factors_test))
my_correction_factors_test$cohort="aControl"
  
confounders_to_test=c("run_day_cat","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling", "host_BMI", "clinical_BowelMovementADayDef")
metabolites_to_test=make.names(c(colnames(q_mtb), colnames(SCFA3)))
flag=1
for (a in metabolites_to_test){
  test_metabolite=a
  for (i in confounders_to_test){
    test_trait=i
    test_factors=confounders_to_test[ !confounders_to_test == i]
    my_f=as.formula(paste(test_metabolite, paste(paste(paste(test_factors, collapse = " + "), "ibd_Diagnosis" , sep="+"), test_trait, sep="+"), sep = " ~ "))
    my_lm=summary(lm(my_f,data = my_correction_factors_test ))
    my_lm_coef=as.data.frame(my_lm$coefficients)
    confounder_result=data.frame(Metabolite=test_metabolite, My_Factor=test_trait,Beta=my_lm_coef[nrow(my_lm_coef),1], Pval=my_lm_coef[nrow(my_lm_coef),4])
    if (flag!=1){
      confounder_result2=rbind(confounder_result2,confounder_result)
    }else{
      confounder_result2=confounder_result
      flag=5
    }
  }
}


write.table(my_corr_factors_results, "~/Desktop/Metabolomics_v2/3.Preliminary_results/technical_confouders_quantitative.txt", sep = "\t", quote = F)
