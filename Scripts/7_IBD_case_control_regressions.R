Case control analyses
===

Testing each metabolite~ disease phenotype association using the several techinicals and host charactersitic covariates( age, sex, BMI, bowel movements a day,day measuring (run day), amount of sample, month sample in a freezer and dietary differences). 

Different tests: 

a) Differences in diet between cases and controls
b) Case-control to assess IBD/CD/UC associated metabolites
c) Case-control to assess IBD/CD/UC presence/absence metabolites (>20 MTB <70%) 


a. Check differences between diet IBD/Controls
---

Fecal metabolites are expected to be (partially) influenced by diet
--
#col 15 = ibd_diagnosis
#col 74:217 = diet

test_diet=cc_pheno6[,c(15,22:24,74:217)]
test_diet$IBD="IBD"
test_diet$IBD[test_diet$ibd_Diagnosis=="Control"]="Control"
test_diet$ibd_Diagnosis=NULL

test_diet$diet_PlanttoAnimal[is.infinite(test_diet$diet_PlanttoAnimal)]=NA

flag=10
for (i in 4:(ncol(test_diet)-1)){
  my_trait=colnames(test_diet)[i]
  my_preds=c("host_Sex","host_Age_sampling", "host_BMI", "IBD")
  my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
  my_lm=summary(lm(my_f,data = test_diet ))
  my_lm_coef=as.data.frame(my_lm$coefficients)
  my_lm_pheno=try(my_lm_coef[grep("IBD", rownames(my_lm_coef)), ])
  my_lm_pheno$diet=my_trait
  my_lm_pheno$phenotype="IBD"
  my_lm_pheno$factor=rownames(my_lm_pheno)
  if (flag!=10){
    diet_diff=rbind(diet_diff,my_lm_pheno)
  }else{
    diet_diff=my_lm_pheno
    flag=5
  }
}

diet_diff=subset(diet_diff, !grepl("diet_group", diet_diff$diet))
diet_diff$FDR=p.adjust(diet_diff$`Pr(>|t|)`,method = "BH")
diet_diff$Bonferroni=p.adjust(diet_diff$`Pr(>|t|)`,method = "bonferroni")
rownames(diet_diff)=diet_diff$diet

b. Metabolites associated to IBD
---


#samples_resections=row.names(subset(phenos_ibd,phenos_ibd$ibd_ResectionAny=="yes"))
#case_control_quantitative2=subset(case_control_quantitative, row.names(case_control_quantitative)%ni% samples_resections)

samples_resections=phenos_ibd[,"ibd_ResectionAny", drop=F]
case_control_quantitative2=merge(samples_resections,case_control_quantitative, by="row.names", all.y=T)
row.names(case_control_quantitative2)=case_control_quantitative2$Row.names
case_control_quantitative2$Row.names=NULL

case_control_quantitative2$ibd_ResectionAny[is.na(case_control_quantitative2$ibd_ResectionAny)]="no"

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling","host_BMI", "clinical_BowelMovementADayDef", rownames(diet_diff)[diet_diff$FDR<0.05], colnames(run_day),"ibd_ResectionAny")



flag=1
for ( i in 1:54) {
  my_pheno=colnames(case_control_quantitative2)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 55:ncol(case_control_quantitative2)){
      my_trait=colnames(case_control_quantitative2)[a]
      my_uni_test=case_control_quantitative2[,c(regressors,my_pheno,my_trait)]
      my_preds=c(regressors,my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
          my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
          my_univariate_results=my_lm_pheno
          flag=5
      }
    }
  }
}

associations_ibd_resec=my_univariate_results
my_univariate_results=subset(my_univariate_results,my_univariate_results$factor!="ibd_DiagnosisIBDU")
my_univariate_results$phenotype="IBD"
my_univariate_results$phenotype[my_univariate_results$factor=="ibd_DiagnosisCD"]="CD"
my_univariate_results$phenotype[my_univariate_results$factor=="ibd_DiagnosisUC"]="UC"

my_cc_quantitative_resec=merge(my_univariate_results,annot, by.x="metabolite", by.y = "for_merg", all.x=T)

my_cc_quantitative_resec$FDR=p.adjust(my_cc_quantitative_resec$`Pr(>|t|)`,method = "BH")
my_cc_quantitative_resec$Bonferroni=p.adjust(my_cc_quantitative_resec$`Pr(>|t|)`,method = "bonferroni")



c. TEST PREVALENCE <20% MTB <70%
---

# LOGISTIC REGRESSION #

test_mtb_presence=merge(cc_pheno6,bi_mtb,by="row.names")
rownames(test_mtb_presence)=test_mtb_presence$Row.names
test_mtb_presence$Row.names=NULL
test_mtb_presence$ibd_FecalCalprotectin=NULL
test_mtb_presence$COMMENT=NULL

case_control_prevalence=merge(case_control_test,bi_mtb,by="row.names")
rownames(case_control_prevalence)=case_control_prevalence$Row.names
case_control_prevalence$Row.names=NULL

colnames(case_control_prevalence)=make.names(colnames(case_control_prevalence))

#case_control_prevalence2=subset(case_control_prevalence, row.names(case_control_prevalence)%ni% samples_resections)
#samples_resections=phenos_ibd[,"ibd_ResectionAny", drop=F]
case_control_prevalence2=merge(samples_resections,case_control_prevalence, by="row.names", all.y=T)
row.names(case_control_prevalence2)=case_control_prevalence2$Row.names
case_control_prevalence2$Row.names=NULL
case_control_prevalence2$ibd_ResectionAny[is.na(case_control_prevalence2$ibd_ResectionAny)]="no"
case_control_prevalence2=merge(day,case_control_prevalence2,by="row.names")
row.names(case_control_prevalence2)=case_control_prevalence2$Row.names
case_control_prevalence2$Row.names=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling","host_BMI", "clinical_BowelMovementADayDef", rownames(diet_diff)[diet_diff$FDR<0.05], "run_day_cat","ibd_ResectionAny")



flag=1
for ( i in 1:35) {
  my_pheno=colnames(case_control_prevalence2)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 36:ncol(case_control_prevalence2)){
      my_trait=colnames(case_control_prevalence2)[a]
      my_uni_test=case_control_prevalence2[,c(regressors,my_pheno,my_trait)]
      my_preds=c(regressors,my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_samples=nrow(my_uni_test)
      #my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm=summary(glm(my_f, family = binomial(link="logit"), data =my_uni_test))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
          my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
          my_univariate_results=my_lm_pheno
          flag=5
      }
    }
  }
}

prev_univariate_ibd2=my_univariate_results

my_prev_results2=subset(prev_univariate_ibd2,prev_univariate_ibd2$factor!="ibd_DiagnosisIBDU")
my_prev_results2$phenotype="IBD"
my_prev_results2$phenotype[my_prev_results2$factor=="ibd_DiagnosisCD"]="CD"
my_prev_results2$phenotype[my_prev_results2$factor=="ibd_DiagnosisUC"]="UC"

my_cc_prev_resec2=merge(my_prev_results2,annot, by.x="metabolite", by.y = "for_merg", all.x=T)

my_cc_prev_resec2$FDR=p.adjust(my_cc_prev_resec2$`Pr(>|z|)`,method = "BH")
my_cc_prev_resec2$Bonferroni=p.adjust(my_cc_prev_resec2$`Pr(>|z|)`,method = "bonferroni")

