Phenotype association to metabolites
===

Testing each metabolite~ host phenotype association using the several techinicals and host charactersitic covariates( age, sex, BMI, bowel movements a day,day measuring (run day), amount of sample, month sample in a freezer) within cohort (controls, CD and UC). 

Tests: 

a. Metabolites (>70% samples) vs phenotypes in controls (linear regression)
b. Metabolites (<70% samples) vs phenotypes in controls (logistic regression)
c. Metabolites (>70% samples) vs phenotypes in CD (linear regression)
d. Metabolites (<70% samples) vs phenotypes in CD (logistic regression)
e. Metabolites (>70% samples) vs phenotypes in UC (linear regression)
f. Metabolites (<70% samples) vs phenotypes in UC (logistic regression)


a. Metabolites (>70% samples) vs phenotypes in controls
===

controls_phenos=read.table("~/Desktop/IBD_Metabolomics/1.Input/Controls_phenos_recoded.txt", sep = "\t", row.names = 1, header = T)

#Remove diet groups (redundant)
controls_phenos=controls_phenos[,!grepl("diet_group",colnames(controls_phenos))]
controls_phenos=merge(run_day,controls_phenos,by="row.names")
row.names(controls_phenos)=controls_phenos$Row.names
controls_phenos$Row.names=NULL

#remove phenotypes with less than 10 cases
controls_phenos$med_anti_epileptics=NULL
controls_phenos$med_alpha_blockers=NULL
controls_phenos$med_anti_androgen_oral_contraceptive=NULL
controls_phenos$med_anti_epileptics=NULL
controls_phenos$med_ca_channel_blocker=NULL
controls_phenos$med_K_saving_diuretic=NULL
controls_phenos$med_melatonine=NULL
controls_phenos$med_metformin=NULL
controls_phenos$med_methylphenidate=NULL
controls_phenos$med_oral_anti_diabetics=NULL
controls_phenos$med_thyrax=NULL
controls_phenos$med_triptans=NULL
controls_phenos$med_vitamin_K_antagonist=NULL
controls_phenos$med_insulin=NULL
controls_phenos$med_opiat=NULL
controls_phenos$med_oral_anti_diabetics=NULL
controls_phenos$med_other_antidepressant=NULL
controls_phenos$med_bisphosphonates=NULL
controls_phenos$med_folic_acid=NULL
controls_phenos$med_mesalazines=NULL
controls_phenos$med_NSAID=NULL
controls_phenos$med_oral_steroid=NULL
controls_phenos$med_paracetamol=NULL
controls_phenos$med_SSRI_antidepressant=NULL
controls_phenos$med_tricyclic_antidepressant=NULL
controls_phenos$med_vitamin_B12=NULL
controls_phenos$med_vitamin_D=NULL
controls_phenos$med_calcium=NULL
controls_phenos$med_benzodiazepine_derivatives_related=NULL

#controls_phenos$med_antibiotics_merged=NULL
#controls_phenos$med_ferrum=NULL
#controls_phenos$med_parasympathicolytic_inhaler=NULL


q_mtbx2=q_mtbx
colnames(q_mtbx2)=make.names(colnames(q_mtbx2))
controls_phenos=merge(controls_phenos,q_mtbx2, by="row.names")
row.names(controls_phenos)=controls_phenos$Row.names
controls_phenos$Row.names=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef", colnames(run_day))
flag=1
for ( i in 1:219) {
  my_pheno=colnames(controls_phenos)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 220:ncol(controls_phenos)){
      my_trait=colnames(controls_phenos)[a]
      my_uni_test=controls_phenos[,c(regressors,my_pheno,my_trait)]
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

associations_controls=my_univariate_results
associations_controls=subset(associations_controls,associations_controls$phenotype!="sp_Shannon_Index")
associations_controls=subset(associations_controls,associations_controls$phenotype!="BioMK_ChromograninA")
associations_controls=subset(associations_controls,associations_controls$phenotype!="BioMK_BetaDefensin2")
associations_controls=subset(associations_controls,associations_controls$phenotype!="BioMK_Calprotectin")
associations_controls$FDR=p.adjust(associations_controls$`Pr(>|t|)`,method = "BH")
associations_controls$Bonferroni=p.adjust(associations_controls$`Pr(>|t|)`,method = "bonferroni")
associations_controls2=merge(annot,associations_controls, by.y="metabolite", by.x = "for_merg", all.y=T)
associations_controls2$SUPER.PATHWAY[is.na(associations_controls2$SUPER.PATHWAY)]="SCFA"
associations_controls2$SUB.PATHWAY[is.na(associations_controls2$SUB.PATHWAY)]="SCFA"

b. Metabolites (<70% samples) vs phenotypes in controls
===


bi_mtb2=bi_mtb
colnames(bi_mtb2)=make.names(colnames(bi_mtb2))
controls_phenos_prev=merge(controls_phenos,bi_mtb2, by="row.names")
row.names(controls_phenos_prev)=controls_phenos_prev$Row.names
controls_phenos_prev$Row.names=NULL


flag=1
for ( i in 1:305) {
  my_pheno=colnames(controls_phenos_prev)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 306:ncol(controls_phenos_prev)){
      my_trait=colnames(controls_phenos_prev)[a]
      my_uni_test=controls_phenos_prev[,c(regressors,my_pheno,my_trait)]
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


prev_univariate_control=my_univariate_results
prev_univariate_control$FDR=p.adjust(prev_univariate_control$`Pr(>|z|)`,method = "BH")
prev_univariate_control$Bonferroni=p.adjust(prev_univariate_control$`Pr(>|z|)`,method = "bonferroni")


Test CD: phenotypes vs prevalent metabolites [quantitative]
---

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef",colnames(run_day))

ibd_test_phenos=read.table("~/Desktop/IBD_metabolomics_2022/1.Input/IBD_phenos_recoded.txt", sep = "\t", row.names = 1, header = T)
ibd_test_phenos=ibd_test_phenos[,!grepl("diet_group",colnames(ibd_test_phenos))]

ibd_test_phenos=merge(run_day,ibd_test_phenos,by="row.names")
row.names(ibd_test_phenos)=ibd_test_phenos$Row.names
ibd_test_phenos$Row.names=NULL

ibd_test_phenos3=merge(ibd_test_phenos,q_mtbx2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL


ibd_test_cd_phenos=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="2")
#ibd_test_uc_phenos=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="1")


ibd_test_cd_phenos$clinical_Diagnosis=NULL

ibd_test_cd_phenos=merge(day,ibd_test_cd_phenos,by="row.names")
row.names(ibd_test_cd_phenos)=ibd_test_cd_phenos$Row.names
ibd_test_cd_phenos$Row.names=NULL

#ibd_test_cd_phenos$blood_IgARF_pos_neg=NULL
#ibd_test_cd_phenos$blood_IgARF_pos_neg=NULL
#ibd_test_cd_phenos$blood_IgMRF_pos_neg=NULL
ibd_test_cd_phenos$med_alpha_blockers=NULL
ibd_test_cd_phenos$med_parasympathicolytic_inhaler=NULL
ibd_test_cd_phenos$run_day_cat=NULL

flag=1
for ( i in 1:229) {
  my_pheno=colnames(ibd_test_cd_phenos)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 230:ncol(ibd_test_cd_phenos)){
      my_trait=colnames(ibd_test_cd_phenos)[a]
      my_uni_test=ibd_test_cd_phenos[,c(regressors,my_pheno,my_trait)]
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

associations_CD_pheno=my_univariate_results
associations_CD_pheno=subset(associations_CD_pheno,associations_CD_pheno$phenotype!="sp_Shannon_Index")
associations_CD_pheno=subset(associations_CD_pheno,associations_CD_pheno$phenotype!="clinical_FecalCalprotectinOver200yesno")
associations_CD_pheno$FDR=p.adjust(associations_CD_pheno$`Pr(>|t|)`,method = "BH")
associations_CD_pheno$Bonferroni=p.adjust(associations_CD_pheno$`Pr(>|t|)`,method = "bonferroni")

# Test CD: phenotypes vs less abundant metabolites [logistic]

ibd_test_phenos3=merge(ibd_test_phenos,bi_mtb2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL
ibd_test_cd_phenos_prev=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="2")
ibd_test_cd_phenos_prev$clinical_Diagnosis=NULL


flag=1
for ( i in 1:195) {
  my_pheno=colnames(ibd_test_cd_phenos_prev)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 196:ncol(ibd_test_cd_phenos_prev)){
      my_trait=colnames(ibd_test_cd_phenos_prev)[a]
      my_uni_test=ibd_test_cd_phenos_prev[,c(regressors,my_pheno,my_trait)]
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

associations_CD_pheno_prev=my_univariate_results


# Test UC: phenotypes vs prevalent metabolites [quantitative]


ibd_test_phenos3=merge(ibd_test_phenos,q_mtbx2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL
ibd_test_uc_phenos=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="1")

ibd_test_uc_phenos$clinical_Diagnosis=NULL
ibd_test_uc_phenos$clinical_NumberOfResectionsIleal=NULL
ibd_test_uc_phenos$clinical_inflammation_ileum=NULL
ibd_test_uc_phenos$clinical_inflammation_colon=NULL
ibd_test_uc_phenos$clinical_ResectionIlealAny=NULL
ibd_test_uc_phenos$clinical_ResectionIlealCecalAny=NULL
ibd_test_uc_phenos$clinical_DrainageAbcesSeton=NULL
ibd_test_uc_phenos$meds_Biologicals=NULL
ibd_test_uc_phenos$clinical_EverHadStomaOrPouch=NULL
ibd_test_uc_phenos$clinical_NumberOfResetionsIleoCecal=NULL


ibd_test_uc_phenos=merge(day,ibd_test_uc_phenos,by="row.names")
row.names(ibd_test_uc_phenos)=ibd_test_uc_phenos$Row.names
ibd_test_uc_phenos$Row.names=NULL
#ibd_test_uc_phenos$blood_IgARF_pos_neg=NULL
#ibd_test_uc_phenos$blood_IgMRF_pos_neg=NULL
#ibd_test_uc_phenos$clinical_Fistulotomy=NULL
#ibd_test_uc_phenos$med_angII_receptor_antagonist=NULL
#ibd_test_uc_phenos$med_anti_histamine=NULL
#ibd_test_uc_phenos$med_antibiotics_merged=NULL
#ibd_test_uc_phenos$med_opiat=NULL
#ibd_test_uc_phenos$med_other_antidepressant=NULL
ibd_test_uc_phenos$med_anti_androgen_oral_contraceptive=NULL
ibd_test_uc_phenos$med_parasympathicolytic_inhaler=NULL
ibd_test_uc_phenos$run_day_cat=NULL
ibd_test_uc_phenos=ibd_test_uc_phenos[,!grepl("diet_group",colnames(ibd_test_uc_phenos))]

flag=1
for ( i in 1:220) {
  my_pheno=colnames(ibd_test_uc_phenos)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 221:ncol(ibd_test_uc_phenos)){
      my_trait=colnames(ibd_test_uc_phenos)[a]
      my_uni_test=ibd_test_uc_phenos[,c(regressors,my_pheno,my_trait)]
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

associations_UC_pheno=my_univariate_results
associations_UC_pheno=subset(associations_UC_pheno,associations_UC_pheno$phenotype!="sp_Shannon_Index")
associations_UC_pheno=subset(associations_UC_pheno,associations_UC_pheno$phenotype!="clinical_FecalCalprotectinOver200yesno")
associations_UC_pheno$FDR=p.adjust(associations_UC_pheno$`Pr(>|t|)`,method = "BH")
associations_UC_pheno$Bonferroni=p.adjust(associations_UC_pheno$`Pr(>|t|)`,method = "bonferroni")

# Test UC: phenotypes vs less abundant metabolites [logistic]


ibd_test_phenos3=merge(ibd_test_phenos,bi_mtb2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL

ibd_test_uc_phenos=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="1")

ibd_test_uc_phenos$clinical_Diagnosis=NULL
ibd_test_uc_phenos$clinical_NumberOfResectionsIleal=NULL
ibd_test_uc_phenos$clinical_inflammation_ileum=NULL
ibd_test_uc_phenos$clinical_ResectionIlealAny=NULL
ibd_test_uc_phenos$clinical_ResectionIlealCecalAny=NULL
ibd_test_uc_phenos$clinical_DrainageAbcesSeton=NULL
#ibd_test_uc_phenos$meds_Biologicals=NULL
ibd_test_uc_phenos$clinical_EverHadStomaOrPouch=NULL
ibd_test_uc_phenos$clinical_NumberOfResetionsIleoCecal=NULL


flag=1
for ( i in 1:187) {
  my_pheno=colnames(ibd_test_uc_phenos)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 188:ncol(ibd_test_uc_phenos)){
      my_trait=colnames(ibd_test_uc_phenos)[a]
      my_uni_test=ibd_test_uc_phenos[,c(regressors,my_pheno,my_trait)]
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

associations_UC_pheno_prev=my_univariate_results
