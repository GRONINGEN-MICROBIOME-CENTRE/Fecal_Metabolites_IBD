Bacterial metabolic functions association to metabolites
===

Testing each metabolite~ host phenotype association using the several techinicals and host charactersitic covariates( age, sex, BMI, bowel movements a day,day measuring (run day), amount of sample, month sample in a freezer) within cohort (controls, CD and UC). 

Tests: 

a. Metabolites (>70% samples) vs metabolic gene clusters in controls (linear regression)
b. Metabolites (>70% samples) vs metabolic gene clusters in CD (linear regression)
c. Metabolites (>70% samples) vs metabolic gene clusters in UC (logistic regression)


d. Metabolites (>70% samples) vs metacyc pathways in controls (linear regression)
e. Metabolites (>70% samples) vs metacyc pathways in CD (linear regression)
f. Metabolites (>70% samples) vs metacyc pathways in UC (logistic regression)


a. Controls
===

#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef")
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef",colnames(run_day))
my_phenos=controls_phenos[,c(regressors)]

controls_mgc=merge(my_phenos,mgc_fil_trans_na, by="row.names")
row.names(controls_mgc)=controls_mgc$Row.names
controls_mgc$Row.names=NULL

controls_mgc=merge(controls_mgc,q_mtbx2, by="row.names")
row.names(controls_mgc)=controls_mgc$Row.names
controls_mgc$Row.names=NULL


flag=1
for ( i in 1:161) {
  my_pheno=colnames(controls_mgc)[i]
  if (my_pheno %in% regressors ) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 162:ncol(controls_mgc)){
      my_trait=colnames(controls_mgc)[a]
      my_uni_test=controls_mgc[,c(regressors,my_pheno,my_trait)]
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

associations_controls_mgc_new=my_univariate_results


b. CD
===

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis",colnames(run_day))

my_phenos=ibd_test_phenos[,c(regressors)]

cd_mgc=merge(my_phenos,mgc_fil_trans, by="row.names")
row.names(cd_mgc)=cd_mgc$Row.names
cd_mgc$Row.names=NULL

cd_mgc=merge(cd_mgc,q_mtbx2, by="row.names")
row.names(cd_mgc)=cd_mgc$Row.names
cd_mgc$Row.names=NULL
cd_mgc=subset(cd_mgc, cd_mgc$clinical_Diagnosis=="2")
cd_mgc$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))


flag=1
for ( i in 1:162) {
  my_pheno=colnames(cd_mgc)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 163:ncol(cd_mgc)){
      my_trait=colnames(cd_mgc)[a]
      my_uni_test=cd_mgc[,c(regressors,my_pheno,my_trait)]
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

associations_cd_mgc_new=my_univariate_results


c. UC
=== 


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis",colnames(run_day))

my_phenos=ibd_test_phenos[,c(regressors)]

uc_mgc=merge(my_phenos,mgc_fil_trans, by="row.names")
row.names(uc_mgc)=uc_mgc$Row.names
uc_mgc$Row.names=NULL

uc_mgc=merge(uc_mgc,q_mtbx2, by="row.names")
row.names(uc_mgc)=uc_mgc$Row.names
uc_mgc$Row.names=NULL
uc_mgc=subset(uc_mgc, uc_mgc$clinical_Diagnosis=="1")
uc_mgc$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))


flag=1
for ( i in 1:162) {
  my_pheno=colnames(uc_mgc)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 163:ncol(uc_mgc)){
      my_trait=colnames(uc_mgc)[a]
      my_uni_test=uc_mgc[,c(regressors,my_pheno,my_trait)]
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

associations_uc_mgc_new=my_univariate_results


#METACYC PATHWAYS [HUMAnN3]


pathways_raw <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/Metacyc/Merged_metacyc_uniref90_filtered_renames.txt", row.names = 1)
pathways_samples=subset(pathways_raw,rownames(pathways_raw) %in% rownames(cc_pheno2))
pathways_samples[pathways_samples==0]=NA
pathways_filt=transform_and_filter_mtb(pathways_samples,samples_row = T,method = "clr",missing_filter = 20)
pathways_filt$UNMAPPED=NULL
pathways_filt$UNINTEGRATED=NULL



d. Controls
===


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef")
#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef",colnames(run_day))
my_phenos=controls_phenos[,c(regressors)]

controls_pwy=merge(my_phenos,pathways_filt, by="row.names")
row.names(controls_pwy)=controls_pwy$Row.names
controls_pwy$Row.names=NULL

controls_pwy=merge(controls_pwy,q_mtbx2, by="row.names")
row.names(controls_pwy)=controls_pwy$Row.names
controls_pwy$Row.names=NULL


flag=1
for ( i in 1:333) {
  my_pheno=colnames(controls_pwy)[i]
  if (my_pheno %in% regressors ) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 334:ncol(controls_pwy)){
      my_trait=colnames(controls_pwy)[a]
      my_uni_test=controls_pwy[,c(regressors,my_pheno,my_trait)]
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

associations_controls_pwy=my_univariate_results


e. CD
===

#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis",colnames(run_day))
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis")
my_phenos=ibd_test_phenos[,c(regressors)]

cd_pwy=merge(my_phenos,pathways_filt, by="row.names")
row.names(cd_pwy)=cd_pwy$Row.names
cd_pwy$Row.names=NULL

cd_pwy=merge(cd_pwy,q_mtbx2, by="row.names")
row.names(cd_pwy)=cd_pwy$Row.names
cd_pwy$Row.names=NULL
cd_pwy=subset(cd_pwy, cd_pwy$clinical_Diagnosis=="2")
cd_pwy$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")


flag=1
for ( i in 1:335) {
  my_pheno=colnames(cd_pwy)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 336:ncol(cd_pwy)){
      my_trait=colnames(cd_pwy)[a]
      my_uni_test=cd_pwy[,c(regressors,my_pheno,my_trait)]
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

associations_cd_pwy=my_univariate_results


f. UC
=== 


#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis",colnames(run_day))

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis")

my_phenos=ibd_test_phenos[,c(regressors)]

uc_pwy=merge(my_phenos,pathways_filt, by="row.names")
row.names(uc_pwy)=uc_pwy$Row.names
uc_pwy$Row.names=NULL

uc_pwy=merge(uc_pwy,q_mtbx2, by="row.names")
row.names(uc_pwy)=uc_pwy$Row.names
uc_pwy$Row.names=NULL
uc_pwy=subset(uc_pwy, uc_pwy$clinical_Diagnosis=="1")
uc_pwy$clinical_Diagnosis=NULL

#regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")

flag=1
for ( i in 1:335) {
  my_pheno=colnames(uc_pwy)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 336:ncol(uc_pwy)){
      my_trait=colnames(uc_pwy)[a]
      my_uni_test=uc_pwy[,c(regressors,my_pheno,my_trait)]
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

associations_uc_pwy=my_univariate_results



