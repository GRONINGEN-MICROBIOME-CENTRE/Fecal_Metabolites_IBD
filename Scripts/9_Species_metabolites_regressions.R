Phenotype association to metabolites
===

Testing each metabolite~ taxa association using the several techinicals and host charactersitic covariates( age, sex, BMI, bowel movements a day,day measuring (run day), amount of sample, month sample in a freezer) and intestinal resections (any, yes/no) within cohort (controls, CD and UC). 

Tests: 

a. Metabolites (>70% samples) vs taxa in controls (linear regression). Taxa coded as presence/absence (1/0)
b. Metabolites (>70% samples) vs taxa in controls (linear regression). Taxa abundance (missing taxa = NA)
c. Metabolites prevalence (Metabolites <70% samples coded as 0/1) vs taxa in controls (logistic regression). Taxa abundance (missing taxa = NA)

d. Metabolites (>70% samples) vs taxa in CD (linear regression). Taxa coded as presence/absence (1/0)
e. Metabolites (>70% samples) vs taxa in CD (linear regression). Taxa abundance (missing taxa = NA)
f. Metabolites prevalence (Metabolites <70% samples coded as 0/1) vs taxa in CD (logistic regression). Taxa abundance (missing taxa = NA)

g. Metabolites (>70% samples) vs taxa in UC (linear regression). Taxa coded as presence/absence (1/0)
h. Metabolites (>70% samples) vs taxa in UC (linear regression). Taxa abundance (missing taxa = NA)
i. Metabolites prevalence (Metabolites <70% samples coded as 0/1) vs taxa in UC (logistic regression). Taxa abundance (missing taxa = NA)



a. Controls: Code bacteria as presence / absence(1/0) value
==

taxa3=taxa2
taxa3[taxa3>0]=1
taxa3=taxa3[row.names(taxa3)%in%row.names(taxa_filt),]
taxa3=taxa3[,colnames(taxa3)%in%colnames(taxa_filt)]
colnames(taxa3)=make.names(colnames(taxa3))

#controls_phenos=merge(run_day,controls_phenos,by="row.names")
#row.names(controls_phenos)=controls_phenos$Row.names
#controls_phenos$Row.names=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef",colnames(run_day))
my_phenos=controls_phenos[,c(regressors)]

controls_species=merge(my_phenos,taxa3, by="row.names")
row.names(controls_species)=controls_species$Row.names
controls_species$Row.names=NULL

q_mtbx2=q_mtbx
colnames(q_mtbx2)=make.names(colnames(q_mtbx2))

controls_species=merge(controls_species,q_mtbx2, by="row.names")
row.names(controls_species)=controls_species$Row.names
controls_species$Row.names=NULL

flag=1
for ( i in 1:136) {
  my_pheno=colnames(controls_species)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 137:ncol(controls_species)){
      my_trait=colnames(controls_species)[a]
      my_uni_test=controls_species[,c(regressors,my_pheno,my_trait)]
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

associations_controls_sp=my_univariate_results
associations_controls$FDR=p.adjust(associations_controls$`Pr(>|t|)`,method = "BH")
associations_controls$Bonferroni=p.adjust(associations_controls$`Pr(>|t|)`,method = "bonferroni")

associations_qp=merge(associations_qp,annot, by.x="metabolite", by.y = "for_merg", all.x=T)


b. Controls: Test associations bacteria-metabolites in controls but without bacterial zeros (bacterial 0 => NAs, CLR transformed values as non-zeros values)
==

taxa_filt2=taxa_filt
colnames(taxa_filt2)=make.names(colnames(taxa_filt2))
table(row.names(taxa3)==row.names(taxa_filt2))
table(colnames(taxa3)==colnames(taxa_filt2))
taxa4=taxa3*taxa_filt2
taxa4[taxa4==0]=NA

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef")
my_phenos=controls_phenos[,c(regressors)]


controls_species_lin=merge(my_phenos,taxa4, by="row.names")
row.names(controls_species_lin)=controls_species_lin$Row.names
controls_species_lin$Row.names=NULL

controls_species_lin=merge(controls_species_lin,q_mtbx2, by="row.names")
row.names(controls_species_lin)=controls_species_lin$Row.names
controls_species_lin$Row.names=NULL


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef")

flag=1
for ( i in 1:116) {
  my_pheno=colnames(controls_species_lin)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 117:ncol(controls_species_lin)){
      my_trait=colnames(controls_species_lin)[a]
      my_uni_test=controls_species_lin[,c(regressors,my_pheno,my_trait)]
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

associations_controls_nz2=my_univariate_results
associations_controls_nz2$FDR=p.adjust(associations_controls_nz2$`Pr(>|t|)`,method = "BH")
associations_controls_nz2$Bonferroni=p.adjust(associations_controls_nz2$`Pr(>|t|)`,method = "bonferroni")

c. Controls: Presence of species is related to presence of metabolites 
===

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef")
my_phenos=controls_phenos[,c(regressors)]

controls_species=merge(my_phenos,my_taxa4, by="row.names")
row.names(controls_species)=controls_species$Row.names
controls_species$Row.names=NULL

controls_species_prev2=merge(controls_species,bi_mtb2, by="row.names")
row.names(controls_species_prev2)=controls_species_prev2$Row.names
controls_species_prev2$Row.names=NULL

flag=1
for ( i in 1:116) {
  my_pheno=colnames(controls_species_prev2)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 117:ncol( controls_species_prev2)){
      my_trait=colnames(controls_species_prev2)[a]
      my_uni_test=controls_species_prev2[,c(regressors,my_pheno,my_trait)]
      my_preds=c(regressors,my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_samples=nrow(my_uni_test)
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

associations_cnt_sp_prev_mtb_prev=my_univariate_results



d. Test CD: Presence of species is associated to metabolites levels
===

#Test only in CD
regressors=c("clinical_Diagnosis","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))

my_phenos=ibd_test_phenos[,c(regressors)]

my_phenos_cd=merge(my_phenos,taxa3, by="row.names")
row.names(my_phenos_cd)=my_phenos_cd$Row.names
my_phenos_cd$Row.names=NULL

ibd_species=merge(my_phenos_cd,q_mtbx2, by="row.names")
row.names(ibd_species)=ibd_species$Row.names
ibd_species$Row.names=NULL

cd_specie=subset(ibd_species, ibd_species$clinical_Diagnosis=="2")
cd_specie$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))

flag=1
for ( i in 1:137) {
  my_pheno=colnames(cd_specie)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 138:ncol(cd_specie)){
      my_trait=colnames(cd_specie)[a]
      my_uni_test=cd_specie[,c(regressors,my_pheno,my_trait)]
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

associations_cd_sp_prev=my_univariate_results
associations_cd_sp_prev$FDR=p.adjust(associations_cd_sp_prev$`Pr(>|t|)`,method = "BH")
associations_cd_sp_prev$Bonferroni=p.adjust(associations_cd_sp_prev$`Pr(>|t|)`,method = "bonferroni")

e. Test CD: Abundance of specie is related to metabolite levels
==


regressors=c("clinical_Diagnosis","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")
my_phenos=ibd_test_phenos[,c(regressors)]

ibd_species_lin=merge(my_phenos,taxa4, by="row.names")
row.names(ibd_species_lin)=ibd_species_lin$Row.names
ibd_species_lin$Row.names=NULL

ibd_species_lin=merge(ibd_species_lin,q_mtbx2, by="row.names")
row.names(ibd_species_lin)=ibd_species_lin$Row.names
ibd_species_lin$Row.names=NULL

cd_specie_lin=subset(ibd_species_lin, ibd_species_lin$clinical_Diagnosis=="2")
cd_specie_lin$clinical_Diagnosis=NULL
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")


flag=1
for ( i in 1:117) {
  my_pheno=colnames(cd_specie_lin)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 118:ncol(cd_specie_lin)){
      my_trait=colnames(cd_specie_lin)[a]
      my_uni_test=cd_specie_lin[,c(regressors,my_pheno,my_trait)]
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

associations_cd_sp_nz=my_univariate_results


f.Test CD: Presence of species is related to presence of metabolites [<70%]
===

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis")

my_phenos=ibd_test_phenos[,c(regressors)]

ibd_species=merge(my_phenos,my_taxa4, by="row.names")
row.names(ibd_species)=ibd_species$Row.names
ibd_species$Row.names=NULL

ibd_test_phenos3=merge(ibd_species,bi_mtb2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL
ibd_test_cd_phenos_prev2=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="2")
ibd_test_cd_phenos_prev2$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")


flag=1
for ( i in 1:117) {
  my_pheno=colnames(ibd_test_cd_phenos_prev2)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 118:ncol( ibd_test_cd_phenos_prev2)){
      my_trait=colnames(ibd_test_cd_phenos_prev2)[a]
      my_uni_test=ibd_test_cd_phenos_prev2[,c(regressors,my_pheno,my_trait)]
      my_preds=c(regressors,my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_samples=nrow(my_uni_test)
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

associations_cd_sp_prev_mtb_prev=my_univariate_results

g. UC: test Species abundance related to metabolite levels
===

uc_specie_lin=subset(ibd_species_lin, ibd_species_lin$clinical_Diagnosis=="1")
uc_specie_lin$clinical_Diagnosis=NULL

flag=1
for ( i in 1:117) {
  my_pheno=colnames(uc_specie_lin)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 118:ncol(uc_specie_lin)){
      my_trait=colnames(uc_specie_lin)[a]
      my_uni_test=uc_specie_lin[,c(regressors,my_pheno,my_trait)]
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

associations_uc_sp_nz=my_univariate_results


h. UC test species presence related to metabolites levels
===
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny",colnames(run_day))

uc_specie=subset(ibd_species, ibd_species$clinical_Diagnosis=="1")
uc_specie$clinical_Diagnosis=NULL

flag=1
for ( i in 1:137) {
  my_pheno=colnames(uc_specie)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 138:ncol(uc_specie)){
      my_trait=colnames(uc_specie)[a]
      my_uni_test=uc_specie[,c(regressors,my_pheno,my_trait)]
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

associations_uc_sp_prev=my_univariate_results


i. Test UC: Presence of species is related to presence of metabolites [<70%]
===


regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny", "clinical_Diagnosis")

my_phenos=ibd_test_phenos[,c(regressors)]

ibd_species=merge(my_phenos,my_taxa4, by="row.names")
row.names(ibd_species)=ibd_species$Row.names
ibd_species$Row.names=NULL

ibd_test_phenos3=merge(ibd_species,bi_mtb2, by="row.names")
row.names(ibd_test_phenos3)=ibd_test_phenos3$Row.names
ibd_test_phenos3$Row.names=NULL
ibd_test_cd_phenos_prev2=subset(ibd_test_phenos3, ibd_test_phenos3$clinical_Diagnosis=="1")
ibd_test_cd_phenos_prev2$clinical_Diagnosis=NULL

regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age", "host_BMI", "clinical_BowelMovementADayDef","clinical_ResectionAny")

flag=1
for ( i in 1:117) {
  my_pheno=colnames(ibd_test_cd_phenos_prev2)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 118:ncol( ibd_test_cd_phenos_prev2)){
      my_trait=colnames(ibd_test_cd_phenos_prev2)[a]
      my_uni_test=ibd_test_cd_phenos_prev2[,c(regressors,my_pheno,my_trait)]
      my_preds=c(regressors,my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_samples=nrow(my_uni_test)
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

associations_uc_sp_prev_mtb_prev=my_univariate_results




