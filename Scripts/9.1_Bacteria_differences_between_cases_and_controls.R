Test microbial differents between cases and controls 
===

a. Taxa
===

case_control_quantitative=read.table("~/Desktop/Metabolomics_v2/2.Input/case_control_input.txt",sep="\t", header = T, row.names = 1)
case_control_bacteria=case_control_quantitative[,1:14]
case_control_bacteria$Shannon_Index=NULL
bugs=taxa_and_mb[,1:109]
colnames(bugs)=make.names(colnames(bugs))
case_control_bacteria=merge(case_control_bacteria, bugs, by="row.names")
rownames(case_control_bacteria)=case_control_bacteria$Row.names
case_control_bacteria$Row.names=NULL
case_control_bacteria=merge(cc_rd,case_control_bacteria, by="row.names")
rownames(case_control_bacteria)=case_control_bacteria$Row.names
case_control_bacteria$Row.names=NULL
case_control_bacteria$LC.COLUMN=NULL
case_control_bacteria$Amount_sample_gram=NULL
case_control_bacteria$metabolon_Month_in_freezer=NULL

flag=1
for ( i in 1:11) {
  my_pheno=colnames(case_control_bacteria)[i]
  if (my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="PF_RD" | my_pheno=="host_BMI" | my_pheno=="host_Sex" | my_pheno=="host_Age") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 12:ncol(case_control_bacteria)){
      my_trait=colnames(case_control_bacteria)[a]
      my_uni_test=case_control_bacteria[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","clinical_BowelMovementADayDef", "host_BMI","Read_Depth",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      #my_uni_test[,8]=as.numeric(as.character(my_uni_test[,8]))
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

associations_ibd_bugs=my_univariate_results
associations_ibd_bugs$FDR=p.adjust(associations_ibd_bugs$`Pr(>|t|)`,method = "BH")
associations_ibd_bugs$Bonferroni=p.adjust(associations_ibd_bugs$`Pr(>|t|)`,method = "bonferroni")



b. Test Biosynthetic gene clusters
===


row.names(cc_rd)=cc_rd$PID
cc_rd$PID=NULL
#cc_rd$IDs_MGS=NULL
my_phenos=merge(my_phenos,cc_rd, by="row.names")
row.names(my_phenos)=my_phenos$Row.names
my_phenos$Row.names=NULL
case_control_mgc=merge(my_phenos,mgc_fil_trans, by="row.names")

rownames(case_control_mgc)=case_control_mgc$Row.names
case_control_mgc$Row.names=NULL
case_control_mgc$LC.COLUMN=NULL
case_control_mgc$Amount_sample_gram=NULL
case_control_mgc$metabolon_Month_in_freezer=NULL
#colnames(case_control_mgc)=make.names(colnames(case_control_mgc))


flag=1
for ( i in 1:11) {
  my_pheno=colnames(case_control_mgc)[i]
  if (my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="PF_RD" | my_pheno=="host_BMI" | my_pheno=="host_Sex" | my_pheno=="host_Age") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 12:ncol(case_control_mgc)){
      my_trait=colnames(case_control_mgc)[a]
      my_uni_test=case_control_mgc[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","PF_RD",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","clinical_BowelMovementADayDef", "host_BMI","PF_RD",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      #my_uni_test[,8]=as.numeric(as.character(my_uni_test[,8]))
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

associations_ibd_mgc2=my_univariate_results