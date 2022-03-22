Calculate metabolic ratios
====
  
  a. Bile-acids
===
  
#Select primary and secondary bile acids
PBA=rownames(annot)[annot$SUB.PATHWAY=="Primary Bile Acid Metabolism"]
SBA=rownames(annot)[annot$SUB.PATHWAY=="Secondary Bile Acid Metabolism"]
BA=all_new_ID_raw[,colnames(all_new_ID_raw) %in% PBA | colnames(all_new_ID_raw) %in% SBA]
PBA2=all_new_ID_raw[,colnames(all_new_ID_raw) %in% PBA]
SBA2=all_new_ID_raw[,colnames(all_new_ID_raw) %in% SBA]

#Replace NA's for a very small value (half of the minimum in the matrix)
PBA2[is.na(PBA2)]=min(SBA2,na.rm = T)/2
SBA2[is.na(SBA2)]=min(SBA2,na.rm = T)/2

# Calculate total PBA and SBA
SBA2$SBA=rowSums(SBA2)
PBA2$PBA=rowSums(PBA2)

#Merge both datasets

BA_norm=rbind(PBA2,SBA2)

#Calculate ratios
#SBA/PBA => DCA/CA || LCA/CDCA || LCA / UDCA
BA_norm$R0_SBA_vs_PBA=scale(log(BA_norm$SBA / BA_norm$PBA))
BA_norm$R1_DCA_vs_CA=scale(log(BA_norm$deoxycholate / BA_norm$cholate))
BA_norm$R2_LCA_vs_CDCA=scale(log(BA_norm$lithocholate / BA_norm$chenodeoxycholate))
BA_norm$R3_UCDA_vs_CDCA=scale(log(BA_norm$ursodeoxycholate / BA_norm$chenodeoxycholate ))

#Conjugated/Unconjugated => (TCA + GCA) / CA || (TCDCA + GCDCA) / CDCA 
BA_norm$R4_Conju_vs_Unconj_CA=scale(log((BA_norm$glycocholate+BA_norm$taurocholate) /BA_norm$cholate))
BA_norm$R5_Conju_vs_Unconj_CDCA=scale(log((BA_norm$glycochenodeoxycholate+BA_norm$taurochenodeoxycholate) /BA_norm$chenodeoxycholate))
BA_norm$R6_Conju_vs_Unconj_LCA=scale(log((BA_norm$glycolithocholate+BA_norm$taurolithocholate) /BA_norm$lithocholate))
BA_norm$R7_Conju_vs_Unconj_UDCA=scale(log((BA_norm$glycoursodeoxycholate+BA_norm$tauroursodeoxycholate) /BA_norm$ursodeoxycholate))
BA_norm$R8_Conju_vs_Unconj_DCA=scale(log((BA_norm$glycodeoxycholate+BA_norm$taurodeoxycholate) /BA_norm$deoxycholate))

# Add phenotypes
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling", "host_BMI", "clinical_BowelMovementADayDef","run_day_cat", rownames(diet_diff)[diet_diff$FDR<0.05])
case_control_test=cc_pheno6[,c(regressors, "IBD", "ibd_Diagnosis")]

BA_p=merge(case_control_test,BA_norm,by="row.names")
row.names(BA_p)=BA_p$Row.names
BA_p$Row.names=NULL
BA_p=subset(BA_p,BA_p$ibd_Diagnosis!="IBDU")
BA_p$ibd_Diagnosis=as.character(BA_p$ibd_Diagnosis)

BA_p=BA_p[,c(1:32,43,76:85)]
BA_p$PBA=scale(log(BA_p$PBA))
BA_p$SBA=scale(log(BA_p$SBA))
BAx=q_mtb[,colnames(q_mtb) %in% PBA | colnames(q_mtb) %in% SBA]

BA_p1=merge(BA_p,BAx, by="row.names")
row.names(BA_p1)=BA_p1$Row.names
BA_p1$Row.names=NULL

#Regress cofounders 
colnames(BA_p1)=make.names(colnames(BA_p1))
for (cc in colnames(BA_p1)) {if (sum(is.na(BA_p1[[cc]])) > 0) {BA_p1[[cc]][is.na(BA_p1[[cc]])] <- median(BA_p1[[cc]],na.rm = T)}}
BA_r=matrix(ncol=28, nrow=nrow(BA_p1))

n=1
for (i in 33:ncol(BA_p1)){
  my_trait=colnames(BA_p1)[i]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = BA_p1)
  BA_r[,n] = my_lm$residuals
  n=n+1
}
BA_r=as.data.frame(BA_r)
rownames(BA_r)=row.names(BA_p1)
colnames(BA_r)=colnames(BA_p1)[33:ncol(BA_p1)]


my_sel=phenos_ibd[,c("ibd_Diagnosis","ibd_ResectionIlealCecalAny", "ibd_ResectionIlealAny", "ibd_DiseaseLocation","ibd_ActiveDisease","ibd_IleocecalValveInSitu")]
my_sel$Resection_ileum="3_CD_intact_ileum"
my_sel$Resection_ileum[my_sel$ibd_ResectionIlealCecalAny=="yes" | my_sel$ibd_ResectionIlealAny=="yes"]="4_CD_resection_in_ileum"
my_sel$Resection_ileum[my_sel$ibd_Diagnosis=="UC" ]="2_UC"

my_sel$Infl=NA
my_sel$Infl[my_sel$ibd_Diagnosis=="UC" & my_sel$ibd_ActiveDisease=="NotActive"]="A_UC_not_active"
my_sel$Infl[my_sel$ibd_Diagnosis=="UC" & my_sel$ibd_ActiveDisease=="Active"]="B_UC_active"
my_sel$Infl[my_sel$ibd_Diagnosis=="CD" & my_sel$ibd_ActiveDisease=="NotActive" & my_sel$ibd_DiseaseLocation=="colon" ]="C_CD_Colon_Remission"
my_sel$Infl[my_sel$ibd_Diagnosis=="CD" & my_sel$ibd_ActiveDisease=="Active" & my_sel$ibd_DiseaseLocation=="colon"]="D_CD_Colon_Active"
my_sel$Infl[my_sel$ibd_Diagnosis=="CD" & my_sel$ibd_ActiveDisease=="NotActive" & my_sel$ibd_DiseaseLocation!="colon"]="E_CD_Ileum_Remission"
my_sel$Infl[my_sel$ibd_Diagnosis=="CD" & my_sel$ibd_ActiveDisease=="Active" & my_sel$ibd_DiseaseLocation!="colon" ]="F_CD_Ileum_Active"

my_sel2=my_sel[,c("Resection_ileum","Infl")]
BA_p2=merge(my_sel2, BA_r, by="row.names", all.y=T)

BA_p2$Resection_ileum[is.na(BA_p2$Resection_ileum)]="1_Controls"
BA_p2$Infl[BA_p2$Resection_ileum=="1_Controls"]="A_Controls"


flag=1
for (a in 4:ncol(BA_p2)){
  x=melt(pairwise.wilcox.test(BA_p2[,a],BA_p2$Resection_ileum, p.adjust.method = "none")$p.value)
  y=melt(pairwise.wilcox.test(BA_p2[,a],BA_p2$Infl, p.adjust.method = "none")$p.value)
  y$pheno="Inflammation"
  x$pheno="Resection"
  z=rbind(x,y)
  z$ratio=colnames(BA_p2)[a]
  if(flag==1){
    associations_BA=z
    flag=7
  }else{
    associations_BA=rbind(associations_BA,z)
  }
}

associations_BA=associations_BA[complete.cases(associations_BA$value),]
associations_BA$FDR=p.adjust(associations_BA$value, method="BH")
associations_BA$Bonferroni=p.adjust(associations_BA$value, method="bonferroni") 
ggplot (BA_p2, aes(Resection_ileum, R8_Conju_vs_Unconj_DCA, fill=Resection_ileum))+ geom_violin(alpha=0.8) + geom_jitter(height = 0, width = 0.1, alpha=0.4) + geom_boxplot(width=0.1, fill="white", outlier.alpha = 0) + theme_bw() + scale_fill_manual(values = c("gray80", "darkolivegreen3", "mediumpurple", "purple4"))


b. Tryptophan metabolites and PUFA ratios
====

case_control_quantitative2$ibd_ResectionAny[is.na(case_control_quantitative2$ibd_ResectionAny)]="no"
regressors=c("LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling","host_BMI", "clinical_BowelMovementADayDef", rownames(diet_diff)[diet_diff$FDR<0.05], colnames(run_day),"ibd_ResectionAny")

phenos_ratios=case_control_quantitative2[, 1:54]

phenos_ratios$CD_vs_UC=phenos_ratios$ibd_Diagnosis

phenos_ratios$CD_vs_UC[phenos_ratios$CD_vs_UC=="A_Control"]=NA
phenos_ratios$CD_vs_UC[phenos_ratios$CD_vs_UC=="IBDU"]=NA
phenos_ratios$ibd_Diagnosis[phenos_ratios$ibd_Diagnosis=="IBDU"]=NA


Try_ratios=all_new_ID_raw[,c("tryptophan", "tryptamine","serotonin", "kynurenine")]


#Calculate ratios

# Kynurenine/try || tryptamine/try || serotonine/try

Try_ratios$R1_KY=scale(log(Try_ratios$kynurenine / Try_ratios$tryptophan))
Try_ratios$R2_TRY=scale(log(Try_ratios$tryptamine / Try_ratios$tryptophan))
Try_ratios$R3_SERO=scale(log(Try_ratios$serotonin / Try_ratios$tryptophan ))

Try_ratios_p=merge(phenos_ratios,Try_ratios,by="row.names")
row.names(Try_ratios_p)=Try_ratios_p$Row.names
Try_ratios_p$Row.names=NULL
Try_ratios_p$tryptophan=NULL
Try_ratios_p$tryptamine=NULL
Try_ratios_p$kynurenine=NULL
Try_ratios_p$serotonin=NULL


# Calculate PUFA ratios

PUFAs=rownames(annot)[annot$SUB.PATHWAY=="Long Chain Polyunsaturated Fatty Acid (n3 and n6)"]
o3=all_new_ID_raw[,c("docosahexaenoate_DHA__22_6n3","docosapentaenoate_DPA__22_5n3","eicosapentaenoate_EPA__20_5n3", "hexadecatrienoate_16_3n3","stearidonate_18_4n3")]
o3$total_o3=rowSums(o3,na.rm = T)
o6=all_new_ID_raw[,c("arachidonate_20_4n6","dihomolinoleate_20_2n6","dihomolinolenate_20_3n3_or_3n6","docosadienoate_22_2n6","hexadecadienoate_16_2n6","linoleate_18_2n6")]
o6$total_o6=rowSums(o6,na.rm = T)

my_PUFA=merge(o3,o6,by="row.names")
row.names(my_PUFA)=my_PUFA$Row.names
my_PUFA$Row.names=NULL
my_PUFA=my_PUFA[,c("total_o3","total_o6")]
my_PUFA=subset(my_PUFA, my_PUFA$total_o3>0)

my_PUFA$ratio_PUFA=scale(log(my_PUFA$total_o6 / my_PUFA$total_o3))
my_PUFA$total_o6=scale(log(my_PUFA$total_o6))
my_PUFA$total_o3=scale(log(my_PUFA$total_o3))


Try_ratios_p=merge(Try_ratios_p,my_PUFA,by="row.names")
row.names(Try_ratios_p)=Try_ratios_p$Row.names
Try_ratios_p$Row.names=NULL



flag=1
for ( i in 1:55) {
  my_pheno=colnames(Try_ratios_p)[i]
  if (my_pheno %in% regressors) {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 56:ncol(Try_ratios_p)){
      my_trait=colnames(Try_ratios_p)[a]
      my_uni_test=Try_ratios_p[,c(regressors,my_pheno,my_trait)]
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

associations_tryp=my_univariate_results
associations_tryp$phenotype="IBD"
associations_tryp$phenotype[associations_tryp$factor=="ibd_DiagnosisCD"]="CD"
associations_tryp$phenotype[associations_tryp$factor=="ibd_DiagnosisUC"]="UC"
associations_tryp$FDR=p.adjust(associations_tryp$`Pr(>|t|)`,method = "BH")
associations_tryp$Bonferroni=p.adjust(associations_tryp$`Pr(>|t|)`,method = "bonferroni")



