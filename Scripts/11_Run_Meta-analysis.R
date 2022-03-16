Run meta-analysis
====

library(meta)


a. Phenotypes
===

Associations_Controls <- read.delim("~/Desktop/IBD_metabolomics_2022/phenos_metabolites_controls.txt")
Associations_CD <- read.delim("~/Desktop/IBD_metabolomics_2022/phenos_metabolites_CD.txt")
Associations_UC <- read.delim("~/Desktop/IBD_metabolomics_2022/phenos_metabolites_UC.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","t_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","t_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Metabolite","Estimate_CNT","SE_CNT","t_value_CNT","pval_CNT","phenotype_CNT", "factor_CNT", "FDR_CNT", "Bonferroni_CNT", "Super_pathway", "Sub_pathway")

Associations_Controls_pheno=Associations_Controls[!grepl("sp_",Associations_Controls$phenotype_CNT),]

Associations_CD=subset(Associations_CD,Associations_CD$phenotype_CD!="clinical_FecalCalprotectinOver200yesno")
Associations_UC=subset(Associations_UC,Associations_UC$phenotype_UC!="clinical_FecalCalprotectinOver200yesno")
Associations_CD=subset(Associations_CD,Associations_CD$phenotype_CD!="sp_Shannon_Index")
Associations_UC=subset(Associations_UC,Associations_UC$phenotype_UC!="sp_Shannon_Index")
Associations_Controls_pheno$phenotype_CNT=gsub("host_smk_current","host_SmokeCurrentSmoker", Associations_Controls_pheno$phenotype_CNT)
Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")
Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")
Associations_Controls_pheno$FDR=p.adjust(Associations_Controls_pheno$pval_CNT, method="BH")
Associations_Controls_pheno$Bonferroni=p.adjust(Associations_Controls_pheno$pval_CNT, method="bonferroni")
Associations_Controls_pheno$factor=paste(Associations_Controls_pheno$Metabolite, Associations_Controls_pheno$phenotype_CNT, sep="_")
pre_table=merge(Associations_Controls_pheno, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#33
table_meta$I2=NA
#34
table_meta$Q=NA
#35
table_meta$pval_Q=NA
#36
table_meta$z_fixed=NA
#37
table_meta$pval_fixed=NA
#38
table_meta$z_random=NA
#39
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,3]) & !is.na(table_meta[i,15]) & !is.na(table_meta[i,24])){
    mytest=metagen(TE=c(table_meta[i,3],table_meta[i,15],table_meta[i,24]), seTE =c(table_meta[i,4],table_meta[i,16],table_meta[i,25]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,33]=mytest$I2
    table_meta[i,34]=mytest$Q
    table_meta[i,35]=mytest$pval.Q
    table_meta[i,36]=mytest$zval.fixed
    table_meta[i,37]=mytest$pval.fixed
    table_meta[i,38]=mytest$zval.random
    table_meta[i,39]=mytest$pval.random
  }
}
table_meta_v2$FDR_meta_fixed=p.adjust(table_meta_v2$pval.fixed, method = "BH")
table_meta_v2$FDR_meta_random=p.adjust(table_meta_v2$pval.random, method = "BH")
write.table(table_meta_v2, "~/Desktop/IBD_Metabolomics/3.Output_statistics/phenos_metabolites/meta-analysis_phenotypes.txt", sep = "\t")


b. Phenotypes less prevalent metabolites
===

Associations_Controls <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/phenos_metabolites/Controls_phenos_prevalence.txt .txt")
Associations_CD <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/phenos_metabolites/CD_phenos_metabolites_prevalence.txt.txt")
Associations_UC <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/phenos_metabolites/UC_phenos_metabolites_prevalence.txt.txt")
colnames(Associations_UC)=c("Estimate_UC","SE_UC","z_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","z_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","z_value_CNT","pval_CNT","metabolite_CNT","phenotype_CNT", "factor_CNT", "FDR_CNT", "Bonferroni_CNT")

Associations_Controls_pheno=Associations_Controls[!grepl("sp_",Associations_Controls$phenotype_CNT),]
Associations_CD=subset(Associations_CD,Associations_CD$phenotype_CD!="clinical_FecalCalprotectinOver200yesno")
Associations_UC=subset(Associations_UC,Associations_UC$phenotype_UC!="clinical_FecalCalprotectinOver200yesno")
Associations_CD=subset(Associations_CD,Associations_CD$phenotype_CD!="sp_Shannon_Index")
Associations_UC=subset(Associations_UC,Associations_UC$phenotype_UC!="sp_Shannon_Index")
Associations_Controls_pheno$phenotype_CNT=gsub("host_smk_current","host_SmokeCurrentSmoker", Associations_Controls_pheno$phenotype_CNT)
Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")
Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")
Associations_Controls_pheno$FDR=p.adjust(Associations_Controls_pheno$pval_CNT, method="BH")
Associations_Controls_pheno$Bonferroni=p.adjust(Associations_Controls_pheno$pval_CNT, method="bonferroni")
Associations_Controls_pheno$factor=paste(Associations_Controls_pheno$metabolite_CNT, Associations_Controls_pheno$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls_pheno, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#31
table_meta$I2=NA
#32
table_meta$Q=NA
#33
table_meta$pval_Q=NA
#34
table_meta$z_fixed=NA
#35
table_meta$pval_fixed=NA
#36
table_meta$z_random=NA
#37
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,13]) & !is.na(table_meta[i,22])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,13],table_meta[i,22]), seTE =c(table_meta[i,3],table_meta[i,14],table_meta[i,23]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,31]=mytest$I2
    table_meta[i,32]=mytest$Q
    table_meta[i,33]=mytest$pval.Q
    table_meta[i,34]=mytest$zval.fixed
    table_meta[i,35]=mytest$pval.fixed
    table_meta[i,36]=mytest$zval.random
    table_meta[i,37]=mytest$pval.random
  }
}

table_meta_prevalence$FDR_meta_fixed=p.adjust(table_meta_prevalence$pval.fixed, method = "BH")
table_meta_prevalence$FDR_meta_random=p.adjust(table_meta_prevalence$pval.random, method = "BH")
write.table(table_meta_prevalence, "~/Desktop/IBD_Metabolomics/3.Output_statistics/phenos_metabolites/meta-analysis_phenotypes_prevalence.txt", sep = "\t")


c. Species (presence/absence)
====

Associations_Controls <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_quant_cnt.txt")
Associations_CD <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_quant_cd.txt")
Associations_UC <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_quant_uc.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","t_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","t_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","t_value_CNT","pval_CNT","metabolite_CNT", "phenotype_CNT","factor_CNT")

Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")

Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")

Associations_Controls$FDR_CNT=p.adjust(Associations_Controls$pval_CNT, method="BH")
Associations_Controls$Bonferroni_CNT=p.adjust(Associations_Controls$pval_CNT, method="bonferroni")
Associations_Controls$factor=paste(Associations_Controls$metabolite_CNT, Associations_Controls$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)


#29
table_meta$I2=NA
table_meta$Q=NA
table_meta$pval_Q=NA
table_meta$z_fixed=NA
table_meta$te_fixed=NA
table_meta$pval_fixed=NA
table_meta$z_random=NA
table_meta$te_random=NA
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,11]) & !is.na(table_meta[i,20])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,11],table_meta[i,20]), seTE =c(table_meta[i,3],table_meta[i,12],table_meta[i,21]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,29]=mytest$I2
    table_meta[i,30]=mytest$Q
    table_meta[i,31]=mytest$pval.Q
    table_meta[i,32]=mytest$zval.fixed
    table_meta[i,33]=mytest$TE.fixed
    table_meta[i,34]=mytest$pval.fixed
    table_meta[i,35]=mytest$zval.random
    table_meta[i,36]=mytest$TE.random
    table_meta[i,37]=mytest$pval.random
  }
}

table_meta_sp_prev=table_meta

table_meta_sp_prev$FDR_meta_fixed=p.adjust(table_meta_sp_prev$pval_fixed, method = "BH")
table_meta_sp_prev$FDR_meta_random=p.adjust(table_meta_sp_prev$pval_random, method = "BH")
write.table(table_meta_sp_prev, "~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/meta-analysis_taxa_prevalence_metabolites", sep = "\t")



d. Species abundance
====

Associations_Controls <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_nz_cnt.txt")
Associations_CD <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_nz_cd.txt")
Associations_UC <- read.delim("~/Desktop/IBD_metabolomics_2022/sp_metabolites_nz_uc.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","t_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","t_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","t_value_CNT","pval_CNT","metabolite_CNT", "phenotype_CNT","factor_CNT")

Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")

Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")

Associations_Controls$FDR_CNT=p.adjust(Associations_Controls$pval_CNT, method="BH")
Associations_Controls$Bonferroni_CNT=p.adjust(Associations_Controls$pval_CNT, method="bonferroni")
Associations_Controls$factor=paste(Associations_Controls$metabolite_CNT, Associations_Controls$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#29
table_meta$I2=NA
table_meta$Q=NA
table_meta$pval_Q=NA
table_meta$z_fixed=NA
table_meta$te_fixed=NA
table_meta$pval_fixed=NA
table_meta$z_random=NA
table_meta$te_random=NA
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,11]) & !is.na(table_meta[i,20])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,11],table_meta[i,20]), seTE =c(table_meta[i,3],table_meta[i,12],table_meta[i,21]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,29]=mytest$I2
    table_meta[i,30]=mytest$Q
    table_meta[i,31]=mytest$pval.Q
    table_meta[i,32]=mytest$zval.fixed
    table_meta[i,33]=mytest$TE.fixed
    table_meta[i,34]=mytest$pval.fixed
    table_meta[i,35]=mytest$zval.random
    table_meta[i,36]=mytest$TE.random
    table_meta[i,37]=mytest$pval.random
  }
}

table_meta_sp_linear=table_meta

table_meta_sp_linear$FDR_meta_fixed=p.adjust(table_meta_sp_linear$pval_fixed, method = "BH")
table_meta_sp_linear$FDR_meta_random=p.adjust(table_meta_sp_linear$pval_random, method = "BH")
write.table(table_meta_sp_linear, "~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/meta-analysis_taxa_abundance_metabolites", sep = "\t")


e. Species prevalence association to metabolite prevalence [less abundant metabolites n<70%]
===

Associations_Controls <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/Controls_taxa_prevalence_metabolites_prevalence.txt")
Associations_CD <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/CD_taxa_prevalence_metabolites_prevalence.txt")
Associations_UC <- read.delim("~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/UC_taxa_prevalence_metabolites_prevalence.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","z_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","z_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","z_value_CNT","pval_CNT","metabolite_CNT", "phenotype_CNT","factor_CNT")

Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")

Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")

Associations_Controls$FDR_CNT=p.adjust(Associations_Controls$pval_CNT, method="BH")
Associations_Controls$Bonferroni_CNT=p.adjust(Associations_Controls$pval_CNT, method="bonferroni")
Associations_Controls$factor=paste(Associations_Controls$metabolite_CNT, Associations_Controls$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#29
table_meta$I2=NA
table_meta$Q=NA
table_meta$pval_Q=NA
table_meta$z_fixed=NA
table_meta$pval_fixed=NA
table_meta$z_random=NA
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,11]) & !is.na(table_meta[i,20])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,11],table_meta[i,20]), seTE =c(table_meta[i,3],table_meta[i,12],table_meta[i,21]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,29]=mytest$I2
    table_meta[i,30]=mytest$Q
    table_meta[i,31]=mytest$pval.Q
    table_meta[i,32]=mytest$zval.fixed
    table_meta[i,33]=mytest$pval.fixed
    table_meta[i,34]=mytest$zval.random
    table_meta[i,35]=mytest$pval.random
  }
}

table_meta_sp_prev_mtb_prev=table_meta

table_meta_sp_prev_mtb_prev$FDR_meta_fixed=p.adjust(table_meta_sp_prev_mtb_prev$pval_fixed, method = "BH")
table_meta_sp_prev_mtb_prev$FDR_meta_random=p.adjust(table_meta_sp_prev_mtb_prev$pval_random, method = "BH")
write.table(table_meta_sp_prev_mtb_prev, "~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/meta-analysis_taxa_prevalence_metabolites_prevalence.txt", sep = "\t")





f. MetaCyc pathways association to metabolite prevalence [less abundant metabolites n<70%]
===

Associations_Controls <- read.delim("~/Desktop/IBD_metabolomics_2022/2.Output/Metabolites_MetaCyc_CNT.txt")
Associations_CD <- read.delim("~/Desktop/IBD_metabolomics_2022/2.Output/Metabolites_MetaCyc_CD.txt")
Associations_UC <- read.delim("~/Desktop/IBD_metabolomics_2022/2.Output/Metabolites_MetaCyc_UC.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","t_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","t_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","t_value_CNT","pval_CNT","metabolite_CNT", "phenotype_CNT","factor_CNT")

Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")

Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")

Associations_Controls$FDR_CNT=p.adjust(Associations_Controls$pval_CNT, method="BH")
Associations_Controls$Bonferroni_CNT=p.adjust(Associations_Controls$pval_CNT, method="bonferroni")
Associations_Controls$factor=paste(Associations_Controls$metabolite_CNT, Associations_Controls$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#29
table_meta$I2=NA
table_meta$Q=NA
table_meta$pval_Q=NA
table_meta$z_fixed=NA
table_meta$pval_fixed=NA
table_meta$z_random=NA
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,11]) & !is.na(table_meta[i,20])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,11],table_meta[i,20]), seTE =c(table_meta[i,3],table_meta[i,12],table_meta[i,21]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,29]=mytest$I2
    table_meta[i,30]=mytest$Q
    table_meta[i,31]=mytest$pval.Q
    table_meta[i,32]=mytest$zval.fixed
    table_meta[i,33]=mytest$pval.fixed
    table_meta[i,34]=mytest$zval.random
    table_meta[i,35]=mytest$pval.random
  }
}

table_meta_sp_prev_mtb_prev=table_meta

table_meta_sp_prev_mtb_prev$FDR_meta_fixed=p.adjust(table_meta_sp_prev_mtb_prev$pval_fixed, method = "BH")
table_meta_sp_prev_mtb_prev$FDR_meta_random=p.adjust(table_meta_sp_prev_mtb_prev$pval_random, method = "BH")
write.table(table_meta_sp_prev_mtb_prev, "~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/meta-analysis_taxa_prevalence_metabolites_prevalence.txt", sep = "\t")






g. Metabolic gene clusters
===


Associations_Controls <- read.delim("~/Desktop/IBD_metabolomics_2022/mgc_metabolites_cnt.txt")
Associations_CD <- read.delim("~/Desktop/IBD_metabolomics_2022/mgc_metabolites_cd.txt")
Associations_UC <- read.delim("~/Desktop/IBD_metabolomics_2022/mgc_metabolites_uc.txt")

colnames(Associations_UC)=c("Estimate_UC","SE_UC","z_value_UC","pval_UC","metabolite_UC","phenotype_UC", "factor_UC")
colnames(Associations_CD)=c("Estimate_CD","SE_CD","z_value_CD","pval_CD","metabolite_CD","phenotype_CD", "factor_CD")
colnames(Associations_Controls)=c("Estimate_CNT","SE_CNT","z_value_CNT","pval_CNT","metabolite_CNT", "phenotype_CNT","factor_CNT")

Associations_UC$FDR_UC=p.adjust(Associations_UC$pval_UC, method="BH")
Associations_UC$Bonferroni_UC=p.adjust(Associations_UC$pval_UC, method="bonferroni")
Associations_UC$factor=paste(Associations_UC$metabolite_UC, Associations_UC$phenotype_UC, sep="_")

Associations_CD$FDR_CD=p.adjust(Associations_CD$pval_CD, method="BH")
Associations_CD$Bonferroni_CD=p.adjust(Associations_CD$pval_CD, method="bonferroni")
Associations_CD$factor=paste(Associations_CD$metabolite_CD, Associations_CD$phenotype_CD, sep="_")

Associations_Controls$FDR_CNT=p.adjust(Associations_Controls$pval_CNT, method="BH")
Associations_Controls$Bonferroni_CNT=p.adjust(Associations_Controls$pval_CNT, method="bonferroni")
Associations_Controls$factor=paste(Associations_Controls$metabolite_CNT, Associations_Controls$phenotype_CNT, sep="_")

pre_table=merge(Associations_Controls, Associations_CD, by="factor", all = T)
table_meta=merge(pre_table, Associations_UC, by="factor", all = T)

#29
table_meta$I2=NA
table_meta$Q=NA
table_meta$pval_Q=NA
table_meta$z_fixed=NA
table_meta$pval_fixed=NA
table_meta$z_random=NA
table_meta$pval_random=NA

for (i in 1:nrow(table_meta)){
  if(!is.na(table_meta[i,2]) & !is.na(table_meta[i,11]) & !is.na(table_meta[i,20])){
    mytest=metagen(TE=c(table_meta[i,2],table_meta[i,11],table_meta[i,20]), seTE =c(table_meta[i,3],table_meta[i,12],table_meta[i,21]), comb.fixed=TRUE, comb.random=TRUE) 
    table_meta[i,29]=mytest$I2
    table_meta[i,30]=mytest$Q
    table_meta[i,31]=mytest$pval.Q
    table_meta[i,32]=mytest$zval.fixed
    table_meta[i,33]=mytest$pval.fixed
    table_meta[i,34]=mytest$zval.random
    table_meta[i,35]=mytest$pval.random
  }
}

table_meta_mgc_mtb=table_meta

table_meta_mgc_mtb$FDR_meta_fixed=p.adjust(table_meta_mgc_mtb$pval_fixed, method = "BH")
table_meta_mgc_mtb$FDR_meta_random=p.adjust(table_meta_mgc_mtb$pval_random, method = "BH")
write.table(table_meta_mgc_mtb, "~/Desktop/IBD_Metabolomics/3.Output_statistics/species_metabolite_models_nov2021/meta-analysis_mgc-metabolites.txt", sep = "\t")





