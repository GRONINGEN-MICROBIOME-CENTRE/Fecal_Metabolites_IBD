PCA analysis on metabolites present in more than 70% of the samples
===

norm_filt=q_mtb
pca_log= prcomp(norm_filt,scale. = T, center = T)
cmp_log=as.data.frame(pca_log$x[,1:5])
cmp_log2=merge(cc_pheno2,cmp_log,by="row.names")
v_exp = pca_log$sdev^2
prop_v_exp= data.frame( PC=1:length(v_exp),perc= 100*(v_exp / sum(v_exp)))

#get PCA loading contributors
get_vars=get_pca_var(pca_log)
pca_contribution=as.data.frame(get_vars$contrib)





Percentage of variance explained per component
---
prop_v_exp2=prop_v_exp[1:20,]
ggplot(prop_v_exp2, aes (PC,perc)) + geom_bar (stat="identity") + geom_text(aes(label=round(perc, 2)), size=3, vjust=-.5) + theme_bw() + xlab ("Principal components") + ylab ("% explained variance")

PC1 & PC2: Diagnosis
---
ggplot(cmp_log2, aes(PC1,PC2, fill=ibd_Diagnosis)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_manual(values =c("purple","black","pink2" ,"darkolivegreen3"))


PC2 & PC3: Diagnosis
---
ggplot(cmp_log2, aes(PC2,PC3, fill=ibd_Diagnosis)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_manual(values =c("purple","black","pink2" ,"darkolivegreen3"))


PC1 & PC2: Color by metabolite abundance 
---
cmp_log3=merge(cmp_log2, q_mtb, by.x="Row.names", by.y="row.names")
ggplot(cmp_log3, aes(PC1,PC2, fill=cholate)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")
ggplot(cmp_log3, aes(PC1,-PC2, fill=azelate_C9_DC)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

__Plot with taxa__


metabolite_pca_and_taxa=merge(cmp_log,taxa_filt, by="row.names")
metabolite_pca_and_taxa=merge(metabolite_pca_and_taxa,q_mtb, by.x="Row.names", by.y="row.names")

ggplot(metabolite_pca_and_taxa, aes(PC1,PC2, fill=Faecalibacterium_prausnitzii)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

ggplot(metabolite_pca_and_taxa, aes(PC1,PC2, fill=Coprococcus_catus)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

ggplot(metabolite_pca_and_taxa, aes(PC1,PC2, fill=Blautia_obeum)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

#calculate correlation PC1 & PC2 vs metabolites and taxa
flag=1
for (i in 7:115){
  if (flag==1){
    corr_pc_taxa=data.frame("BUG"=colnames(metabolite_pca_and_taxa)[i], "PC1_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$estimate, "PC1_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$p.value, "PC2_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$estimate, "PC2_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$p.value)
    flag=50
  } else{
    my_corr=data.frame("BUG"=colnames(metabolite_pca_and_taxa)[i], "PC1_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$estimate, "PC1_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$p.value, "PC2_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$estimate, "PC2_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$p.value)
    corr_pc_taxa=rbind(corr_pc_taxa,my_corr)
  }
}


flag=1
for (i in 116:ncol(metabolite_pca_and_taxa)){
  if (flag==1){
    corr_pc_metabolites=data.frame("BUG"=colnames(metabolite_pca_and_taxa)[i], "PC1_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$estimate, "PC1_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$p.value, "PC2_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$estimate, "PC2_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$p.value)
    flag=50
  } else{
    my_corr=data.frame("BUG"=colnames(metabolite_pca_and_taxa)[i], "PC1_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$estimate, "PC1_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC1)$p.value, "PC2_rho" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$estimate, "PC2_p" = cor.test(metabolite_pca_and_taxa[,i],metabolite_pca_and_taxa$PC2)$p.value)
    corr_pc_metabolites=rbind(corr_pc_metabolites,my_corr)
  }
}




5.PCA & PCoA analysis on bacteria present in more than 20% of the samples
===

__Calculate PCA based on filtered clr relative abundaces__


print(paste("Number of species used for PCoA analysis:", ncol(taxa_filt)))

#Euclidean distances on clr transformed values = Aitchison distance
my_dist=vegdist(taxa_filt, method = "euclidean")
#PCoA analyses
mypcoa=cmdscale(my_dist, k = 5)
colnames(mypcoa)=c("PC1","PC2","PC3","PC4","PC5")
mypcoa2=cmdscale(my_dist, k = 5, eig = T)
round(mypcoa2$eig*100/sum(mypcoa2$eig),2)[1:5]
cmp_clr3=merge(cc_pheno2,mypcoa,by="row.names")

ggplot(cmp_clr3, aes(PC1,PC2, fill=ibd_Diagnosis)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_manual(values =c("purple","black","pink2" ,"darkolivegreen3"))

6.Calculate dysbiotic score
===

#Calculate dysbiotic score: adapted from Lloyd-Price et al. Nature 2019 (https://bitbucket.org/biobakery/hmp2_analysis/src/master/common/disease_activity.r)

my_dist_matrix=as.matrix(my_dist)
#ref_set=rownames(cc_pheno2)[cc_pheno2$ibd_Diagnosis=="Control"]
#Take non-IBD controls as a reference group to calculate the median distance
ref_set=my_dist_matrix[,rownames(cc_pheno2)[cc_pheno2$ibd_Diagnosis=="Control"]]
dysbiosis_score=as.data.frame(ref_set)
dysbiosis_score_na=dysbiosis_score
#To avoid calculating the distance of a sample including the sample itself replace 0 distances for NAs
dysbiosis_score_na[dysbiosis_score_na==0]=NA
#Calculate median distance of each sample relative to the non-IBD controls. 
dysbiosis_score_na$row_median = apply(dysbiosis_score, 1, median,na.rm=TRUE)
dysbiosis_score2=dysbiosis_score_na[,"row_median", drop=F]

#Merge score with phenotypes
dysbiosis_score_pheno=merge(cc_pheno2,dysbiosis_score2, by="row.names")
dysbiosis_score_pheno$IBD=dysbiosis_score_pheno$ibd_Diagnosis
dysbiosis_score_pheno$IBD[dysbiosis_score_pheno$ibd_Diagnosis!="Control"]="IBD"

#Define dysbiosis as the 95th quantile of the median distance of non-IBD samples 
dysbiosis_threshold=quantile(dysbiosis_score_pheno$row_median[dysbiosis_score_pheno$IBD=="Control"], 0.95)
dysbiosis_score_pheno$dysbiotic="no"
dysbiosis_score_pheno$dysbiotic[dysbiosis_score_pheno$row_median>=dysbiosis_threshold]="yes"

row.names(dysbiosis_score_pheno)=dysbiosis_score_pheno$Row.names
dysbiosis_score_pheno$Row.names=NULL

#Merge with PCoA coordinates and plot
cmp_clr3.d=merge(dysbiosis_score_pheno,mypcoa,by="row.names")
ggplot(cmp_clr3.d, aes(PC1,PC2, fill=IBD)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_manual(values =c("gray80","red"))
ggplot(cmp_clr3.d, aes(PC1,PC2, fill=dysbiotic)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_manual(values =c("gray80","purple"))


7. Clustering metabolites
===


Method 1: Regress disease and (potential) technical confouders and perform spearman correlation
---

__Regress batch & diseases effects__

Here we use all metabolites present in more than 20% of the samples. 

We extract residuals from the following linear model (lm) : 

__lm(Metabolite ~ Day_of_measurment  + Amount_of_input_material + LC Column + Months_sample_in_freezer + ibd_Diagnosis (Control/CD/UC) )__



#Regress cofounders
rownames(cc_pheno4)=cc_pheno4$Row.names
cc_pheno4$Row.names=NULL
cc_pheno5=merge(my_batch_id,cc_pheno4, by="row.names")
#my_regress_phenos=cc_pheno5[,c("Row.names","run_day_cat", "Amount_sample_gram", "LC.COLUMN" , "metabolon_Month_in_freezer" , "ibd_Diagnosis")]
my_regress_phenos2=read.table("~/Desktop/IBD_metabolomics_2022/1.Input/mini_phenotypes.txt", sep = "\t", header = T, row.names = 1)
my_regress_phenos3=merge(my_regress_phenos2,run_day, by="row.names")
rownames(my_regress_phenos3)=my_regress_phenos3$Row.names
my_regress_phenos3$Row.names=NULL

rgs_mb1=merge(my_regress_phenos3,q_mtb, by="row.names")
rownames(rgs_mb1)=rgs_mb1$Row.names
rgs_mb1$Row.names=NULL
colnames(rgs_mb1)=make.names(colnames(rgs_mb1))
regressors=colnames(my_regress_phenos3)

data_clust3=matrix(nrow=nrow(my_regress_phenos3), ncol = ncol(q_mtb))
n=1
for (i in 30:ncol(rgs_mb1)){
  my_trait=colnames(rgs_mb1)[i]
  my_uni_test=rgs_mb1[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  data_clust3[,n] = my_lm$residuals
 n=n+1
}
data_clust3=as.data.frame(data_clust3)
rownames(data_clust3)=rgs_mb1$Row.names
colnames(data_clust3)=colnames(rgs_mb1)[30:ncol(rgs_mb1)]



#Calculate clusters
my_corr=cor(data_clust3, method = "spearman")
my_dist=as.dist(1-my_corr)
metabolite_clusters=hclust(my_dist, method = "ward.D" )
metabolite_cluster_def=cutree(metabolite_clusters, h=1)
data_clust_id=as.data.frame(t(data_clust3))
data_clust_id=data_clust_id %>% mutate (cluster=metabolite_cluster_def) 
data_clust_id=data_clust_id[,c("cluster"), drop=F]

#Explore regions of the dendogram
my_dendo=as.dendrogram(metabolite_clusters)
plot(cut(my_dendo,h=10)$lower[[1]], cex=0.2)
plot(cut(my_dendo,h=10)$lower[[5]], cex=0.2,horiz=T)

#Figure

annot2=annot
annot2$SUB.PATHWAY=NULL
annot2$SUPER.PATHWAY=gsub(" ","_", annot2$SUPER.PATHWAY)

my_path_annot= list(
  SUPER.PATHWAY = c(Amino_Acid="firebrick3", Carbohydrate = "steelblue2", Cofactors_and_Vitamins = "gold3", Energy = "limegreen", Lipid = "grey77", Nucleotide = "pink1", Peptide = "salmon", SCFA="turquoise", Unknown_compount="black", Xenobiotics="azure"))

pheatmap(my_corr, fontsize_col = 3, fontsize_row  = 3, annotation_row = annot2, clustering_method = "ward.D", scale = "none", show_colnames = F, annotation_colors = my_path_annot)




