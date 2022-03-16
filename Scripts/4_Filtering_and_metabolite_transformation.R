3. Filter metabolites by missigness
===

#Split per platform

my_batch_mtb$PLATFORM= make.names(my_batch_mtb$PLATFORM)
colnames(all_new_ID_raw)=gsub(",","_",colnames(all_new_ID_raw))
summary_short2$names=gsub(",","_",summary_short2$Metabolite)


id_LCMS_NEG=subset(my_batch_mtb, my_batch_mtb$PLATFORM=="LC.MS.Neg")
id_LSMS_Pol=subset(my_batch_mtb, my_batch_mtb$PLATFORM=="LC.MS.Polar")
id_LCMS_PosE=subset(my_batch_mtb, my_batch_mtb$PLATFORM=="LC.MS.Pos.Early")
id_LCMS_PosL=subset(my_batch_mtb, my_batch_mtb$PLATFORM=="LC.MS.Pos.Late")


all_new_ID_raw_zeros=all_new_ID_raw
#all_new_ID_raw_zeros[,is.na(all_new_ID_raw_zeros)]=0

all_LCMS_NEG=all_new_ID_raw_zeros[,colnames(all_new_ID_raw_zeros) %in% row.names(id_LCMS_NEG)]
all_LSMS_Pol=all_new_ID_raw_zeros[,colnames(all_new_ID_raw_zeros) %in% row.names(id_LSMS_Pol)]
all_LCMS_PosE=all_new_ID_raw_zeros[,colnames(all_new_ID_raw_zeros) %in% row.names(id_LCMS_PosE)]
all_LCMS_PosL=all_new_ID_raw_zeros[,colnames(all_new_ID_raw_zeros) %in% row.names(id_LCMS_PosL)]


#Perform CLR transformation

tr_LCMS_NEG=transform_and_filter_mtb(all_LCMS_NEG,samples_row = T,method = "clr",missing_filter = -1)
tr_LCMS_Pol=transform_and_filter_mtb(all_LSMS_Pol,samples_row = T,method = "clr",missing_filter = -1)
tr_LCMS_PosE=transform_and_filter_mtb(all_LCMS_PosE,samples_row = T,method = "clr",missing_filter = -1)
tr_LCMS_PosL=transform_and_filter_mtb(all_LCMS_PosL,samples_row = T,method = "clr",missing_filter = -1)


#Merge

all_mtb_clr=bind_cols(tr_LCMS_NEG,tr_LCMS_Pol,tr_LCMS_PosE,tr_LCMS_PosL)

__We remove metabolites that are present in less than 20% of the samples in the 3 cohorts (CD,UC,Controls)__


to_remove=subset(summary_short2,prevalence_Control<20 & prevalence_CD<20 & prevalence_UC<20)
q_mtb_v2=all_new_ID_raw[,colnames(all_new_ID_raw) %ni%to_remove$names]
print(paste("Number of metabolites removed:",length(to_remove$names)))


__Metabolites present in more than 70% of the samples will be analyzed quantitatively__



#to_keep=subset(summary_short2,prevalence_Control>70 & prevalence_CD>70 & prevalence_UC>70)
to_keep=subset(summary_short2,prevalence_all>70)
q_mtb_t=all_mtb_clr[,colnames(all_mtb_clr) %in%to_keep$names]
print(paste("Number of metabolites for quantitative analysis:",length(to_keep$names)))

q_mtb=impute_missing_values(q_mtb_t,samples_row = T, method = "knn", missing_filter = 0)
row.names(q_mtb)=row.names(q_mtb_t)

__Metabolites that are present between 20% to 70% of the samples will be subset and transform to binary feature (present/absent)__


rest=c(to_remove$names,to_keep$names)
bi_mtb=all_new_ID_raw[,colnames(all_new_ID_raw) %ni% rest]
bi_mtb[!is.na(bi_mtb)]=1
bi_mtb[is.na(bi_mtb)]=0

print(paste("Number of metabolites for presence/absence analysis:",length(colnames(bi_mtb))))




__Filter bacteria (detection in >20% of the samples) and normalization__


#Subset samples from the taxonomic table for which we also have metabolomics 
unk=t(my_taxa[1,,drop=F])
taxa2=merge(unk,taxa,by="row.names")
row.names(taxa2)=taxa2$Row.names
taxa2$Row.names=NULL
taxa2=subset(taxa2,rownames(taxa2) %in% rownames(cc_pheno2))
#taxa_gn2=taxa_gn[,colnames(taxa_gn) %in% rownames(cc_pheno2)]

print (paste("Number of taxa detected:",ncol(taxa2)))

#Filter 20% and to clr transformation
taxa_filt=transform_and_filter_taxa(taxa2,samples_row = T,method = "clr",missing_filter = 20)
taxa_filt$UNKNOWN=NULL
#taxa_gn_filt=transform_and_filter_taxa(t(taxa_gn2),samples_row = T,method = "clr",missing_filter = 20)

print (paste("Number of taxa after filtering",ncol(taxa_filt)))


#Merge taxa and phenotypes
cc_pheno3=merge(cc_pheno2,taxa_filt,by="row.names")
rownames(cc_pheno3)=cc_pheno3$Row.names
cc_pheno3$Row.names=NULL

__Merge taxa and metabolites__


print("Missing metagenomic data:")
setdiff( rownames(q_mtb), rownames(taxa_filt))
taxa_and_mb=merge(taxa_filt,q_mtb,by="row.names")
rownames(taxa_and_mb)=taxa_and_mb$Row.names
taxa_and_mb$Row.names=NULL



__Merge phenotypes to SCFA__


SCFA2=SCFA_new_ID
SCFA2[,3:10]=sapply(SCFA2[,c(3:10)], function(x) log(as.numeric(as.character(x))) )
SCFA_with_phenos=merge(SCFA_new_ID,cc_pheno3, by="row.names")
cc_pheno4=merge(SCFA2,cc_pheno3, by="row.names")

SCFAsum=SCFA_new_ID[,3:10, drop=F]
SCFAsum=sapply(SCFAsum, function(x) as.numeric(as.character(x)) )
row.names(SCFAsum)=row.names(SCFA_new_ID)
SCFAsum=SCFAsum[row.names(SCFAsum) %in% row.names(q_mtb),]
summary_SCFA=summary_stats(SCFAsum, select, missing_vals = "NA")




__Plot SCFA per diagnosis__


plot_SCFA=SCFA_with_phenos[,c("Row.names","ibd_Diagnosis","Methylbutyric_acid", "acetic_acid", "butyric_acid", "hexanoic_acid", "isobutyric_acid","isovaleric_acid" ,"propionic_acid", "valeric_acid")]

plot_SCFA$diag=as.character(plot_SCFA$ibd_Diagnosis)
plot_SCFA$diag[plot_SCFA$ibd_Diagnosis=="Control"]="A_Control"
plot_SCFA=subset(plot_SCFA, plot_SCFA$ibd_Diagnosis!="IBDU")


ggplot (plot_SCFA, aes (diag,log(as.numeric(plot_SCFA$valeric_acid)), fill=diag)) +geom_violin(alpha=0.8) + geom_jitter(height = 0, width = 0.1, alpha=0.4) + geom_boxplot(width=0.1,position=position_dodge(1), fill="white")+ theme_bw() +scale_fill_manual(values=c("grey80","purple","darkolivegreen3"))  + ylab ("log(Valeric acid µg/g)") + xlab ("Cohort") + scale_x_discrete(labels=c("(n=255)","(n=239)", "(n=176)"))


