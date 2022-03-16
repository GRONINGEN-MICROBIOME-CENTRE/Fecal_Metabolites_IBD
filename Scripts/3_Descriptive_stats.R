Summary statistics
===


Annotation by METABOLON
---

Describe the categories of metabolites
---


met_2=as.data.frame(table(annot$SUPER.PATHWAY))
met_1=as.data.frame(table(annot$SUB.PATHWAY))

ggplot(met_2, aes(x = reorder(Var1, -Freq),Freq)) + geom_bar(stat = "identity", color="black") + theme_classic() + theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 12)) + coord_flip() + xlab("") + ylab("Number of metabolites")




ggplot(met_1, aes(x = reorder(Var1, -Freq),Freq)) + geom_bar(stat = "identity", color="black") + theme_classic() + theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 5)) + coord_flip() + xlab("") + ylab("Number of metabolites")



Remove participants with a pouch or stoma
---

__Rational: Patients with pouch or stoma do not capture the colonic metabolite/microbiota profile__


remove_pouch_stoma=rownames(phenos_ibd)[phenos_ibd$ibd_CurrentStomaOrPouch=="yes"]
#Keep data of the participants with stoma
all_new_ID_raw_with_stoma=all_new_ID_raw
#Remove (n=68)
all_new_ID_raw=all_new_ID_raw[row.names(all_new_ID_raw)%ni%remove_pouch_stoma,]
cc_pheno=cc_pheno[row.names(cc_pheno)%ni%remove_pouch_stoma,]
print (paste("Number of samples excluded due to pouch or stoma:", length(remove_pouch_stoma)))
print (paste("Number of samples left:", length(rownames(all_new_ID_raw))))



Summary stats per metabolite per cohort
---



#Raw data
select=cc_pheno[,1, drop=F]
select$ibd_Diagnosis=as.factor(select$ibd_Diagnosis)
summary_met=summary_stats(all_new_ID_raw,select, missing_vals = "NA")
summary_met=as.data.frame(summary_met)



__Check prevalence per metabolite and per phenotype (CD/UC/non-IBD)__


summary_short=summary_met[,c("Metabolite", "Non_NAs")]
summary_short2=summary_met[,c("Metabolite","Non_NAs","Non_NAs_Control","Non_NAs_CD", "Non_NAs_UC")]
summary_short2$prevalence_all=((as.numeric(as.character(summary_short2$Non_NAs))/nrow(all_new_ID_raw))*100)
summary_short2$prevalence_Control=((as.numeric(as.character(summary_short2$Non_NAs_Control))/length(select[select$ibd_Diagnosis=="Control",])*100))
summary_short2$prevalence_CD=((as.numeric(as.character(summary_short2$Non_NAs_CD))/length(select[select$ibd_Diagnosis=="CD",])*100))
summary_short2$prevalence_UC=((as.numeric(as.character(summary_short2$Non_NAs_UC))/length(select[select$ibd_Diagnosis=="UC",])*100))
summary_short2$Non_NAs=NULL
summary_short2$Non_NAs_Control=NULL
summary_short2$Non_NAs_CD=NULL
summary_short2$Non_NAs_UC=NULL
summary_short2[,-1] <-round(summary_short2[,-1],0)
summary_short3=melt(summary_short2)


__Presence of each metabolite per sample__


ggplot(summary_short, aes(reorder(Metabolite,-as.numeric(as.character(Non_NAs))),as.numeric(as.character(Non_NAs))))+ geom_bar(color="darkblue",stat = "identity") + theme_classic() + ylab("Number of samples (n=682)") + xlab("Metabolites") +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())



__Percentage of prevalence of metabolites per cohort__


ggplot(summary_short3, aes(value, fill=variable)) + geom_histogram(position = "dodge", bins = 10, binwidth = 7, color="black") + theme_bw() + scale_fill_manual(values = c("grey79","white","blue2", "firebrick2")) + ylab("Number of metabolites") + xlab("% presence in samples")



__Prevalence/Missigness per metabolite per cohort with annotations__


summary_short4=summary_short2
rownames(summary_short4)=summary_short4$Metabolite
summary_short4$Metabolite=NULL
summary_short4=merge(annot,summary_short4, by="row.names")
rownames(summary_short4)=summary_short4$Row.names
summary_short4$Row.names=NULL
colnames(summary_short4)=c("Super pathway","Annotation", "Detected in cohort (%)", "Detected in controls (%)", "Detected in CD (%)", "Dectected in UC (%)")
datatable(summary_short4)
#write.table(summary_short4, "~/Desktop/Metabolomics_v2/3.Preliminary_results/Metabolites_prevalence_short_2020.txt", sep = "\t", quote = F)



__Number of metabolites per samples__



all_new_ID_raw2=as.data.frame(t(all_new_ID_raw))
summary_sample=matrix(nrow=ncol(all_new_ID_raw2), ncol=3)

for (i in 1:ncol(all_new_ID_raw2)){
  summary_sample[i,2]=sum(is.na(all_new_ID_raw2[,i]))
  summary_sample[i,3]=sum(!is.na(all_new_ID_raw2[,i]))
}

summary_sample[,1]=colnames(all_new_ID_raw2)
summary_sample=as.data.frame(summary_sample)
colnames(summary_sample)=c("Metabolite", "NAs","Non_NAs")
summary_sample$perc_detected=((as.numeric(as.character(summary_sample$Non_NAs))/nrow(all_new_ID_raw2))*100)
rownames(summary_sample)=summary_sample$Metabolite
summary_sample$Metabolite=NULL
summary_sample2=merge(summary_sample,cc_pheno,by="row.names")
#write.table(summary_sample, "~/Desktop/Metabolomics_v2/3.Preliminary_results/metabolites_per_samples.txt", sep = "\t", quote = F)


ggplot(summary_sample2, aes(perc_detected, fill=ibd_Diagnosis)) + geom_histogram(position = "dodge", bins = 100, binwidth = 7, color="black") + theme_bw() + scale_fill_manual(values = c("firebrick2","white","grey42","blue2")) + ylab("Number of samples") + xlab("% detected metabolites")



__Summary statistics microbiome__

1. Focus only at species taxa level


2. Filter out taxa present in less than 20% of all samples



#Select one taxa level from rownames and clean

# annot_sp$Genus=gsub(".*g__","",row.names(annot_sp))
# annot_sp$Genus=sapply(strsplit(annot_sp$Genus, "|", fixed=TRUE), head, 1)
# annot_sp$Family=gsub(".*f__","",row.names(annot_sp))
# annot_sp$Family=sapply(strsplit(annot_sp$Family, "|", fixed=TRUE), head, 1)


#Remove the percentage of unaligned reads
taxa=my_taxa
preFilt=data.frame(colnames(taxa))
preFilt$In_MGS="Yes"

#Remove samples that failed in the microbiome prediction: IBDFEC0535, IBDFEC1024,IBDFEC0515
taxa=taxa[, colSums(taxa != 0) > 0]

#Keep track of the samples that failed microbiome prediction 
preFilt[preFilt$colnames.taxa.%ni%colnames(taxa),]$In_MGS="Failed"
rownames(preFilt)=preFilt$colnames.taxa.
preFilt$colnames.taxa.=NULL
tract_id=merge(tract_id,preFilt, by="row.names", all=T)
tract_id=tract_id[complete.cases(tract_id$In_metabolomics),]
row.names(tract_id)=tract_id$Row.names
tract_id$Row.names=NULL

#Select genus level

taxa_gn=taxa[unlist(lapply(rownames(taxa),function(x) length(unlist(strsplit(x, "|", fixed=TRUE)))==6)),]
#Remove other taxa levels from names
rownames(taxa_gn)=sapply(strsplit(rownames(taxa_gn), "|", fixed=TRUE), tail, 1)

#Clean species names

taxa_sp=taxa[grep("s__", rownames(taxa)), ]
taxa_sp=taxa_sp[!grepl("t__", rownames(taxa_sp)), ]
row.names(taxa_sp)=gsub(".*s__","",row.names(taxa_sp))

taxa=as.data.frame(t(taxa_sp))

my_shannon=as.data.frame(diversity(taxa, index="shannon"))
colnames(my_shannon)=c("Shannon_Index")
my_simpson=as.data.frame(diversity(taxa, index="simpson"))
colnames(my_simpson)=c("Simpson_Index")
richness_microbiota=merge(my_shannon,my_simpson,by="row.names")
rownames(richness_microbiota)=richness_microbiota$Row.names
richness_microbiota$Row.names=NULL

#Remove samples with low richeness
norich=rownames(richness_microbiota)[richness_microbiota$Shannon_Index==0] 
tract_id[rownames(tract_id)==norich,]$In_MGS="Failed"

richness_microbiota2=subset(richness_microbiota,rownames(richness_microbiota) %in% rownames(summary_sample))
richness_microbiota2=richness_microbiota2[richness_microbiota2$Shannon_Index>0,]
cc_pheno2=merge(richness_microbiota2,cc_pheno,by="row.names")
rownames(cc_pheno2)=cc_pheno2$Row.names
cc_pheno2$Row.names=NULL

print (paste("Number of samples excluded due failed metagenomics:", length(tract_id[tract_id$In_MGS!="yes",])))
print (paste("Number of samples left:", length(rownames(cc_pheno2))))


__Microbial richness indexes (Shannon & Simpson)__

cc_pheno2.1=cc_pheno2
cc_pheno2.1$ibd_Diagnosis=as.character(cc_pheno2.1$ibd_Diagnosis)
cc_pheno2.1$ibd_Diagnosis[cc_pheno2.1$ibd_Diagnosis=="Control"]="A_Control"
ggplot(cc_pheno2.1, aes(ibd_Diagnosis,Shannon_Index, fill=ibd_Diagnosis)) + theme_bw() + geom_violin(alpha=0.8) + geom_jitter(height = 0, width = 0.1, alpha=0.4) + geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2, fill="white")   + theme_bw()  + scale_fill_manual(values =c("grey80", "purple","pink2" ,"darkolivegreen3"))



4.Summary of SCFA 
===

__Show SCFA names__

unique(SCFA$Analyte)


__Load data and replace values under the detection levels for missing values and restructure the data__

Original report contain 1 row per sample per measurment, we want to restructure in a way that all samples are in rows and each column is a SCFA.  



#Duplicate original values (before remplacing for NA)

SCFA$ori_Result=SCFA$Result

#Replace values below level quantification (BLOQ) to NA
SCFA$Result[SCFA$Comment.1=="BLOQ"] = "N/Q"
SCFA$Result[SCFA$Result=="N/Q"] = NA
```

__Summary of the values under the detection levels (TRUE column)__


table(SCFA$Analyte,is.na(SCFA$Result))



#Remove unecessary columns (sample type, group/project name, detection limits and comments)

SCFAsub=SCFA[,c("Unique.Sample.ID...as.labeled.on.tube.","Sample.Amount..gram.","Comment","Analyte","Result")]

#Reshape the SCFA table
SCFA_table <- dcast(SCFAsub, ...~Analyte)

#Chance sample names to match metadata

rownames(SCFA_table)=SCFA_table$Unique.Sample.ID...as.labeled.on.tube.
SCFA_table$Unique.Sample.ID...as.labeled.on.tube.=NULL

SCFA_new_ID=merge(IDs,SCFA_table, by="row.names")
SCFA_new_ID$Row.names=NULL
row.names(SCFA_new_ID)=SCFA_new_ID$V5
SCFA_new_ID$V5=NULL
colnames(SCFA_new_ID)=c("Amount_sample_gram", "Box", "Methylbutyric_acid", "acetic_acid", "butyric_acid", "hexanoic_acid", "isobutyric_acid","isovaleric_acid" ,"propionic_acid", "valeric_acid")




