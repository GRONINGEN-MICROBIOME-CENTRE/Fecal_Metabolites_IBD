1. Import data
===

a) Metabolites: Here we use both raw values provided by Metabolon (metabolites AUC).


b) Short Chain Fatty Acids (SFCA): Together with the untargetted metabolite assessment, SCFA concentrations were measured in the same fecal samples (μg/g)


c) Host phenotypes: Information about the host (age,sex,BMI, diet, etc.). This consist of 3 files: phenotypes shared and used for the case-control analysis, phenotypes specific for the population controls & phenotypes specific for cases (disease location, activity, etc.)


d) Experiment variables: Different technical variables from the untargetted metabolic experiments. We will use this information as a potential counfounding factors. 


e) Metagenomics taxa: MetaPhlans profiles (v2.9)


f) Pathway annotation provided by Metabolon


g) IDs for each data modality: Allows to combine different datasets / information 




# a) Raw values (AUC of peaks)
all_metabolites_raw <- read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/all_metabolites_raw2.txt", row.names=1)

# b) SCFA

SCFA <- read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/SCFA.txt")

# c) Phenotypes

phenos_ibd=read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/phenos_IBD_clean_v2.txt",row.names = 1)
phenos_lld=read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/phenos_LLD_clean_v2.txt", row.names=1)
cc_pheno=read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/phenos_case_control_v1.txt", row.names=1)

cc_stoma=cc_pheno[,c("ibd_CurrentStomaOrPouch"),drop=F]
montreal= read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/phenos_montreal.txt", row.names=1)

# d) Experiment variables 

my_batch=read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/Metabolon_batches_v2.txt", row.names=1)
my_batch_mtb=read.delim("~/Desktop/IBD_Metabolomics_2022/13.Try_new_norm/1.Input/mtb_plat.txt", row.names=1)

# e) Metagenomic taxa (add metagenomic genes/pathways)

#MetaPhlAn v3
my_taxa = read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/Taxa_mpa3.txt", row.names=1, check.names = F)

# f) Metabolites annotation
annot <- read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/Metabolon_annotation_v2.txt", row.names=1)
annot$SUPER.PATHWAY.1=NULL


# g) ID's index (translate from different data modalities metabolomics <-> metagenomic <-> phenotypes)

IDs_LLD <- read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("~/Desktop/IBD_Metabolomics_2022/1.Input/IDs_IBD.txt", header=FALSE)


Check phenotypes files (count NAs) and inpute missing values
---


na_count = sapply(cc_pheno, function(y) sum(length(which(is.na(y)))))
na_count=data.frame(na_count)
cc_pheno$host_BMI[is.na(cc_pheno$host_BMI)]=median(cc_pheno$host_BM, na.rm = T)
cc_pheno$med_ACE_inhibitor[is.na(cc_pheno$med_ACE_inhibitor)]="Non_user"
cc_pheno$med_PPI[is.na(cc_pheno$med_PPI)]="Non_user"

diet=colnames(cc_pheno)[grep("diet",colnames(cc_pheno))]

for (a in diet){
  cc_pheno[,a][is.na(cc_pheno[,a])]=median(cc_pheno[,a], na.rm = T)
}

cc_pheno$host_smk_current[is.na(cc_pheno$host_smk_current)]="no"

na_count = sapply(cc_pheno, function(y) sum(length(which(is.na(y)))))
na_count=data.frame(na_count)



Adjust all ids
---

Translate all ids to the same id (easier to deal with multiple tables)



#Change Metabolon ids to UMCG ids to connect later to phenotypes 
all_raw=as.data.frame(t(all_metabolites_raw))

# Select IDs
IDs_IBD=IDs_IBD[,c(1,5)]
IDs_LLD$V2=NULL

# [IMPORTANT] Table with matching ids between UMCG and Metabolon
IDs=data.frame(rbind(as.matrix(IDs_IBD), as.matrix(IDs_LLD)))
rownames(IDs)=IDs$V1
IDs$V1=NULL

# Merge to replace the ids Metabolon => UMCG

all_new_ID_raw=merge(IDs,all_raw, by="row.names")
rownames(all_new_ID_raw)=all_new_ID_raw$V5
all_new_ID_raw$Row.names=NULL
all_new_ID_raw$V5=NULL

#Repeat for the batch information: Metabolon ID => UMCG ID

my_batch_id=merge(IDs,my_batch, by="row.names")
rownames(my_batch_id)=my_batch_id$V5
my_batch_id$Row.names=NULL
my_batch_id$V5=NULL


tract_id=data.frame(row.names(all_new_ID_raw))
rownames(tract_id)=tract_id$row.names.all_new_ID_raw.
tract_id$In_metabolomics="Yes"
tract_id$row.names.all_new_ID_raw.=NULL
tract_id=merge(tract_id,cc_stoma, by="row.names")
rownames(tract_id)=tract_id$Row.names
tract_id$Row.names=NULL
