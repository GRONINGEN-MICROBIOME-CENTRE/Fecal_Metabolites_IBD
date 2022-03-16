Process, filter, normalize metabolic gene clusters
====


setwd("~/Desktop/IBD_Metabolomics_2022/1.Input/BGC_new/")

file_list <- list.files()
flag=1
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (flag==1){
    bgc <- read.table(file, header=TRUE,  check.names = F, sep=",", row.names = 1, stringsAsFactors = F)
    flag=5
  }
  
  # if the merged dataset does exist, append to it
  else {
    temp_bgc <-read.table(file, header=TRUE,  check.names = F, sep=",", row.names = 1, stringsAsFactors = F)
    bgc<-merge(bgc, temp_bgc, by="row.names", all.x=T)
    row.names(bgc)=bgc$Row.names
    bgc$Row.names=NULL
    rm(temp_bgc)
  }
  
}

cc_rd=read.table("~/Desktop/IBD_Metabolomics/1.Input/transform_ID.txt", header = T, row.names = "IDs_MGS")


#Get coverage of the core regions
bgc_coverage=select(bgc, contains(".corecov"))
colnames(bgc_coverage)=gsub(".corecov","",colnames(bgc_coverage))
bgc_coverage2=merge(cc_rd,t(bgc_coverage), by="row.names")
rownames(bgc_coverage2)=bgc_coverage2$PID2
bgc_coverage2$PID2=NULL
bgc_coverage2$PID=NULL
bgc_coverage2$Row.names=NULL
bgc_coverage2$IDs_MGS=NULL
bgc_coverage2$PF_RD=NULL
bgc_coverage3=subset(bgc_coverage2, row.names(bgc_coverage2) %in% row.names(cc_pheno3))


#Set a filtering

bgc_coverage4=bgc_coverage3
bgc_mean_coverage=data.frame(BCG=colnames(bgc_coverage4), Mean=colMeans(bgc_coverage4))
ggplot(bgc_mean_coverage, aes(reorder(BCG,-Mean), Mean )) + geom_bar(stat = "identity", fill="red4") + theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab ("Mean coverage")  + xlab ("Biosynthetic genes clusters")

select_coverage=bgc_mean_coverage$BCG[bgc_mean_coverage$Mean>0.05]


#Collapse RPKM to pathways
library(stringr)

bgc_rpkm=select(bgc, contains(".RPKM"))
colnames(bgc_rpkm)=gsub(".RPKM","",colnames(bgc_rpkm))


bgc_rpkm2=merge(cc_rd,t(bgc_rpkm), by="row.names")
rownames(bgc_rpkm2)=bgc_rpkm2$PID2
bgc_rpkm2$PID2=NULL
bgc_rpkm2$PID=NULL
bgc_rpkm2$Row.names=NULL
bgc_rpkm2$IDs_MGS=NULL
bgc_rpkm2$PF_RD=NULL
bgc_rpkm3=subset(bgc_rpkm2, row.names(bgc_rpkm2) %in% row.names(cc_pheno3))

#Subset MGC based on coverage

bgc_rpkm4=bgc_rpkm3[,colnames(bgc_rpkm3)%in%select_coverage]
bgc_rpkm4=as.data.frame(t(bgc_rpkm4))
bgc_rpkm4$type=str_split_fixed(row.names(bgc_rpkm4),"--", 3)[,2]
bgc_rpkm4= bgc_rpkm4 %>% group_by(type) %>% summarise_all(sum)
bgc_rpkm4=as.data.frame(bgc_rpkm4)
row.names(bgc_rpkm4)=make.names(bgc_rpkm4$type)
bgc_rpkm4$type=NULL
bgc_rpkm4=as.data.frame(t(bgc_rpkm4))


#Transform 
mgc_fil_trans=transform_and_filter_taxa(bgc_rpkm4, samples_row = T, method = "clr",missing_filter = 20)