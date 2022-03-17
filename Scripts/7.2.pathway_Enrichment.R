library(ggplot2)
library(factoextra)
library(vegan)
library(ggsci)
library(crayon)

my_cc_categories=my_cc_quantitative
my_cc_categories$SUPER.PATHWAY[is.na(my_cc_categories$SUPER.PATHWAY)]="SCFA"
my_cc_categories$SUB.PATHWAY[is.na(my_cc_categories$SUB.PATHWAY)]="SCFA"

my_cc_categories$level=NA
my_cc_categories$level[my_cc_categories$Estimate>0]="UpRegulate"
my_cc_categories$level[my_cc_categories$Estimate<0]="DownRegulate"
my_cc_categories$level[my_cc_categories$FDR>0.05]="Notsignificant"

my_cc_categories_IBD=my_cc_categories[my_cc_categories$phenotype=="IBD",]
my_cc_categories_CD=my_cc_categories[my_cc_categories$phenotype=="CD",]
my_cc_categories_UC=my_cc_categories[my_cc_categories$phenotype=="UC",]

sub_pathway_count=data.frame(table(my_cc_categories_CD$SUB.PATHWAY))

# extract those with number >5 (the minimun sample size for test)
#sub_pathway_count=sub_pathway_count[sub_pathway_count$Freq>=5,]


enrich_two_compare_2 = foreach(i=1:nrow(sub_pathway_count),.combine = rbind) %do%  { # compare t statistics between CD and UC
  
  tmp.name=as.character(sub_pathway_count$Var1[i])
  
  tmp.order.CD=my_cc_categories_CD[order(my_cc_categories_CD$t.value,decreasing = T),] # rank
  tmp.order.CD=tmp.order.CD[tmp.order.CD$SUB.PATHWAY==tmp.name,]
  
  tmp.order.UC=my_cc_categories_UC[order(my_cc_categories_UC$t.value,decreasing = T),] # rank
  tmp.order.UC=tmp.order.UC[tmp.order.UC$SUB.PATHWAY==tmp.name,]
  
  #tmp.order.CD=tmp.order.CD[order(tmp.order.CD$metabolite),]
  #tmp.order.UC=tmp.order.UC[order(tmp.order.UC$metabolite),]
  
  test=(wilcox.test(tmp.order.CD$t.value, tmp.order.UC$t.value, paired = TRUE))
  
  return.string=data.frame(sub.pathway=tmp.name,total_metabolites=nrow(tmp.order.CD), UpinCD=nrow(tmp.order.CD[tmp.order.CD$level=="UpRegulate",]),DowninCD=nrow(tmp.order.CD[tmp.order.CD$level=="DownRegulate",]),UpinUC=nrow(tmp.order.UC[tmp.order.UC$level=="UpRegulate",]), DowninUC=nrow(tmp.order.UC[tmp.order.UC$level=="DownRegulate",]), Mean_CD=mean(tmp.order.CD$t.value),Mean_UC=mean(tmp.order.UC$t.value),
                           diff=mean(tmp.order.CD$t.value)-mean(tmp.order.UC$t.value),statistics=test$statistic,Pvalue=test$p.value)
  
}
enrich_two_compare_2$FDR=p.adjust(enrich_two_compare$Pvalue,method="BH")


enrich_CD = foreach(i=1:nrow(sub_pathway_count),.combine = rbind) %do%  {
  
  tmp.name=as.character(sub_pathway_count$Var1[i])
  
  tmp.order.all=my_cc_categories_CD[order(my_cc_categories_CD$t.value,decreasing = T),] # rank
  tmp.order.all$rank=1:nrow(tmp.order.all)
  tmp.order.all=tmp.order.all[tmp.order.all$SUB.PATHWAY==tmp.name,]
  
  tmp.order.part=my_cc_categories_CD[my_cc_categories_CD$SUB.PATHWAY==tmp.name,]
  tmp.order.part=tmp.order.part[order(tmp.order.part$t.value,decreasing = T),] # rank
  tmp.order.part$rank=1:nrow(tmp.order.part)
  
  tmp.order.all=tmp.order.all[order(tmp.order.all$metabolite),]
  tmp.order.part=tmp.order.part[order(tmp.order.part$metabolite),]
  
  n_up=as.numeric(length(tmp.order.part$metabolite[tmp.order.part$FDR<0.05 & tmp.order.part$level=="UpRegulate"]))
  n_down=as.numeric(length(tmp.order.part$metabolite[tmp.order.part$FDR<0.05 & tmp.order.part$level=="DownRegulate"]))
  
  test=(wilcox.test(tmp.order.part$rank, tmp.order.all$rank, paired = TRUE))
  
  return.string=data.frame(sub.pathway=tmp.name,Median_t=median(tmp.order.part$t.value),N_up=n_up,N_down=n_down,Pvalue=test$p.value)
  
}
enrich_CD$FDR=p.adjust(enrich_CD$Pvalue,method = "BH")

# UC vs. control
enrich_UC = foreach(i=1:nrow(sub_pathway_count),.combine = rbind) %do%  {
  
  tmp.name=as.character(sub_pathway_count$Var1[i])
  
  tmp.order.all=my_cc_categories_UC[order(my_cc_categories_UC$t.value,decreasing = T),] # rank
  tmp.order.all$rank=1:nrow(tmp.order.all)
  tmp.order.all=tmp.order.all[tmp.order.all$SUB.PATHWAY==tmp.name,]
  
  tmp.order.part=my_cc_categories_UC[my_cc_categories_UC$SUB.PATHWAY==tmp.name,]
  tmp.order.part=tmp.order.part[order(tmp.order.part$t.value,decreasing = T),] # rank
  tmp.order.part$rank=1:nrow(tmp.order.part)
  
  tmp.order.all=tmp.order.all[order(tmp.order.all$metabolite),]
  tmp.order.part=tmp.order.part[order(tmp.order.part$metabolite),]
  
  n_up=as.numeric(length(tmp.order.part$metabolite[tmp.order.part$FDR<0.05 & tmp.order.part$level=="UpRegulate"]))
  n_down=as.numeric(length(tmp.order.part$metabolite[tmp.order.part$FDR<0.05 & tmp.order.part$level=="DownRegulate"]))
  
  test=(wilcox.test(tmp.order.part$rank, tmp.order.all$rank, paired = TRUE))
  
  return.string=data.frame(sub.pathway=tmp.name,Median_t=median(tmp.order.part$t.value),N_up=n_up,N_down=n_down,Pvalue=test$p.value)
  
}
enrich_UC$FDR=p.adjust(enrich_UC$Pvalue,method = "BH")
