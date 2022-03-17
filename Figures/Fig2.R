library(tidyverse)
meta_phenos <- read.delim("~/Desktop/IBD_metabolomics_2022/6.Supplementary_tables/meta_phenos.txt")
 cnt_sig=subset(meta_phenos,meta_phenos$FDR_CNT<0.05)

cd_sig=subset(meta_phenos,meta_phenos$FDR_CD<0.05)
UC_sig=subset(meta_phenos,meta_phenos$FDR_UC<0.05)
meta_sig=subset(meta_phenos,meta_phenos$FDR_meta_random<0.05)

cnt_sig2=as.data.frame(table(cnt_sig$phenotype))
cnt_sig2$Cohort="1_Control"
cnt_sig2$Order=1
cd_sig2=as.data.frame(table(cd_sig$phenotype))
cd_sig2$Cohort="3_CD"
cd_sig2$Order=1
UC_sig2=as.data.frame(table(UC_sig$phenotype))
UC_sig2$Cohort="2_UC"
UC_sig2$Order=1
meta_sig2=as.data.frame(table(meta_sig$phenotype))
meta_sig2$Cohort="4_Meta-analysis"
meta_sig2$Order=meta_sig2$Freq

my_meta=bind_rows(cnt_sig2, UC_sig2,cd_sig2,meta_sig2)
my_meta=bind_rows(cnt_sig2, UC_sig2,cd_sig2,meta_sig2)


ggplot(my_meta,aes(reorder(phenotype, Order),Freq)) + geom_bar(stat = "identity") + coord_flip() + facet_wrap(variable~., nrow = 1) + theme_bw()


#coffee
ggplot(taxa_bugs_pheno, aes(scale(diet_coffee.Res), taxa_bugs_pheno$`5_acetylamino_6_amino_3_methyluracil`, color=ibd_Diagnosis))  + 
geom_point(alpha=0.9, size=3) + geom_smooth(method = "lm")+ 
scale_color_manual(values = c("black","mediumpurple","darkolivegreen3")) + theme_bw()   + facet_wrap(ibd_Diagnosis~.) 

#Statin

ggplot(taxa_bugs_pheno, aes(med_statin, taxa_bugs_pheno$`7_ketocholesterol`, fill=med_statin))  + geom_violin(alpha=0.8,) + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, )  + geom_smooth(method = "lm") + theme_bw()   + facet_wrap(ibd_Diagnosis~.) + scale_fill_manual(values = c("grey80","firebrick3"))


         Non_user User
  Control      240   15
  IBD_UC       156   18
  IBD_zCD      228   10


#Mesalazine
taxa_bugs_pheno_IBD=subset(taxa_bugs_pheno,taxa_bugs_pheno$ibd_Diagnosis!="Control")
ggplot(taxa_bugs_pheno_IBD, aes(as.factor(med_mesalazines),  gentisate, fill=as.factor(med_mesalazines))) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + facet_wrap(ibd_Diagnosis~.) + theme_bw() + scale_fill_manual(values = c("grey80","firebrick3")) + ylab("clr(gentisate)")
table(taxa_bugs_pheno_IBD$ibd_Diagnosis,taxa_bugs_pheno_IBD$med_mesalazines)


          Non_user User
  IBD_UC        48  126
  IBD_zCD      214   24
         