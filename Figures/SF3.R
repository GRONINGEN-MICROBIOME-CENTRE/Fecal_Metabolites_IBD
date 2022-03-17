#Mesalazine
taxa_bugs_pheno_IBD=subset(taxa_bugs_pheno,taxa_bugs_pheno$ibd_Diagnosis!="Control")
ggplot(taxa_bugs_pheno_IBD, aes(as.factor(med_mesalazines),  gentisate, fill=as.factor(med_mesalazines))) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + facet_wrap(ibd_Diagnosis~.) + theme_bw() + scale_fill_manual(values = c("grey80","firebrick3")) + ylab("clr(gentisate)")
table(taxa_bugs_pheno_IBD$ibd_Diagnosis,taxa_bugs_pheno_IBD$med_mesalazines)


          Non_user User
  IBD_UC        48  126
  IBD_zCD      214   24


ggplot(ibd_test_phenos3, aes(as.factor(clinical_ResectionColonicAny),  pyridoxamine, fill=as.factor(clinical_ResectionColonicAny))) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + facet_wrap(clinical_Diagnosis~.) + theme_bw() + scale_fill_manual(values = c("grey80","firebrick3")) + ylab("clr(pyridoxamine)")

table(ibd_test_phenos3$clinical_ResectionColonicAny,ibd_test_phenos3$clinical_Diagnosis)
   
      1   2
  1 177 213
  2   9  25


my_phenoscd2=my_phenoscd[complete.cases(my_phenoscd$ibd_montreal_B),]
  ggplot(my_phenoscd2, aes(ibd_montreal_B,  butyric_acid, fill=as.factor(ibd_montreal_B))) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + theme_bw() + scale_fill_manual(values = c("red","salmon","firebrick3")) + ylab("scaled (log(butyrate))")

table(my_phenoscd$ibd_montreal_B)

 B1  B2  B3 
129  76  28 


ggplot (BA_p2, aes(Resection_ileum, R1_DCA_vs_CA, fill=Resection_ileum))+ geom_violin(alpha=0.8) + geom_jitter(height = 0, width = 0.1, alpha=0.4) + geom_boxplot(width=0.1, fill="white", outlier.alpha = 0) + theme_bw() + scale_fill_manual(values = c("gray80", "darkolivegreen3", "mediumpurple", "purple4"))