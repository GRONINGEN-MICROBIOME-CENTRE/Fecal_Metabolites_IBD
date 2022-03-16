Evaluating the effect of technical factors and host age, sex, BMI, Bowel movements per day
===


rownames(cc_pheno5)=cc_pheno5$Row.names
cc_pheno5$Row.names=NULL
cc_pheno6=cc_pheno5
table(cc_pheno6$run_day_cat, cc_pheno6$COMMENT)
cc_pheno6=subset(cc_pheno6,select = -c(Preparation_batch,Preparation_day,Preparation_order,CLIENT.IDENTIFIER,Box,Preparation_day_cat, RUN.DAY,SPL_AMNT_GRAM))

#Correction factors: day measurment was done, LC column, ammount sample per gram, month sample in a freezer, sex and age. 
my_correction_factors=cc_pheno6[,c("run_day_cat","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","host_Sex","host_Age_sampling", "host_BMI", "clinical_BowelMovementADayDef")]
my_correction_factors_test=merge(my_correction_factors,q_mtb,by="row.names")
rownames(my_correction_factors_test)=my_correction_factors_test$Row.names
my_correction_factors_test$Row.names=NULL


my_corr_factors_results=matrix(ncol = 5, nrow = 6832)
b=1
for ( a in 1:ncol(my_correction_factors)){
  for (x in (ncol(my_correction_factors)+1):ncol(my_correction_factors_test)){
    if(is.numeric(my_correction_factors[,a])){
      my_corr_factors_results[b,1]=colnames(my_correction_factors_test)[a]
      my_corr_factors_results[b,2]=colnames(my_correction_factors_test)[x]
      my_test=cor.test(my_correction_factors_test[,a],my_correction_factors_test[,x], method = "spearman")
      my_corr_factors_results[b,3]=my_test$estimate
      my_corr_factors_results[b,4]=my_test$p.value
      my_corr_factors_results[b,5]="spearman"
      b=b+1
    } else{
      my_corr_factors_results[b,1]=colnames(my_correction_factors_test)[a]
      my_corr_factors_results[b,2]=colnames(my_correction_factors_test)[x]
      my_test=kruskal.test(my_correction_factors_test[,x]~my_correction_factors_test[,a])
      my_corr_factors_results[b,3]=my_test$statistic
      my_corr_factors_results[b,4]=my_test$p.value
      my_corr_factors_results[b,5]="krustal"
      b=b+1
    }
  }
}

my_corr_factors_results=as.data.frame(my_corr_factors_results)
colnames(my_corr_factors_results)=c("Factor", "Metabolite","Estimate","p.value", "test")
my_corr_factors_results$p.value=as.numeric(as.character(my_corr_factors_results$p.value))
my_corr_factors_results$FDR=p.adjust(my_corr_factors_results$p.value,method = "fdr")
my_corr_factors_results$Bonferroni=p.adjust(my_corr_factors_results$p.value,method = "bonferroni")

x=as.data.frame(table(my_corr_factors_results$Factor[my_corr_factors_results$Bonferroni<0.05]))
y=as.data.frame(table(my_corr_factors_results$Factor[my_corr_factors_results$FDR<0.05]))
associations_factors=as.data.frame(cbind(x,y))
associations_factors[,3]=NULL
colnames(associations_factors)=c("Name","Bonferroni","FDR")
plot_asso=melt(associations_factors)


write.table(my_corr_factors_results, "~/Desktop/Metabolomics_v2/3.Preliminary_results/technical_confouders_quantitative.txt", sep = "\t", quote = F)