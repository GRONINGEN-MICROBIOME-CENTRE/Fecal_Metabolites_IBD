## Metabolite prediction

lasso <- read.delim("~/Desktop/IBD_metabolomics_2022/5.Figures/Fig5/lasso.txt")
mlasso=melt(lasso)

ggplot(mlasso, aes(Metabolite,variable, fill=variable, alpha=value)) + geom_tile( colour="white") + theme_bw() + scale_fill_manual(values = c("grey80","yellow","lightblue","cyan","blue3","salmon","red3","purple")) 


## IBD prediction

my_pred_ratios=case_control_quantitative[,c("IBD","ibd_Diagnosis")]
my_pred_ratios_mtb=all_new_ID_raw[,c("lactosyl_N_palmitoyl_sphingosine_d18_1_16_0","L_urobilin")]

my_pred_ratios=merge(my_pred_ratios,my_pred_ratios_mtb, by="row.names")

my_pred_ratios$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0[is.na(my_pred_ratios$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0)]=min(my_pred_ratios$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0, na.rm=T)/2


my_pred_ratios$L_urobilin[is.na(my_pred_ratios$L_urobilin)]=min(my_pred_ratios$L_urobilin, na.rm=T)/2


my_pred_ratios_corr=case_control_quantitative[,c("IBD","ibd_Diagnosis","lactosyl_N_palmitoyl_sphingosine_d18_1_16_0","L_urobilin")]
my_pred_ratios_corr=subset(my_pred_ratios_corr,my_pred_ratios_corr$ibd_Diagnosis!="IBDU")
my_pred_ratios_corr$ibd_Diagnosis=gsub("CD","xCD",my_pred_ratios_corr$ibd_Diagnosis)

ggplot(my_pred_ratios, aes(IBD,Ratio, fill=IBD)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + theme_bw() + scale_fill_manual(values = c("grey80","orange2")) 



ggplot(my_pred_ratios_corr, aes(ibd_Diagnosis,lactosyl_N_palmitoyl_sphingosine_d18_1_16_0, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + theme_bw()  + ylab("clr(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + xlab("") + scale_fill_manual(values = c("grey80","darkolivegreen3","mediumpurple"))


ggplot(my_pred_ratios_corr, aes(ibd_Diagnosis,L_urobilin, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + theme_bw()  + ylab("clr(L-urobilin)") + xlab("") + scale_fill_manual(values = c("grey80","darkolivegreen3","mediumpurple"))

lasso <- read.delim("~/Desktop/IBD_metabolomics_2022/5.Figures/Fig5/lasso.txt")
mlasso=melt(lasso)
ggplot(mlasso, aes(variable,value, fill=variable)) + geom_jitter(alpha=0.9, shape=21, colour="black", size=3) + geom_boxplot(width=0.8, position = position_dodge(0.5), colour="black", alpha=0.7, outlier.alpha = 0 ) + theme_bw() + scale_fill_manual(values = c("grey80","yellow","lightblue","cyan","blue3","salmon","red3","purple")) +theme(legend.position = "none")