library(ggrepel)
library (ggplot2)

Fig.1A
---
ggplot(cmp_log2, aes(-PC1,PC2, fill=ibd_Diagnosis)) + geom_point(size=3, pch=21) + theme_bw() + scale_fill_manual(values =c("purple","black","pink2" ,"darkolivegreen3"))
Fig 1B
---
ggplot(cmp_log3, aes(-PC1,PC2, fill=cholate)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

Fig 1C
---
ggplot(cmp_log3, aes(-PC1,PC2, fill=phenylalanylalanine)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")

Fig 1D
---

ggplot(cmp_log3, aes(-PC1,PC2, fill=suberate_C8_DC)) + geom_point(size=3,pch=21) + theme_bw() + scale_fill_viridis(direction = -1, option = "magma")


Fig 1E
---
cc_asso=corrected_for_resection_quantitative
colnames(cc_asso)=make.names(colnames(cc_asso))
cc_asso$SUPER.PATHWAY[is.na(cc_asso$SUPER.PATHWAY)]="SCFA"
cc_asso$SUB.PATHWAY[is.na(cc_asso$SUB.PATHWAY)]="SCFA"

cd_asso=subset(cc_asso,cc_asso$phenotype=="CD")
cd_num=as.data.frame(table(cd_asso$SUB.PATHWAY))

cd_median=aggregate(t.value ~ SUB.PATHWAY, cd_asso, median)

cd_path_median=merge(cd_num, cd_median, by.x="Var1", by.y="SUB.PATHWAY")
cc_asso2=subset(cc_asso, cc_asso$phenotype!="IBD")
my_plot=merge(cc_asso2,cd_path_median, by.x="SUB.PATHWAY", by.y = "Var1", all.x = T)

my_plot$names=paste(my_plot$SUB.PATHWAY, my_plot$Freq, sep = " (n=")
my_plot$names=paste(my_plot$names, ")", sep = "")
my_plot_lite=subset(my_plot, my_plot$Freq>5)

#Vertical
ggplot(my_plot_lite, aes(reorder(names, t.value.y),t.value.x, fill=phenotype)) + geom_boxplot( outlier.size = 0.1) +   geom_hline(aes(yintercept = 0))+ coord_flip() + facet_grid(.~phenotype) + theme_bw() + scale_fill_manual(values = c("purple", "darkolivegreen3")) + ylim(-10,10)
#Horitzontal
ggplot(my_plot_lite, aes(reorder(names, -t.value.y),t.value.x, fill=phenotype)) + geom_boxplot( outlier.size = 0.1) +   geom_hline(aes(yintercept = 0))+ facet_wrap(.~phenotype,ncol = 1) + theme_bw() + scale_fill_manual(values = c("purple", "darkolivegreen3")) + theme(axis.text.x = element_text(angle = 75, hjust = 1, size=6)) + xlab("") + ylab (" t-statistic")


