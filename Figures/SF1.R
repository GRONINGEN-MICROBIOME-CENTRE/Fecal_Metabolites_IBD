library(ggrepel)
library (ggplot2)


Suppl Fig 1A
---

for_volcano2=subset(my_cc_quantitative_resec, my_cc_quantitative_resec$phenotype!="IBD")
for_volcano2$SUPER.PATHWAY[is.na(for_volcano2$SUPER.PATHWAY)]="SCFA"
for_volcano2$SUB.PATHWAY[is.na(for_volcano2$SUB.PATHWAY)]="SCFA"
for_volcano2$my_label=NA
for_volcano2$my_label[-log10(for_volcano2$FDR)>4] <- for_volcano2$metabolite[ -log10(for_volcano2$FDR)>4]

ggplot(for_volcano2, aes(Estimate,-log10(FDR), fill=SUPER.PATHWAY, label=my_label)) + geom_point(shape=21, size=2) + theme_bw() + geom_hline(yintercept = -log(0.05), col="red") + xlim (-2.1,2.1) + facet_wrap(.~phenotype, ncol = 1) + scale_fill_manual(values  = c("firebrick3",  "steelblue2",  "gold3", "limegreen","grey77",  "pink1", "salmon", "turquoise", "black", "azure")) + geom_text_repel(size=2)


Fig 1E
---

ggplot(taxa_bugs_pheno, aes(ibd_Diagnosis,L_urobilin, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.4, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.2, color="black" ) + theme_bw() + scale_fill_manual(values =c("grey80" ,"darkolivegreen3","purple")) + xlab ("") + ylab("clr(L-urobilin)")
ggplot(taxa_bugs_pheno, aes(ibd_Diagnosis,diaminopimelate, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.4, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.2, color="black" ) + theme_bw() + scale_fill_manual(values =c("grey80" ,"darkolivegreen3","purple")) + xlab ("") + ylab("clr(diaminopimelate)")
ggplot(Try_ratios_p2, aes(ibd_Diagnosis,R2_TRY, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.4, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.2, color="black" ) + theme_bw() + scale_fill_manual(values =c("grey80" ,"darkolivegreen3","purple")) + xlab ("") + ylab("Ratio Tryptamine/Tryptophan")
ggplot(Try_ratios_p2, aes(ibd_Diagnosis,Ratio_PUFA, fill=ibd_Diagnosis)) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.4, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.2, color="black" ) + theme_bw() + scale_fill_manual(values =c("grey80" ,"darkolivegreen3","purple")) + xlab ("") + ylab("Ratio o6/o3 PUFA")

