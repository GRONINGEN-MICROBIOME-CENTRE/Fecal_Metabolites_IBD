Biplot (mmvec output)
---

all_ord=read_qza("~/Desktop/IBD_metabolomics_2022/2.Output/7.mmvec/biplot.qza")
info <- read.delim("~/Desktop/IBD_metabolomics_2022/2.Output/CC_corrected_for_resection_quantitative.txt")
info$SIG="No"
info$SIG[info$FDR<0.05 & info$Estimate>0]="Yes"
#info$SIG[info$FDR<0.05 & info$Estimate<0]="YesDown"
info_fig=subset(info, info$phenotype=="IBD")


colnames(info_fig)[1]="SampleID"


info_fig$SampleID=gsub("glutamate__gamma_methyl_ester","glutamate._gamma_methyl_ester", info_fig$SampleID)
info_fig$SampleID=gsub("N_N_dimethyl_5_aminovalerate","N.N.dimethyl_5_aminovalerate", info_fig$SampleID)
info_fig$SampleID=gsub("N_N_dimethylalanine","N.N_dimethylalanine", info_fig$SampleID)
info_fig$SampleID=gsub("N_N_N_trimethyl_5_aminovalerate","N.N.N_trimethyl_5_aminovalerate", info_fig$SampleID)
info_fig$SampleID=gsub("N_N_N_trimethyl_alanylproline_betaine_TMAP","N.N.N_trimethyl_alanylproline_betaine_TMAP", info_fig$SampleID)
info_fig$SampleID=gsub("N1_N12_diacetylspermine","N1.N12_diacetylspermine", info_fig$SampleID)
info_fig$SampleID=gsub("N.N.dimethyl_5_aminovalerate","N.N_dimethyl_5_aminovalerate", info_fig$SampleID)
info_fig$SampleID=gsub("N2_N5_diacetylornithine","N2.N5_diacetylornithine", info_fig$SampleID)
info_fig$SampleID=gsub("N6.N6_dimethyllysine","N6_N6_dimethyllysine", info_fig$SampleID)
info_fig$SampleID=gsub("N6_N6_N6_trimethyllysine","N6.N6.N6_trimethyllysine", info_fig$SampleID)
info_fig$SampleID=gsub("uridine_2_3_cyclic_monophosphate","uridine_2.3_cyclic_monophosphate", info_fig$SampleID)
info_fig$SampleID=gsub("N6_N6_dimethyllysine","N6.N6_dimethyllysine", info_fig$SampleID)
info_fig$SampleID=gsub("X1_3_7_trimethylurate","X1.3.7_trimethylurate", info_fig$SampleID)
info_fig$SampleID=gsub("X1_3_dimethylurate","X1.3_dimethylurate", info_fig$SampleID)
info_fig$SampleID=gsub("X1_3_propanediol","X1.3_propanediol", info_fig$SampleID)
info_fig$SampleID=gsub("X1_7_dimethylurate","X1.7_dimethylurate", info_fig$SampleID)
info_fig$SampleID=gsub("X12_13_DiHOME","X12.13_DiHOME", info_fig$SampleID)
info_fig$SampleID=gsub("X2_3_dihydroxyisovalerate","X2.3_dihydroxyisovalerate", info_fig$SampleID)
info_fig$SampleID=gsub("X2_3_dimethylsuccinate","X2.3_dimethylsuccinate", info_fig$SampleID)
info_fig$SampleID=gsub("X2_4_6_trihydroxybenzoate","X2.4.6_trihydroxybenzoate", info_fig$SampleID)
info_fig$SampleID=gsub("X2_4_dihydroxybutyrate","X2.4_dihydroxybutyrate", info_fig$SampleID)
info_fig$SampleID=gsub("X2_8_quinolinediol","X2.8_quinolinediol", info_fig$SampleID)
info_fig$SampleID=gsub("X2R_3R_dihydroxybutyrate","X2R.3R_dihydroxybutyrate", info_fig$SampleID)
info_fig$SampleID=gsub("X2_6_dihydroxybenzoic_acid","X2.6_dihydroxybenzoic_acid", info_fig$SampleID)
info_fig$SampleID=gsub("X2S_3R_dihydroxybutyrate","X2S.3R_dihydroxybutyrate", info_fig$SampleID)
info_fig$SampleID=gsub("X3_7_dimethylurate","X3.7_dimethylurate", info_fig$SampleID)
info_fig$SampleID=gsub("X5_6_dihydrouridine","X5.6_dihydrouridine", info_fig$SampleID)
info_fig$SampleID=gsub("X9_10_DiHOME","X9.10_DiHOME", info_fig$SampleID)


baseplot_cnt = ggplot() + 
     geom_point(
         data=all_ord$data$Vectors%>%
             left_join(info_fig),
         aes(x=PC1, y=PC2, fill=SIG), alpha=1, 
         shape=21, size = 3
     ) + scale_fill_manual(values=c("white", "orange1"))

 
 baseplot_cnt +
     theme_bw() +
     xlab(paste(round(100*all_ord$data$ProportionExplained[1],2),"%")) +
     ylab(paste(round(100*all_ord$data$ProportionExplained[2],2),"%")) +
     ggtitle("Biplot mmvec") +
     geom_segment(size=2, alpha=0.6, data=all_ord$data$Species %>% 
                      mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                      top_n(10, a) %>% #keep 10 furthest away points
                      mutate(PC1=PC1*0.5, PC2=PC2*0.5), 
                  aes(x=0, xend=PC1, y=0, yend=PC2,color=FeatureID),
                  arrow = arrow(length = unit(0.3,"cm"))
     ) + scale_color_manual(values =c ("black","red3" ,"red","steelblue2","blue3","magenta" ,"gray40" ,"gray80","violetred" ,"purple"))


ggplot(taxa_bugs_pheno, aes(Ruminococcus_gnavus, tryptamine, color=ibd_Diagnosis))  + geom_point(alpha=0.9, size=3) + geom_smooth(method='lm')+ scale_color_manual(values = c("black","mediumpurple","darkolivegreen3")) + theme_bw()   + facet_wrap(ibd_Diagnosis~.) 


ggplot(taxa_bugs_pheno, aes(as.factor(Bilophila_wadsworthia),  taurine, fill=as.factor(Bilophila_wadsworthia))) + geom_violin(alpha=0.8, color="black") + geom_jitter(position = position_dodge(0.9),alpha=0.2, color="black") + geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black" ) + facet_wrap(ibd_Diagnosis~.) + theme_bw() + scale_fill_manual(values = c("grey80","firebrick3"))


 ggplot(mgc_for_plot, aes(mgc_for_plot$Entryname.EUT_pathway, mgc_for_plot$oleoyl_ethanolamide, color=ibd_CD_and_UC))  +     geom_point(alpha=0.9, size=3) + geom_smooth(method='lm')+  scale_color_manual(values = c("black","darkolivegreen3","mediumpurple")) + theme_bw()   + facet_wrap(ibd_CD_and_UC~., scales = "free_x") 


ggplot(mgc_for_plot, aes(mgc_for_plot$Entryname.bai_operon, mgc_for_plot$cholate, color=ibd_CD_and_UC))  +    geom_point(alpha=0.9, size=3) + geom_smooth(method='lm')+    scale_color_manual(values = c("black","darkolivegreen3","mediumpurple")) + theme_bw()   + facet_wrap(ibd_CD_and_UC~., scales = "free") 


ggplot(pwy_for_plot, aes(pwy_for_plot$HISDEG.PWY..L.histidine.degradation.I, histidine, color=ibd_CD_and_UC))  +  geom_point(alpha=0.9, size=3) + geom_smooth(method='lm') +  scale_color_manual(values = c("black","darkolivegreen3","mediumpurple")) + theme_bw()   + facet_wrap(ibd_CD_and_UC~., scales = "free_x")

















