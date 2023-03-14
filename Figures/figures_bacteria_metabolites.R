
#GNAVUS
ggplot(models_for_interac_lin, aes(Ruminococcus_gnavus, tryptamine, color=cohort))  + 
  geom_point(alpha=0.5, size=3) + 
  geom_smooth(method='lm' ,color="black")+ 
  scale_color_manual(values = c("black","purple","darkolivegreen3")) + theme_bw() + 
  xlab ("clr(Ruminococcus gnavus)") + ylab ("clr(tryptamine)")


#HISTIDINE
ggplot(models_for_interac_pwy, aes(models_for_interac_pwy$HISDEG.PWY..L.histidine.degradation.I, histidine, color=cohort))  + 
  geom_point(alpha=0.5, size=3) + geom_smooth(method='lm' ,color="black")+ 
  scale_color_manual(values = c("black","purple","darkolivegreen3")) + theme_bw() + 
  xlab ("clr(Histidine degradation pathway (HISDEG-PWY) 470 / 199") + ylab ("clr(histidine)")

#Bilophila

ggplot(models_for_interac, aes(as.factor(Bilophila_wadsworthia), taurine, color=cohort))  + 
  geom_jitter(size=3, alpha=0.5, position = position_jitter(width = 0.3, height = 0.3)) +
  geom_boxplot(width=0.1, position = position_dodge(0.9),alpha=0.8, color="black", outlier.alpha = 0)+ 
  scale_fill_manual(values = c("grey80","firebrick3")) + scale_color_manual(values = c("black","purple","darkolivegreen3")) + theme_bw()

#Bai gene 

ggplot(models_for_interac_mgc, aes(Entryname.bai_operon, lithocholate, color=cohort))  + 
  geom_point(alpha=0.5, size=3) + 
  geom_smooth(method='lm' ,color="black")+ 
  scale_color_manual(values = c("black","purple","darkolivegreen3")) + theme_bw() + 
  xlab ("clr(Bai operon)") + ylab ("clr(lithocholate)")

ggplot(models_for_interac_mgc, aes(Entryname.bai_operon, lithocholate, color=cohort, fill=dysbiotic))  + 
  geom_point(size=3, alpha=0.5) + geom_smooth(method='lm' ,color="black")+ 
  scale_fill_manual(values = c("grey80","firebrick3")) + theme_bw() + 
  scale_color_manual(values = c("black","purple","darkolivegreen3")) +
  facet_wrap(~dysbiotic, labeller = labeller( dysbiotic=c("no"="eubiosis", "yes"="dysbiosis"))) + 
  xlab ("clr(Bai operon)") + ylab ("clr(lithocholate)")


ggplot(models_for_interac_mgc, aes(EUT_pathway, palmitoyl_ethanolamide, color=cohort))  + 
  geom_point(alpha=0.5, size=3) + 
  geom_smooth(method='lm' ,color="black")+ 
  scale_color_manual(values = c("black","purple","darkolivegreen3")) + theme_bw() + 
  xlab ("clr(EUT operon)") + ylab ("clr(palmitoyl-ethanolamide)")
