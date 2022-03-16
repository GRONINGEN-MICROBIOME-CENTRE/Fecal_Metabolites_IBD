my_cc_categories=my_cc_quantitative
my_cc_categories$SUPER.PATHWAY[is.na(my_cc_categories$SUPER.PATHWAY)]="SCFA"
my_cc_categories$SUB.PATHWAY[is.na(my_cc_categories$SUB.PATHWAY)]="SCFA"

my_cc_categories_IBD=my_cc_categories[my_cc_categories$phenotype=="IBD",]
my_cc_categories_CD=my_cc_categories[my_cc_categories$phenotype=="CD",]
my_cc_categories_UC=my_cc_categories[my_cc_categories$phenotype=="UC",]

CD_sig=subset(my_cc_categories_CD,my_cc_categories_CD$FDR<0.05)
UC_sig=subset(my_cc_categories_UC,my_cc_categories_UC$FDR<0.05)

CD_sig_up=subset(CD_sig,CD_sig$Estimate>0)
UC_sig_up=subset(UC_sig,UC_sig$Estimate>0)

my_cc_categories_resec=my_cc_quantitative_resec
my_cc_categories_resec$SUPER.PATHWAY[is.na(my_cc_categories_resec$SUPER.PATHWAY)]="SCFA"
my_cc_categories_resec$SUB.PATHWAY[is.na(my_cc_categories_resec$SUB.PATHWAY)]="SCFA"

my_cc_categories_IBD_r=my_cc_categories_resec[my_cc_categories_resec$phenotype=="IBD",]
my_cc_categories_CD_r=my_cc_categories_resec[my_cc_categories_resec$phenotype=="CD",]
my_cc_categories_UC_r=my_cc_categories_resec[my_cc_categories_resec$phenotype=="UC",]

CD_sig_r=subset(my_cc_categories_CD_r,my_cc_categories_CD_r$FDR<0.05)
UC_sig_r=subset(my_cc_categories_UC_r,my_cc_categories_UC_r$FDR<0.05)

CD_sig_up_r=subset(CD_sig_r,CD_sig_r$Estimate>0)
UC_sig_up_r=subset(UC_sig_r,UC_sig_r$Estimate>0)

CD_sig2=subset(CD_sig,CD_sig$metabolite %in% CD_sig_r$metabolite)
CD_sig_r2=subset(CD_sig_r,CD_sig_r$metabolite %in% CD_sig$metabolite)

table(sign(CD_sig2$Estimate) == sign(CD_sig_r2$Estimate))


UC_sig2=subset(UC_sig,UC_sig$metabolite %in% UC_sig_r$metabolite)
UC_sig_r2=subset(UC_sig_r,UC_sig_r$metabolite %in% UC_sig$metabolite)

table(sign(UC_sig2$Estimate) == sign(UC_sig_r2$Estimate))
