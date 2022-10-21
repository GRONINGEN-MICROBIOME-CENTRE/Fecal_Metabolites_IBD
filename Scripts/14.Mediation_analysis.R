#Mediation analysis

#This code is still a bit messy, I will clean it an implemented as a function soon (sorry!)

#For commodity I import 
#load("~/Desktop/IBD_metabolomics_2022/3.Workspace/Lasso.RData")

#install.packages("regmed")
library(regmed)
library(tidyverse)
library(readxl)
library(mediation)
#Test controls diet-microbiome-metabolite mediation
controls_data=subset(d3,d3$IBD==1)


#Import results phenotype-metabolite regression
phenos_cnt=read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_LLD_clean_v2.txt", row.names = 1)
phenos_cnt=subset(phenos_cnt, row.names(phenos_cnt)%in%row.names(controls_data))
phenos_cnt$clinical_Strict_IBS[phenos_cnt$clinical_Strict_IBS=="no"]=0
phenos_cnt$clinical_Strict_IBS[phenos_cnt$clinical_Strict_IBS=="yes"]=1
s12 <- read_excel("~/Desktop/IBD_metabolomics_2022/7.Manuscript_v3_coauthors/Submission_Gut/Suppl_Tables_1_22.xlsx", sheet = "S12")
colnames(s12)=s12[5,]
s12=s12[-c(1:5),]

s12$FDR_CNT=as.numeric(as.character(s12$FDR_CNT))
s12$FDR_CD=as.numeric(as.character(s12$FDR_CD))
s12$FDR_UC=as.numeric(as.character(s12$FDR_UC))

#Subset significant association between phenotypes and metabolites 
cnt_phenos_sig=subset(s12,s12$FDR_CNT<0.05)
cd_phenos_sig=subset(s12,s12$FDR_CD<0.05)
uc_phenos_sig=subset(s12,s12$FDR_UC<0.05)
cnt_phenos=unique(s12$phenotype[s12$FDR_CNT<0.05])
cnt_metabolites=unique(s12$metabolite[s12$FDR_CNT<0.05])



#################################################

## In Controls - 1 exposure, multiple mediator, 1 outcome
#################################################



#Select exposure 
cnt_phenos_recoded <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_LLD_clean_v2.txt", row.names = 1)
pheno_sel_cnt=unique(cnt_phenos_sig$phenotype)

#Rename
pheno_sel_cnt=gsub("Calprot200","ibd_FecalCalprotectinOver200yesno" , pheno_sel_cnt)
pheno_sel_cnt=gsub("ChrA","blood_ChromograninA" , pheno_sel_cnt)
pheno_sel_cnt=gsub("host_SmokeCurrentSmoker","host_smk_current" , pheno_sel_cnt)
pheno_sel_cnt=gsub("HBD2","blood_BetaDefensin2", pheno_sel_cnt)
pheno_sel_cnt=gsub("antrop_SBP","host_SBP" , pheno_sel_cnt)
pheno_sel_cnt=gsub("Biochem_Glucose","blood_Glucose" , pheno_sel_cnt)
pheno_sel_cnt=gsub("Biochem_Insulin","blood_Insulin" , pheno_sel_cnt)
pheno_sel_cnt=gsub("Biochem_LDL","blood_LDL" , pheno_sel_cnt)
pheno_sel_cnt=gsub("Biochem_TG","blood_TG" , pheno_sel_cnt)
pheno_sel_cnt=gsub("BlCells_Lympho","blood_Lympho" , pheno_sel_cnt)
cnt_phenos_sig$phenotype=gsub("Calprot200","ibd_FecalCalprotectinOver200yesno" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("ChrA","blood_ChromograninA" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("host_SmokeCurrentSmoker","host_smk_current" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("HBD2","blood_BetaDefensin2", cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("antrop_SBP","host_SBP" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("Biochem_Glucose","blood_Glucose" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("Biochem_Insulin","blood_Insulin" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("Biochem_LDL","blood_LDL" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("Biochem_TG","blood_TG" , cnt_phenos_sig$phenotype)
cnt_phenos_sig$phenotype=gsub("BlCells_Lympho","blood_Lympho" , cnt_phenos_sig$phenotype)
#Recode
control_exposure=cnt_phenos_recoded[,colnames(cnt_phenos_recoded)%in%pheno_sel_cnt]
control_exposure$ibd_FecalCalprotectinOver200yesno[is.na(control_exposure$ibd_FecalCalprotectinOver200yesno)]=0
control_exposure$ibd_FecalCalprotectinOver200yesno[control_exposure$ibd_FecalCalprotectinOver200yesno=="no"]=0
control_exposure$ibd_FecalCalprotectinOver200yesno[control_exposure$ibd_FecalCalprotectinOver200yesno=="yes"]=1
control_exposure$clinical_Strict_IBS[is.na(control_exposure$clinical_Strict_IBS)]=0
control_exposure$clinical_Strict_IBS[control_exposure$clinical_Strict_IBS=="no"]=0
control_exposure$clinical_Strict_IBS[control_exposure$clinical_Strict_IBS=="yes"]=1
control_exposure=mutate_all(control_exposure, function(x) as.numeric(as.character(x)))
#diet=colnames(control_exposure)[grep("diet",colnames(control_exposure))]
for (a in 1:ncol(control_exposure)){
  control_exposure[,a][is.na(control_exposure[,a])]=median(control_exposure[,a], na.rm = T)
}
control_exposure=as.data.frame(control_exposure)


#Select mediator: bacterial species (n=109)
control_mediator=controls_data[,grep("sp_",colnames(controls_data))]
control_mediator=as.data.frame(mutate_all(control_mediator, function(x) as.numeric(as.character(x))))
#Select outcome
control_outcome=controls_data[,colnames(controls_data)%in%cnt_metabolites]
control_outcome=as.data.frame(mutate_all(control_outcome, function(x) as.numeric(as.character(x))))
control_outcome=as.data.frame(scale (control_outcome))

control_exposure=control_exposure[match(row.names(control_outcome),row.names(control_exposure)),] 

#Regress confounders before mediation analysis
phenos_regression <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/phenos_regression.txt", row.names=1)
row.names(phenos_regression)=phenos_regression$PID2
phenos_regression$PID2=NULL

control_outcome2=merge(phenos_regression, control_outcome, by="row.names")
control_outcome2$PF_RD=NULL
control_outcome2$clinical_BowelMovementADayDef[is.na(control_outcome2$clinical_BowelMovementADayDef)]=median(control_outcome2$clinical_BowelMovementADayDef, na.rm = T)

regressors=colnames(control_outcome2)[2:26]

control_outcome3=matrix(nrow=nrow(control_outcome2), ncol = ncol(control_outcome))
n=1
for (i in 27:ncol(control_outcome2)){
  my_trait=colnames(control_outcome2)[i]
  my_uni_test=control_outcome2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  control_outcome3[,n] = my_lm$residuals
  n=n+1
}
control_outcome3=as.data.frame(control_outcome3)
rownames(control_outcome3)=control_outcome2$Row.names
colnames(control_outcome3)=colnames(control_outcome)

control_outcome3=as.data.frame(scale(control_outcome3))


#Regress phenotyes on bacteria

control_mediator2=merge(phenos_regression, control_mediator, by="row.names")
regressors=c("PF_RD","clinical_BowelMovementADayDef","host_Age_sampling","host_BMI","host_Sex" )

control_mediator2$clinical_BowelMovementADayDef[is.na(control_mediator2$clinical_BowelMovementADayDef)]=median(control_mediator2$clinical_BowelMovementADayDef, na.rm = T)


control_mediator3=matrix(nrow=nrow(control_mediator), ncol = ncol(control_mediator))
n=1
for (i in 28:ncol(control_mediator2)){
  my_trait=colnames(control_mediator2)[i]
  my_uni_test=control_mediator2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  control_mediator3[,n] = my_lm$residuals
  n=n+1
}
control_mediator3=as.data.frame(control_mediator3)
rownames(control_mediator3)=control_mediator2$Row.names
colnames(control_mediator3)=colnames(control_mediator)

control_mediator3=as.data.frame(scale(control_mediator3))


####

control_exposure[,18:39]=scale(control_exposure[,18:39])
control_outcome_test=control_outcome3[,unique(cnt_phenos_sig$metabolite)]
control_exposure_test=control_exposure

### Test influence of diet on bacteria and select taxa 

flag=1
for ( i in 1:ncol(control_exposure_test)) {
  my_pheno=colnames(control_exposure_test)[i]
  print (paste("Testing:",my_pheno))
  for (a in 1:ncol(control_mediator3)){
    my_trait=colnames(control_mediator3)[a]
    my_data=as.data.frame(cbind(control_exposure_test[,my_pheno],control_mediator3[,my_trait]))
    colnames(my_data)=c(my_pheno, my_trait)
    my_f=as.formula(paste(my_trait, my_pheno, sep = " ~ "))
    my_lm=summary(lm(my_f,data = my_data ))
    #Extract regression results
    my_lm_coef=as.data.frame(my_lm$coefficients)
    #Subset the results of the phenotype of interest
    my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
    my_lm_pheno$trait=my_trait
    my_lm_pheno$phenotype=my_pheno
    my_lm_pheno$test_group=gsub(my_pheno,"", rownames(my_lm_pheno))
    #Merge results (if this is not the first iteration, then merge with previous results)
    if (flag!=1){
      my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
    }else{
      my_univariate_results=my_lm_pheno
      flag=5
    } 
  }
}
  
# Select bacteria (mediator) associated with at least 1 dietary element (exposure) (here I use a relaxed threshold p<0.05, we could consider FDR<0.05 for a more strict analysis)
selected=subset(my_univariate_results, my_univariate_results$`Pr(>|t|)`<0.05)


cnt_phenos_sig=cnt_phenos_sig[-214,]


flag=9
lambda.grid <- c(0.4,0.3,0.2,0.1,0.08,0.04,0.03,0.02,0.01)
best_lambdas_cnt=matrix(ncol=2, nrow=nrow(cnt_phenos_sig))
results_mediation_cnt=data.frame(Species=colnames(control_mediator), Cohort="Control")
sum_fit_report_final=data.frame(Exposure= NA, Mediator=NA, Outcome=NA, alpha=NA, beta=NA, alphaxbeta=NA,delta=NA , pvalue=NA)

for ( i in 1:nrow(cnt_phenos_sig)){
  #Define exposure 
  my_exposure=as.character(cnt_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(cnt_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=control_mediator3[,select_mediators$trait, drop=F]
  #If there is only one mediator then we can simply use the mediate function
  if (ncol(my_mediator)==1){
    #Get mediator name
    name_mediator=colnames(my_mediator)
    #Create a dataframe containing the elements to test
    my_dataset=data.frame(Exposure=control_exposure[,my_exposure], Outcome=control_outcome3[,my_outcome],Mediator=my_mediator)
    # Give the right names (not necessary but nice for direct network plotting with regmed)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    # First we model the impact of the exposure to the mediator (alpha)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    my_alpha=as.numeric(fit.mediator$coefficients[2])
    # Next we model the outcome (dependent var) in relation with the exposure (independent var) while controlling for the mediator
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    #extract impact of mediator to the outcome (beta)
    my_beta=as.numeric(fit.out$coefficients[3])
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    # Run mediation analysis 
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    # Extract direct effect exposure to outcome (delta)
    my_delta1=summary(my_mediation$model.y)
    my_delta1.1=my_delta1$coefficients[2,1]
    #If there is a significant mediation effect, then report the results
    if(my_mediation$d0.p<0.05){
      sum_fit=data.frame(a=my_alpha, b=my_beta, x=my_alpha*my_beta)
      sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta,delta=my_delta1.1, pvalue=my_mediation$d0.p)
      sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta")
      row.names(sum_fit)=name_mediator
    }else{
      my_mess=sprintf("No significant mediation: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
      print(my_mess)
    }
  # If we have more than one (potential) mediator, then we run regmed (regularized mediationn analysis) 
  }else {
    # Warning: Standarization set off because was performed before
    my_fit<- regmed.grid(control_exposure[,my_exposure],my_mediator,control_outcome3[,my_outcome],lambda.grid, frac.lasso = 0.8, x.std = F, med.std = F,print.iter=T, max.outer=2000)
    # We extract the best model
    fit.best <- regmed.grid.bestfit(my_fit)
    # Extract estimated delta 
    my_delta=fit.best$delta
    # Make an overview of the best lambda for each model
    best_fit=my_fit$grid.data$lambda[my_fit$grid.data$bic==fit.best$bic]
    #my_mess=sprintf("%s -> %s best lambda  %s", my_exposure, my_outcome, best_fit)
    #print(my_mess)
    best_lambdas_cnt[i,1]=sprintf("%s -> %s", my_exposure, my_outcome)
    best_lambdas_cnt[i,2]=best_fit[1]
    #Summarize results of regmed
    sum_fit=data.frame(a=fit.best$alpha, b=fit.best$beta, x=fit.best$alpha*fit.best$beta, y=my_delta)
    colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta", paste(my_exposure,my_outcome, sep = "->"))
    results_mediation_cnt=merge(results_mediation_cnt,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    # If there is mediation (effect exposure->mediator x mediator -> outcome != 0) then we test each relation with with mediate package
    # This is not necessary but can be useful for: 
    # a) consistency with methods used in case of one exposure - one mediator - one outcome
    # b) to get a p-value of each mediated relation
    if (sum(sum_fit$alphaxbeta)!=0){
      #Subset mediators
      sum_fit2=subset(sum_fit, sum_fit$alphaxbeta!=0)
      for (x in 1:nrow(sum_fit2)){
        # This is basically the same as done in the case of 1 exposure 1 mediator 1 outcome (check above for commented code)
        name_mediator=row.names(sum_fit2)[x]
        my_dataset=data.frame(Exposure=control_exposure[,my_exposure], Outcome=control_outcome3[,my_outcome],Mediator=my_mediator[,name_mediator])
        colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
        f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
        fit.mediator=lm(f1, my_dataset)
        my_alpha=as.numeric(fit.mediator$coefficients[2])
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        my_beta=as.numeric(fit.out$coefficients[3])
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        my_delta1=summary(my_mediation$model.y)
        my_delta1.1=my_delta1$coefficients[2,1]
        sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta,delta=my_delta1.1,  pvalue=my_mediation$d0.p)
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}




#################################################

## In UC - 1 exposure, multiple mediator, 1 outcome
#################################################


#Select UC samples
IBD_phenos_recoded <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_IBD_clean_v2.txt", row.names = 1)
IBD_data2=subset(d3,d3$IBD==2)
uc_samples=row.names(IBD_phenos_recoded)[IBD_phenos_recoded$ibd_Diagnosis=="UC" & IBD_phenos_recoded$ibd_Stoma=="No" ]

#Rename to match files
pheno_sel_uc=uc_phenos_sig$phenotype
pheno_sel_uc=c(pheno_sel_uc,"ibd_FecalCalprotectinOver200yesno.x")
pheno_sel_uc=gsub("Calprot200","ibd_FecalCalprotectinOver200yesno.x" , pheno_sel_uc)
pheno_sel_uc=gsub("clinical_IleocecalValveInSitu","ibd_IleocecalValveInSitu" , pheno_sel_uc)
pheno_sel_uc=gsub("clinical_NumberOfResectionColonic","ibd_NumberOfResectionColonic", pheno_sel_uc)
pheno_sel_uc=gsub("clinical_ResectionAny","ibd_ResectionAny" , pheno_sel_uc)
pheno_sel_uc=gsub("clinical_ResectionColonicAny","ibd_ResectionColonicAny" , pheno_sel_uc)

#Re-code some variables
uc_exposure=IBD_phenos_recoded[,colnames(IBD_phenos_recoded)%in%pheno_sel_uc]
uc_exposure$med_BileAcids[is.na(uc_exposure$med_BileAcids)]=0
uc_exposure$med_BileAcids[uc_exposure$med_BileAcids=="Non_user"]=0
uc_exposure$med_BileAcids[uc_exposure$med_BileAcids=="User"]=1
uc_exposure$med_mesalazines[is.na(uc_exposure$med_mesalazines)]=0
uc_exposure$med_mesalazines[uc_exposure$med_mesalazines=="Non_user"]=0
uc_exposure$med_mesalazines[uc_exposure$med_mesalazines=="User"]=1
uc_exposure$meds_Immunosuppressants[is.na(uc_exposure$meds_Immunosuppressants)]=0
uc_exposure$meds_Immunosuppressants[uc_exposure$meds_Immunosuppressants=="Non_user"]=0
uc_exposure$meds_Immunosuppressants[uc_exposure$meds_Immunosuppressants=="User"]=1
uc_exposure$ibd_FecalCalprotectinOver200yesno.x[is.na(uc_exposure$ibd_FecalCalprotectinOver200yesno.x)]=0
uc_exposure$ibd_FecalCalprotectinOver200yesno.x[uc_exposure$ibd_FecalCalprotectinOver200yesno.x=="no"]=0
uc_exposure$ibd_FecalCalprotectinOver200yesno.x[uc_exposure$ibd_FecalCalprotectinOver200yesno.x=="yes"]=1
uc_exposure$ibd_IleocecalValveInSitu[is.na(uc_exposure$ibd_IleocecalValveInSitu)]=0
uc_exposure$ibd_IleocecalValveInSitu[uc_exposure$ibd_IleocecalValveInSitu=="no"]=0
uc_exposure$ibd_IleocecalValveInSitu[uc_exposure$ibd_IleocecalValveInSitu=="yes"]=1
uc_exposure$ibd_ResectionAny[is.na(uc_exposure$ibd_ResectionAny)]=0
uc_exposure$ibd_ResectionAny[uc_exposure$ibd_ResectionAny=="no"]=0
uc_exposure$ibd_ResectionAny[uc_exposure$ibd_ResectionAny=="yes"]=1
uc_exposure$ibd_ResectionColonicAny[is.na(uc_exposure$ibd_ResectionColonicAny)]=0
uc_exposure$ibd_ResectionColonicAny[uc_exposure$ibd_ResectionColonicAny=="no"]=0
uc_exposure$ibd_ResectionColonicAny[uc_exposure$ibd_ResectionColonicAny=="yes"]=1
uc_exposure$ibd_NumberOfResectionColonic[is.na(uc_exposure$ibd_NumberOfResectionColonic)]=0
uc_exposure=mutate_all(uc_exposure, function(x) as.numeric(as.character(x)))
diet=colnames(uc_exposure)[grep("diet",colnames(uc_exposure))]
for (a in diet){
  uc_exposure[,a][is.na(uc_exposure[,a])]=median(uc_exposure[,a], na.rm = T)
}

uc_exposure=as.data.frame(uc_exposure)

###

uc_mediator=IBD_data2[row.names(IBD_data2)%in%uc_samples,grep("sp_",colnames(IBD_data2))]
uc_outcome=IBD_data2[row.names(IBD_data2)%in%uc_samples,444:ncol(IBD_data2)]
uc_exposure=subset(uc_exposure,row.names(uc_outcome)%in%row.names(uc_exposure))
uc_exposure=uc_exposure[match(row.names(uc_outcome),row.names(uc_exposure)),] 

#Regress phenotypes on metabolites

phenos_regression <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/phenos_regression.txt", row.names=1)
row.names(phenos_regression)=phenos_regression$PID2
phenos_regression$PID2=NULL

uc_outcome2=merge(phenos_regression, uc_outcome, by="row.names")
uc_outcome2$PF_RD=NULL
uc_outcome2$valeric_acid[is.na(uc_outcome2$valeric_acid)]=min(uc_outcome2$valeric_acid, na.rm = T)/2
uc_outcome2$acetic_acid[is.na(uc_outcome2$acetic_acid)]=min(uc_outcome2$acetic_acid, na.rm = T)/2
uc_outcome2$butyric_acid[is.na(uc_outcome2$butyric_acid)]=min(uc_outcome2$butyric_acid, na.rm = T)/2
uc_outcome2$Methylbutyric_acid[is.na(uc_outcome2$Methylbutyric_acid)]=min(uc_outcome2$Methylbutyric_acid, na.rm = T)/2
uc_outcome2$hexanoic_acid[is.na(uc_outcome2$hexanoic_acid)]=min(uc_outcome2$hexanoic_acid, na.rm = T)/2
uc_outcome2$isobutyric_acid[is.na(uc_outcome2$isobutyric_acid)]=min(uc_outcome2$isobutyric_acid, na.rm = T)/2
uc_outcome2$isovaleric_acid[is.na(uc_outcome2$isovaleric_acid)]=min(uc_outcome2$isovaleric_acid, na.rm = T)/2
uc_outcome2$propionic_acid[is.na(uc_outcome2$propionic_acid)]=min(uc_outcome2$propionic_acid, na.rm = T)/2

regressors=colnames(uc_outcome2)[2:26]

uc_outcome3=matrix(nrow=nrow(uc_outcome2), ncol = ncol(uc_outcome))
n=1
for (i in 27:ncol(uc_outcome2)){
  my_trait=colnames(uc_outcome2)[i]
  my_uni_test=uc_outcome2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  uc_outcome3[,n] = my_lm$residuals
  n=n+1
}
uc_outcome3=as.data.frame(uc_outcome3)
rownames(uc_outcome3)=uc_outcome2$Row.names
colnames(uc_outcome3)=colnames(uc_outcome)

uc_outcome3=as.data.frame(scale(uc_outcome3))


#Regress phenotyes on bacteria

uc_mediator2=merge(phenos_regression, uc_mediator, by="row.names")
regressors=c("PF_RD","clinical_BowelMovementADayDef","host_Age_sampling","host_BMI","host_Sex" )

uc_mediator3=matrix(nrow=nrow(uc_mediator2), ncol = ncol(uc_mediator))
n=1
for (i in 28:ncol(uc_mediator2)){
  my_trait=colnames(uc_mediator2)[i]
  my_uni_test=uc_mediator2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  uc_mediator3[,n] = my_lm$residuals
  n=n+1
}
uc_mediator3=as.data.frame(uc_mediator3)
rownames(uc_mediator3)=uc_mediator2$Row.names
colnames(uc_mediator3)=colnames(uc_mediator)

uc_mediator3=as.data.frame(scale(uc_mediator3))


####

uc_exposure[,1:8]=scale(uc_exposure[,1:8])
uc_phenos_sig$phenotype[uc_phenos_sig$phenotype=="Calprot200"]="ibd_FecalCalprotectinOver200yesno.x"
uc_phenos_sig$phenotype[uc_phenos_sig$phenotype=="clinical_IleocecalValveInSitu"]="ibd_IleocecalValveInSitu"
uc_phenos_sig$phenotype[uc_phenos_sig$phenotype=="clinical_NumberOfResectionColonic"]="ibd_NumberOfResectionColonic"
uc_phenos_sig$phenotype[uc_phenos_sig$phenotype=="clinical_ResectionAny"]="ibd_ResectionAny"
uc_phenos_sig$phenotype[uc_phenos_sig$phenotype=="clinical_ResectionColonicAny"]="ibd_ResectionColonicAny"

### Test influence of exposure on mediator (phenotype on bacteria)


uc_outcome_test=uc_outcome3[,unique(uc_phenos_sig$metabolite)]
uc_exposure_test=uc_exposure


flag=1
for ( i in 1:ncol(uc_exposure_test)) {
  my_pheno=colnames(uc_exposure_test)[i]
  print (paste("Testing:",my_pheno))
  for (a in 1:ncol(uc_mediator3)){
    my_trait=colnames(uc_mediator3)[a]
    my_data=as.data.frame(cbind(uc_exposure_test[,my_pheno],uc_mediator3[,my_trait]))
    colnames(my_data)=c(my_pheno, my_trait)
    my_f=as.formula(paste(my_trait, my_pheno, sep = " ~ "))
    my_lm=summary(lm(my_f,data = my_data ))
    #Extract regression results
    my_lm_coef=as.data.frame(my_lm$coefficients)
    #Subset the results of the phenotype of interest
    my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
    my_lm_pheno$trait=my_trait
    my_lm_pheno$phenotype=my_pheno
    my_lm_pheno$test_group=gsub(my_pheno,"", rownames(my_lm_pheno))
    #Merge results (if this is not the first iteration, then merge with previous results)
    if (flag!=1){
      my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
    }else{
      my_univariate_results=my_lm_pheno
      flag=5
    } 
  }
}

# Select bacteria associated with at least 1 dietary element (here I use a relaxed threshold p<0.05, we could consider FDR<0.05 for a more strict analysis)

selected=subset(my_univariate_results, my_univariate_results$`Pr(>|t|)`<0.05)


flag=9
lambda.grid <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01)
best_lambdas_uc=matrix(ncol=2, nrow=nrow(uc_phenos_sig))
results_mediation_uc=data.frame(Species=colnames(uc_mediator), Cohort="UC")
UC_sum_fit_report_final=data.frame(Exposure= NA, Mediator=NA, Outcome=NA, alpha=NA, beta=NA, alphaxbeta=NA, pvalue=NA)

for ( i in 1:nrow(uc_phenos_sig)){
  my_exposure=as.character(uc_phenos_sig[i,1])
  my_outcome=as.character(uc_phenos_sig[i,2])
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=uc_mediator3[,select_mediators$trait, drop=F]
  if (ncol(my_mediator)==1){
    name_mediator=colnames(my_mediator)
    my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=uc_outcome3[,my_outcome],Mediator=my_mediator)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    my_alpha=as.numeric(fit.mediator$coefficients[2])
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    my_beta=as.numeric(fit.out$coefficients[3])
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    if(my_mediation$d0.p<0.05){
      sum_fit=data.frame(a=my_alpha, b=my_beta, x=my_alpha*my_beta)
      sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta, pvalue=my_mediation$d0.p)
      sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta")
      row.names(sum_fit)=name_mediator
    }else{
      my_mess=sprintf("No significant mediation: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
      print(my_mess)
    }
  }else {
    #standarization set off because was performed before
    my_fit<- regmed.grid(uc_exposure[,my_exposure],my_mediator,uc_outcome3[,my_outcome],lambda.grid, frac.lasso = 0.8, x.std = F, med.std = F)
    fit.best <- regmed.grid.bestfit(my_fit)
    best_fit=my_fit$grid.data$lambda[my_fit$grid.data$bic==fit.best$bic]
    my_mess=sprintf("%s -> %s best lambda  %s", my_exposure, my_outcome, best_fit)
    print(my_mess)
    best_lambdas_uc[i,1]=sprintf("%s -> %s", my_exposure, my_outcome)
    best_lambdas_uc[i,2]=best_fit[1]
    sum_fit=data.frame(a=fit.best$alpha, b=fit.best$beta, x=fit.best$alpha*fit.best$beta)
    colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta")
    results_mediation_uc=merge(results_mediation_uc,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    if (sum(sum_fit$alphaxbeta)!=0){
      sum_fit2=subset(sum_fit, sum_fit$alphaxbeta!=0)
      for (x in 1:nrow(sum_fit2)){
        name_mediator=row.names(sum_fit2)[x]
        my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=uc_outcome3[,my_outcome],Mediator=my_mediator[,name_mediator])
        colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
        f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
        fit.mediator=lm(f1, my_dataset)
        my_alpha=as.numeric(fit.mediator$coefficients[2])
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        my_beta=as.numeric(fit.out$coefficients[3])
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta, pvalue=my_mediation$d0.p)
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}


#########################################

## RUN CD ##

#########################################



#Select CD samples
IBD_phenos_recoded <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_IBD_clean_v2.txt", row.names = 1)
stomas=row.names(IBD_phenos_recoded)[IBD_phenos_recoded$ibd_Stoma=="Yes" ]

IBD_phenos_recoded <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/IBD_phenos_recoded.txt", row.names = 1)
IBD_data2=subset(d3,d3$IBD==2)
cd_samples=row.names(IBD_phenos_recoded)[IBD_phenos_recoded$clinical_Diagnosis=="2" ]

#Rename to match files
pheno_sel_cd=cd_phenos_sig$phenotype
cd_exposure=IBD_phenos_recoded[,colnames(IBD_phenos_recoded)%in%pheno_sel_cd]


#Re-code some variables

for (a in 1:ncol(cd_exposure)){
  cd_exposure[,a][is.na(cd_exposure[,a])]=median(cd_exposure[,a], na.rm = T)
}

cd_exposure=as.data.frame(cd_exposure)

###

cd_mediator=IBD_data2[row.names(IBD_data2)%in%cd_samples,grep("sp_",colnames(IBD_data2))]
cd_outcome=IBD_data2[row.names(IBD_data2)%in%cd_samples,444:ncol(IBD_data2)]
cd_exposure=subset(cd_exposure,row.names(cd_outcome)%in%row.names(cd_exposure))
cd_exposure=cd_exposure[match(row.names(cd_outcome),row.names(cd_exposure)),] 

#Regress phenotypes on metabolites

phenos_regression <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/phenos_regression.txt", row.names=1)
row.names(phenos_regression)=phenos_regression$PID2
phenos_regression$PID2=NULL

cd_outcome2=merge(phenos_regression, cd_outcome, by="row.names")
cd_outcome2$PF_RD=NULL
cd_outcome2$valeric_acid[is.na(cd_outcome2$valeric_acid)]=min(cd_outcome2$valeric_acid, na.rm = T)/2
cd_outcome2$acetic_acid[is.na(cd_outcome2$acetic_acid)]=min(cd_outcome2$acetic_acid, na.rm = T)/2
cd_outcome2$butyric_acid[is.na(cd_outcome2$butyric_acid)]=min(cd_outcome2$butyric_acid, na.rm = T)/2
cd_outcome2$Methylbutyric_acid[is.na(cd_outcome2$Methylbutyric_acid)]=min(cd_outcome2$Methylbutyric_acid, na.rm = T)/2
cd_outcome2$hexanoic_acid[is.na(cd_outcome2$hexanoic_acid)]=min(cd_outcome2$hexanoic_acid, na.rm = T)/2
cd_outcome2$isobutyric_acid[is.na(cd_outcome2$isobutyric_acid)]=min(cd_outcome2$isobutyric_acid, na.rm = T)/2
cd_outcome2$isovaleric_acid[is.na(cd_outcome2$isovaleric_acid)]=min(cd_outcome2$isovaleric_acid, na.rm = T)/2
cd_outcome2$propionic_acid[is.na(cd_outcome2$propionic_acid)]=min(cd_outcome2$propionic_acid, na.rm = T)/2

regressors=colnames(cd_outcome2)[2:26]

cd_outcome3=matrix(nrow=nrow(cd_outcome2), ncol = ncol(cd_outcome))
n=1
for (i in 27:ncol(cd_outcome2)){
  my_trait=colnames(cd_outcome2)[i]
  my_uni_test=cd_outcome2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  cd_outcome3[,n] = my_lm$residuals
  n=n+1
}
cd_outcome3=as.data.frame(cd_outcome3)
rownames(cd_outcome3)=cd_outcome2$Row.names
colnames(cd_outcome3)=colnames(cd_outcome)

cd_outcome3=as.data.frame(scale(cd_outcome3))


#Regress phenotyes on bacteria

cd_mediator2=merge(phenos_regression, cd_mediator, by="row.names")
regressors=c("PF_RD","clinical_BowelMovementADayDef","host_Age_sampling","host_BMI","host_Sex" )

cd_mediator3=matrix(nrow=nrow(cd_mediator2), ncol = ncol(cd_mediator))
n=1
for (i in 28:ncol(cd_mediator2)){
  my_trait=colnames(cd_mediator2)[i]
  my_uni_test=cd_mediator2[,c(regressors,my_trait)]
  my_f=as.formula(paste(my_trait, paste(regressors, collapse = " + "), sep = " ~ "))
  my_lm=lm(my_f,data = my_uni_test )
  cd_mediator3[,n] = my_lm$residuals
  n=n+1
}
cd_mediator3=as.data.frame(cd_mediator3)
rownames(cd_mediator3)=cd_mediator2$Row.names
colnames(cd_mediator3)=colnames(cd_mediator)

cd_mediator3=as.data.frame(scale(cd_mediator3))


####

cd_exposure$clinical_DiseaseDurationYears.x=scale(cd_exposure$clinical_DiseaseDurationYears.x)[,1]
cd_exposure$ChrA=scale(cd_exposure$ChrA)[,1]
cd_exposure$HBD2=scale(cd_exposure$HBD2)[,1]
cd_exposure[,37:84]=scale(cd_exposure[,37:84])


### Test influence of exposure on mediator (phenotype on bacteria)


cd_outcome_test=cd_outcome3[,unique(cd_phenos_sig$metabolite)]
cd_exposure_test=cd_exposure


flag=1
for ( i in 1:ncol(cd_exposure_test)) {
  my_pheno=colnames(cd_exposure_test)[i]
  print (paste("Testing:",my_pheno))
  for (a in 1:ncol(cd_mediator3)){
    my_trait=colnames(cd_mediator3)[a]
    my_data=as.data.frame(cbind(cd_exposure_test[,my_pheno],cd_mediator3[,my_trait]))
    colnames(my_data)=c(my_pheno, my_trait)
    my_f=as.formula(paste(my_trait, my_pheno, sep = " ~ "))
    my_lm=summary(lm(my_f,data = my_data ))
    #Extract regression results
    my_lm_coef=as.data.frame(my_lm$coefficients)
    #Subset the results of the phenotype of interest
    my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
    my_lm_pheno$trait=my_trait
    my_lm_pheno$phenotype=my_pheno
    my_lm_pheno$test_group=gsub(my_pheno,"", rownames(my_lm_pheno))
    #Merge results (if this is not the first iteration, then merge with previous results)
    if (flag!=1){
      my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
    }else{
      my_univariate_results=my_lm_pheno
      flag=5
    } 
  }
}

# Select bacteria associated with at least 1 dietary element (here I use a relaxed threshold p<0.05, we could consider FDR<0.05 for a more strict analysis)

selected=subset(my_univariate_results, my_univariate_results$`Pr(>|t|)`<0.05)


flag=9
#lambda.grid <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01)
lambda.grid <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
best_lambdas_cd=matrix(ncol=2, nrow=nrow(cd_phenos_sig))
results_mediation_cd=data.frame(Species=colnames(cd_mediator), Cohort="CD")
CD_sum_fit_report_final=data.frame(Exposure= NA, Mediator=NA, Outcome=NA, alpha=NA, beta=NA, alphaxbeta=NA, pvalue=NA)

for ( i in 1:nrow(cd_phenos_sig)){
  my_exposure=as.character(cd_phenos_sig[i,1])
  my_outcome=as.character(cd_phenos_sig[i,2])
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=cd_mediator3[,select_mediators$trait, drop=F]
  if (ncol(my_mediator)==0){
    print ("No association between exposure and mediators")
  } else if (ncol(my_mediator)==1){
    name_mediator=colnames(my_mediator)
    my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=cd_outcome3[,my_outcome],Mediator=my_mediator)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    my_alpha=as.numeric(fit.mediator$coefficients[2])
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    my_beta=as.numeric(fit.out$coefficients[3])
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    if(my_mediation$d0.p<0.05){
      sum_fit=data.frame(a=my_alpha, b=my_beta, x=my_alpha*my_beta)
      sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta, pvalue=my_mediation$d0.p)
      sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta")
      row.names(sum_fit)=name_mediator
    }else{
      my_mess=sprintf("No significant mediation: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
      print(my_mess)
    }
  }else {
    #standarization set off because was performed before
    my_fit<- regmed.grid(cd_exposure[,my_exposure],my_mediator,cd_outcome3[,my_outcome],lambda.grid, frac.lasso = 0.8, x.std = F, med.std = F)
    fit.best <- regmed.grid.bestfit(my_fit)
    best_fit=my_fit$grid.data$lambda[my_fit$grid.data$bic==fit.best$bic]
    my_mess=sprintf("%s -> %s best lambda  %s", my_exposure, my_outcome, best_fit)
    print(my_mess)
    best_lambdas_cd[i,1]=sprintf("%s -> %s", my_exposure, my_outcome)
    best_lambdas_cd[i,2]=best_fit[1]
    sum_fit=data.frame(a=fit.best$alpha, b=fit.best$beta, x=fit.best$alpha*fit.best$beta)
    colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta")
    results_mediation_cd=merge(results_mediation_cd,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    if (sum(sum_fit$alphaxbeta)!=0){
      sum_fit2=subset(sum_fit, sum_fit$alphaxbeta!=0)
      for (x in 1:nrow(sum_fit2)){
        name_mediator=row.names(sum_fit2)[x]
        my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=cd_outcome3[,my_outcome],Mediator=my_mediator[,name_mediator])
        colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
        f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
        fit.mediator=lm(f1, my_dataset)
        my_alpha=as.numeric(fit.mediator$coefficients[2])
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        my_beta=as.numeric(fit.out$coefficients[3])
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        sum_fit_report=data.frame(Exposure= my_exposure, Mediator=name_mediator, Outcome=my_outcome, alpha=my_alpha, beta=my_beta, alphaxbeta=my_alpha*my_beta, pvalue=my_mediation$d0.p)
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}

