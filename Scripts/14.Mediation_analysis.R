#Mediation analysis

#This code is still a bit messy, I will clean it an implemented as a function soon (sorry!)

load("~/Desktop/IBD_metabolomics_2022/3.Workspace/Lasso.RData")

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

#Select exposure 
cnt_phenos_recoded <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_LLD_clean_v2.txt", row.names = 1)
pheno_sel_cnt=unique(cnt_phenos_sig$phenotype)

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
control_exposure$blood_ChromograninA=(scale(log(control_exposure$blood_ChromograninA)))[,1]
control_exposure$host_SBP=(scale(log(control_exposure$host_SBP)))[,1]
control_exposure$blood_Lympho=(scale(log(control_exposure$blood_Lympho)))[,1]
control_exposure$blood_Glucose=(scale(log(control_exposure$blood_Glucose)))[,1]
control_exposure$blood_LDL=(scale(log(control_exposure$blood_LDL)))[,1]
control_exposure$blood_TG=(scale(log(control_exposure$blood_TG)))[,1]
control_exposure$blood_Insulin=(scale(log(control_exposure$blood_Insulin)))[,1]

###########################

#Run mediation in Controls

###########################

#
# Function copied from: https://stackoverflow.com/questions/41582486/how-to-convert-r-mediation-summary-to-data-frame
#
extract_mediation_summary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat}



flag=9
lambda.grid <- c(0.4,0.3,0.2,0.1,0.08,0.04,0.03,0.02,0.01)
best_lambdas_cnt=matrix(ncol=2, nrow=nrow(cnt_phenos_sig))

output_log=data.frame(Test=NA,Message=NA)
results_mediation_cnt=data.frame(F1=NA,F2=NA,type=NA,coeff=NA)

sum_fit_report_final=data.frame(Exposure= NA, 
                                Mediator=NA, 
                                Outcome=NA, 
                                total_est=NA,
                                total_p=NA,
                                impact_on_mediator_est=NA,
                                impact_on_mediator_p=NA,
                                Mediated_est=NA,
                                Mediated_p=NA,
                                ACME_estimate=NA,
                                ACME_p_value=NA, 
                                ADE_estimate=NA, 
                                ADE_p_value=NA, 
                                Total_effect_estimate=NA, 
                                Total_effect_p_value=NA,
                                Prop_mediated_estimate=NA, 
                                Prop_mediated_p_value=NA )


my_univariate_results$fdr=p.adjust(my_univariate_results$`Pr(>|t|)`, method = "BH")
selected=subset(my_univariate_results, my_univariate_results$fdr<0.1)

for ( i in 1:nrow(cnt_phenos_sig)){
  #Define exposure 
  my_exposure=as.character(cnt_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(cnt_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=control_mediator3[,select_mediators$trait, drop=F]
  #If there is only one mediator then we can simply use the mediate function
  my_exposure=as.character(cnt_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(cnt_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=control_mediator3[,select_mediators$trait, drop=F]
  if (ncol(my_mediator)==0){
    #my_mess=paste(name_exposure,name_outcome, sep = "->")
    my_mess2="No mediator"
    print (my_mess2)
  }else if (ncol(my_mediator)==1){
    #Get mediator name
    name_mediator=colnames(my_mediator)
    #Create a dataframe containing the elements to test
    my_dataset=data.frame(Exposure=control_exposure[,my_exposure], Outcome=control_outcome3[,my_outcome],Mediator=my_mediator)
    # Give the right names (not necessary but nice for direct network plotting with regmed)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    # First we model the impact of the exposure to the mediator (alpha)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    # Next we model the outcome (dependent var) in relation with the exposure (independent var) while ucling for the mediator
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    #extract impact of mediator to the outcome (beta)
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    # Run mediation analysis 
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    my_mediation_df=extract_mediation_summary(my_mediation)
    
    # Calculate each step of the mediation
    
    # Step 1: total effect
    f0=as.formula(paste( my_outcome,my_exposure, sep = " ~ "))
    fit.total=summary(lm(f0, my_dataset))
    my_total=fit.total$coefficients[2,1]
    my_total_p=fit.total$coefficients[2,4]
    
    #Step 2: Efect of indepedent variable on the mediator 
    my_alpha1=summary(fit.mediator)
    my_alpha=my_alpha1$coefficients[2,1]
    my_alpha_p=my_alpha1$coefficients[2,4]
    
    #Step 3: Effect mediator on the dependent variable (ucling for the independent)
    my_delta1=summary(fit.out)
    my_delta=my_delta1$coefficients[3,1]
    my_delta_p=my_delta1$coefficients[3,4]
    
    sum_fit_report=data.frame(Exposure= my_exposure,
                              Mediator=name_mediator, 
                              Outcome=my_outcome, 
                              total_est=my_total,
                              total_p=my_total_p,
                              impact_on_mediator_est=my_alpha,
                              impact_on_mediator_p=my_alpha_p,
                              Mediated_est=my_delta,
                              Mediated_p=my_delta_p,
                              ACME_estimate=my_mediation_df[1,1], 
                              ACME_p_value=my_mediation_df[1,4],
                              ADE_estimate=my_mediation_df[2,1], 
                              ADE_p_value=my_mediation_df[2,4],
                              Total_effect_estimate=my_mediation_df[3,1], 
                              Total_effect_p_value=my_mediation_df[3,4],
                              Prop_mediated_estimate=my_mediation_df[4,1], 
                              Prop_mediated_p_value=my_mediation_df[4,4]
    )
    sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
    
    
    # Do the same but now: Exposure -> Metabolite -> Bacteria 
    # Our previous mediator is now the outcome and viceversa
    
    #Create a dataframe containing the elements to test
    #my_dataset=data.frame(Exposure=control_exposure[,my_exposure], Outcome=my_mediator,Mediator=control_outcome3[,my_outcome])
    # Give the right names (not necessary but nice for direct network plotting with regmed)
    #colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
    # First we model the impact of the exposure to the mediator (alpha)
    #f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
    #fit.mediator=lm(f1, my_dataset)
    
    # Next we model the outcome (dependent var) in relation with the exposure (independent var) while ucling for the mediator
    #f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
    #fit.out=lm(f2, my_dataset)
    
    #my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
    #print(my_mess)
    # Run mediation analysis 
    #my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
    # Extract direct effect exposure to outcome (delta)
    
    # Step 1: total effect
    #f0=as.formula(paste(my_exposure, name_mediator, sep = " ~ "))
    #fit.total=summary(lm(f0, my_dataset))
    #my_total=fit.total$coefficients[2,1]
    #my_total_p=fit.total$coefficients[2,4]
    
    #Step 2: Efect of indepedent variable on the mediator 
    #my_alpha1=summary(fit.mediator)
    #my_alpha=my_alpha1$coefficients[2,1]
    #my_alpha_p=my_alpha1$coefficients[2,4]
    
    #Step 3: Effect mediator on the dependent variable (ucling for the independent)
    #my_delta1=summary(fit.out)
    #my_delta=my_delta1$coefficients[3,1]
    #my_delta_p=my_delta1$coefficients[3,4]
    
    #my_mediation_df=extract_mediation_summary(my_mediation)
    #sum_fit_report=data.frame(Exposure= my_exposure, 
    #                          Mediator=my_outcome, 
    #                          Outcome=name_mediator, 
    #                         total_est=my_total,
    #                         total_p=my_total_p,
    #                          impact_on_mediator_est=my_alpha,
    #                          impact_on_mediator_p=my_alpha_p,
    #                          Mediated_est=my_delta,
    #                          Mediated_p=my_delta_p,
    #                          ACME_estimate=my_mediation_df[1,1], 
    #                          ACME_p_value=my_mediation_df[1,4],
    #                          ADE_estimate=my_mediation_df[2,1], 
    #                          ADE_p_value=my_mediation_df[2,4],
    #                          Total_effect_estimate=my_mediation_df[3,1], 
    #                          Total_effect_p_value=my_mediation_df[3,4],
    #                          Prop_mediated_estimate=my_mediation_df[4,1], 
    #                          Prop_mediated_p_value=my_mediation_df[4,4])
    
    #sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
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
    #results_mediation_uc_=merge(results_mediation_uc_,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    results_mediation_alpha=data.frame(F1=my_exposure, F2=row.names(fit.best$alpha), type="alpha", coeff=as.numeric(fit.best$alpha))
    results_mediation_beta=data.frame(F1=row.names(fit.best$beta), F2=my_outcome, type="beta", coeff=as.numeric(fit.best$beta))
    results_mediation_ab=data.frame(F1=my_exposure, F2=my_outcome, type="alphaxbeta", coeff=as.numeric(fit.best$alpha*fit.best$beta))
    results_mediation_delta=data.frame(F1=my_exposure, F2=my_outcome, type="delta", coeff=fit.best$delta)
    sum_fit3=bind_rows(results_mediation_alpha, results_mediation_beta, results_mediation_ab,results_mediation_delta)
    results_mediation_cnt=rbind(results_mediation_cnt,sum_fit3)
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
        
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        
        # Step 1: total effect
        f0=as.formula(paste( my_outcome,my_exposure, sep = " ~ "))
        fit.total=summary(lm(f0, my_dataset))
        my_total=fit.total$coefficients[2,1]
        my_total_p=fit.total$coefficients[2,4]
        
        #Step 2: Efect of indepedent variable on the mediator 
        my_alpha1=summary(fit.mediator)
        my_alpha=my_alpha1$coefficients[2,1]
        my_alpha_p=my_alpha1$coefficients[2,4]
        
        #Step 3: Effect mediator on the dependent variable (ucling for the independent)
        my_delta1=summary(fit.out)
        my_delta=my_delta1$coefficients[3,1]
        my_delta_p=my_delta1$coefficients[3,4]
        
        
        my_mediation_df=extract_mediation_summary(my_mediation)
        sum_fit_report=data.frame(Exposure= my_exposure, 
                                  Mediator=name_mediator, 
                                  Outcome=my_outcome, 
                                  total_est=my_total,
                                  total_p=my_total_p,
                                  impact_on_mediator_est=my_alpha,
                                  impact_on_mediator_p=my_alpha_p,
                                  Mediated_est=my_delta,
                                  Mediated_p=my_delta_p,
                                  ACME_estimate=my_mediation_df[1,1], 
                                  ACME_p_value=my_mediation_df[1,4],
                                  ADE_estimate=my_mediation_df[2,1], 
                                  ADE_p_value=my_mediation_df[2,4],
                                  Total_effect_estimate=my_mediation_df[3,1], 
                                  Total_effect_p_value=my_mediation_df[3,4],
                                  Prop_mediated_estimate=my_mediation_df[4,1], 
                                  Prop_mediated_p_value=my_mediation_df[4,4])
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
        
        # Do the same but now: Exposure -> Metabolite -> Bacteria 
        # Our previous mediator is now the outcome and viceversa
        
        #my_dataset=data.frame(Exposure=control_exposure[,my_exposure], Outcome=my_mediator[,name_mediator],Mediator=control_outcome3[,my_outcome])
        #colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
        #f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
        #fit.mediator=lm(f1, my_dataset)
        
        #f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
        #fit.out=lm(f2, my_dataset)
        
        #my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
        #print(my_mess)
        #my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
        
        # Step 1: total effect
        #f0=as.formula(paste(my_exposure, name_mediator, sep = " ~ "))
        #fit.total=summary(lm(f0, my_dataset))
        #my_total=fit.total$coefficients[2,1]
        #my_total_p=fit.total$coefficients[2,4]
        
        #Step 2: Efect of indepedent variable on the mediator 
        #my_alpha1=summary(fit.mediator)
        #my_alpha=my_alpha1$coefficients[2,1]
        #my_alpha_p=my_alpha1$coefficients[2,4]
        
        #Step 3: Effect mediator on the dependent variable (ucling for the independent)
        #my_delta1=summary(fit.out)
        #my_delta=my_delta1$coefficients[3,1]
        #my_delta_p=my_delta1$coefficients[3,4]
        
        #my_mediation_df=extract_mediation_summary(my_mediation)
        #sum_fit_report=data.frame(Exposure= my_exposure, 
        #                          Mediator=my_outcome,
        #                          Outcome=name_mediator, 
        #                          total_est=my_total,
        #                          total_p=my_total_p,
        #                          impact_on_mediator_est=my_alpha,
        #                          impact_on_mediator_p=my_alpha_p,
        #                          Mediated_est=my_delta,
        #                          Mediated_p=my_delta_p,
        #                          ACME_estimate=my_mediation_df[1,1], 
        #                          ACME_p_value=my_mediation_df[1,4],
        #                          ADE_estimate=my_mediation_df[2,1], 
        #                          ADE_p_value=my_mediation_df[2,4],
        #                          Total_effect_estimate=my_mediation_df[3,1], 
        #                          Total_effect_p_value=my_mediation_df[3,4],
        #                          Prop_mediated_estimate=my_mediation_df[4,1], 
        #                          Prop_mediated_p_value=my_mediation_df[4,4])
        #sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}

sum_fit_report_final$total_fdr=p.adjust(sum_fit_report_final$total_p, method = "BH")
sum_fit_report_final$impact_on_mediator_fdr=p.adjust(sum_fit_report_final$impact_on_mediator_p, method = "BH")
sum_fit_report_final$Mediated_fdr=p.adjust(sum_fit_report_final$Mediated_p, method = "BH")
sum_fit_report_final$ACME_fdr=p.adjust(sum_fit_report_final$ACME_p_value, method = "BH")
sum_fit_report_final$ADE_fdr=p.adjust(sum_fit_report_final$ADE_p_value, method = "BH")

###########################

#Run mediation in UC

###########################

flag=9
lambda.grid <- c(0.4,0.3,0.2,0.1,0.08,0.04,0.03,0.02,0.01)
best_lambdas_uc=matrix(ncol=2, nrow=nrow(uc_phenos_sig))

output_log=data.frame(Test=NA,Message=NA)
results_mediation_uc=data.frame(F1=NA,F2=NA,type=NA,coeff=NA)

sum_fit_report_final=data.frame(Exposure= NA, 
                                Mediator=NA, 
                                Outcome=NA, 
                                total_est=NA,
                                total_p=NA,
                                impact_on_mediator_est=NA,
                                impact_on_mediator_p=NA,
                                Mediated_est=NA,
                                Mediated_p=NA,
                                ACME_estimate=NA,
                                ACME_p_value=NA, 
                                ADE_estimate=NA, 
                                ADE_p_value=NA, 
                                Total_effect_estimate=NA, 
                                Total_effect_p_value=NA,
                                Prop_mediated_estimate=NA, 
                                Prop_mediated_p_value=NA )


my_univariate_results$fdr=p.adjust(my_univariate_results$`Pr(>|t|)`, method = "BH")
selected=subset(my_univariate_results, my_univariate_results$fdr<0.1)

for ( i in 1:nrow(uc_phenos_sig)){
  #Define exposure 
  my_exposure=as.character(uc_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(uc_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=uc_mediator3[,select_mediators$trait, drop=F]
  #If there is only one mediator then we can simply use the mediate function
  my_exposure=as.character(uc_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(uc_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=uc_mediator3[,select_mediators$trait, drop=F]
  if (ncol(my_mediator)==0){
    #my_mess=paste(name_exposure,name_outcome, sep = "->")
    my_mess2="No mediator"
    print (my_mess2)
  }else if (ncol(my_mediator)==1){
    #Get mediator name
    name_mediator=colnames(my_mediator)
    #Create a dataframe containing the elements to test
    my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=uc_outcome3[,my_outcome],Mediator=my_mediator)
    # Give the right names (not necessary but nice for direct network plotting with regmed)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    # First we model the impact of the exposure to the mediator (alpha)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    # Next we model the outcome (dependent var) in relation with the exposure (independent var) while ucling for the mediator
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    #extract impact of mediator to the outcome (beta)
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    # Run mediation analysis 
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    my_mediation_df=extract_mediation_summary(my_mediation)
    
    # Calculate each step of the mediation
    
    # Step 1: total effect
    f0=as.formula(paste( my_outcome,my_exposure, sep = " ~ "))
    fit.total=summary(lm(f0, my_dataset))
    my_total=fit.total$coefficients[2,1]
    my_total_p=fit.total$coefficients[2,4]
    
    #Step 2: Efect of indepedent variable on the mediator 
    my_alpha1=summary(fit.mediator)
    my_alpha=my_alpha1$coefficients[2,1]
    my_alpha_p=my_alpha1$coefficients[2,4]
    
    #Step 3: Effect mediator on the dependent variable (ucling for the independent)
    my_delta1=summary(fit.out)
    my_delta=my_delta1$coefficients[3,1]
    my_delta_p=my_delta1$coefficients[3,4]
    
    sum_fit_report=data.frame(Exposure= my_exposure,
                              Mediator=name_mediator, 
                              Outcome=my_outcome, 
                              total_est=my_total,
                              total_p=my_total_p,
                              impact_on_mediator_est=my_alpha,
                              impact_on_mediator_p=my_alpha_p,
                              Mediated_est=my_delta,
                              Mediated_p=my_delta_p,
                              ACME_estimate=my_mediation_df[1,1], 
                              ACME_p_value=my_mediation_df[1,4],
                              ADE_estimate=my_mediation_df[2,1], 
                              ADE_p_value=my_mediation_df[2,4],
                              Total_effect_estimate=my_mediation_df[3,1], 
                              Total_effect_p_value=my_mediation_df[3,4],
                              Prop_mediated_estimate=my_mediation_df[4,1], 
                              Prop_mediated_p_value=my_mediation_df[4,4]
    )
    sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
    
    
    # # Do the same but now: Exposure -> Metabolite -> Bacteria 
    # # Our previous mediator is now the outcome and viceversa
    
    # #Create a dataframe containing the elements to test
    # #my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=my_mediator,Mediator=uc_outcome3[,my_outcome])
    # # Give the right names (not necessary but nice for direct network plotting with regmed)
    # #colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
    # # First we model the impact of the exposure to the mediator (alpha)
    # #f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
    # #fit.mediator=lm(f1, my_dataset)
    
    # # Next we model the outcome (dependent var) in relation with the exposure (independent var) while ucling for the mediator
    # f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
    # fit.out=lm(f2, my_dataset)
    
    # my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
    # print(my_mess)
    # # Run mediation analysis 
    # my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
    # # Extract direct effect exposure to outcome (delta)
    
    # # Step 1: total effect
    # f0=as.formula(paste(my_exposure, name_mediator, sep = " ~ "))
    # fit.total=summary(lm(f0, my_dataset))
    # my_total=fit.total$coefficients[2,1]
    # my_total_p=fit.total$coefficients[2,4]
    
    # #Step 2: Efect of indepedent variable on the mediator 
    # my_alpha1=summary(fit.mediator)
    # my_alpha=my_alpha1$coefficients[2,1]
    # my_alpha_p=my_alpha1$coefficients[2,4]
    
    # #Step 3: Effect mediator on the dependent variable (ucling for the independent)
    # my_delta1=summary(fit.out)
    # my_delta=my_delta1$coefficients[3,1]
    # my_delta_p=my_delta1$coefficients[3,4]
    
    # my_mediation_df=extract_mediation_summary(my_mediation)
    # sum_fit_report=data.frame(Exposure= my_exposure, 
    #                           Mediator=my_outcome, 
    #                           Outcome=name_mediator, 
    #                           total_est=my_total,
    #                           total_p=my_total_p,
    #                           impact_on_mediator_est=my_alpha,
    #                           impact_on_mediator_p=my_alpha_p,
    #                           Mediated_est=my_delta,
    #                           Mediated_p=my_delta_p,
    #                           ACME_estimate=my_mediation_df[1,1], 
    #                           ACME_p_value=my_mediation_df[1,4],
    #                           ADE_estimate=my_mediation_df[2,1], 
    #                           ADE_p_value=my_mediation_df[2,4],
    #                           Total_effect_estimate=my_mediation_df[3,1], 
    #                           Total_effect_p_value=my_mediation_df[3,4],
    #                           Prop_mediated_estimate=my_mediation_df[4,1], 
    #                           Prop_mediated_p_value=my_mediation_df[4,4])
    
    # sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
  }else {
    # Warning: Standarization set off because was performed before
    my_fit<- regmed.grid(uc_exposure[,my_exposure],my_mediator,uc_outcome3[,my_outcome],lambda.grid, frac.lasso = 0.8, x.std = F, med.std = F,print.iter=T, max.outer=2000)
    # We extract the best model
    fit.best <- regmed.grid.bestfit(my_fit)
    # Extract estimated delta 
    my_delta=fit.best$delta
    # Make an overview of the best lambda for each model
    best_fit=my_fit$grid.data$lambda[my_fit$grid.data$bic==fit.best$bic]
    #my_mess=sprintf("%s -> %s best lambda  %s", my_exposure, my_outcome, best_fit)
    #print(my_mess)
    best_lambdas_uc[i,1]=sprintf("%s -> %s", my_exposure, my_outcome)
    best_lambdas_uc[i,2]=best_fit[1]
    #Summarize results of regmed
    sum_fit=data.frame(a=fit.best$alpha, b=fit.best$beta, x=fit.best$alpha*fit.best$beta, y=my_delta)
    colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta", paste(my_exposure,my_outcome, sep = "->"))
    #results_mediation_uc_=merge(results_mediation_uc_,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    results_mediation_alpha=data.frame(F1=my_exposure, F2=row.names(fit.best$alpha), type="alpha", coeff=as.numeric(fit.best$alpha))
    results_mediation_beta=data.frame(F1=row.names(fit.best$beta), F2=my_outcome, type="beta", coeff=as.numeric(fit.best$beta))
    results_mediation_ab=data.frame(F1=my_exposure, F2=my_outcome, type="alphaxbeta", coeff=as.numeric(fit.best$alpha*fit.best$beta))
    results_mediation_delta=data.frame(F1=my_exposure, F2=my_outcome, type="delta", coeff=fit.best$delta)
    sum_fit3=bind_rows(results_mediation_alpha, results_mediation_beta, results_mediation_ab,results_mediation_delta)
    results_mediation_uc=rbind(results_mediation_uc,sum_fit3)
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
        my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=uc_outcome3[,my_outcome],Mediator=my_mediator[,name_mediator])
        colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
        f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
        fit.mediator=lm(f1, my_dataset)
        
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        
        # Step 1: total effect
        f0=as.formula(paste(my_outcome, my_exposure,sep = " ~ "))
        fit.total=summary(lm(f0, my_dataset))
        my_total=fit.total$coefficients[2,1]
        my_total_p=fit.total$coefficients[2,4]
        
        #Step 2: Efect of indepedent variable on the mediator 
        my_alpha1=summary(fit.mediator)
        my_alpha=my_alpha1$coefficients[2,1]
        my_alpha_p=my_alpha1$coefficients[2,4]
        
        #Step 3: Effect mediator on the dependent variable (ucling for the independent)
        my_delta1=summary(fit.out)
        my_delta=my_delta1$coefficients[3,1]
        my_delta_p=my_delta1$coefficients[3,4]
        
        
        my_mediation_df=extract_mediation_summary(my_mediation)
        sum_fit_report=data.frame(Exposure= my_exposure, 
                                  Mediator=name_mediator, 
                                  Outcome=my_outcome, 
                                  total_est=my_total,
                                  total_p=my_total_p,
                                  impact_on_mediator_est=my_alpha,
                                  impact_on_mediator_p=my_alpha_p,
                                  Mediated_est=my_delta,
                                  Mediated_p=my_delta_p,
                                  ACME_estimate=my_mediation_df[1,1], 
                                  ACME_p_value=my_mediation_df[1,4],
                                  ADE_estimate=my_mediation_df[2,1], 
                                  ADE_p_value=my_mediation_df[2,4],
                                  Total_effect_estimate=my_mediation_df[3,1], 
                                  Total_effect_p_value=my_mediation_df[3,4],
                                  Prop_mediated_estimate=my_mediation_df[4,1], 
                                  Prop_mediated_p_value=my_mediation_df[4,4])
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
       
      #   # Do the same but now: Exposure -> Metabolite -> Bacteria 
      #   # Our previous mediator is now the outcome and viceversa
        
      #   my_dataset=data.frame(Exposure=uc_exposure[,my_exposure], Outcome=my_mediator[,name_mediator],Mediator=uc_outcome3[,my_outcome])
      #   colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
      #   f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
      #   fit.mediator=lm(f1, my_dataset)
        
      #   f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
      #   fit.out=lm(f2, my_dataset)
        
      #   my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
      #   print(my_mess)
      #   my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
        
      #   # Step 1: total effect
      #   f0=as.formula(paste(my_exposure, name_mediator, sep = " ~ "))
      #   fit.total=summary(lm(f0, my_dataset))
      #   my_total=fit.total$coefficients[2,1]
      #   my_total_p=fit.total$coefficients[2,4]
        
      #   #Step 2: Efect of indepedent variable on the mediator 
      #   my_alpha1=summary(fit.mediator)
      #   my_alpha=my_alpha1$coefficients[2,1]
      #   my_alpha_p=my_alpha1$coefficients[2,4]
        
      #   #Step 3: Effect mediator on the dependent variable (ucling for the independent)
      #   my_delta1=summary(fit.out)
      #   my_delta=my_delta1$coefficients[3,1]
      #   my_delta_p=my_delta1$coefficients[3,4]
        
      #   my_mediation_df=extract_mediation_summary(my_mediation)
      #   sum_fit_report=data.frame(Exposure= my_exposure, 
      #                             Mediator=my_outcome,
      #                             Outcome=name_mediator, 
      #                             total_est=my_total,
      #                             total_p=my_total_p,
      #                             impact_on_mediator_est=my_alpha,
      #                             impact_on_mediator_p=my_alpha_p,
      #                             Mediated_est=my_delta,
      #                             Mediated_p=my_delta_p,
      #                             ACME_estimate=my_mediation_df[1,1], 
      #                             ACME_p_value=my_mediation_df[1,4],
      #                             ADE_estimate=my_mediation_df[2,1], 
      #                             ADE_p_value=my_mediation_df[2,4],
      #                             Total_effect_estimate=my_mediation_df[3,1], 
      #                             Total_effect_p_value=my_mediation_df[3,4],
      #                             Prop_mediated_estimate=my_mediation_df[4,1], 
      #                             Prop_mediated_p_value=my_mediation_df[4,4])
      #   sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}

sum_fit_report_final$total_fdr=p.adjust(sum_fit_report_final$total_p, method = "BH")
sum_fit_report_final$impact_on_mediator_fdr=p.adjust(sum_fit_report_final$impact_on_mediator_p, method = "BH")
sum_fit_report_final$Mediated_fdr=p.adjust(sum_fit_report_final$Mediated_p, method = "BH")
sum_fit_report_final$ACME_fdr=p.adjust(sum_fit_report_final$ACME_p_value, method = "BH")
sum_fit_report_final$ADE_fdr=p.adjust(sum_fit_report_final$ADE_p_value, method = "BH")


###########################

#Run mediation in CD

###########################

flag=9
lambda.grid <- c(0.4,0.3,0.2,0.1,0.08,0.04,0.03,0.02,0.01)
best_lambdas_cd=matrix(ncol=2, nrow=nrow(cd_phenos_sig))

output_log=data.frame(Test=NA,Message=NA)
results_mediation_cd=data.frame(F1=NA,F2=NA,type=NA,coeff=NA)

sum_fit_report_final=data.frame(Exposure= NA, 
                                Mediator=NA, 
                                Outcome=NA, 
                                total_est=NA,
                                total_p=NA,
                                impact_on_mediator_est=NA,
                                impact_on_mediator_p=NA,
                                Mediated_est=NA,
                                Mediated_p=NA,
                                ACME_estimate=NA,
                                ACME_p_value=NA, 
                                ADE_estimate=NA, 
                                ADE_p_value=NA, 
                                Total_effect_estimate=NA, 
                                Total_effect_p_value=NA,
                                Prop_mediated_estimate=NA, 
                                Prop_mediated_p_value=NA )


my_univariate_results$fdr=p.adjust(my_univariate_results$`Pr(>|t|)`, method = "BH")
selected=subset(my_univariate_results, my_univariate_results$fdr<0.1)

for ( i in 1:nrow(cd_phenos_sig)){
  #Define exposure 
  my_exposure=as.character(cd_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(cd_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=cd_mediator3[,select_mediators$trait, drop=F]
  #If there is only one mediator then we can simply use the mediate function
  my_exposure=as.character(cd_phenos_sig[i,1])
  #Define outcome
  my_outcome=as.character(cd_phenos_sig[i,2])
  #Define set of mediators that are significantly associated with the exposure: if there is no relation between exposure and the mediator is unlikely that there's a mediation effect
  select_mediators=subset(selected,selected$phenotype==my_exposure)
  my_mediator=cd_mediator3[,select_mediators$trait, drop=F]
  if (ncol(my_mediator)==0){
    #my_mess=paste(name_exposure,name_outcome, sep = "->")
    my_mess2="No mediator"
    print (my_mess2)
  }else if (ncol(my_mediator)==1){
    #Get mediator name
    name_mediator=colnames(my_mediator)
    #Create a dataframe containing the elements to test
    my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=cd_outcome3[,my_outcome],Mediator=my_mediator)
    # Give the right names (not necessary but nice for direct network plotting with regmed)
    colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
    # First we model the impact of the exposure to the mediator (alpha)
    f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
    fit.mediator=lm(f1, my_dataset)
    # Next we model the outcome (dependent var) in relation with the exposure (independent var) while cdling for the mediator
    f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
    fit.out=lm(f2, my_dataset)
    #extract impact of mediator to the outcome (beta)
    my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
    print(my_mess)
    # Run mediation analysis 
    my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
    my_mediation_df=extract_mediation_summary(my_mediation)
    
    # Calculate each step of the mediation
    
    # Step 1: total effect
    f0=as.formula(paste(my_exposure, my_outcome, sep = " ~ "))
    fit.total=summary(lm(f0, my_dataset))
    my_total=fit.total$coefficients[2,1]
    my_total_p=fit.total$coefficients[2,4]
    
    #Step 2: Efect of indepedent variable on the mediator 
    my_alpha1=summary(fit.mediator)
    my_alpha=my_alpha1$coefficients[2,1]
    my_alpha_p=my_alpha1$coefficients[2,4]
    
    #Step 3: Effect mediator on the dependent variable (cdling for the independent)
    my_delta1=summary(fit.out)
    my_delta=my_delta1$coefficients[3,1]
    my_delta_p=my_delta1$coefficients[3,4]
    
    sum_fit_report=data.frame(Exposure= my_exposure,
                              Mediator=name_mediator, 
                              Outcome=my_outcome, 
                              total_est=my_total,
                              total_p=my_total_p,
                              impact_on_mediator_est=my_alpha,
                              impact_on_mediator_p=my_alpha_p,
                              Mediated_est=my_delta,
                              Mediated_p=my_delta_p,
                              ACME_estimate=my_mediation_df[1,1], 
                              ACME_p_value=my_mediation_df[1,4],
                              ADE_estimate=my_mediation_df[2,1], 
                              ADE_p_value=my_mediation_df[2,4],
                              Total_effect_estimate=my_mediation_df[3,1], 
                              Total_effect_p_value=my_mediation_df[3,4],
                              Prop_mediated_estimate=my_mediation_df[4,1], 
                              Prop_mediated_p_value=my_mediation_df[4,4]
    )
    sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
    
    
    # Do the same but now: Exposure -> Metabolite -> Bacteria 
    # Our previous mediator is now the outcome and viceversa
    
    # #Create a dataframe containing the elements to test
    # my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=my_mediator,Mediator=cd_outcome3[,my_outcome])
    # # Give the right names (not necessary but nice for direct network plotting with regmed)
    # colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
    # # First we model the impact of the exposure to the mediator (alpha)
    # f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
    # fit.mediator=lm(f1, my_dataset)
    
    # # Next we model the outcome (dependent var) in relation with the exposure (independent var) while cdling for the mediator
    # f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
    # fit.out=lm(f2, my_dataset)
    
    # my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
    # print(my_mess)
    # # Run mediation analysis 
    # my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
    # # Extract direct effect exposure to outcome (delta)
    
    # # Step 1: total effect
    # f0=as.formula(paste(my_exposure, name_mediator, sep = " ~ "))
    # fit.total=summary(lm(f0, my_dataset))
    # my_total=fit.total$coefficients[2,1]
    # my_total_p=fit.total$coefficients[2,4]
    
    # #Step 2: Efect of indepedent variable on the mediator 
    # my_alpha1=summary(fit.mediator)
    # my_alpha=my_alpha1$coefficients[2,1]
    # my_alpha_p=my_alpha1$coefficients[2,4]
    
    # #Step 3: Effect mediator on the dependent variable (cdling for the independent)
    # my_delta1=summary(fit.out)
    # my_delta=my_delta1$coefficients[3,1]
    # my_delta_p=my_delta1$coefficients[3,4]
    
    # my_mediation_df=extract_mediation_summary(my_mediation)
    # sum_fit_report=data.frame(Exposure= my_exposure, 
    #                           Mediator=my_outcome, 
    #                           Outcome=name_mediator, 
    #                           total_est=my_total,
    #                           total_p=my_total_p,
    #                           impact_on_mediator_est=my_alpha,
    #                           impact_on_mediator_p=my_alpha_p,
    #                           Mediated_est=my_delta,
    #                           Mediated_p=my_delta_p,
    #                           ACME_estimate=my_mediation_df[1,1], 
    #                           ACME_p_value=my_mediation_df[1,4],
    #                           ADE_estimate=my_mediation_df[2,1], 
    #                           ADE_p_value=my_mediation_df[2,4],
    #                           Total_effect_estimate=my_mediation_df[3,1], 
    #                           Total_effect_p_value=my_mediation_df[3,4],
    #                           Prop_mediated_estimate=my_mediation_df[4,1], 
    #                           Prop_mediated_p_value=my_mediation_df[4,4])
    
    # sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
  }else {
    # Warning: Standarization set off because was performed before
    my_fit<- regmed.grid(cd_exposure[,my_exposure],my_mediator,cd_outcome3[,my_outcome],lambda.grid, frac.lasso = 0.8, x.std = F, med.std = F,print.iter=T, max.outer=2000)
    # We extract the best model
    fit.best <- regmed.grid.bestfit(my_fit)
    # Extract estimated delta 
    my_delta=fit.best$delta
    # Make an overview of the best lambda for each model
    best_fit=my_fit$grid.data$lambda[my_fit$grid.data$bic==fit.best$bic]
    #my_mess=sprintf("%s -> %s best lambda  %s", my_exposure, my_outcome, best_fit)
    #print(my_mess)
    best_lambdas_cd[i,1]=sprintf("%s -> %s", my_exposure, my_outcome)
    best_lambdas_cd[i,2]=best_fit[1]
    #Summarize results of regmed
    sum_fit=data.frame(a=fit.best$alpha, b=fit.best$beta, x=fit.best$alpha*fit.best$beta, y=my_delta)
    colnames(sum_fit)=c(paste(my_exposure,"species", sep = "->"), paste("species", my_outcome,sep = "->"), "alphaxbeta", paste(my_exposure,my_outcome, sep = "->"))
    #results_mediation_cd_=merge(results_mediation_cd_,sum_fit, all.x =T, by.x="Species", by.y = "row.names")
    results_mediation_alpha=data.frame(F1=my_exposure, F2=row.names(fit.best$alpha), type="alpha", coeff=as.numeric(fit.best$alpha))
    results_mediation_beta=data.frame(F1=row.names(fit.best$beta), F2=my_outcome, type="beta", coeff=as.numeric(fit.best$beta))
    results_mediation_ab=data.frame(F1=my_exposure, F2=my_outcome, type="alphaxbeta", coeff=as.numeric(fit.best$alpha*fit.best$beta))
    results_mediation_delta=data.frame(F1=my_exposure, F2=my_outcome, type="delta", coeff=fit.best$delta)
    sum_fit3=bind_rows(results_mediation_alpha, results_mediation_beta, results_mediation_ab,results_mediation_delta)
    results_mediation_cd=rbind(results_mediation_cd,sum_fit3)
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
        my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=cd_outcome3[,my_outcome],Mediator=my_mediator[,name_mediator])
        colnames(my_dataset)=c(my_exposure,my_outcome,name_mediator)
        f1=as.formula(paste(name_mediator, my_exposure, sep = " ~ "))
        fit.mediator=lm(f1, my_dataset)
        
        f2=as.formula(paste(my_outcome, paste(my_exposure,name_mediator, sep = "+"), sep = " ~ "))
        fit.out=lm(f2, my_dataset)
        
        my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure,name_mediator, my_outcome)
        print(my_mess)
        my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=name_mediator, boot=T))
        
        # Step 1: total effect
        f0=as.formula(paste(my_outcome,my_exposure, sep = " ~ "))
        fit.total=summary(lm(f0, my_dataset))
        my_total=fit.total$coefficients[2,1]
        my_total_p=fit.total$coefficients[2,4]
        
        #Step 2: Efect of indepedent variable on the mediator 
        my_alpha1=summary(fit.mediator)
        my_alpha=my_alpha1$coefficients[2,1]
        my_alpha_p=my_alpha1$coefficients[2,4]
        
        #Step 3: Effect mediator on the dependent variable (cdling for the independent)
        my_delta1=summary(fit.out)
        my_delta=my_delta1$coefficients[3,1]
        my_delta_p=my_delta1$coefficients[3,4]
        
        
        my_mediation_df=extract_mediation_summary(my_mediation)
        sum_fit_report=data.frame(Exposure= my_exposure, 
                                  Mediator=name_mediator, 
                                  Outcome=my_outcome, 
                                  total_est=my_total,
                                  total_p=my_total_p,
                                  impact_on_mediator_est=my_alpha,
                                  impact_on_mediator_p=my_alpha_p,
                                  Mediated_est=my_delta,
                                  Mediated_p=my_delta_p,
                                  ACME_estimate=my_mediation_df[1,1], 
                                  ACME_p_value=my_mediation_df[1,4],
                                  ADE_estimate=my_mediation_df[2,1], 
                                  ADE_p_value=my_mediation_df[2,4],
                                  Total_effect_estimate=my_mediation_df[3,1], 
                                  Total_effect_p_value=my_mediation_df[3,4],
                                  Prop_mediated_estimate=my_mediation_df[4,1], 
                                  Prop_mediated_p_value=my_mediation_df[4,4])
        sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
        # # Do the same but now: Exposure -> Metabolite -> Bacteria 
        # # Our previous mediator is now the outcome and viceversa
        
        # my_dataset=data.frame(Exposure=cd_exposure[,my_exposure], Outcome=my_mediator[,name_mediator],Mediator=cd_outcome3[,my_outcome])
        # colnames(my_dataset)=c(my_exposure,name_mediator,my_outcome)
        # f1=as.formula(paste(my_outcome, my_exposure, sep = " ~ "))
        # fit.mediator=lm(f1, my_dataset)
        
        # f2=as.formula(paste(name_mediator, paste(my_exposure,my_outcome, sep = "+"), sep = " ~ "))
        # fit.out=lm(f2, my_dataset)
        
        # my_mess=sprintf("Testing 1 mediator: %s -> %s -> %s", my_exposure, my_outcome,name_mediator)
        # print(my_mess)
        # my_mediation=summary(mediate(fit.mediator, fit.out, treat=my_exposure, mediator=my_outcome, boot=T))
        
        # # Step 1: total effect
        # f0=as.formula(paste( name_mediator,my_exposure, sep = " ~ "))
        # fit.total=summary(lm(f0, my_dataset))
        # my_total=fit.total$coefficients[2,1]
        # my_total_p=fit.total$coefficients[2,4]
        
        # #Step 2: Efect of indepedent variable on the mediator 
        # my_alpha1=summary(fit.mediator)
        # my_alpha=my_alpha1$coefficients[2,1]
        # my_alpha_p=my_alpha1$coefficients[2,4]
        
        # #Step 3: Effect mediator on the dependent variable (cdling for the independent)
        # my_delta1=summary(fit.out)
        # my_delta=my_delta1$coefficients[3,1]
        # my_delta_p=my_delta1$coefficients[3,4]
        
        # my_mediation_df=extract_mediation_summary(my_mediation)
        # sum_fit_report=data.frame(Exposure= my_exposure, 
        #                           Mediator=my_outcome,
        #                           Outcome=name_mediator, 
        #                           total_est=my_total,
        #                           total_p=my_total_p,
        #                           impact_on_mediator_est=my_alpha,
        #                           impact_on_mediator_p=my_alpha_p,
        #                           Mediated_est=my_delta,
        #                           Mediated_p=my_delta_p,
        #                           ACME_estimate=my_mediation_df[1,1], 
        #                           ACME_p_value=my_mediation_df[1,4],
        #                           ADE_estimate=my_mediation_df[2,1], 
        #                           ADE_p_value=my_mediation_df[2,4],
        #                           Total_effect_estimate=my_mediation_df[3,1], 
        #                           Total_effect_p_value=my_mediation_df[3,4],
        #                           Prop_mediated_estimate=my_mediation_df[4,1], 
        #                           Prop_mediated_p_value=my_mediation_df[4,4])
        # sum_fit_report_final=rbind(sum_fit_report_final, sum_fit_report)
      }
    }
  }
}

sum_fit_report_final$total_fdr=p.adjust(sum_fit_report_final$total_p, method = "BH")
sum_fit_report_final$impact_on_mediator_fdr=p.adjust(sum_fit_report_final$impact_on_mediator_p, method = "BH")
sum_fit_report_final$Mediated_fdr=p.adjust(sum_fit_report_final$Mediated_p, method = "BH")
sum_fit_report_final$ACME_fdr=p.adjust(sum_fit_report_final$ACME_p_value, method = "BH")
sum_fit_report_final$ADE_fdr=p.adjust(sum_fit_report_final$ADE_p_value, method = "BH")

