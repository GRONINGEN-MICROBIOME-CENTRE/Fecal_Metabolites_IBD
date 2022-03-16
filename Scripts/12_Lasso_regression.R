library(tidyverse)
library (dplyr)
library(caret)
library(glmnet)
library(reshape2)
library(ggplot2)
#library(glmnetUtils, lib.loc="../Rlib/")
library(glmnetUtils)
#withr::with_libpaths(new="../Rlib/", install_github("hongooi73/glmnetUtils"))

library(doParallel)
registerDoParallel(5)


#Import data
setwd("~/Desktop/IBD_metabolomics_2022/1.Input/lasso/")
phenos <- read.delim("./CC_phenos_recoded_v2.txt")
phenos$host_SmokeCurrentSmoker=NULL

#Replace missing for median values
phenos$host_BMI[is.na(phenos$host_BMI)]=median(phenos$host_BMI, na.rm = T)
phenos$clinical_ChrA[is.na(phenos$clinical_ChrA)]=median(phenos$clinical_ChrA, na.rm = T)
phenos$clinical_HBD2[is.na(phenos$clinical_HBD2)]=median(phenos$clinical_HBD2, na.rm = T)
phenos$clinical_BowelMovementADayDef[is.na(phenos$clinical_BowelMovementADayDef)]=median(phenos$clinical_BowelMovementADayDef, na.rm = T)

mtb <- read.delim("./mtb_scfa.txt", sep="\t")
mgc <- read.delim("./mgc.txt",sep="\t")
taxa <- read.delim("./taxa.txt",sep="\t")
colnames(taxa)=paste("sp_",colnames(taxa), sep ="")

d1=merge(phenos,taxa, by.x="ID", by.y="row.names")
d2=merge(d1,mgc, by.x="ID", by.y="row.names")
d3=merge(d2,mtb, by.x="ID", by.y="row.names")
row.names(d3)=d3$ID
d3$ID=NULL
d3$sp_Shannon_Index=NULL

list_metabolites=colnames(mtb)

Metabolite_iteration_v3(Input = d3, Summary= list_metabolites) -> Output_model

var_expl=as.data.frame(Output_model[[1]])

var_expl_selected=as.data.frame(Output_model[[2]])
var_expl_selected_filt=subset(var_expl_selected,var_expl_selected$Model=="Complete")
var_expl$Metabolite=NULL
colnames(var_expl)=c("0_Null", "6_Complete", "4_Taxa", "1_Diet", "2_Medication", "3_IBD", "5_Metabolic_gene_cluster", "2.1_Biomarkers","Metabolite")
ve_plot=melt(var_expl)
ve_plot(variable)=as.character(ve_plot(variable))

#Plot
ggplot(ve_plot, aes(variable, value, fill=variable))  + theme_bw() + geom_jitter(shape=21, size=2) + geom_boxplot(outlier.alpha = 0) + scale_fill_manual(values = c("grey80","yellow" ,"lightblue","cyan","blue1","salmon","red3", "purple"))

var_expl_selected_filt=subset(var_expl_selected,var_expl_selected$Model=="Complete")
var_expl_selected_filt2=melt(var_expl_selected_filt,id.vars = c("Model", "Metabolite"))
var_expl_selected_filt2$value=as.numeric(as.character(var_expl_selected_filt2$value))
var_expl_selected_filt2=subset(var_expl_selected_filt2, var_expl_selected_filt2$value!=0)

var_expl_sd=as.data.frame(Output_model[[3]])

var_expl_mean=as.data.frame(Output_model[[4]])

 var_expl_mean2=var_expl_mean
 var_expl_mean2[var_expl_mean2<0]=0
 write.table(var_expl, "~/Desktop/IBD_metabolomics_2022/lasso_var_expl.txt", sep = "\t", quote = F)
 write.table(var_expl_selected_filt2, "~/Desktop/IBD_metabolomics_2022/lasso_var_expl_model_complete.txt", sep = "\t", quote = F)
 write.table(var_expl_mean2, "~/Desktop/IBD_metabolomics_2022/lasso_var_expl_mean.txt", sep = "\t", quote = F)
 write.table(var_expl_sd, "~/Desktop/IBD_metabolomics_2022/lasso_var_expl_sd.txt", sep = "\t", quote = F)


Fit_lasso_cv = function(Dependent, Regressors,my_folds=5){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors2
  if (dim(Regressors)[[2]] < 2){
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    Param_best = tibble(Variable = colnames(Regressors), Beta= Beta)
  } else{
    #Create folds
    set.seed(2021)
    row_folds=createFolds(1:nrow(Regressors2),k=my_folds)
    Best_error = NULL
    my_var_iter = numeric(0)
    for (i in names(row_folds)){
      
      #print (paste("Running fold" ,i))
      #Split training / testing
      select_rows=row_folds[[i]]
      test_reg=Regressors2[select_rows,]
      train_reg=Regressors2[-select_rows,]
      #Lasso with 10 fold 
      cvfit=cv.glmnet(Dependent ~ ., train_reg, alpha = 1, nfolds = 10, type.measure="mse",standardize=T, parallel=TRUE)
      #Predict in the test dataset and estimate the R2
      lasso_predicted = predict(cvfit, s=cvfit$lambda.min, newdata=data.matrix(test_reg[,-ncol(test_reg)]))
      tmp_cv=data.frame(real=test_reg$Dependent,pred=lasso_predicted)
      colnames(tmp_cv)=c("real", "pred")
      RSSM = sum((tmp_cv$real - tmp_cv$pred)^2)/dim(tmp_cv)[1]
      if (length(Best_error) == 0){ 
          Best_error = RSSM
          Best_model = cvfit
      }else if  (RSSM < Best_error){
        Best_error = RSSM
        Best_model = cvfit
      }
      #Keep r2 per iteration
      #lasso_pred_iter = predict(cvfit, s=cvfit$lambda.min, newdata=test_reg$Dependent)
      Exp_var_iter = R2_calc(test_reg$Dependent, lasso_predicted)
      my_var_iter=c(my_var_iter,Exp_var_iter)
    }
    cvfit = Best_model
    Param_best <- coef(cvfit, s = "lambda.min")
    as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
    colnames(Param_best) = c("Variable","Beta")
    lasso_predicted = predict(cvfit, s=cvfit$lambda.min, newdata=data.matrix(Regressors2))
    Explained_variance = R2_calc(Regressors2$Dependent, lasso_predicted)
    #Calculate the mean R2 and its SD per iteration
    Mean_var_iter=mean(my_var_iter)
    sd_var_iter=sd(my_var_iter)
    if (Explained_variance < 0){ Explained_variance = 0 }
    #print(Explained_variance)
  }
  return(list(Param_best, Explained_variance, Mean_var_iter, sd_var_iter))
}


R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}



Metabolite_iteration_v3 = function(Input, Summary, my_folds){

  ##there are some duplicated records with a _1 at the end, remove them from Input and from summary (in summary instead of _ ther is a .)
  #colnames(Input)[grepl("[a-z]_1$",colnames(Input))] -> Remove
  #Input %>% select(-Remove) -> Input ; Summary %>% select(-one_of(str_replace(Remove,"_","."))) -> Summary
  #List of all metabolites to iterate 
  Metabolites <- Summary
  #Divide phenotypes in categories
  Microbes <- colnames(Input)[grepl("sp_",colnames(Input))]
  Diet <- colnames(Input)[grepl("diet_",colnames(Input))]
  Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex", "Amount_sample_gram") #Possibly include others?
  Medication <- colnames(Input)[grepl("med_",colnames(Input))]
  IBD <- c("IBD")
  genetics=colnames(Input)[grepl("Entryname",colnames(Input))]
  biomarkers=c("clinical_HBD2", "clinical_ChrA", "clinical_Calprot200")
  
  #Output dataframe
  Variability_explained = tibble()
  Variability_explained_mean = tibble()
  Variability_explained_sd = tibble()
  #Name of the models
  Name_models <- c("Null","Complete", "Microbes", "Diet", "Medication", "Disease", "Genetics","Biomarkers")
  #Make all variables in numberic/character, currently the function does not accept highly multifactorial variables. 2 level factors become numeric.
  #select(Input, -one_of(c("Row.names", "run_day_cat"))) -> Variables_Input
  #apply(Variables_Input,2, FUN=Make_numeric) %>% as_tibble() -> Variables_Input #Check Make_numeric function
  #This are the column names that are used to save the Betas, so if you need the beta of a specific varaible, should be in Variables_col vector
  Input %>% select(c(Microbes, Diet, Covariates, Medication, IBD, genetics, biomarkers)) -> Variables_col 
  colnames(Variables_col) -> Variables_col
  hmm=0
  All_model_info ={}
    for (Metabolite in Metabolites){
    #foreach(i=1:length(Metabolites), .combine=rbind) %dopar% {
    #Metabolite=Metabolites[i]
    Logit = F
    hmm=hmm+1
    tested=paste(hmm,length(Metabolites), sep="/")
    #Get the metabolite of interest and their associations 
    #Variables_col  -> Selected_phenotypes
    Selected_phenotypes=c(Microbes, Diet, Medication, IBD, genetics, biomarkers)
    print(paste(tested,Metabolite, sep=" "))
    #From the Input after transforming it to numeric select the Metabolite (dependent), phenotypes assocaited and Covariates. Remove all records with NA.
    Input %>% select(one_of(c(Metabolite, Selected_phenotypes,Covariates ))) %>% drop_na() -> Input_model
    #Make a vector out of dependent
    Dependent <- as.vector(as_vector(Input_model[,1]))
    #If dependent is a character, then do logistic
    if (class(Dependent[0]) == "character"){ Logit = T}
    #Prepare the different inputs for each model
    Variables_complete <- Input_model[,2:dim(Input_model)[2]]
    Variables_IBD  <- select(Input_model, one_of(c(Covariates, IBD)))
    Variables_microbiome <- select(Input_model, one_of(c(Microbes, Covariates)))
    Variables_diet <- select(Input_model, one_of(c(Diet, Covariates)))
    Variables_medication <- select(Input_model, one_of(c(Medication, Covariates)))
    Variables_null <- select(Input_model, one_of(Covariates))
    Variables_genetics <- select(Input_model, one_of(c(genetics, Covariates)))
    Variables_biomarkers <- select(Input_model, one_of(c(biomarkers, Covariates)))
    
    Models <- list( Variables_null, Variables_complete,Variables_microbiome, Variables_diet,Variables_medication, Variables_IBD, Variables_genetics, Variables_biomarkers)
    #Vector of R2s for output 1
    Variability_model = c()
    Variability_mean = c()
    Variability_sd = c()
    #Data.frame of variables for output 2
    tibble(Variable = Variables_col) -> Variables
    #For each model, fit lasso (normal or logistic) and save R2 and variables
    for (N in seq(1:length(Name_models)) ){
      Input_model <- Models[[N]] ; Name <- Name_models[[N]]
      #If 0 significant features add as 0 all beta and R2 and go to next
      if (dim(Input_model)[2] < 1){ 
          Variability_model = c(Variability_model,0)
          Model_summary = tibble(Variable=Variables, Beta=NA) %>% t() %>% as_tibble() %>% `colnames<-`(Metabolites)
          Model_summary %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
          next 
        }

      if (Logit == F){ Fit_lasso_cv(Dependent = Dependent, Regressors = Input_model) -> Lasso_results
      }else{ Fit_logistic_lasso(Dependent = Dependent, Regressors = Input_model) -> Lasso_results }
      #Add the betas to the Data.frame of features (only features included in that data.frame are going to get the Beta saved)   
      left_join(Variables,Lasso_results[[1]],by = "Variable") %>% t() %>% as_tibble() %>% `colnames<-`(Variables$Variable) -> Model_summary
      Model_summary[2,] %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
      All_model_info = rbind(All_model_info, Model_summary)
      #Save R2
      Variability_model = c(Variability_model,as.numeric(Lasso_results[[2]]))
      #Mean R2 per iteration
      Variability_mean = c(Variability_mean,as.numeric(Lasso_results[[3]]))
      #SD R2 per iteration
      Variability_sd = c(Variability_sd,as.numeric(Lasso_results[[4]]))
    }
    #Make R2s into a data.frame with each model per column
    as_tibble(matrix(Variability_model,nrow = 1,ncol = 9)) %>% mutate(V10 = Metabolite) -> Variability_model
    colnames(Variability_model) = c(Name_models, "Metabolite")
    rbind(Variability_explained, Variability_model) -> Variability_explained


    as_tibble(matrix(Variability_mean,nrow = 1,ncol = 9)) %>% mutate(V10 = Metabolite) -> Variability_mean
    colnames(Variability_mean) = c(Name_models, "Metabolite")
    rbind(Variability_explained_mean, Variability_mean) -> Variability_explained_mean


    as_tibble(matrix(Variability_sd,nrow = 1,ncol = 9)) %>% mutate(V10 = Metabolite) -> Variability_sd
    colnames(Variability_sd) = c(Name_models, "Metabolite")
    rbind(Variability_explained_sd, Variability_sd) -> Variability_explained_sd
    #Variability_explained

  }
  return(list(Variability_explained, All_model_info, Variability_explained_sd, Variability_explained_mean))
}





