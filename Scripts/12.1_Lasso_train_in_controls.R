library(tidyverse)
library (dplyr)
library(caret)
library(glmnet)
#library(glmnetUtils, lib.loc="../Rlib/")
library(glmnetUtils)
#withr::with_libpaths(new="../Rlib/", install_github("hongooi73/glmnetUtils"))

library(doParallel)
registerDoParallel(5)


#Import data

phenos <- read.delim("~/Desktop/IBD_Metabolomics/7.Intermediate_files/CC_phenos_recoded_v2.txt")
phenos$host_SmokeCurrentSmoker=NULL

#Replace missing for median values
phenos$host_BMI[is.na(phenos$host_BMI)]=median(phenos$host_BMI, na.rm = T)
phenos$clinical_ChrA[is.na(phenos$clinical_ChrA)]=median(phenos$clinical_ChrA, na.rm = T)
phenos$clinical_HBD2[is.na(phenos$clinical_HBD2)]=median(phenos$clinical_HBD2, na.rm = T)
phenos$clinical_BowelMovementADayDef[is.na(phenos$clinical_BowelMovementADayDef)]=median(phenos$clinical_BowelMovementADayDef, na.rm = T)

mtb <- read.delim("~/Desktop/IBD_Metabolomics/7.Intermediate_files/mtb.txt", sep=" ")
mgc <- read.delim("~/Desktop/IBD_Metabolomics/7.Intermediate_files/mgc.txt",sep=" ")
taxa <- read.delim("~/Desktop/IBD_Metabolomics/7.Intermediate_files/taxa.txt",sep=" ")

d1=merge(phenos,taxa, by.x="ID", by.y="row.names")
d2=merge(d1,mgc, by.x="ID", by.y="row.names")
d3=merge(d2,mtb, by.x="ID", by.y="row.names")
row.names(d3)=d3$ID
d3$ID=NULL
d3$sp_Shannon_Index=NULL

list_metabolites=colnames(mtb)
#list_metabolites=colnames(mtb)[1:10]

Metabolite_iteration_v3(Input = d3, Summary= list_metabolites) -> Output_model


Fit_lasso_cv = function(Dependent, Regressors,my_folds=5){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors_all
  Regressors2=subset(Regressors_all, Regressors_all$IBD==1)
  Regressors2$IBD=NULL
  Regressors3=subset(Regressors_all, Regressors_all$IBD==2)
  Regressors3$IBD=NULL
  if (dim(Regressors)[[2]] < 2){
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    Param_best = tibble(Variable = colnames(Regressors), Beta= Beta)
  } else{
    #Create folds
    set.seed(2021)
    cvfit=cv.glmnet(Dependent ~ ., Regressors2, alpha = 1, nfolds = 10, type.measure="mse",standardize=T, parallel=TRUE)
    Param_best <- coef(cvfit, s = "lambda.min")
    as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
    colnames(Param_best) = c("Variable","Beta")
    lasso_predicted = predict(cvfit, s=cvfit$lambda.min, newdata=data.matrix(Regressors3))
    Explained_variance = R2_calc(Regressors3$Dependent, lasso_predicted)
    if (Explained_variance < 0){ Explained_variance = 0 }
  }

  return(list(Param_best, Explained_variance))
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
  Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex", "Amount_sample_gram","IBD") #Possibly include others?
  Medication <- colnames(Input)[grepl("med_",colnames(Input))]
  #IBD <- c("IBD")
  genetics=colnames(Input)[grepl("Entryname",colnames(Input))]
  biomarkers=c("clinical_HBD2", "clinical_ChrA", "clinical_Calprot200")
  
  #Output dataframe
  Variability_explained = tibble()
  Variability_explained_mean = tibble()
  Variability_explained_sd = tibble()
  #Name of the models
  Name_models <- c("Null","Complete", "Microbes", "Diet", "Medication", "Genetics","Biomarkers")
  #Make all variables in numberic/character, currently the function does not accept highly multifactorial variables. 2 level factors become numeric.
  #select(Input, -one_of(c("Row.names", "run_day_cat"))) -> Variables_Input
  #apply(Variables_Input,2, FUN=Make_numeric) %>% as_tibble() -> Variables_Input #Check Make_numeric function
  #This are the column names that are used to save the Betas, so if you need the beta of a specific varaible, should be in Variables_col vector
  Input %>% select(c(Microbes, Diet, Covariates, Medication, genetics, biomarkers)) -> Variables_col 
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
    Selected_phenotypes=c(Microbes, Diet, Medication, genetics, biomarkers)
    print(paste(tested,Metabolite, sep=" "))
    #From the Input after transforming it to numeric select the Metabolite (dependent), phenotypes assocaited and Covariates. Remove all records with NA.
    Input %>% select(one_of(c(Metabolite, Selected_phenotypes,Covariates ))) %>% drop_na() -> Input_model
    #Make a vector out of dependent
    Dependent <- as.vector(as_vector(Input_model[,1]))
    #If dependent is a character, then do logistic
    if (class(Dependent[0]) == "character"){ Logit = T}
    #Prepare the different inputs for each model
    Variables_complete <- Input_model[,2:dim(Input_model)[2]]
    #Variables_IBD  <- select(Input_model, one_of(c(Covariates, IBD)))
    Variables_microbiome <- select(Input_model, one_of(c(Microbes, Covariates)))
    Variables_diet <- select(Input_model, one_of(c(Diet, Covariates)))
    Variables_medication <- select(Input_model, one_of(c(Medication, Covariates)))
    Variables_null <- select(Input_model, one_of(Covariates))
    Variables_genetics <- select(Input_model, one_of(c(genetics, Covariates)))
    Variables_biomarkers <- select(Input_model, one_of(c(biomarkers, Covariates)))
    
    Models <- list( Variables_null, Variables_complete,Variables_microbiome, Variables_diet,Variables_medication, Variables_genetics, Variables_biomarkers)
    #Vector of R2s for output 1
    Variability_model = c()
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
    }
    #Make R2s into a data.frame with each model per column
    as_tibble(matrix(Variability_model,nrow = 1,ncol = 8)) %>% mutate(V9 = Metabolite) -> Variability_model
    colnames(Variability_model) = c(Name_models, "Metabolite")
    rbind(Variability_explained, Variability_model) -> Variability_explained

  }
  return(list(Variability_explained, All_model_info))
}





