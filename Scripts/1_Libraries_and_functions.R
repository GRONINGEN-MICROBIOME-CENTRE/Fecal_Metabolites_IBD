Load libraries

library(corrplot)
library(plotly)
library(pheatmap)
library(NbClust)
library(reshape2)
library(ggrepel)
library(dplyr)
library(eulerr)
library(vegan)
library(ape)
library(ggridges)
library (UpSetR)
library(ggplot2)
library(psych)
library(ggord)
library(data.table)
library(mice)
library(caret)
library(factoextra)
library(Rtsne)
library(heatmaply)
library(ggpubr)
library(ggbiplot)
library(coin)
library(DT)
library(pvclust)
library(glmnet)
library(compositions)
library(RColorBrewer)
library(iClusterPlus)
library (ggvegan)
library (VIM)
library(factoextra)


Create functions


####################################
#Inverse rank transformation function.
####################################

invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

###############################
#negate in (to exclude columns)
###############################

`%ni%` <- Negate(`%in%`)

#####################################
#Imputation and filtering metabolites
#####################################

impute_missing_values=function(x, samples_row=T, method="knn", missing_filter=0){
  #if samples are in columns transpose
  if (!samples_row){
    x=as.data.frame(t(x))
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>1){
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 1") 
  }
  col_ori=ncol(x)
  col_keep=colnames(x)[colSums(!is.na(x))/nrow(x)>missing_filter]
  x=x[,col_keep]
  my_num_removed=col_ori-ncol(x)
  warning (paste(my_num_removed, "columns removed due to many missing values"))
  if (method=="zero"){
    x[is.na(x)]=0
  }
  else if (method=="min"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), min(a, na.rm = TRUE), a))
  }
  else if (method=="hfmin"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), min(a, na.rm = TRUE)/2, a))
  }
  else if (method=="mean"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), mean(a, na.rm = TRUE), a))
  }
  else if (method=="median"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), median(a, na.rm = TRUE), a))
  }
  else if (method=="knn"){
    x=kNN(x,k=10,metric = "euclidean", imp_var = F)
  }
  else if (method=="none"){
    x=x
  }
  else{
    stop("\n Hey! \n Method parameters for imputing missing values per column are: \n min -> min (col), hfmin -> min(col)/2, mean -> mean (col), median -> median(col), mice -> using mice package [Multivariate Imputation By Chained Equations] defaults \n")
  }
  return(as.data.frame(x))
}


##############################
# Summary stats
##############################

summary_stats=function(feature_matrix,pheno, missing_vals=0){
  a=1
  num_var=length(levels(pheno[,1]))
  print (paste("Number of phenotypes detected for summary statistics:", num_var))
  summary_met=matrix(nrow=ncol(feature_matrix),ncol=(9*(num_var+1))+1)
  for (i in 1:ncol(feature_matrix)){
    #print (colnames(feature_matrix)[i])
    if (missing_vals==0){
      pos=2
      summary_met[i,pos]=sum(feature_matrix[,i]!=0)
      pos=pos+1
      summary_met[i,pos]=sum(feature_matrix[,i]==0)
    }
    else {  
      pos=2
      summary_met[i,pos]=sum(is.na(feature_matrix[,i]))
      pos=pos+1
      summary_met[i,pos]=sum(!is.na(feature_matrix[,i]))
    }
    if (sum(is.na(feature_matrix[,i]) == nrow(feature_matrix))){
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
      pos=pos+1
      summary_met[i,pos]="NA"
    }
    else{
      pos=pos+1
      summary_met[i,pos]=min(feature_matrix[,i], na.rm = T)
      pos=pos+1
      summary_met[i,pos]=max(feature_matrix[,i], na.rm = T)
      pos=pos+1
      summary_met[i,pos]=mean(feature_matrix[,i], na.rm = T)
      pos=pos+1
      summary_met[i,pos]=median(feature_matrix[,i], na.rm = T)
      pos=pos+1
      summary_met[i,pos]=sd(feature_matrix[,i], na.rm = T)
      pos=pos+1
      summary_met[i,pos]=length(boxplot(feature_matrix[,i], plot=FALSE, range = 3)$out)
      if (sum(!is.na(feature_matrix[,i])) < 10){
        pos=pos+1
        summary_met[i,pos]="NA"
        } else {
        nor=shapiro.test(feature_matrix[,i])
        pos=pos+1
        summary_met[i,pos]=nor$p.value>=0.05
      }
    }
    for (x in levels (pheno[,1])){
      id_list=rownames(pheno)[pheno[,1]==x]
      tmp_all=feature_matrix[rownames(feature_matrix)%in%id_list,]
      pos=pos+1
      summary_met[i,pos]=sum(is.na(tmp_all[,i]))
      pos=pos+1
      summary_met[i,pos]=sum(!is.na(tmp_all[,i]))
      if (sum(is.na(tmp_all[,i]) == nrow(tmp_all))){
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
        pos=pos+1
        summary_met[i,pos]="NA"
      }else{
        pos=pos+1
        summary_met[i,pos]=min(tmp_all[,i], na.rm = T)
        pos=pos+1 
        summary_met[i,pos]=max(tmp_all[,i], na.rm = T)
        pos=pos+1 
        summary_met[i,pos]=mean(tmp_all[,i], na.rm = T)
        pos=pos+1
        summary_met[i,pos]=median(tmp_all[,i], na.rm = T)
        pos=pos+1 
        summary_met[i,pos]=sd(tmp_all[,i], na.rm = T)
        pos=pos+1 
        summary_met[i,pos]=length(boxplot(tmp_all[,i], plot=FALSE, range = 3)$out)
        if (sum(!is.na(tmp_all[,i])) < 10){
          pos=pos+1 
          summary_met[i,pos]="NA"
          } else {
          #nor=NA
          nor=shapiro.test(tmp_all[,i])
          pos=pos+1 
          summary_met[i,pos]=nor$p.value>=0.05
          }
        }
      }
    }
  summary_met[,1]=colnames(feature_matrix)
  my_colnames=c("Metabolite", "NAs","Non_NAs", "Min", "Max", "Mean", "Median","SD", "Outliers_x3IQR","Normal_distrib")
  for (o in levels (pheno[,1])){
    my_colnames=c(my_colnames,paste("NAs",o, sep="_"))
    my_colnames=c(my_colnames,paste("Non_NAs",o,sep="_"))
    my_colnames=c(my_colnames,paste("Min",o,sep="_"))
    my_colnames=c(my_colnames,paste("Max",o,sep="_"))
    my_colnames=c(my_colnames,paste("Mean",o,sep="_"))
    my_colnames=c(my_colnames,paste("Median",o,sep="_"))
    my_colnames=c(my_colnames,paste("SD",o,sep="_"))
    my_colnames=c(my_colnames,paste("Outliers_x3IQR",o,sep="_"))
    my_colnames=c(my_colnames,paste("Normal_distrib",o,sep="_"))
  }
  colnames(summary_met)=my_colnames
  return(summary_met)
}  



###############################################
# Filtering and transforming taxa             #
###############################################

transform_and_filter_mtb=function(x, samples_row=T, method="asin", missing_filter=0){
  #x[x=="NA"]=0
  #x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
  
    x=as.data.frame(t(x))
    print("transposing matrix!")
  
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
  
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100") 
  
  }
  
  #x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  x_filt=x[,((colSums(!is.na(x)) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  
  if (method=="asin"){
    print ("ASIN")
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))

  } else if (method=="log"){
    print ("LOG10")
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0])/2
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)

  }else if (method=="clr"){
    print ("CLR")
    #Adapted from Alexander Kurilshikov 
    #x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    #my_min=min(x[x>0], na.rm=T)/2
    #x=x+my_min
    #Calculate geometric mean
    gm_mean = function(x, na.rm=TRUE){
      exp(mean(log(x),na.rm=T))
    }
    Gmean_core = apply(x, 1, gm_mean)
    x = cbind(Gmean_core,x)
    d <- t(apply(x, 1, function(b) {
      log(b/b[1])[-1]
    }))
    x=d
    x_filt=x[,colnames(x) %in%colnames(x_filt)]
  }
  return(as.data.frame(x_filt))
}


transform_and_filter_taxa=function(x, samples_row=T, method="asin", missing_filter=0){
  x[x=="NA"]=0
  x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
  
    x=as.data.frame(t(x))
  
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
  
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100") 
  
  }
  
  x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  #x_filt=x[,((colSums(is.na(x)) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  
  if (method=="asin"){
    print ("ASIN")
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))

  } else if (method=="log"){
    print ("LOG10")
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0])/2
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)

  }else if (method=="clr"){
    print ("CLR")
    #Adapted from Alexander Kurilshikov 
    #x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    my_min=min(x[x>0], na.rm=T)/2
    x=x+my_min
    #Calculate geometric mean
    gm_mean = function(x, na.rm=TRUE){
      exp(mean(log(x)))
    }
    Gmean_core = apply(x, 1, gm_mean)
    x = cbind(Gmean_core,x)
    d <- t(apply(x, 1, function(b) {
      log(b/b[1])[-1]
    }))
    x=d
    x_filt=x[,colnames(x) %in%colnames(x_filt)]
  }
  return(as.data.frame(x_filt))
}

