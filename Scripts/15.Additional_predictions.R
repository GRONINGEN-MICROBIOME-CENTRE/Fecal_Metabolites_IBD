#Load table containing participants without flares after samples collection
load("~/Desktop/IBD_metabolomics_2022/3.Workspace/CoDaCore_pred.RData")
timetoflare <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/timetoflare.txt")
coltime <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/timeproduction.txt")
IDs_IBD = read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/IDs_IBD.txt", header=FALSE)
tt2=merge(IDs_IBD,timetoflare,by.x="V2", by.y="UMCGIBDResearchIDorLLDeepID")
tt3=merge(tt2,coltime,by.x="V5", by.y="AA_ID")

#Select participants that didn't had a flare in the past 1 years
no_flare=subset(tt3, tt3$TimeEndPreviousExacerbation>365)
no_flare_phenos=subset(ibd_phenos2, row.names(ibd_phenos2) %in% no_flare$V5)
no_flare_phenos=subset(no_flare_phenos, no_flare_phenos$clinical_Calprot200==1)
no_flare_2 <- read.delim("~/Desktop/IBD_metabolomics_2022/1.Input/phenos_IBD_clean_v2.txt")
no_flare_2 = no_flare_2[,c("AA_ID", "ibd_SSCAI", "ibd_ActiveDisease", "ibd_HarveyBradshaw")]
no_flare_3=subset(no_flare_2, no_flare_2$AA_ID %in%  row.names(no_flare_phenos))
no_flare_3=subset(no_flare_3,no_flare_3$ibd_ActiveDisease!="Active")

my_ratio_IBD_model=getLogRatios(IBD_model_noresc)
train_correct_IBD$IBD2="aControl"
train_correct_IBD$IBD2[train_correct_IBD$IBD==2]="bIBD"
train_correct_IBD$IBD2[train_correct_IBD$IBD==1]="bIBD"
train_correct_IBD$IBD2[row.names(train_correct_IBD)%in%no_flare_3$AA_ID]="cNo_Flare"

ggplot(train_correct_IBD, aes(IBD2, metabolite, fill=IBD2)) + geom_boxplot() + theme_classic() + ylab ("Metabolite ratio") + xlab("") + scale_fill_manual(values = c("white",  "firebrick3","salmon"))

#### Model with calprotectin

my_phenos=mtb_cc2[,c("host_Age","host_BMI","host_Sex", "clinical_Calprot200","IBD")]
my_phenos$IBD[my_phenos$IBD=="Control"]=0
my_phenos$IBD[my_phenos$IBD=="IBD"]=1
my_phenos$IBD=as.numeric(as.character(my_phenos$IBD))
my_pseudocount=min(all_metabolites_raw, na.rm = T)/2
my_ratio_mtb=all_metabolites_raw[,c("lactosyl_N_palmitoyl_sphingosine_d18_1_16_0", "L_urobilin")]
my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0[is.na(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0)]=min(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0, na.rm = T) /2
my_ratio_mtb$L_urobilin[is.na(my_ratio_mtb$L_urobilin)]=min(my_ratio_mtb$L_urobilin, na.rm = T) /2
my_ratio_mtb$ratio=log(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0/my_ratio_mtb$L_urobilin)
my_phenos2=merge(my_phenos, my_ratio_mtb, by="row.names")
row.names(my_phenos2)=my_phenos2$Row.names
my_phenos2$Row.names=NULL
my_phenos2$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0=NULL
my_phenos2$L_urobilin=NULL

my_test_cnt=sample(row.names(my_phenos[my_phenos$IBD==0,]), size = 62)
my_test_cnt=c(my_test_cnt, no_flare_3$AA_ID )

aucs=c()
library(Hmisc)
for (i in 1:100){
  my_test_cnt=sample(row.names(my_phenos[my_phenos$IBD==0,]), size = 62)
  my_test_cnt=c(my_test_cnt, no_flare_3$AA_ID )
  my_train_fl=my_phenos2[row.names(my_phenos2) %nin%my_test_cnt,]
  my_test_fl=my_phenos2[row.names(my_phenos2) %in%my_test_cnt,]
  model_fl=glm(IBD ~ ., data = my_train_fl, family = "binomial")

  test_model_fl = predict(model_fl, newdata = my_test_fl, type = "response")
  test_roc = roc(my_test_fl$IBD ~ test_model_fl, plot = TRUE, print.auc = TRUE)
  aucs=c(aucs,as.numeric(test_roc$auc))
}

my_test_fl_no_ratio=my_test_fl
my_test_fl_no_ratio$ratio=NULL
withratio <- glm(IBD ~ ., data=my_test_fl, family='binomial')
withoutratio <- glm(IBD ~ ., data=my_test_fl_no_ratio, family='binomial')
library(lmtest)
lrtest(withratio, withoutratio)

####

#### Model without calprotectin

my_phenos=mtb_cc2[,c("host_Age","host_BMI","host_Sex","IBD")]
my_phenos$IBD[my_phenos$IBD=="Control"]=0
my_phenos$IBD[my_phenos$IBD=="IBD"]=1
my_phenos$IBD=as.numeric(as.character(my_phenos$IBD))
my_pseudocount=min(all_metabolites_raw, na.rm = T)/2
my_ratio_mtb=all_metabolites_raw[,c("lactosyl_N_palmitoyl_sphingosine_d18_1_16_0", "L_urobilin")]
my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0[is.na(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0)]=min(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0, na.rm = T) /2
my_ratio_mtb$L_urobilin[is.na(my_ratio_mtb$L_urobilin)]=min(my_ratio_mtb$L_urobilin, na.rm = T) /2
my_ratio_mtb$ratio=log(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0/my_ratio_mtb$L_urobilin)
my_phenos2=merge(my_phenos, my_ratio_mtb, by="row.names")
row.names(my_phenos2)=my_phenos2$Row.names
my_phenos2$Row.names=NULL
my_phenos2$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0=NULL
my_phenos2$L_urobilin=NULL

my_test_cnt=sample(row.names(my_phenos[my_phenos$IBD==0,]), size = 62)
my_test_cnt=c(my_test_cnt, row.names(no_flare_phenos) )

library(Hmisc)
my_train_fl=my_phenos2[row.names(my_phenos2) %nin%my_test_cnt,]
my_test_fl=my_phenos2[row.names(my_phenos2) %in%my_test_cnt,]
model_fl=glm(IBD ~ ., data = my_train_fl, family = "binomial")

test_model_fl = predict(model_fl, newdata = my_test_fl, type = "response")
test_roc = roc(my_test_fl$IBD ~ test_model_fl, plot = TRUE, print.auc = TRUE)

my_test_fl_no_ratio=my_test_fl
my_test_fl_no_ratio$ratio=NULL
withratio <- glm(IBD ~ ., data=my_test_fl, family='binomial')
withoutratio <- glm(IBD ~ ., data=my_test_fl_no_ratio, family='binomial')
library(lmtest)
lrtest(withratio, withoutratio)



## Check other data

## Download two datasets available in https://www.metabolomicsworkbench.org/ 
# Blood metabolomics: Scoville, E.A., Allaman, M.M., Brown, C.T. et al. Alterations in lipid, amino acid, and energy metabolism distinguish Crohn’s disease from ulcerative colitis and control subjects by serum metabolomic profiling. Metabolomics 14, 17 (2018). https://doi.org/10.1007/s11306-017-1311-y
# Faecal metabolomics: HMP2: doi: 10.21228/M82T15

hmp2_ratio=read.table("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/prediction/faecal_ratios.txt")
blood_ratios <- read.delim("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/prediction/blood_ratios.txt")

hmp2_metadata_metabolites <- read_excel("~/Desktop/IBD_metabolomics_2022/7.2.Rebuttal-Gut/hmp2_metadata_metabolites.xlsx")
hmp2_metadata=hmp2_metadata_metabolites[,c(2,3,7,8,110)]

hmp2_ratios2=merge(hmp2_ratios, hmp2_metadata, by.x = "Samples", by.y = "Sample name")
hmp2_ratios2$urobilin[is.na(hmp2_ratios2$urobilin)]=min(hmp2_ratios2$urobilin, na.rm = T)/2
hmp2_ratios2$ratio=log(hmp2_ratios2$C16.0.Ceramide..d18.1./hmp2_ratios2$urobilin)

ggplot(hmp2_ratios2, aes(as.factor(`Visit number`), ratio, fill=IBD)) + geom_boxplot() + theme_classic() + ylab ("Metabolite ratio") + xlab("") + scale_fill_manual(values = c("white",  "firebrick3","salmon"))
wilcox.test(hmp2_ratios2$ratio~hmp2_ratios2$IBD)

# Subset visit with at least 5 controls

table(hmp2_ratios2$`Visit number`,hmp2_ratios2$IBD)
visits_to_keep=c(1,4,5,9,14,15,17,19,21,24)

hmp2_ratios3=subset(hmp2_ratios2, hmp2_ratios2$`Visit number` %in% visits_to_keep)

for (i in unique(hmp2_ratios3$`Visit number`)){
  print(i)
  my_sub=subset(hmp2_ratios3,hmp2_ratios3$`Visit number`==i)
  my_test=wilcox.test(my_sub$ratio~my_sub$IBD)
  print(my_test$p.value)
}

ggplot(hmp2_ratios3, aes(as.factor(`Visit number`), ratio, fill=IBD)) +  theme_classic() +  ylab ("ratio:  Ceramide (18:1/16:0) / Urobilin") + xlab("Week") + geom_line(aes(group=`Subject name`, color=IBD), position = position_dodge(0.2), alpha=0.3) +geom_point(aes(fill=IBD,group=`Subject name`),size=2,shape=21, position = position_dodge(0.2)) + scale_fill_manual(values = c("white",  "firebrick3","salmon")) + scale_colour_manual(values = c("black",  "firebrick3"))
ggplot(hmp2_ratios3, aes(as.factor(`Visit number`), ratio, fill=IBD)) + geom_boxplot() + theme_classic() + ylab ("ratio:  Ceramide (18:1/16:0) / Urobilin") + xlab("Week") + scale_fill_manual(values = c("white",  "firebrick3","salmon")) +geom_point(aes(fill=IBD,group=`Subject name`),size=2,shape=21, position = position_dodge(0.2)) 
ggplot(hmp2_ratios3, aes(as.factor(`Visit number`), ratio, fill=IBD)) + geom_boxplot() + theme_classic() + ylab ("ratio:  Ceramide (18:1/16:0) / Urobilin") + xlab("Week") + scale_fill_manual(values = c("white",  "firebrick3","salmon")) +geom_point(aes(fill=IBD),size=2,shape=21, position = position_dodge(0.5)) +  stat_compare_means( aes(group = IBD),label = "p.format", size=3) 


#Blood
blood_ratios$I.urobilinogen[is.na(blood_ratios$I.urobilinogen)]=min(blood_ratios$I.urobilinogen, na.rm = T)/2
blood_ratios$ratio=log(blood_ratios$lactosyl.N.palmitoyl.sphingosine/blood_ratios$I.urobilinogen)

ggplot(blood_ratios, aes(IBD, ratio, fill=IBD)) + geom_boxplot() + theme_classic() + ylab ("ratio: lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)/urobilinogen") + xlab("") + scale_fill_manual(values = c("white",  "firebrick3","salmon")) + stat_compare_means()
wilcox.test(blood_ratios$ratio~blood_ratios$IBD)




#### Test with and without smoking
host_SmokeCurrentSmoker

my_phenos=cc_pheno[,c("host_Age","host_BMI","host_Sex", "clinical_Calprot200","IBD","host_SmokeCurrentSmoker")]
my_phenos$IBD[my_phenos$IBD==1]=0
my_phenos$IBD[my_phenos$IBD==2]=1
my_phenos$IBD=as.numeric(as.character(my_phenos$IBD))
my_pseudocount=min(all_metabolites_raw, na.rm = T)/2
my_ratio_mtb=all_metabolites_raw[,c("lactosyl_N_palmitoyl_sphingosine_d18_1_16_0", "L_urobilin")]
my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0[is.na(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0)]=min(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0, na.rm = T) /2
my_ratio_mtb$L_urobilin[is.na(my_ratio_mtb$L_urobilin)]=min(my_ratio_mtb$L_urobilin, na.rm = T) /2
my_ratio_mtb$ratio=log(my_ratio_mtb$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0/my_ratio_mtb$L_urobilin)
my_phenos2=merge(my_phenos, my_ratio_mtb, by="row.names")
row.names(my_phenos2)=my_phenos2$Row.names
my_phenos2$Row.names=NULL
my_phenos2$lactosyl_N_palmitoyl_sphingosine_d18_1_16_0=NULL
my_phenos2$L_urobilin=NULL

aucs=c()
ps=c()

for (i in 1:100){
  my_test_cnt=sample(row.names(my_phenos2), size = 170)
  my_train_fl=my_phenos2[row.names(my_phenos2) %nin%my_test_cnt,]
  my_test_fl=my_phenos2[row.names(my_phenos2) %in%my_test_cnt,]
  model_fl=glm(IBD ~ ., data = my_train_fl, family = "binomial")
  
  test_model_fl = predict(model_fl, newdata = my_test_fl, type = "response")
  test_roc = roc(my_test_fl$IBD ~ test_model_fl, plot = TRUE, print.auc = TRUE)
  model_with=data.frame(Test=paste("test",i,sep="_"), Sensitivity=test_roc$sensitivities, Specificity=test_roc$specificities, Model="Age + Sex + BMI + Calprotectin + Metabolite ratio + Smoking")
  
  my_train_fl2=my_train_fl
  my_train_fl2$host_SmokeCurrentSmoker=NULL
  my_test_fl2=my_test_fl
  my_test_fl2$host_SmokeCurrentSmoker=NULL
  
  model_fl2=glm(IBD ~ ., data = my_train_fl2, family = "binomial")
  test_model_fl2 = predict(model_fl2, newdata = my_test_fl2, type = "response")
  test_roc2 = roc(my_test_fl2$IBD ~ test_model_fl2, plot = TRUE, print.auc = TRUE)
  model_with2=data.frame(Test=paste("test",i,sep="_"), Sensitivity=test_roc2$sensitivities, Specificity=test_roc2$specificities, Model="Age + Sex + BMI + Calprotectin + Metabolite ratio")
  
  my_result=rbind(model_with, model_with2)
  
  if (i==1){
    my_result_final=my_result
  } else{
    my_result_final=rbind(my_result_final,my_result)
  }
  
  #my_test_fl_no_ratio=my_test_fl
  #my_test_fl_no_ratio$host_SmokeCurrentSmoker=NULL
  #withsmk <- glm(IBD ~ ., data=my_test_fl, family='binomial')
  #withoutsmk <- glm(IBD ~ ., data=my_test_fl_no_ratio, family='binomial')
  #comp=lrtest(withsmk, withoutsmk)
  #aucs=c(aucs,as.numeric(test_roc$auc))
  #ps=c(ps,as.numeric(comp$`Pr(>Chisq)`[2]))
}

ci.list <- lapply(list(test_roc2, test_roc), ci.se, specificities = seq(0, 1, l = 25))
dat.ci.list <- lapply(ci.list, function(ciobj)  data.frame(x = as.numeric(rownames(ciobj)),lower = ciobj[, 1], upper = ciobj[, 3]))
p <- ggroc(list(test_roc2, test_roc)) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal()

for(i in 1:2) {
  p <- p + geom_ribbon(data = dat.ci.list[[i]],aes(x = x, ymin = lower, ymax = upper),fill = i + 1,alpha = 0.2,inherit.aes = F) 
} 

#In CD
cd_sample=row.names(ibd_phenos2)[ibd_phenos2$CD_or_UC=="CD"]
aucs=c()
ps=c()

for (i in 1:100){
  my_test_cnt=sample(row.names(my_phenos2), size = 170)
  my_train_fl=my_phenos2[row.names(my_phenos2) %nin%my_test_cnt,]
  my_test_fl=my_phenos2[row.names(my_phenos2) %in%my_test_cnt,]
  model_fl=glm(IBD ~ ., data = my_train_fl, family = "binomial")
  
  test_model_fl = predict(model_fl, newdata = my_test_fl, type = "response")
  test_roc = roc(my_test_fl$IBD ~ test_model_fl, plot = TRUE, print.auc = TRUE)
  my_test_fl_no_ratio=my_test_fl
  my_test_fl_no_ratio$host_SmokeCurrentSmoker=NULL
  withsmk <- glm(IBD ~ ., data=my_test_fl, family='binomial')
  withoutsmk <- glm(IBD ~ ., data=my_test_fl_no_ratio, family='binomial')
  comp=lrtest(withsmk, withoutsmk)
  aucs=c(aucs,as.numeric(test_roc$auc))
  ps=c(ps,as.numeric(comp$`Pr(>Chisq)`[2]))
}


my_test_fl_no_ratio=my_test_fl
my_test_fl_no_ratio$ratio=NULL
withratio <- glm(IBD ~ ., data=my_test_fl, family='binomial')
withoutratio <- glm(IBD ~ ., data=my_test_fl_no_ratio, family='binomial')
library(lmtest)
lrtest(withratio, withoutratio)



