

rare_metabolites=all_new_ID_raw[,to_remove$names]
rare_metabolites[is.na(rare_metabolites)]=0
rare_metabolites[rare_metabolites>0]=1
meds_IBD=ibd_test_phenos[,grepl("med",colnames(ibd_test_phenos))]

test_rare_meds=merge(meds_IBD,rare_metabolites, by="row.names")


m_res=matrix(nrow=24175,ncol = 5)
v=1
for (i in 2:50){
  for (a in 51:ncol(test_rare_meds)){
    m_res[v,1]=colnames(test_rare_meds)[i]
    m_res[v,2]=colnames(test_rare_meds)[a]
    m_res[v,3]=sum(test_rare_meds[,a] !=0) #N.samples metabolite detected
    m_res[v,4]=sum(test_rare_meds[,i] ==2) #N. users
    if(m_res[v,3]>3){
      m_res[v,5]=chisq.test(test_rare_meds[,i],as.factor(test_rare_meds[,a]))$p.value
      v=v+1
    } else {
      m_res[v,5]=1
      v=v+1
    }
  }
}

m_res=m_res[complete.cases(m_res[,1]),]
m_res=as.data.frame(m_res)
m_res$bnf=p.adjust(as.numeric(as.character(m_res[,5])), method = "bonferroni")
m_res$FDR=p.adjust(as.numeric(as.character(m_res[,5])), method = "BH")
colnames(m_res)=c("Drug", "Metabolite", "Number_samples_metabolite_detected","Number of medication users","P-value","BNF","FDR")



meds_cnt=phenos_lld[,grepl("med",colnames(phenos_lld))]

test_rare_meds=merge(meds_cnt,rare_metabolites, by="row.names")


m_res=matrix(nrow=24175,ncol = 5)
v=1
for (i in 2:46){
  for (a in 47:ncol(test_rare_meds)){
    m_res[v,1]=colnames(test_rare_meds)[i]
    m_res[v,2]=colnames(test_rare_meds)[a]
    m_res[v,3]=sum(test_rare_meds[,a] !=0) #N.samples metabolite detected
    m_res[v,4]=sum(test_rare_meds[,i] ==1) #N. users
    if(m_res[v,3]>3 && m_res[v,4] >1 ){
      m_res[v,5]=chisq.test(test_rare_meds[,i],as.factor(test_rare_meds[,a]))$p.value
      v=v+1
    } else {
      m_res[v,5]=1
      v=v+1
    }
  }
}

m_res=m_res[complete.cases(m_res[,1]),]
m_res=as.data.frame(m_res)
m_res$bnf=p.adjust(as.numeric(as.character(m_res[,5])), method = "bonferroni")
m_res$FDR=p.adjust(as.numeric(as.character(m_res[,5])), method = "BH")
colnames(m_res)=c("Drug", "Metabolite", "Number_samples_metabolite_detected","Number of medication users","P-value","BNF","FDR")