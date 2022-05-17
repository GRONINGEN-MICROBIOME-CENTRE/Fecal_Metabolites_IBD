
# this script is to estimate
1) the relation between sample size and detection power while giving a grid search in variance expalained by SNP (0.01~0.1)
2) the relation bewteen variance explanation and detection power while giving a grid search in present rate (0.1~1)

library(pwr)

##################################
# sample size vs. power
##################################

N <- seq(10,150,by=5) # 0.9 means the 90% present rate of metabolite
H2 <- seq(0.001,0.1,by=0.0001) # set a range of h2
alpha <- 0.05/length(PCA_feature) # genome-wide/study-wide significance

threshold <- qchisq(alpha, df = 1, lower.tail = FALSE) 

paraComb <- expand.grid(N, H2)
ncp <- with(paraComb, Var1 * Var2)

present90=setNames(cbind(paraComb,
               sapply(ncp, 
                      function(ncp) pchisq(threshold, 
                                           df = 1, lower.tail = FALSE, 
                                           ncp = ncp))
),
c("N", "H2", "power"))

ggplot(data=present90, aes(x=N, y=power, group=H2,color=H2)) +
  geom_line(size=0.5)+scale_color_gradient(low="#FDB366", high="#A50026",limits=c(0.1,0.4),oob = scales::squish)+
  theme_classic()+
  geom_hline(yintercept = 0.8)+
  geom_segment(aes(x = 50, y = 0, xend = 50, yend = 0.8),color="black")
ggsave("Plot/Power.calculation.pdf",width = 5,height = 5)

##################################
# variation vs. sample size
##################################

N <- c(100*0.9,500*0.9,1000*0.9,2000)
H2 <- seq(0.001,0.1,by=0.0001)
alpha <- 0.05/length(PCA_feature)

calculate=matrix(nrow = length(H2),ncol = 4)
calculate=as.data.frame(calculate)
colnames(calculate)=c("SampleSize","H2","Power","PresenRate")
for(i in 1:length(H2)){
  
  H2.tmp=H2[i]
  aa=pwr.t.test(n = NULL , d = H2.tmp, sig.level = alpha, power = 0.8, type = "one.sample")
  calculate$H2[i]=H2.tmp
  calculate$SampleSize[i]=aa$n
  calculate$Power[i]=0.8
  calculate$PresenRate[i]="100%"
  
}

calculate1=calculate
calculate2=calculate
calculate3=calculate
calculate4=calculate
calculate5=calculate
calculate1$SampleSize=calculate1$SampleSize/0.8
calculate1$PresenRate="80%"
calculate2$SampleSize=calculate2$SampleSize/0.6
calculate2$PresenRate="60%"
calculate3$SampleSize=calculate3$SampleSize/0.4
calculate3$PresenRate="40%"
calculate4$SampleSize=calculate4$SampleSize/0.2
calculate4$PresenRate="20%"
calculate5$SampleSize=calculate5$SampleSize/0.1
calculate5$PresenRate="10%"

calculate=rbind(calculate,calculate1,calculate2,calculate3,calculate4,calculate5)
calculate$PresenRate=factor(calculate$PresenRate,levels = c("10%","20%","40%","60%","80%","100%"))
ggplot(data=calculate, aes(x=H2, y=SampleSize, color=PresenRate)) +
  geom_line(size=1)+scale_color_nejm()+
  theme_classic()+
ggsave("Plot/Power.calculation.present.pdf",width = 10,height = 5)
