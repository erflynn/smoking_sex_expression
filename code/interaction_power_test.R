# https://davidbaranger.com/2019/08/06/interaction-analyses-power-part-1/
# non-parametric methods
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-186#Sec9
library(pwr)
pwr.r.test(r=0.2,sig.level = 0.05,power = 0.8) 
# effect size



# possibilty:
# two step:
# - main effect in one study
# - interaction effect in follow up studies


# microarray power test
# sd matched pairs
# number of genes = num_non_de + num_de
# effect size
# FDR

library("ssize.fdr")
d<-1 ##effect size
s<-0.6   ##standard deviation
a<-0.05 ##false discovery rate to be controlled
pwr<-0.8 ##desired power
p0<-c(0.5,0.9,0.95, 0.99) ##proportions of non-differentially expressed genes
N<-30 ##maximum sample size for calculations
os<-ssize.oneSamp(delta=d,sigma=s,fdr=a,power=pwr,
                  pi0=p0,maxN=N,side="two-sided")
os$ssize ##first sample sizes to reach desired power

# y ~ sex + smoking + beta*(smoking*sex)

# y ~ beta_m*smoking (in males)
# y ~ beta_f*smoking (in females)
# beta_m - beta_f <-- cohen's d


###
d<-0.1 ##interaction effect size
s<-0.6   #residual std deviation of linear model x sqrt(3/4)
# sqrt(3)/2: assumes balanced
a<-0.05 ##false discovery rate to be controlled
pwr<-0.8 ##desired power
p0<-c(0.5,0.9,0.95, 0.99) ##proportions of non-differentially expressed genes
N<-1000 ##maximum sample size for calculations
os<-ssize.oneSamp(delta=d,sigma=s,fdr=a,power=pwr,
                  pi0=p0,maxN=N,side="two-sided")
os$ssize ##first sample sizes to reach desired power



# FDR of 0.05 with 0.95 NDE genes
# 1000 genes, identify 50 genes, 0.05*50=2.5 FPs

# FDR of 0.05 with 0.99 NDE genes
# 2000 genes, identify 20 genes, 0.05*20=1 FP
# 20,000 genes, identify 200 genes, 10 FPs
# power 0.8 =  80% prob of finding an effect if it exists

d<-0.7 ##effect size
s<-0.6   ##standard deviation
a<-0.05 ##false discovery rate to be controlled
pwr<-0.9 ##desired power
p0<-c(0.9,0.95, 0.99) ##proportions of non-differentially expressed genes
N<-30 ##maximum sample size for calculations
os<-ssize.oneSamp(delta=d,sigma=s,fdr=a,power=pwr,pi0=p0,maxN=N,side="two-sided")
os$ssize ##first sample sizes to reach desired power



# effect size 0.5 to 1.5 vs sample_size vs proportions
# log2 diff in means?
list_es <- seq(0.05, 0.6,0.05)
dfs <- do.call(rbind, lapply(list_es, function(d){
  os<-ssize.oneSamp(delta=d,
                    sigma=s*sqrt(3)/2,fdr=a,power=0.8,
                    pi0=p0,maxN=1000,side="two-sided")
  df <- data.frame(os$ssize)
  df$d <- d
  return(df)
})) %>%
  mutate(power=0.8)

ggplot(dfs %>% mutate(prop_non_de=factor(pi0)), aes(x=d, y=ssize, col=prop_non_de))+
  geom_point()+
  geom_line(lty=2)+
  theme_bw()+
  xlab("effect size")+
  ylab("sample size")+
  ggtitle("power = 0.8, sd=0.6")


dfs2 <- data.frame(do.call(rbind, lapply(list_es, function(d){
  os<-ssize.oneSamp(delta=d,sigma=s,fdr=a,power=0.9,pi0=p0,maxN=50,side="two-sided")
  df <- data.frame(os$ssize)
  df$d <- d
  return(df)
})) %>%
  mutate(power=0.9))

ggplot(data.frame(dfs %>% bind_rows(dfs2)) %>% mutate(prop_non_de=factor(pi0)), 
       aes(x=d, y=ssize, col=prop_non_de))+
  geom_point()+
  geom_line(aes(lty=factor(power)))+
  theme_bw()+
  xlab("effect size")+
  ylab("sample size")+
  geom_hline(yintercept=16, col="gray")


#### SET UP WITH INTERACTION TEST? ####
# assume ordinal (we have higher power for disordinal but they are rare)

df <- tibble(sex=c(0,0,1,1), smok=c(0,1,0,1))
design <- model.matrix(~sex + smok + smok:sex+0, data=df)
b<-c(0.5, 0.6, 0.1 )		
df<-	
p0.F<-c(0.9,0.95,0.995)		##proportions of non-differentially expressed genes

ftv<-ssize.F(X=design,beta=b,dn=df,
                 s=1,fdr=0.05,power=0.8,
                 pi0=p0,maxN=100) # a,b inverse variance

ftv$ssize	##first sample sizes to reach desired power
ftv$power	##calculated power for each sample size
ftv$crit.vals




