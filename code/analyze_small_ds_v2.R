library(tidyverse)
library(limma)
library(variancePartition)
library(MetaIntegrator)
source("code/00_utils.R")
# load

my_studies <- c('GSE103174',
                'GSE13896','GSE16149','GSE17913','GSE18723','GSE19027',
                'GSE20189','GSE2125','GSE21862','GSE31210','GSE32539',
                'GSE42057','GSE42743','GSE4302','GSE44456','GSE46699',
                'GSE56768','GSE7895','GSE87072', 'GSE8987','GSE994')
my_studies2 <- c("GSE14633", "GSE5056", "GSE27002",  "E-MTAB-5278", "E-MTAB-5279" ,
                 "GSE20681" , "GSE23323" , "GSE23515" , "GSE30272"  , "GSE32504"  ,  "E-MTAB-3604")



args <- commandArgs(trailingOnly = T)
idx <- as.numeric(args[[1]])
print(idx)
if (idx < 22 ){
  ds <- my_studies[[idx]] 
  orig = T
} else{
  ds <- my_studies2[[idx]]
  orig = F
}


# load the data
if (!orig){
  study_file <- sprintf("data/small_gses2/%s.RData", ds)
  pDat0 <- read_csv(sprintf("data/pdata_filt2_v2/%s.csv", tolower(ds)))
  load(study_file) # my_gse
}


# old way
if (orig){
  study_file <- sprintf("data/small_studies/%s.RData", ds)
  pDat0 <- read_csv(sprintf("data/pdata_filt/%s.csv", tolower(ds))) %>% rename(sex=sex_lab)
  my_gse <- gse$originalData[[1]] # note - manually removed extra
}


# prepare the data
eDat0 <- my_gse$expr
platform <- my_gse$platform

pDat <- pDat0 %>% 
  select(sample_acc, sex, smok) %>%
  rename(smoking=smok) %>%
  filter(sex %in% c("male", "female"),
         smoking %in% c("S", "NS"))

eDat <- eDat0[,pDat$sample_acc]


# ------  smoking model ------- #
if (length(unique(pDat$sex))==1){
  ss <- T
  # only one sex
  
  design <- model.matrix(~smoking, data=pDat) 
  
} else {
  ss <- F
  design <- model.matrix(~ sex + smoking + sex*smoking, data=pDat) 
  
}
fit <- lmFit(eDat, design)
fit <- eBayes(fit)
df_smok <- data.frame(topTable(fit, coef="smokingS",  number=nrow(fit), confint=T)) 


# add gene annotations
# way 1:
df_smok1 <- add_annot_mi(df_smok, platform)
print(sprintf("%s unique genes present from %s probes", length(unique(df_smok1$gene)),
              nrow(df_smok1)))

# run MA
smok_ma <- run_clean_ma(df_smok1 %>% 
                           rename(ID=probes, geneSymbol=gene) %>% 
                           mutate(chromosome=""))
save(df_smok1, smok_ma, file=sprintf("data/results/%s_de_s.RData", ds))

print("saved smoking")
if (!ss){
  df_sex <- data.frame(topTable(fit, coef="sexmale",  number=nrow(fit), confint=T)) 
  df_int <- data.frame(topTable(fit2, coef="sexmale:smokingS",  number=nrow(fit2), confint=T))
 
  df_sex1 <- add_annot_mi(df_sex, platform)
  df_int1 <- add_annot_mi(df_int, platform)
  
  int_ma <- run_clean_ma(df_int1 %>% 
                             rename(ID=probes, geneSymbol=gene) %>% 
                             mutate(chromosome=""))
  
  sex_ma <- run_clean_ma(df_sex1 %>% 
                           rename(ID=probes, geneSymbol=gene) %>% 
                           mutate(chromosome=""))
  
  save(df_smok1, df_int1, df_sex1,
       smok_ma, int_ma, sex_ma, file=sprintf("data/results/%s_de_ma.RData", ds))
}
print("running PVCA")
# -----  RUN PVCA  ----- #
set.seed(1014)
rs <- rowSums(is.na(eDat))
eDat_complete <- eDat[rs==0,]
print(sprintf("Removed %s rows out of %s because of NAs", nrow(eDat)-nrow(eDat_complete),
              nrow(eDat)))
pcs_l <- varPartPC(eDat_complete)
expr_pcs_t <- t(pcs_l$pcs)
props <- pcs_l$props

if (ss){
  my_model <- as.formula("~smoking")
} else {
  my_model <- as.formula("~ sex + smoking + sex:smoking")
  
}
(est_var <- get_pvca(expr_pcs_t, pDat, my_model, props))

save(est_var,
     file=sprintf("data/results/%s_pvca.RData", ds))

rand_var <- rand_pvca(100, expr_pcs_t, pDat, my_model, props)
save(rand_var, 
     file=sprintf("data/results/%s_rand.RData", ds))


#######
grab_study <- function(ds){
  print(ds)
  smok_f <- sprintf("data/results/%s_de_s.RData", ds)
  if (file.exists(smok_f)){
    load(smok_f) # df_smok1.1, smok_ma1
    gene_smok <- smok_ma
    probe_smok <- df_smok1
  } else {
    gene_smok <- NA
    probe_smok <- NA
  }
  var_f <- sprintf("data/results/%s_pvca.RData", ds)
  if (file.exists(var_f)){
    load(var_f) # est_var
  } else {
    est_var <- NA
  }
  r_var <- sprintf("data/results/%s_rand.RData", ds)
  if (file.exists(r_var)){
    load(r_var) # rand_var
  } else {
    rand_var <- NA
  }
  
  full_f <- sprintf("data/results/%s_de_ma.RData", ds)
  if (file.exists(full_f)){
    load(full_f)

    gene_smok <- smok_ma
    probe_smok <- df_smok1
    gene_int <- int_ma
    probe_int <- df_int1
    probe_sex <- df_sex1
    gene_sex <- sex_ma
  } else {

    gene_int  <- NA
    probe_int <- NA
    probe_sex <- NA
    gene_sex <- NA
  }
  my_study = list("study"=ds, "var"=est_var, "rand"=rand_var,
                  "probe_smok"=probe_smok, "gene_smok"=gene_smok,
                  "probe_sex"=probe_sex, "gene_sex"=gene_sex,
                  "probe_int"=probe_int, "gene_int"=gene_int)
  return(my_study)
}

my_studies <- c('GSE103174',
                'GSE13896','GSE16149','GSE17913','GSE18723','GSE19027',
                'GSE20189','GSE2125','GSE21862','GSE31210','GSE32539',
                'GSE42057','GSE42743','GSE4302','GSE44456','GSE46699',
                'GSE56768','GSE7895','GSE87072', 'GSE8987','GSE994',
                "GSE14633", "GSE5056", "GSE27002",  "E-MTAB-5278", "E-MTAB-5279" ,
                "GSE20681" , "GSE23323" , "GSE23515" , "GSE30272"  , 
                "GSE32504"  ,  "E-MTAB-3604")
                
                
all_studies <- lapply(my_studies, grab_study)
save(all_studies, file="data/all_studies_v2.RData")

