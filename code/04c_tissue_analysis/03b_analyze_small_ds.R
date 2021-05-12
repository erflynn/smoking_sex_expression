
library(tidyverse)
library(limma)
library(variancePartition)
library(MetaIntegrator)
source("code/00_utils.R")
# load


args <- commandArgs(trailingOnly = T)
ds <- args[[1]]
print(ds)
study_file <- sprintf("data/small_studies/%s.RData", ds)
if (!file.exists(study_file)){
  gse <- getGEOData(ds )
  save(gse, file=study_file)
} else {
  load(study_file) # gse
}


pDat0 <- read_csv(sprintf("data/pdata_filt/%s.csv", tolower(ds)))

geo_obj <- gse$originalData[[1]] # note - manually removed extra
eDat0 <- geo_obj$expr
platform <- geo_obj$platform

pDat <- pDat0 %>% 
  select(sample_acc, sex_lab, smok) %>%
  rename(sex=sex_lab, smoking=smok) %>%
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
  design <- model.matrix(~ sex + smoking, data=pDat) 
  
}
fit <- lmFit(eDat, design)
fit <- eBayes(fit)
df_smok1 <- data.frame(topTable(fit, coef="smokingS",  number=nrow(fit), confint=T)) 



# add gene annotations

df_smok1.1 <- add_annot_mi(df_smok1, platform)

print(sprintf("%s unique genes present from %s probes", length(unique(df_smok1.1$gene)),
              nrow(df_smok1.1)))# 13041
smok_ma1 <- run_clean_ma(df_smok1.1 %>% 
                           rename(ID=probes, geneSymbol=gene) %>% 
                           mutate(chromosome=""))
save(df_smok1.1, smok_ma1, file=sprintf("data/results/%s_de_s.RData", ds))


# --> pDat, eDat
# variance
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
#rand_var <- rand_pvca(100, expr_pcs_t, pDat, my_model, props)

#plot_pvca_rand(est_var, rand_var)
save(est_var,# rand_var, 
     file=sprintf("data/results/%s_pvca.RData", ds))


if (!ss){
# smoking*sex model
design2 <- model.matrix(~ sex + smoking + smoking*sex, data=pDat) 
fit2 <- lmFit(eDat, design2)
fit2 <- eBayes(fit2)
df_int <- data.frame(topTable(fit2, coef="sexmale:smokingS",  number=nrow(fit2), confint=T))
df_smok2 <- data.frame(topTable(fit2, coef="smokingS",  number=nrow(fit2), confint = T))
save(df_smok1, df_int, df_smok2, file=sprintf("data/results/%s_de.RData", ds))

df_smok2.1 <- add_annot_mi(df_smok2, platfomr)
df_int1 <- add_annot_mi(df_int, platform)

# TODO: do any inconsistently multi-map?

# how many unique 


# ------ Meta-analysis ------ #


smok_ma2 <- run_clean_ma(df_smok2.1%>% 
                                  rename(ID=probes, geneSymbol=gene) %>% 
                                  mutate(chromosome=""))

smok_int <- run_clean_ma(df_int1 %>% 
                           rename(ID=probes, geneSymbol=gene) %>% 
                           mutate(chromosome=""))

save(df_smok1.1, df_int1, df_smok2.1,
     smok_ma1, smok_ma2, smok_int, file=sprintf("data/results/%s_de_ma.RData", ds))
}


### pull the results


grab_study <- function(ds){
  print(ds)
  smok_f <- sprintf("data/results/%s_de_s.RData", ds)
  if (file.exists(smok_f)){
    load(smok_f) # df_smok1.1, smok_ma1
    gene_smok <- smok_ma1
    probe_smok <- df_smok1.1
  } else {
    df_smok1.1 <- NA
    smok_ma1 <- NA
  }
  var_f <- sprintf("data/results/%s_pvca.RData", ds)
  if (file.exists(var_f)){
    load(var_f) # est_var
  } else {
    est_var <- NA
  }
  full_f <- sprintf("data/results/%s_de_ma.RData", ds)
  if (file.exists(full_f)){
    load(full_f)
    gene_smok <- smok_ma2
    probe_smok <- df_smok2.1
    gene_int <- smok_int
    probe_int <- df_int1
  } else {
    df_smok1.1 <- NA
    df_int1 <- NA
    df_smok2.1 <- NA
    smok_ma1 <- NA 
    smok_ma2 <- NA 
    smok_int <- NA
    gene_int  <- NA
    probe_int <- NA
  }
  my_study = list("study"=ds, "var"=est_var, "probe_smok"=probe_smok, "gene_smok"=gene_smok,
                  "probe_int"=probe_int, "gene_int"=gene_int)
  return(my_study)
}

my_studies <- c('GSE103174',
                'GSE13896','GSE16149','GSE17913','GSE18723','GSE19027',
                'GSE20189','GSE2125','GSE21862','GSE31210','GSE32539',
                'GSE42057','GSE42743','GSE4302','GSE44456','GSE46699',
                'GSE56768','GSE7895','GSE87072', 'GSE8987','GSE994')
all_studies <- lapply(my_studies, grab_study)
save(all_studies, file="data/all_studies.RData")

