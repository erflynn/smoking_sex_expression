library('tidyverse')
library('glmnet')
library('mice')
library('limma')

run_imp <- function(expr, pdat, my_vars, ...){
  expr_pdat <- data.frame(cbind(t(expr[my_vars,]), pdat))
  return(mice(expr_pdat, ...) )
}

run_tt <- function(imp_data, expr, coef_str, my_model, ...){
  design2 <- model.matrix(as.formula(my_model), data=imp_data)
  fit <- lmFit(expr, design2)
  fit <- eBayes(fit)
  tt <- data.frame(topTable(fit=fit, coef=coef_str, ...))
  tt$probes <- rownames(tt)
  return(tt)
}

extract_tt <- function(pdat_imp, expr, list_pvars, coef_str, ...){
  #my_model <- paste("~", paste(list_pvars, collapse=" + "))
  my_model <- "~ ."
  comb_res <- lapply(1:pdat_imp$m, function(x) { 
    df <- run_tt(imp_data=complete(pdat_imp, x), 
                 expr=expr, coef_str=coef_str, my_model=my_model, ...); 
    return(df$probes)
  })
  return(unique(unlist(comb_res)))
}


extract_var <- function(pdat_imp, expr, list_pvars, coef_str,  out_var){
  # select top variables
 # my_mod <-  #as.formula(paste("y", paste(list_pvars, collapse="+"), sep="~") )
  y = unlist(expr[out_var,]) 
  mipo <-pool(with(data=pdat_imp, 
                   lm(as.formula(paste("y", paste(colnames(pdat_imp[[1]]), collapse="+"), sep="~")))))
  my_summ <- as.tibble(summary(mipo))
  comb_summ <- my_summ %>% 
    left_join(as.tibble(mipo$pooled) %>% dplyr::select(-estimate, -df), by="term")
  comb_summ %>%
    filter(term==coef_str) %>%
    dplyr::select(-term) %>%
    mutate(probe:={{out_var}})
}

extract_var2 <- function(pdat_imp, y, my_mod,  out_var){
  # select top variables
  mipo <-pool(with(data=pdat_imp, 
                   lm(as.formula(paste("y", paste(colnames(pdat_imp[[1]]), collapse="+"), sep="~")))))
  return(summary(mipo, "all") %>% mutate(probe=out_var))
}


extract_top_var <- function(pdat_imp, expr, list_pvars, coef_str, my_l){
  data.frame(data.table::rbindlist(lapply(my_l, function(x) 
    extract_var(pdat_imp, expr, list_pvars, coef_str, x) )
  ))
}


extract_all_var <- function(pdat_imp, expr, list_pvars){
  my_mod <- "y ~ ." #paste("y", paste(colnames(pdat_imp[[1]]), collapse="+"), sep="~") 
  lapply(rownames(expr), function(x) 
    extract_var2(pdat_imp, unlist(expr[x,]), "", x) )
}

load("data/ae_full_exp.RData") # expDat5, pDat5


# use the completed data to select variables for the remainder
pDat5.1 <- pDat5 %>% 
  mutate(age=as.numeric(age),
         pack_years=as.numeric(pack_years),
         pack_years=ifelse(smoking=="NS", 0, pack_years)) %>%
  mutate(submission_date=as.factor(submission_date)) %>%
  select(geo_accession, sex, age, smoking, race_ethnicity, pack_years, submission_date)

complete_pDat <- pDat5.1 %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years)) # 244

exp_rot <- t(expDat5[,complete_pDat$geo_accession])


select_probes <- function(pDat, exp, var, covars, max_probes, multinomial=F, ...){
  
  Y <- pDat %>% pull(var)
  mod_mat <- model.matrix(~.-1, pDat %>% dplyr::select({{covars}}))
  X <- as.matrix(cbind(exp, mod_mat))
  train_idx <- 1:length(Y) #sample(1:length(Y), floor(0.8*length(Y)))
  Y_train <- Y[train_idx]
  X_train <- X[train_idx,]

  cvfit = cv.glmnet(X_train, 
                    Y_train, ...)

  coef_l <- coef(cvfit, s="lambda.1se")
  if (multinomial==T){
    coef_l <- do.call(rbind, coef_l)
  }
  mat_coef <- coef_l %>% as.matrix()
  nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
  probes <- nonzero_coef[str_detect(names(nonzero_coef), "_at")]
  print(nonzero_coef)
  if (length(probes) > max_probes){
    probes = sort(desc(abs(probes)))[1:max_probes]
  }
  return(names(probes))
}

set.seed(414)

age_probes <- select_probes( pDat=complete_pDat %>%
                               mutate(age=scale(age)), 
                             exp=exp_rot, 
                             var="age", 
                             covars=c("sex", "smoking", 
                                      "race_ethnicity", "pack_years", "submission_date"), 
                             alpha=0.5, 
                             family="gaussian")
length(age_probes) # 18
smokers <- complete_pDat %>% 
  filter(smoking=="S")
pkyrs_probes <- select_probes( pDat=smokers %>%
                                 mutate(pack_years=scale(pack_years)), 
                               exp=exp_rot[smokers$geo_accession,],
                               var="pack_years", 
                               covars=c("sex",  "race_ethnicity", "age", 
                                        "submission_date"),  
                               alpha=0.5, 
                               family="gaussian")
# none
race_eth2 <- complete_pDat %>%
  filter(race_ethnicity %in% c("hispanic", "black", "white")) %>%
  mutate(race_ethnicity=as.factor(race_ethnicity))
re_probes <- select_probes( pDat=race_eth2, exp=exp_rot[race_eth2$geo_accession,], 
                            var="race_ethnicity", 
                            covars=c("sex", "smoking", "pack_years", 
                                     "age", "submission_date"), 
                            multinomial=T,
                            alpha=0.5, 
                            family="multinomial",
                            type.multinomial="grouped")

length(re_probes) # 312
save(age_probes, re_probes, file="data/ae_imp_probes.RData")

# ------ RUN IMPUTATION ------ #
set.seed(415)
pdat_imp <- run_imp(expr=expDat5, pdat=pDat5.1, 
                    my_vars=c(age_probes, re_probes), m=5, defaultMethod="pmm")
list_pvars <- c("age", "sex", "smoking", "submission_date", "smoking:sex")
my_l <- extract_tt(pdat_imp=pdat_imp, expr=expDat5, list_pvars=list_pvars, 
                   coef_str="smokingS", number=10)

print("about to run")

start_time = Sys.time()
df_out <- extract_all_var(pdat_imp, expDat5, list_pvars)
end_time=Sys.time()
end_time-start_time
save(df_out, file="data/ae_imputed.RData")


### FOR AFTER
library('tidyverse')
library('mice')
load("data/ae_imputed.RData")

start_time = Sys.time()
out_l <- lapply(names(df_out), function(my_probe){
  mipo <- df_out[[my_probe]]
  summary(mipo) %>% 
  as_tibble() %>%
  left_join(mipo$pooled %>% as_tibble() %>% dplyr::select(-estimate, -df), by="term") %>%
  filter(str_detect(term, "smok") | str_detect(term, "sex")) %>%
  mutate(probe=my_probe)})
end_time=Sys.time()
end_time-start_time
df <- data.table::rbindlist(out_l)
save(df, file="data/ae_imputed_df.RData")


## meta analysis
df2 <- df %>% as_tibble() %>% filter(term=="sexmale:smokingS")
df3 <- df2 %>% rename(probes=probe) %>% 
  add_gene() %>% arrange(p.value) %>%
  mutate(p.adj=p.adjust(p.value, "fdr")) %>%
  mutate(CI.L=estimate - 1.96*std.error)
ma_imp <- run_clean_ma(df3 %>% 
                         rename(logFC=estimate, ID=probes, geneSymbol=gene, P.Value=p.value) %>% 
                         mutate(chromosome=""))

df2_s <- df %>% as_tibble() %>% filter(term=="smokingS")
df3_s <- df2_s %>% rename(probes=probe) %>% 
  add_gene() %>% arrange(p.value) %>%
  mutate(p.adj=p.adjust(p.value, "fdr")) %>%
  mutate(CI.L=estimate - 1.96*std.error) # 74 have p==0 :/ 

ma_imp_smok <- run_clean_ma(df3_s %>% 
                         rename(logFC=estimate, ID=probes, geneSymbol=gene, P.Value=p.value) %>% 
                         mutate(chromosome=""))
save(ma_imp, ma_imp_smok, file="data/imp_ae_res.RData")
#  estimate  std.error 
load(file="data/results/ae_smok.RData") # --> ma_smok, ma_sex, smok_ci_smok, smok_ci_sex,
head(ma_smok )
disc_s_g <- ma_smok %>% filter(adj.p < 0.05)  
imp_s_g <- ma_imp_smok %>% filter(adj.p < 0.05)
disc_s_g %>% inner_join(ma_imp_smok, by="gene") %>%
  ggplot(aes(logFC.x, logFC.y))+
  geom_point()

imp_s_g %>% inner_join(ma_smok, by="gene") %>%
  ggplot(aes(logFC.x, logFC.y))+
  geom_point()

imp_s_g %>% inner_join(ma_smok, by="gene") %>%
  ggplot(aes(log10(p.x), log10(p.y)))+
  geom_point()

disc_s_g %>% inner_join(ma_imp_smok, by="gene") %>%
  ggplot(aes(log10(p.x), log10(p.y)))+
  geom_point()

ma_imp_smok %>%
  filter(p < 10**-100) %>%
  View()
# ... there are 12 probes with overinflated p-values

load("data/ae_imp_probes.RData") # age_probes, re_probes
df3_s %>% filter(probes %in% age_probes |
                   probes %in% re_probes)
df3_s %>% filter(p.adj ==0)


ma_imp %>% filter(adj.p < 0.05) %>%
  inner_join(comb_ma_dat, by="gene") %>%
  ggplot(aes(log10(p.x), log10(p.y)))+
  geom_point()# 42 

ma_imp %>% filter(adj.p < 0.05) %>%
  inner_join(comb_ma_dat, by="gene") %>%
  ggplot(aes(logFC.x, logFC.y))+
  geom_point()# 42 



comb_ma_dat %>% filter(adj.p < 0.05) %>%
  inner_join(ma_imp, by="gene") %>%
  ggplot(aes(log10(p.x), log10(p.y)))+
  geom_point()# 72

comb_ma_dat %>% filter(adj.p < 0.05) %>%
  inner_join(ma_imp, by="gene") %>%
  ggplot(aes(logFC.x, logFC.y))+
  geom_point()# 72 
