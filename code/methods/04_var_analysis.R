
library(tidyverse)
library('variancePartition')

varPartPC <- function(expr, cutoff=0.8){
  expr_pcs <- prcomp(t(expr))
  last_pc <- max(which(summary(expr_pcs)$importance["Cumulative Proportion",] < cutoff))
  
  expr_pcs2 <- expr_pcs$x[,1:last_pc]
  props <- summary(expr_pcs)$importance["Proportion of Variance", 1:last_pc]
  return(list("pcs"=expr_pcs2, "props"=props))
}

get_pvca <- function(expr_pcs_t, phe, form, props, cutoff=0.8){
  varPart <- fitExtractVarPartModel(expr_pcs_t, form, phe)
  return( colSums(varPart*props)/cutoff)
}

rand_pvca1 <- function(expr_pcs_t, phe, form, props, cutoff=0.8){
  phe2 <- phe
  phe2$sex <- sample(phe$sex, nrow(phe2), replace=F)
  phe2$smoking <- sample(phe$smoking, nrow(phe2), replace=F)
  return(get_pvca(expr_pcs_t, phe2, form, props, cutoff))
}

rand_pvca <- function(n, expr_pcs_t, phe, form, props, cutoff=0.8){
  res <- lapply(1:n,  function(i) 
    rand_pvca1(expr_pcs_t, phe, form, props, cutoff))
  df <- do.call(cbind, lapply(res[!is.na(res)], data.frame)) 
  data.frame(t(df) )  %>% 
    as_tibble() %>%
    mutate(n=1:n()) %>%
    pivot_longer(-n, names_to="covariate", values_to="variance") %>%
    select(-n) %>%
    mutate(src="random") %>%
    mutate(covariate=str_replace_all(covariate, "\\.", ":"))
}

plot_pvca_rand <- function(est_var, rand_var){
  tibble("covariate"=names(est_var),
         "variance"=unlist(est_var)) %>%
    mutate(src="estimated") %>%
    bind_rows(rand_var) %>%
    filter(covariate!="Residuals") %>%
    ggplot(aes(x=covariate, y=variance, col=src))+
    geom_boxplot(position=position_dodge(0), alpha=0.5)+
    theme_bw()+
    ylab("fraction of variance")+
    xlab("")+
    scale_color_manual(values=c("red", "black"))+
    theme(legend.position = "None")
}


# ---- apply to AE data ---- #
load("data/ae_full_exp.RData") # pDat5, expDat5
complete_pdat <- pDat5 %>% 
  select(geo_accession, sex, age, smoking, pack_years, 
         race_ethnicity, submission_date) %>% 
  mutate(across(c(age, pack_years), as.numeric)) %>%
  mutate(pack_years=ifelse(smoking=="NS", 0, pack_years)) %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years))
complete_exp <- expDat5[,complete_pdat$geo_accession]

pcs_l <- varPartPC(complete_exp)
expr_pcs_t <- t(pcs_l$pcs)
props <- pcs_l$props

my_model <- as.formula("~ age + pack_years + race_ethnicity + sex + smoking + sex:smoking")
(est_var <- get_pvca(expr_pcs_t, complete_pdat, my_model, props))
set.seed(1011)
rand_var <- rand_pvca(1000, expr_pcs_t, complete_pdat, my_model, props)
plot_pvca_rand(est_var, rand_var)
ggsave("figures/paper_figs/s_fig_pvca_ae_int.png")

est_var[["sex:smoking"]] # 0.0066
rand_var %>% 
  filter(covariate=="sex:smoking") %>%
  filter(variance > est_var[["sex:smoking"]]) %>%
  nrow() # 86/1000 = 0.086

my_model2 <- as.formula("~ age + race_ethnicity+ pack_years + sex + smoking")
(est_var2 <- get_pvca(expr_pcs_t, complete_pdat, my_model2, props))
set.seed(1012)

rand_var2 <- rand_pvca(1000, expr_pcs_t, complete_pdat, my_model2, props)
plot_pvca_rand(est_var2, rand_var2)
ggsave("figures/paper_figs/s_fig_pvca_ae_smok.png")

save(est_var, rand_var, est_var2, rand_var2, file="data/results/ae_pvca.RData")

# ---- apply to blood data ---- #
load("data/blood_in_data.RData") # --> eDat3, pDat3


pcs_l <- varPartPC(eDat3)
expr_pcs_t <- t(pcs_l$pcs)
props <- pcs_l$props

my_model <- as.formula("~ age + sex + smoking + sex:smoking")
(est_var <- get_pvca(expr_pcs_t, pDat3, my_model, props))
set.seed(1011)
rand_var <- rand_pvca(29, expr_pcs_t, pDat3, my_model, props)
plot_pvca_rand(est_var, rand_var)
ggsave("figures/paper_figs/s_fig_pvca_blood_int.png")

est_var[["sex:smoking"]] 
rand_var %>% 
  filter(covariate=="sex:smoking") %>%
  filter(variance > est_var[["sex:smoking"]]) %>%
  nrow() 

# my_model2 <- as.formula("~ age + sex + smoking")
# (est_var2 <- get_pvca(expr_pcs_t, pDat3, my_model2, props))
# set.seed(1012)
# 
# rand_var2 <- rand_pvca(1000, expr_pcs_t, pDat3, my_model2, props)
# plot_pvca_rand(est_var2, rand_var2)
# ggsave("figures/paper_figs/s_fig_pvca_blood_smok.png")

save(est_var, rand_var, file="data/results/blood_pvca.RData")
load("data/results/ae_pvca.RData")
