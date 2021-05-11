# --- Leveraging Mice to Perform Imputation of Missing Covariates in Expression Data --- #

# Setup:
# 1. impute age, sex 
# 2. looking at smoking-related effects
# 3. smokers are on average older, more male smokers
#
# look at dataset with HDL!
#
# Evaluation:
# 1. comparison to results without imputation
# 2. comparison to expanded dataset vs complete cases
# 3. vs basic mean imputation with no gex covars
#
# TODO: figure out how to nest nicely...


library(mice)
library(tidyverse)
library(limma)
library(data.table)

# create an object:
#  expr, pdat, list_pvars, coef, my_vars


# id variables for imputation
id_vars <- function(expr, pdat, list_pvars){
  
}


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
  my_model <- paste("~", paste(list_pvars, collapse=" + "))
  comb_res <- lapply(1:pdat_imp$m, function(x) { 
    df <- run_tt(imp_data=complete(pdat_imp, x), 
                 expr=expr, coef_str=coef_str, my_model=my_model, ...); 
    return(df$probes)
  })
  return(unique(unlist(comb_res)))
}


extract_var <- function(pdat_imp, expr, list_pvars, coef_str,  out_var){
  # select top variables
  y = unlist(expr[out_var,]) 
  mipo <-pool(with(data=pdat_imp, 
                   lm(as.formula(paste("y", paste(list_pvars, collapse="+"), sep="~") ))))
  my_summ <- as.tibble(summary(mipo))
  comb_summ <- my_summ %>% 
    left_join(as.tibble(mipo$pooled) %>% dplyr::select(-estimate, -df), by="term")
  comb_summ %>%
    filter(term==coef_str) %>%
    dplyr::select(-term) %>%
    mutate(probe:={{out_var}})
}


extract_top_var <- function(pdat_imp, expr, list_pvars, coef_str, my_l){
  data.frame(rbindlist(lapply(my_l, function(x) 
                       extract_var(pdat_imp, expr, list_pvars, coef_str, x) )
  ))
}


extract_all_var <- function(pdat_imp, expr, list_pvars, coef_str){
  data.frame(rbindlist(lapply(rownames(expr), function(x) 
    extract_var(pdat_imp, expr, list_pvars, coef_str, x) )
  ))
}


# example code
load("data/ae_filtered_pheno.RData") # expDat3, pDat5
probes_age <-  c("205009_at", "205337_at", "208476_s_at", "229159_at", "233813_at")
exp4 <- expDat3[,pDat5$gsm]
pDat6 <- pDat5
set.seed(411)
train_idx <- sample(1:length(nrow(pDat6)), floor(0.6*length(pDat6)))
pDat6$age[-train_idx] <- NA

pdat_imp <- run_imp(expr=exp4, pdat=pDat5, 
                    my_vars=probes_age, m=4, defaultMethod="pmm")
list_pvars <- c("age", "sex", "smoking", "submission_date")
my_l <- extract_tt(pdat_imp=pdat_imp, expr=exp4, list_pvars=list_pvars, 
                   coef_str="sexm", number=5)

df_out <- extract_top_var(pdat_imp, exp4, list_pvars, "sexm", my_l)
rownames(df_out) <- df_out$probe
df_out %>% add_gene()

# NOTE: need to do a removal of variables used in the analysis!

