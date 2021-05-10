library('tidyverse')
library('glmnet')
library('mice')
library('limma')
library('groupdata2')

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


extract_var2 <- function(pdat_imp, y, my_mod,  out_var){
  # select top variables
  mipo <-pool(with(data=pdat_imp, 
                   lm(as.formula(paste("y", paste(colnames(pdat_imp[[1]]), collapse="+"), sep="~")))))
  return(summary(mipo, "all") %>% mutate(probe=out_var))
}



extract_all_var <- function(pdat_imp, expr, list_pvars){
  my_mod <- "y ~ ." #paste("y", paste(colnames(pdat_imp[[1]]), collapse="+"), sep="~") 
  lapply(rownames(expr), function(x) 
    extract_var2(pdat_imp, unlist(expr[x,]), "", x) )
}


run_fold <- function(fold_list, my.fold){
  nfolds = length(unique(fold_list$partition))
  train_data <- fold_list %>% filter(!partition %in% my.fold)
  valid_data <- fold_list %>% filter(partition %in% my.fold)
  
  # redo the training fold
  if (my.fold!=nfolds){
    train_data <- train_data %>% mutate(partition=case_when(
      partition==nfolds ~ my.fold,
      TRUE ~ partition
    ))  
  }
  
  X_train <- as.matrix(train_data %>% select(-class, -partition))
  Y_train <- train_data$class
  train_folds <- train_data$partition

  X_valid <- as.matrix(valid_data %>% select(-class, -partition))
  Y_valid <- valid_data$class


  test_params <- function(my.alpha){
    cvfit = cv.glmnet(X_train, Y_train, 
                      family="gaussian",
                      foldid=train_folds,
                      alpha=my.alpha, 
                      standardize=FALSE,
                      trace=FALSE)
    
    train_assess <- assess.glmnet(cvfit, newx=X_train, newy=Y_train, s="lambda.1se")
    valid_assess <- assess.glmnet(cvfit, newx=X_valid, newy=Y_valid, s="lambda.1se")
    my.lambda <- cvfit$lambda.1se
    
    ta <- unlist(train_assess) 
    va <- unlist(valid_assess) 
    train_valid_assess <- bind_rows(ta, va)
    train_valid_assess$grp <- c("train", "valid")
    train_valid_assess$alpha <- my.alpha
    train_valid_assess$lambda <- my.lambda
    
    # turn this into a dataframe
    # fold, my.alpha, lambda, train_assess, valid_assess
    return(train_valid_assess)
  }

  #res1 <- lapply(seq(0,1,0.1), function(my.alpha) test_params(my.alpha) )
  res1 <- test_params(0.5)
  return(res1)
}


select_probes <- function(pDat, exp, var, covars, max_probes, multinomial=F, ...){
  
  Y <- pDat %>% pull(var)
  mod_mat <- model.matrix(~.-1, pDat %>% dplyr::select({{covars}}))
  X <- as.matrix(cbind(exp, mod_mat))
  #print(nrow(X))
  
  # balanced(?) split in 5 folds
  indat <- data.frame(cbind(Y, exp, mod_mat)) %>% rename(class=V1)
  indat2 <- fold(indat, k=5) %>% rename(partition=.folds) %>%
    ungroup() %>%
    mutate(partition=as.integer(partition)) 
  res2 <- do.call(rbind, lapply(1:5, function(x) run_fold(indat2, x)))
  # train on 4 folds
  Y_train <- #Y[train_idx]
  X_train <- #X[train_idx,]
  
  # UPDATE TO NESTED CV + TRY FOR ALPHA?
  cvfit = cv.glmnet(X_train, 
                    Y_train, 
                    standardize=F, 
                    ...)
  
  # test on last fold
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



test_impute <- function(run_id, complete_pDat, complete_eDat, fraction_miss,  
                        vars_miss, max_probes, coef_str){
  
  nsamples <- nrow(complete_pDat)
  train_idx <- sample(1:nsamples, floor((1-fraction_miss)*nsamples))
  complete_pDat0 <- complete_pDat[train_idx,]
  complete_pDat1 <- complete_pDat[-train_idx,]
  complete_eDat0 <- complete_eDat[,train_idx]
  complete_eDat1 <- complete_eDat[,-train_idx]
  
  # select probes using complete0
  exp_rot <- t(complete_eDat0)

  age_probes <- select_probes( pDat=complete_pDat0 %>%
                                 mutate(age=scale(age)), 
                               exp=exp_rot, 
                               var="age", 
                               covars=c("sex", "smoking"), 
                               max_probes=max_probes)
  
  # add missingness to complete1
  list_probes <- age_probes
  complete_pDat1$age <- NA
  
  miss_pDat <- complete_pDat0 %>% bind_rows(complete_pDat1)
  miss_eDat <- complete_eDat[,miss_pDat$geo_accession]
  
  # run DE w mice
  pdat_imp <- run_imp(expr=miss_eDat, pdat=miss_pDat %>% select(-geo_accession), 
                      my_vars=list_probes, m=5, defaultMethod="pmm")
  start_time = Sys.time()
  df_out <- extract_all_var(pdat_imp, miss_eDat, c(""))
  end_time=Sys.time()
  print(end_time-start_time)
  save(df_out, file=sprintf("data/de_imp_run_%s_%s.RData", run_id, fraction_miss))
  
  imp_out <- data.table::rbindlist(lapply(df_out, function(x) x %>% filter(str_detect(term,"smok"))))
  
  # add pval + adj-pval
  tt_imp <- imp_out %>% filter(term==coef_str) %>%
    mutate(P.Value=2*pt(-abs(statistic),df=df)) %>%
    arrange(desc(P.Value)) %>%
    mutate(adj.P.Val=p.adjust(P.Value, "fdr")) %>%
    rename(probes=probe, logFC=estimate)
  
  # comparisons  
  tt_true <- run_tt(complete_pDat %>% select(-geo_accession), 
                    complete_eDat, coef_str, my_model="~ .",  
                    number=nrow(complete_eDat))
  
  tt_sm <- run_tt(complete_pDat0 %>% select(-geo_accession), 
                  complete_eDat0, coef_str, my_model="~ .",
                  number=nrow(complete_eDat))
  
  # cutoffs
  calc_out <- function(pcut){
    imp_cut <- tt_imp %>% filter(adj.P.Val < pcut)
    true_cut <- tt_true %>% filter(adj.P.Val < pcut)
    sm_cut <- tt_sm %>% filter(adj.P.Val < pcut)
    
    # counts
    num_true = nrow(true_cut)
    num_imp = nrow(imp_cut)
    num_sm = nrow(sm_cut)
    true_imp_overlap = length(intersect(true_cut$probes, imp_cut$probes))
    true_sm_overlap = length(intersect(true_cut$probes, sm_cut$probes))
    imp_sm_overlap = length(intersect(imp_cut$probes, sm_cut$probes))
    
    # calculate stats
    true_imp_df = true_cut %>% select(probes, logFC, P.Value) %>% 
      inner_join(imp_cut  %>% select(probes, logFC, P.Value), by="probes") 
    true_sm_df = true_cut  %>% select(probes, logFC, P.Value) %>% 
      inner_join(sm_cut  %>% select(probes, logFC, P.Value), by="probes") 
    
    est_cor_imp = cor.test(true_imp_df$logFC.x, true_imp_df$logFC.y)
    est_cor_sm = cor.test(true_sm_df$logFC.x, true_sm_df$logFC.y)
    p_cor_impS =cor.test(true_imp_df$P.Value.x, true_imp_df$P.Value.y, method="spearman")
    p_cor_impP =cor.test(log10(true_imp_df$P.Value.x), log10(true_imp_df$P.Value.y))
    p_cor_smS =cor.test(true_sm_df$P.Value.x, true_sm_df$P.Value.y, method="spearman")
    p_cor_smP =cor.test(log10(true_sm_df$P.Value.x), log10(true_sm_df$P.Value.y))
    return(list(
      "run_id"= run_id,
      "frac_miss"= fraction_miss,
      "cutoff" = pcut,
      "num_true" = num_true,
      "num_imp" = num_imp,
      "num_sm" = num_sm,
      "true_imp_overlap" = true_imp_overlap,
      "imp_sm_overlap" = true_sm_overlap,
      "imp_sm_overlap" = imp_sm_overlap,
      "est_cor_imp"= est_cor_imp$estimate[[1]],
      "est_cor_imp_p" =est_cor_imp$p.value,
      "est_cor_sm"= est_cor_sm$estimate[[1]],
      "est_cor_sm_p" =est_cor_sm$p.value,
      "p_cor_impS"= p_cor_impS$estimate[[1]],
      "p_cor_impS_p" =p_cor_impS$p.value,
      "p_cor_impP"= p_cor_impP$estimate[[1]],
      "p_cor_impP_p" =p_cor_impP$p.value,
      "p_cor_smS"= p_cor_smS$estimate[[1]],
      "p_cor_smS_p" =p_cor_smS$p.value,
      "p_cor_smP"= p_cor_smP$estimate[[1]],
      "p_cor_smP_p" =p_cor_smP$p.value
    ))
  }
  out1 <- calc_out(0.05)
  out2 <- calc_out(0.01)
  return(rbind(out1, out2))
}

args <- commandArgs(trailingOnly=T)
run_idx <- as.numeric( args[1])
fraction_miss <- as.numeric(args[2])


load("data/ae_full_exp.RData") # expDat5, pDat5.1
complete_pDat <- pDat5.1 %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years)) # 244
complete_eDat <- expDat5[,complete_pDat$geo_accession]

run_id <- sprintf("%s_age", run_idx)
df <- test_impute(run_id, complete_pDat, complete_eDat, fraction_miss,  
                  c("age"), max_probes, "smokingS")
save(df, file=sprintf("data/imp_test_%s.RData", run_id, fraction_miss))
