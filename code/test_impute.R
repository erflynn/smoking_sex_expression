

#####
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

# do with correlation/chisq
select_probes2 <- function(pDat, eDat, max_probes, var){
  
  design2 <- model.matrix(~ .-1, data=pDat)
  fit <- lmFit(eDat, design2)
  fit <- eBayes(fit)
  coefs <- colnames(coef(fit))
  list_probes <- lapply(coefs[str_detect(coefs, var)], 
         function(coef_str){
           tt <- data.frame(topTable(fit=fit, coef=coef_str, number=max_probes))
           return(rownames(tt))
         })
  return(unique(unlist(list_probes)))
}
  #pDat_cov <- complete_pDat0$age
  # if (var_type=="numeric"){
  #   cor_list <- apply(complete_eDat0, 1, function(x) cor(x, pDat_cov))
  #   cor_srt <- sort(abs(cor_list), decreasing=T)
  # } 
  # if (var_type=="categorical"){ # point-biserial
  #   
  # }

  #return(names(cor_srt[1:max_probes]))
#}



####

# # train/test split
# set.seed(424)
# nsamples <- nrow(complete_pDat)
# train_idx <- sample(1:nsamples, floor(0.7*nsamples))
# complete_pDat0 <- complete_pDat[train_idx,]
# complete_pDat1 <- complete_pDat[-train_idx,]
# complete_eDat0 <- complete_eDat[,train_idx]
# complete_eDat1 <- complete_eDat[,-train_idx]
# 
# # select probes using complete0
# exp_rot <- t(complete_eDat0)
# age_probes <- select_probes( pDat=complete_pDat0 %>%
#                                mutate(age=scale(age)), 
#                              exp=exp_rot, 
#                              var="age", 
#                              covars=c("sex", "smoking", 
#                                       "race_ethnicity", "pack_years", "submission_date"), 
#                              alpha=0.5, 
#                              family="gaussian")
# length(age_probes) # 10
# 
# 
# 
# 
# # add missingness to complete_pDat1
# complete_pDat1$age <- NA
# miss_pDat <- complete_pDat0 %>% bind_rows(complete_pDat1)
# miss_eDat <- complete_eDat[,miss_pDat$geo_accession]
# 
# # run DE w mice
# pdat_imp <- run_imp(expr=miss_eDat, pdat=miss_pDat %>% select(-geo_accession), 
#                     my_vars=c(age_probes), m=5, defaultMethod="pmm")
# comb_res <- lapply(1:pdat_imp$m, function(x) { 
#   df <- run_tt(imp_data=complete(pdat_imp, x), 
#                expr=miss_eDat, coef_str="smokingS", my_model="~ ."); 
#   return(df)
# })
# head(comb_res[[1]])
# tt_out <- extract_top_var(pdat_imp, miss_eDat, c(""), "smokingS", my_l) 
# 
# 
# # fix the p-val calculation + adjust
# calc_p <- function(t,n){ 2*pt(-abs(t),df=n-1)}
# tt_out$p <- sapply(tt_out$statistic, function(x) calc_p(x, 205))
# head(tt_out)
# # summary(mipo, conf.int=T, type="all") # -- this is useful
# 
# # compare to running DE with whole complete dataset
# tt_true <- run_tt(complete_pDat %>% select(-geo_accession), 
#                   complete_eDat, coef_str="smokingS", my_model="~ .", confint=TRUE, number=nrow(complete_eDat))
# head(tt_true)
# 
# # and only the complete cases
# tt_sm <- run_tt(complete_pDat0 %>% select(-geo_accession), complete_eDat0, 
#                 coef_str="smokingS", my_model="~ .", confint=TRUE)
# head(tt_sm)
# 
# # compare effect sizes
# # compare set of genes
# imp_probes <- tt_out$probe
# true_probes <- head(tt_true$probes, 10)
# sm_probes <- head(tt_sm$probes, 10)
# intersect(imp_probes, true_probes) # 10 
# length(intersect(sm_probes, true_probes)) # 8
# 
# joint_tt <- tt_true %>% select(probes, logFC, P.Value) %>% 
#   inner_join(tt_out %>% select(probe, estimate, p), by=c("probes"="probe"))
# cor.test(joint_tt$logFC, joint_tt$estimate) # 0.994
# cor.test(joint_tt$p, joint_tt$P.Value, method="spearman") # 0.745
# cor.test(log10(joint_tt$p), log10(joint_tt$P.Value)) # 0.708
# 
# 
# ggplot(joint_tt,
#        aes(x=logFC, y=estimate))+
#   geom_point()
# 
# tt_sm_compare <- tt_true %>% 
#   select(probes, P.Value) %>% 
#   inner_join(tt_sm %>% select(probes, P.Value), by="probes")
# cor.test(tt_sm_compare$P.Value.x, tt_sm_compare$P.Value.y, method="spearman") # 0.642
# cor.test(log10(tt_sm_compare$P.Value.x), log10(tt_sm_compare$P.Value.y)) # 0.807
# 
# ggplot(tt_sm_compare,
#        aes(x=logFC.x, logFC.y))+
#   geom_point()
# 


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
  list_probes <- c()
  # add missingness to complete1
  num_miss <- length(vars_miss)
  
  if ("age" %in% vars_miss){
    age_probes <- select_probes( pDat=complete_pDat0 %>%
                                   mutate(age=scale(age)), 
                                 exp=exp_rot, 
                                 var="age", 
                                 covars=c("sex", "smoking", 
                                          "race_ethnicity", "pack_years", "submission_date"), 
                                 max_probes=15,
                                 alpha=0.5, 
                                 family="gaussian")
    list_probes <- c(list_probes, age_probes)
    miss_idx <- sample(1:nrow(complete_pDat1), floor(1/num_miss))
    complete_pDat1$age[miss_idx] <- NA
  }
  if ("pack_years" %in% vars_miss){
    smokers <- complete_pDat0 %>% 
      filter(smoking=="S")
    pkyrs_probes <- select_probes( pDat=smokers %>%
                                     mutate(pack_years=scale(pack_years)), 
                                   exp=exp_rot[smokers$geo_accession,],
                                   var="pack_years", 
                                   covars=c("sex",  "race_ethnicity", "age", 
                                            "submission_date"),  
                                   max_probes=15,
                                   alpha=0.5, 
                                   family="gaussian")
    list_probes <- c(list_probes, pkyrs_probes)
    miss_idx <- sample(1:nrow(complete_pDat1), floor(1/num_miss))
    complete_pDat1$pack_years[miss_idx] <- NA
  }

  if ("race_ethnicity" %in% vars_miss){
    race_eth2 <- complete_pDat0 %>%
      filter(race_ethnicity %in% c("hispanic", "black", "white")) %>%
      mutate(race_ethnicity=as.factor(race_ethnicity))
    re_probes <- select_probes( pDat=race_eth2, 
                                exp=exp_rot[race_eth2$geo_accession,], 
                                var="race_ethnicity", 
                                covars=c("sex", "smoking", "pack_years", 
                                         "age", "submission_date"), 
                                max_probes=15,
                                multinomial=T,
                                alpha=0.5, 
                                family="multinomial",
                                type.multinomial="grouped")
    list_probes <- c(list_probes, re_probes)
    miss_idx <- sample(1:nrow(complete_pDat1), floor(1/num_miss))
    complete_pDat1$race_ethnicity[miss_idx] <- NA
  }

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



# # simplify submission date
# adj_date = pDat5 %>%
#   group_by(submission_date) %>%
#   count() %>%
#   ungroup() %>%
#   mutate(submission_date2=lubridate::ymd(submission_date)) %>%
#   arrange(submission_date2) %>%
#   mutate(n_prev=lag(n),
#          n_next=lead(n)) %>%
#   mutate(
#     prevs=case_when(
#       n_prev > 9 ~ lag(submission_date2),
#       lag(n,2) > 9 ~ lag(submission_date2,2),
#       lag(n,3) > 9 ~ lag(submission_date2, 3)
#     ),
#     nexts=case_when(
#       n_next > 9 ~ lead(submission_date2),
#       lead(n,2) > 9 ~ lead(submission_date2,2),
#       lead(n,3) > 9 ~ lead(submission_date2, 3)
#     )) %>%
#   mutate( prev_d=as.numeric(submission_date2 - prevs),
#           next_d=as.numeric(nexts-submission_date2)) %>%
#   mutate(new_date=case_when(
#     n > 9  ~ submission_date2,
#     is.na(prev_d)  ~ nexts,
#     is.na(next_d) ~ prevs,
#     prev_d < next_d  ~ prevs,
#     next_d < prev_d  ~ nexts
#   ))
# date_convert <- adj_date%>%
#   select(submission_date, new_date)
# 
# pDat5.1 <- pDat5 %>% 
#   left_join(date_convert, by="submission_date") %>%
#   select(-submission_date) %>%
#   rename(submission_date=new_date) %>%
#   mutate(age=as.numeric(age),
#          pack_years=as.numeric(pack_years),
#          pack_years=ifelse(smoking=="NS", 0, pack_years)) %>%
#   mutate(submission_date=as.factor(submission_date)) %>%
#   group_by(race_ethnicity) %>%
#   mutate(n=n()) %>%
#   ungroup() %>%
#   mutate(race_ethnicity=ifelse(n<5, "other", race_ethnicity)) %>%
#   select(geo_accession, sex, age, smoking, race_ethnicity, pack_years, submission_date) # 353
# save(pDat5.1, expDat5, pDat5, file="data/ae_full_exp.RData") 
### --- cleaning finished above --- #

args <- commandArgs(trailingOnly=T)
run_idx <- as.numeric( args[1])
load("data/ae_full_exp.RData") # expDat5, pDat5.1
complete_pDat <- pDat5.1 %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years)) # 244
complete_eDat <- expDat5[,complete_pDat$geo_accession]

pDat5.1 %>% filter(is.na(age) | is.na(race_ethnicity)) %>%
  fct_summ()

pDat5.1 %>% filter(!is.na(race_ethnicity) & is.na(pack_years))

fraction_miss <- 0.3
run_id <- sprintf("%s_age", run_idx)
df <- test_impute(run_id, complete_pDat, complete_eDat, fraction_miss,  
                        c("age"), 10, "smokingS")
save(df, file=sprintf("data/imp_test_%s.RData", run_id, fraction_miss))

run_id2 <- sprintf("%s_age_re", run_idx)
df2 <- test_impute(run_id2, complete_pDat, complete_eDat, fraction_miss,  
                  c("age", "race_ethnicity"), 10, "smokingS")
save(df2, file=sprintf("data/imp_test_%s.RData", run_id2, fraction_miss))


run_id3 <- sprintf("%s_age_re_pk", run_idx)
df3 <- test_impute(run_id3, complete_pDat, complete_eDat, fraction_miss,  
                   c("age", "race_ethnicity"), 10, "smokingS")
save(df3, file=sprintf("data/imp_test_%s.RData", run_id3, fraction_miss))

# TO SHOW THAT THIS WORKS:
# - fix issues: 
#    [x]smaller number of probes per category:
#      order by coefficient (other: chisq, IG, correlation?)
#    [x] multiple different values missing
#   [x] too many submission dates -- group together small ones? (if <5 or 10 --> other?)
#      -OR cluster?
#      - group together w nearest date: if < =5 --> add to nearest date
#   [x] possibly recode RE so that we can have ppl coded w multiple (or put 5 ppl in other?)
#    [x] pvalue extraction, calculate adjusted
#
# - vary range of missing values 0.1 --> 0.4 (what is it in the actual dataset?)
#     - in this case: 69.1% of the data are complete cases, 29.2(RE), 30.3(age), 35.9(PY of smok)
# - stick in a randomized loop (mb 10 random runs?)

# submission date clustering?
# try a test example
# get it to work on the server


# all_f = lapply(head(my_f), function(f) {load.Rdata(f, "my_df"); return(data.frame(my_df))})
# f_df = do.call(rbind, all_f)
# f_df2 = data.frame(apply(f_df, c(1,2), unlist))
# f_df2 %>% write_csv("summary_imp.csv")
