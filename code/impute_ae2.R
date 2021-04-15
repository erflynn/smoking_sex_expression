# impute AE

load("data/ae_full_exp.RData") # expDat4, pDat4
library('glmnet')

# use the completed data to select variables for the remainder
pDat5 <- pDat4 %>% 
  mutate(age=as.numeric(age),
         pack_years=as.numeric(pack_years),
         pack_years=ifelse(is.na(pack_years) & smoking=="NS", 0, pack_years)) %>%
  mutate(submission_date=as.factor(submission_date))

complete_pDat <- pDat5 %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years))

# - age
# - pack_years
# - race_ethnicity
fct_summ(complete_pDat)

# glmnet
exp_rot <- t(expDat4[,complete_pDat$geo_accession])

select_probes <- function(pDat, exp, var, covars, ...){

  Y <- pDat %>% pull(var)
  mod_mat <- model.matrix(~.-1, pDat %>% 
                            select({{covars}}))
  X <- as.matrix(cbind(exp, mod_mat))
  #train_idx <- sample(1:length(Y), floor(0.8*length(Y)))
  Y_train <- Y[train_idx]
  X_train <- X[train_idx,]
  Y_test <- Y[-train_idx]
  X_test <- X[-train_idx,]
  cvfit = cv.glmnet(X_train, 
                    Y_train, ...)
  #train_assess <- assess.glmnet(cvfit, newx=X_train, 
  #                              newy=Y_train, s="lambda.1se")
  #valid_pred <- predict(cvfit, newx=X_test, 
  #                      type="response", s="lambda.1se")
  #train_pred <- predict(cvfit, newx=X_train, 
  #                      type="response", s="lambda.1se")
  #print(cor(valid_pred, Y_test)) # 0.32
  my.lambda <- "lambda.1se"
  coef_l <- coef(cvfit, s=my.lambda)
  if (length(coef_l) > 1){
    cv_coef <- do.call(rbind, coef_l)
  }
  mat_coef <- cv_coef %>% as.matrix()
  nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
  print(nonzero_coef)
  my_probes <- names(nonzero_coef)[str_detect(names(nonzero_coef), "_at")]
  return(my_probes)
}

set.seed(414)

age_probes <- select_probes( pDat=complete_pDat %>%
                               mutate(age=scale(age)), 
                             exp=exp_rot, 
                             var="age", 
                             covars=c("sex", "smoking", "race_ethnicity", "pack_years", "submission_date"),  
                             alpha=0.5, 
                             family="gaussian")

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
race_eth2 <- complete_pDat %>%
  filter(race_ethnicity %in% c("hispanic", "black", "white")) %>%
  mutate(race_ethnicity=as.factor(race_ethnicity))
re_probes <- select_probes( pDat=race_eth2, exp=exp_rot[race_eth2$geo_accession,], 
                            var="race_ethnicity", 
                               covars=c("sex", "smoking", "pack_years", 
                                        "age", "submission_date"),  
                               alpha=0.5, 
                               family="multinomial",
                              type.multinomial="grouped")
 
# var selection with stepwise

# or just select top 5?

# ------ RUN IMPUTATION ------ #

pdat_imp <- run_imp(expr=expDat4, pdat=pDat5, 
                    my_vars=c(age_probes, re_probes), m=4, defaultMethod="pmm")
list_pvars <- c("age", "sex", "smoking", "submission_date", "smoking:sex")
my_l <- extract_tt(pdat_imp=pdat_imp, expr=expDat4, list_pvars=list_pvars, 
                   coef_str="smokingS", number=10)

df_out <- extract_top_var(pdat_imp, expDat4, list_pvars, "smokingS", my_l)
rownames(df_out) <- df_out$probe
df_out %>% add_gene() 

# ----- COMPARE RESULTS WITH pDat_complete vs imputed dat ----- #
mod <-  paste("~", paste(list_pvars, collapse=" + "))
tt_complete <- run_tt(complete_pDat, expDat4[,complete_pDat$geo_accession], "smokingS", mod)
tt1 <- tt_complete %>% add_gene()
tt1


