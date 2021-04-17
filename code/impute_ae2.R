# impute AE

load("data/ae_full_exp.RData") # expDat5, pDat5
library('glmnet')

# PC plot
library('bigpca')
start_time = Sys.time(); 
ae_bp1 <- big.PCA(expDat5, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time # 40s

ae_pcs2 <- data.frame(ae_bp1$PCs)
ae_pcs2$submission_date <- pDat5$submission_date
ae_pcs2$geo_accession <- pDat5$geo_accession


ggplot(ae_pcs2, 
       aes(x=PC1, y=PC2, col=submission_date))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/ae_date_pc.png")
ggplot(blood_pcs2, 
       aes(x=PC1, y=PC2, col=hdl))+
  geom_point(alpha=0.5)+
  theme_bw()


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

# 
# # try stepwise regression
# library(caret)
# library(leaps)
# library(MASS)
# 
# 
# cor.dat <- apply(exp_rot, 2, function(x) cor(x, complete_pDat$age))
# ordered_cor <- cor.dat[order(-abs(cor.dat))]
# keep_probes <- names(head(ordered_cor, 200))
# 
# probe_pdata <- data.frame(cbind(complete_pDat %>% 
#                                   dplyr::select(-geo_accession) , 
#                                 exp_rot[,keep_probes]))
# ss <- step(lm(age ~ 1,  data=probe_pdata),
#   scope=as.formula(paste("age ~", 
#               paste(colnames(probe_pdata), collapse="+"))),
#   direction="both",
#   trace=-1
# )
# summary(ss)
# sum(str_detect(names(coef(ss)), "_at")) # 127



# - age
# - pack_years
# - race_ethnicity
fct_summ(complete_pDat)

# glmnet

select_probes <- function(pDat, exp, var, covars, multinomial=F, ...){
  
  Y <- pDat %>% pull(var)
  mod_mat <- model.matrix(~.-1, pDat %>% dplyr::select({{covars}}))
  X <- as.matrix(cbind(exp, mod_mat))
  train_idx <- 1:length(Y) #sample(1:length(Y), floor(0.8*length(Y)))
  Y_train <- Y[train_idx]
  X_train <- X[train_idx,]
  #Y_test <- Y[-train_idx]
  #X_test <- X[-train_idx,]
  cvfit = cv.glmnet(X_train, 
                    Y_train, ...)
  #train_assess <- assess.glmnet(cvfit, newx=X_train, 
  #                              newy=Y_train, s="lambda.1se")
  #valid_pred <- predict(cvfit, newx=X_test, 
  #                      type="response", s="lambda.1se")
  #train_pred <- predict(cvfit, newx=X_train, 
  #                      type="response", s="lambda.1se")
  #print(cor(valid_pred, Y_test)) # 0.32
  coef_l <- coef(cvfit, s="lambda.1se")
  if (multinomial==T){
    coef_l <- do.call(rbind, coef_l)
  }
  mat_coef <- coef_l %>% as.matrix()
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


# ------ RUN IMPUTATION ------ #
set.seed(415)
pdat_imp <- run_imp(expr=expDat5, pdat=pDat5.1, 
                    my_vars=c(age_probes, re_probes), m=5, defaultMethod="pmm")
list_pvars <- c("age", "sex", "smoking", "submission_date", "smoking:sex")
my_l <- extract_tt(pdat_imp=pdat_imp, expr=expDat5, list_pvars=list_pvars, 
                   coef_str="smokingS", number=10)

df_out <- extract_top_var(pdat_imp, expDat5, list_pvars, "smokingS", my_l)
rownames(df_out) <- df_out$probe
df_out %>% add_gene() 

my_l <- extract_tt(pdat_imp=pdat_imp, expr=expDat5, list_pvars=list_pvars, 
                   coef_str="sexmale:smokingS", number=100)
my_l2 <- setdiff(my_l, c(age_probes, re_probes)) 
length(intersect(my_l, c(age_probes, re_probes)))
start_time = Sys.time();
df_out <- extract_top_var2(pdat_imp, expDat5, list_pvars, "sexmale:smokingS", my_l2)
end_time=Sys.time()
end_time-start_time
df_all <- extract_all_var(pdat_imp, expDat5, list_pvars, "sexmale:smokingS")


df_out2 <- df_out %>% 
  rename(probes=probe) %>% add_gene() # TODO: add adjusted p


# ----- COMPARE RESULTS WITH pDat_complete vs imputed dat ----- #
mod <-  paste("~", paste(list_pvars, collapse=" + "))
tt_complete <- run_tt(complete_pDat, 
                      expDat5[,complete_pDat$geo_accession], 
                      "sexmale:smokingS", mod, number=100)
tt1 <- tt_complete %>% add_gene()
tt1

design <- model.matrix(~ sex + smoking + smoking*sex + 
                         age+race_ethnicity+pack_years+submission_date, data=complete_pDat) 
fit <- lmFit(expDat5[,complete_pDat$geo_accession], design)
fit <- eBayes(fit)
smok_ci_int_p <- topTable(fit, coef="sexmale:smokingS", number=nrow(fit), 
                          confint = TRUE)
smok_ci_int_p <- data.frame(smok_ci_int_p) %>% add_gene()
comb_ma_dat <- run_clean_ma(smok_ci_int_p %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome=""))
head(smok_ci_int_p)
comb_ma_dat %>% filter(adj.p < 0.05) # 72....

save(smok_ci_int_p, comb_ma_dat, file="data/results/ae_int.RData")

smok_ci_smok <- topTable(fit, coef="smokingS", number=nrow(fit), 
                          confint = TRUE)
smok_ci_smok <- data.frame(smok_ci_smok) %>% add_gene()
ma_smok <- run_clean_ma(smok_ci_smok %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome=""))
head(smok_ci_smok)
ma_smok %>% filter(adj.p < 0.05) # 72....

smok_ci_sex <- topTable(fit, coef="sexmale", number=nrow(fit), 
                         confint = TRUE)
smok_ci_sex <- data.frame(smok_ci_sex) %>% add_gene()
ma_sex <- run_clean_ma(smok_ci_sex %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome=""))
head(smok_ci_sex)
ma_sex %>% filter(adj.p < 0.05) # 72....
save(ma_smok, ma_sex, smok_ci_smok, smok_ci_sex, file="data/results/ae_smok.RData")



# ---- write out the results ---- #

# ---- summarize with meta-analysis ---- #


