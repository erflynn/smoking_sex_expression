# code for imputing missing phenotype values in AE data


# ------ 0. prep data for MI ------ #

# load pheno + expr data
load("data/ae_only_eset.RData") # --> ae_only, ae_only_meta

# TODO: check on sex labels?

# 539 samples --> 532
# all have expr_sex labels, 7 have mismatched expr + metadata
#ae_only_meta %>% 
#  filter(expr_sex!=metadata_sex & !is.na(metadata_sex)) 

ae_only_meta2 <- ae_only_meta %>%
  filter(expr_sex==metadata_sex | is.na(metadata_sex))
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
dates <- dbGetQuery(con, sprintf("SELECT gsm, submission_date FROM gsm WHERE gsm IN ('%s')", 
                                 paste(ae_only_meta2$geo_accession, collapse="','")))
dbDisconnect(con)
ae_only_meta3 <- ae_only_meta2 %>% left_join(dates, by=c("geo_accession"="gsm"))

# add submission date (all of them have)

# remove trachea --> 487
ae_only_meta4 <- ae_only_meta3 %>% 
  filter(str_detect(tissue, "airway")) %>%
  mutate(tissue=as.factor(as.character(tissue)))

ae_expr <- ae_only[,ae_only_meta4$geo_accession]

# how many have identical phe data?
ae_only_meta4 %>% fct_summ()
extended_id <- ae_only_meta4 %>% 
  separate_rows(id, sep=";") 
repeated_ids <- extended_id %>% 
  group_by(id) %>% 
  count() %>% 
  arrange(desc(n)) %>%
  filter(n >1)


repeated_ids %>%
  group_by(n) %>%
  count() # 316 have only 1


single_id <- ae_only_meta4 %>% 
  filter(is.na(id)) %>% # 487
  bind_rows(extended_id %>% anti_join(repeated_ids)) %>% # 320
  group_by(geo_accession) %>%
  summarize(across(everything(), ~paste(., collapse=";"))) # 238

# look at DGM vs not DGM
mult_condensed_dgm <- extended_id %>% 
  semi_join(repeated_ids) %>%
  filter(str_detect(id, "DGM")) %>%
  group_by(id) %>%
  summarize(across(everything(), ~paste(unique(.[!is.na(.)]), collapse=";")))

mult_condensed_dgm %>%
  dplyr::select(metadata_sex, race_ethnicity, smok, age, pack_years,
                expr_sex, submission_date) %>%
  View()
# these are identical where the data is available. 

# these are not identical but are more similar than expected by chance...
# ... BUT there's a lot missing. so just going to use DGM id
mult_condensed_non_dgm <- extended_id %>% 
  semi_join(repeated_ids) %>%
  filter(!str_detect(id, "DGM")) %>%
  group_by(id) %>%
  summarize(across(everything(), ~paste(unique(.[!is.na(.)]), collapse=";"))) 

apply(mult_condensed_non_dgm %>% 
        dplyr::select(metadata_sex, race_ethnicity, expr_sex, smok, age, pack_years, expr_sex), 2, 
      function(x) sum(str_detect(x, ";")))
nrow(mult_condensed_non_dgm)  

# what is the pairwise correlated in expression if repeated vs not repeated
cor2 <- coop::pcor(ae_expr)
list_crs <- cor2[upper.tri(cor2)]
stopifnot(length(list_crs)==choose(nrow(cor2), 2))

ggplot(tibble("cor"=list_crs), aes(x=cor))+
  geom_histogram(bins=60)+
  theme_bw()+
  xlab("pairwise correlation")+
  ylab("")


ggplot(tibble("cor"=list_crs), aes(x=cor))+
  geom_histogram(bins=60)+
  theme_bw()+
  xlab("pairwise correlation")+
  ylab("")+
  xlim(0.975, 1)

dgm_l <- mult_condensed_dgm %>% 
  dplyr::select(id, geo_accession) %>%
  separate_rows(geo_accession, sep=";") %>%
  group_by(id) %>%
  summarize(geo_accession=list(geo_accession))
dgm_l$pairs <- lapply(dgm_l$geo_accession, function(x) combn(x, 2, simplify=FALSE))
dgm_l$cors <- lapply(dgm_l$pairs, function(x) lapply(x, function(y) cor2[y[1], y[2]]) )
length(unlist(dgm_l$cors)) # 318

ggplot(tibble("cor"=unlist(dgm_l$cors)), 
       aes(x=cor))+
  geom_histogram(bins=60)+
  theme_bw()+
  xlab("pairwise correlation")+
  ylab("") 

# 0.990 seems like a good cutoff - dedup everything above this?
length(list_crs[list_crs > 0.990]) # 1588 pairs

# doesn't appear to be
non_dgm_l <- mult_condensed_non_dgm %>% 
  dplyr::select(id, geo_accession) %>%
  separate_rows(geo_accession, sep=";") %>%
  group_by(id) %>%
  summarize(geo_accession=list(geo_accession))
non_dgm_l$pairs <- lapply(non_dgm_l$geo_accession, function(x) combn(x, 2, simplify=FALSE))
non_dgm_l$cors <- lapply(non_dgm_l$pairs, function(x) lapply(x, function(y) cor2[y[1], y[2]]) )
summary(unlist(non_dgm_l$cors))
ggplot(tibble("cor"=unlist(non_dgm_l$cors)), 
       aes(x=cor))+
  geom_histogram(bins=60)+
  theme_bw()+
  xlab("pairwise correlation")+
  ylab("") 

to.remove <- which(cor2>0.990 & cor2!=1)


# filter + collapse by just dgm_id
dup_dgm <- ae_only_meta4 %>%
  mutate(dgm_id=str_extract(id, "DGM-[0-9]+")) %>%
  filter(!is.na(dgm_id)) %>%
  group_by(dgm_id) %>%
  summarize(list_gses=list(geo_accession), n=length(geo_accession)) %>%
  filter(n>1)

condense_dgm <- ae_only_meta4 %>%
  mutate(dgm_id=str_extract(id, "DGM-[0-9]+")) %>%
  semi_join(dup_dgm, by="dgm_id") %>%
  group_by(dgm_id) %>%
  arrange(dgm_id, submission_date) %>%
  summarize(across(everything(), ~paste(unique(.[!is.na(.)]), collapse=";"))) %>%
  mutate(across(c(age, pack_years), as.numeric)) %>%
  mutate(geo_accession=str_extract(geo_accession, "GSM[0-9]+"))
# selected the first submission date
ae_dgm_condensed <- ae_only_meta4 %>%
  mutate(dgm_id=str_extract(id, "DGM-[0-9]+")) %>%
  dplyr::select(-pred) %>%
  anti_join(dup_dgm, by="dgm_id") %>%
  bind_rows(condense_dgm %>% dplyr::select(-pred)) # 487 --> 314
ae_expr2 <- ae_expr[,ae_dgm_condensed$geo_accession]  
ae_pheno2 <- ae_dgm_condensed %>%
  mutate(pack_years=ifelse(smok=="NS", 0, pack_years))

# filter to remove > 0.99 pairwise correlation? -- not for now
# cor2f <- coop::pcor(ae_expr2)
# list_crs_f <- cor2f[upper.tri(cor2f)]
# ggplot(tibble("cor"=list_crs_f), 
#        aes(x=cor))+
#   geom_histogram(bins=60)+
#   theme_bw()+
#   xlab("pairwise correlation")+
#   ylab("") +
#   geom_vline(xintercept=0.990)
# 
#high_cor = data.frame(which(cor2f>0.990 & cor2f != 1, arr.ind=T))
#high_cor %>% filter(row > col) %>%
#   mutate(row_name=rownames(cor2f)[row], col_name=colnames(cor2f)[col])
# # --> 243 pairs

save(ae_expr2, ae_pheno2, file="data/large_ae_dedup.RData") # --> 314

# how many values are missing?
apply(ae_pheno2, 2, function(x) sum(is.na(x)))


# why is this not that much more than the other data 
# despite having a lot missing?
load("data/ae_filtered_pheno.RData") # expDat3, pDat5

setdiff(pDat5$gsm, ae_only_meta$geo_accession) # 87... hmm
# -- are these the ones we don't have expr_sex values for?
load("data/eset_comb.RData") # eset_comb, dup_split
comb_metadata <- read_csv("data/comb_sample_epithelium_metadata.csv")

setdiff(pDat5$gsm, comb_metadata$geo_accession)
pDat5 %>% anti_join(comb_metadata, by=c("gsm"="geo_accession"))

# TODO fix this

# ------ 1. run MI ------- #

# method 1: standard MI


# method 2: gene-based MI
# variable selection: ~15?


#ae_expr2, ae_pheno2
# predict missing race_ethnicity
re <- ae_pheno2 %>% 
  filter(!is.na(race_ethnicity)) %>%
  mutate(race_ethnicity=ifelse(race_ethnicity=="black/hispanic", "black",
                               race_ethnicity))
exp_re <- ae_expr2[,re$geo_accession]
design1 <- model.matrix(~ race_ethnicity,data=re)
fit1 <- lmFit(exp_re, design1)
fit1 <- eBayes(fit1)
topTable(fit1, coef="race_ethnicitywhite")
topTable(fit1, coef="race_ethnicityblack")
topTable(fit1, coef="race_ethnicityhispanic")
topTable(fit1, coef="race_ethnicityasian") # only three people here...


# predict missing age
phe_age <- ae_pheno2 %>% 
  filter(!is.na(age)) 
exp_age <- ae_expr2[,phe_age$geo_accession]
design2 <- model.matrix(~ age,data=phe_age)
fit2 <- lmFit(exp_age, design2)
fit2 <- eBayes(fit2)
rownames(topTable(fit2))
formula_str <- paste("age ~", paste(lapply(rownames(topTable(fit2)), 
                                           function(x) paste("X", x, sep="")), collapse=" + "))
lm( as.formula(formula_str), 
   data=cbind(data.frame(t(exp_age), "age"=phe_age$age)))

glmnet()

# try doing this with the complete data, leaving out the vars!
load("data/ae_filtered_pheno.RData") # expDat3, pDat5
Y <- pDat5$age
exp_rot <- t(expDat3[,pDat5$gsm])
mod_mat <- model.matrix(~.-1, pDat5 %>% 
  dplyr::select(sex, smoking, 
                race_ethnicity, 
                pack_years, 
                submission_date) %>%
  mutate(submission_date=as.factor(submission_date))) 
X <- as.matrix(cbind(exp_rot, mod_mat))
# normalize the data?
#Y_std <- (Y-mean(Y))/sd(Y)

set.seed(411)
train_idx <- sample(1:length(Y), floor(0.6*length(Y)))
Y_train <- Y_std[train_idx]
X_train <- X[train_idx,]
Y_test <- Y_std[-train_idx]
X_test <- X[-train_idx,]
cvfit = cv.glmnet(X[train_idx,], Y_train, alpha=0.5, 
                  family="gaussian")
train_assess <- assess.glmnet(cvfit, newx=X_train, 
                              newy=Y_train, s="lambda.1se")
valid_pred <- predict(cvfit, newx=X_test, 
                      type="response", s="lambda.1se")
train_pred <- predict(cvfit, newx=X_train, 
                      type="response", s="lambda.1se")
plot(train_pred, Y_train) # hmmm

plot(valid_pred, Y_test) # hmmm
valid_assess <- assess.glmnet(cvfit, newx=X_test, 
                              newy=Y_test, s="lambda.1se")

my.lambda <- "lambda.1se"
mat_coef <- coef(cvfit, s=my.lambda) %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
nonzero_coef # --> then we exclude these gene expression probes from later
probes_age <- names(nonzero_coef)[str_detect(names(nonzero_coef), "_at")]


exp4 <- expDat3[,pDat5$gsm]
pDat_probe <- data.frame(cbind(t(exp4[probes_age,]), pDat5))
pDat_probe$age[-train_idx] <- NA
pdat_imp <- mice(pDat_probe,m=5, method="pmm") 
tail(complete(pdat_imp, 3))
pdat_comp <- complete(pdat_imp, "long", include = TRUE)
pdat_comp$age.NA <- cci(pDat_probe$age)
ggplot(pdat_comp, 
       aes(x = .imp, y = age, color = age.NA)) + 
  geom_jitter(show.legend = FALSE, 
              width = .1)

design2 <- model.matrix(~ age+sex+smoking+smoking*sex+submission_date+
                         race_ethnicity,
             data=pDat5)
fit <- lmFit(expDat3[,pDat5$gsm], design2)
fit <- eBayes(fit)
colnames(topTable(fit))
tt <- data.frame(topTable(fit, coef="sexm:smokingS"))
tt %>% add_gene()


design2 <- model.matrix(~ age+sex+smoking+smoking*sex+submission_date+
                          race_ethnicity,
                        data=complete(pdat_imp, 1))
fit <- lmFit(expDat3[,pDat5$gsm], design2)
fit <- eBayes(fit)
colnames(topTable(fit))
tt <- data.frame(topTable(fit, coef="sexm:smokingS"))
tt %>% add_gene()





# create a mira by creating the fit objects?

mipo <-pool(with(data=pdat_imp, 
     lm(expDat3["228790_at",pDat5$gsm] ~ 
          sex*smoking + sex+smoking+age+submission_date+race_ethnicity))) 
mipo$pooled

#mipo <- pool(with(data = imp, lm(bmi ~ hyp + chl)))

comb_res %>% filter(gene=="FAM110B") %>% dplyr::select(logFC, CI.L, CI.R)
  ggplot(aes(x=logFC, y=log(P.Value)))+
  geom_point()

  
# code for running a differential expression analysis 
# with multiple imputation
# input: expr_dat
  
id_vars <- function(df){
  
}


run_imp <- function(expr, pdat, my_vars, ...){
  expr_pdat <- data.frame(cbind(t(expr[my_vars,]), pdat))
  return(mice(expr_pdat, ...) )
}

pdat_imp <- run_imp(expr=exp4, pdat=pDat5, my_vars=probes_age, m=4, defaultMethod="pmm")

run_tt <- function(imp_data, expr, coef_str, my_model){
  design2 <- model.matrix(as.formula(my_model), data=imp_data)
  fit <- lmFit(expr, design2)
  fit <- eBayes(fit)
  tt <- data.frame(topTable(fit, coef=coef_str, confint=T))
  tt$probes <- rownames(tt)
  return(tt)
}

extract_tt <- function(pdat_imp, expr, list_pvars, coef_str){
  my_model <- paste("~", paste(list_pvars, collapse=" + "))
   comb_res <- lapply(1:pdat_imp$m, function(x) { 
        df <- run_tt(complete(pdat_imp, x), expr, coef_str, my_model); 
        return(df$probes)
      })
   return(unique(unlist(comb_res)))
}
list_pvars <- c("age", "sex", "smoking", "submission_date")
my_l <- extract_tt(pdat_imp, expr, list_pvars, "smokingS")

extract_var <- function(pdat_imp, expr, list_pvars, coef_str,  out_var){
    # select top variables
    y = unlist(expr[out_var,]) 
    mipo <-pool(with(data=pdat_imp, 
                     lm(as.formula(paste("y", paste(list_pvars,collapse="+"), sep="~") ))))
    my_summ <- as.tibble(summary(mipo))
    comb_summ <- my_summ %>% 
      left_join(as.tibble(mipo$pooled) %>% dplyr::select(-estimate, -df), by="term")
    comb_summ %>%
      filter(term==coef_str) %>%
      dplyr::select(-term) %>%
      mutate(probe:={{out_var}})
}

extract_top_var <- function(pdat_imp, expr, list_pvars, coef_str, my_l){
  data.frame(do.call(rbind,
          lapply(my_l, function(x) 
            extract_var(pdat_imp, expr, list_pvars, coef_str, x) )
  ))
}

#extract_var(pdat_imp, exp4, list_pvars, "smokingS", "228790_at")
extract_top_var(pdat_imp, exp4, list_pvars, "smokingS", my_l)$p.value



# select the most correlated predictors
design_age <- model.matrix(~ age+sex+smoking+submission_date+race_ethnicity,
                        data=pDat5[train_idx,])
fit_age <- lmFit(expDat3[,pDat5$gsm][,train_idx], design_age)
fit_age <- eBayes(fit_age)
topTable(fit_age, coef="age")
adj_data <- lm(age ~ sex+smoking+submission_date+race_ethnicity, data=pDat5[train_idx,])
resid <- adj_data$residuals
cor(resid, expDat3[,pDat5$gsm][,train_idx]["229159_at",]) # -0.311
# predict missing pack-years
