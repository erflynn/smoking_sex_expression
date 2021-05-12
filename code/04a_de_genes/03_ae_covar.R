#
# Exploratory analysis of AE present covariates
# Looks at PCs +
#
#
# TODO:
# --> REDO PC analysis after excluding the basal cell study
#    - add percentage for PCs
#  - look at PC variance decomposition, how many PCs do we want to look at?
#  - check for duplicates again with dgm_id
#  - rescue missing expr_sex labels
#  - replication: redo analysis at probe level & check concordance with Yang et al


library(tidyverse)
library(limma)
source("code/00_utils.R")
load("data/ae_only_eset.RData") # ae_only, ae_only_meta

# filter for complete-cases
complete_ae <- ae_only_meta %>% # 539!
  filter(!is.na(race_ethnicity), !is.na(age)) %>%
  mutate(pack_years=case_when(
    smok=="NS" ~ 0,
    !is.na(pack_years) ~ pack_years
  )) %>%
  filter(!is.na(pack_years))
# 237

# TODO: consider checking for duplicates again??

# TODO: rescue missing expr_sex labels?
#  note - these arent in the data so we need to filter differently

# remove mismatched
mismatched <- complete_ae %>% 
  filter(metadata_sex!="unlabeled" & expr_sex!=metadata_sex) %>%
  dplyr::select(geo_accession, contains("sex"), pred)
mismatched %>% nrow()

complete_ae2 <- complete_ae %>%
  anti_join(mismatched, by="geo_accession")
#--> 230
table(complete_ae2$race_ethnicity)

# look at packyears breakdown by group
complete_ae2 %>% 
  dplyr::select(smok, expr_sex, pack_years) %>%
  filter(smok=="S") %>%
  ggplot(aes(x=expr_sex, y=pack_years))+
  geom_violin(trim=FALSE)+
  #stat_summary(fun.data=mean_sdl, mult=1, 
  #               geom="pointrange")+
  geom_boxplot(width=0.1)+
  theme_bw()+
  xlab("sex")+
  ylab("pack years of smokers")
ggsave("figures/pack_years_violin.png")

pkyrs_stat <- complete_ae2 %>% 
  dplyr::select(smok, expr_sex, pack_years) %>%
  filter(smok=="S") 
f_pkyrs <- pkyrs_stat %>%
  filter(expr_sex=="female") %>%
  pull(pack_years)
m_pkyrs <- pkyrs_stat %>%
  filter(expr_sex=="male") %>%
  pull(pack_years)

# run a statistical test -- these are the same thing
#  this is n.s.
t.test(f_pkyrs, m_pkyrs) 
summary(aov(pack_years ~ expr_sex, data = pkyrs_stat))


# add in other covariates:
#  - check on DGM ID for duplication
#  - extract date, platform data
#  - look for other covar?
library('GEOmetadb')
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020
gsm_l <- paste(complete_ae2$geo_accession, collapse="','")
gsm_add <- dbGetQuery(con, sprintf("SELECT * from gsm WHERE gsm IN ('%s')" , gsm_l))
colnames(gsm_add)
gsm_add2 <- gsm_add %>% dplyr::select(gsm, gpl, status, 
                                      submission_date,
                                      characteristics_ch1, title, 
                                      source_name_ch1,
                                      description)
dbDisconnect(con)
# create cov table
ae_input_meta <- complete_ae2 %>% 
  dplyr::select(geo_accession, expr_sex, smok, race_ethnicity, age, pack_years,
                tissue) %>%
  dplyr::rename(sex=expr_sex, smoking=smok) %>%
  mutate(tissue=factor(as.character(tissue))) %>%
  left_join(gsm_add2 %>% dplyr::select(gsm, gpl, submission_date), 
            by=c("geo_accession"="gsm"))
ae_input_expr <- ae_only[,ae_input_meta$geo_accession]
table(ae_input_meta$gpl)

# VISUALIZATION - in PC space
# use bigpca b/c it is faster
library(bigpca)
start_time = Sys.time(); 
ae_pc <- big.PCA(as.matrix(ae_input_expr), pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

# put together with the metadata
pcs2 <- data.frame(ae_pc$PCs)
pcs2$sample_acc <- rownames(pcs2) 
pcs2.2 <- pcs2 %>% left_join(ae_input_meta, by=c("sample_acc"="geo_accession"))




# --- tissue --- #
# tracheal epithelium looks different
# airway epithelium blends in w SAE
# large airway has an interesting distribution
# -- what is going on with data to the L??
plotPC3(pcs2.2, tissue)

# look at outlier
outlier_pc1 <- pcs2.2 %>% filter(PC1 < -0.2)

gsm_add %>% filter(submission_date=="2010-01-04")
# <-- all this submission date
# - actually all basal cell culture
# fix this tissue annotation
ae_input_meta3 <- ae_input_meta %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=case_when(
    submission_date=="2010-01-04" ~ "large airway basal culture",
    tissue=="airway epithelium" ~ "small airway epithelium",
    TRUE ~ tissue))

pcs2.3 <- pcs2 %>% left_join(ae_input_meta3, by=c("sample_acc"="geo_accession"))

pcs2.4 <- pcs2.3 %>% filter(tissue!="large airway basal culture")
plotPC3(pcs2.3, tissue)

# platform - skip, actually all GPL570
ae_input_meta4 <- ae_input_meta3 %>% dplyr::select(-gpl)
# ---  date of sample --- #
# definitely correlated by date!
plotPC3(pcs2.3, submission_date)
library(lubridate)
pcs2.5 <- pcs2.4 %>% mutate(submission_date=ymd(submission_date),
                  submission_year=year(submission_date),
                  submission_month=month(submission_date)) %>%
  unite(year_month, c(submission_year, submission_month), remove=FALSE, sep="-")
plotPC3(pcs2.5, submission_year)
plotPC3(pcs2.5, as.factor(submission_year))
plotPC3(pcs2.5, year_month)
length(unique(pcs2.4$submission_date))

# race_ethnicity
plotPC3(pcs2.5, race_ethnicity)

# sex
plotPC3(pcs2.5, sex)

# pack_years
plotPC3(pcs2.5 %>% filter(smoking=="S"), pack_years)

# ---  look at correlations --- #
summary(aov(PC1 ~ tissue, data=pcs2.5))
summary(aov(PC2 ~ tissue, data=pcs2.5))
summary(aov(PC3 ~ tissue, data=pcs2.5))

summary(aov(PC1 ~ race_ethnicity, data=pcs2.5))
summary(aov(PC2 ~ race_ethnicity, data=pcs2.5))
summary(aov(PC3 ~ race_ethnicity, data=pcs2.5))

summary(aov(PC1 ~ submission_date, data=pcs2.5))
summary(aov(PC2 ~ submission_date, data=pcs2.5))
summary(aov(PC3 ~ submission_date, data=pcs2.5))

cor.test(pcs2.5$PC1, pcs2.5$age)
cor.test(pcs2.5$PC2, pcs2.5$age)
cor.test(pcs2.5$PC3, pcs2.5$age)

pcs_s <- pcs2.5 %>% filter(smoking=="S")
cor.test(pcs_s$PC1, pcs_s$pack_years)
cor.test(pcs_s$PC2, pcs_s$pack_years)
cor.test(pcs_s$PC3, pcs_s$pack_years)



# update matrix
ae_input_meta5 <- ae_input_meta4 %>% 
  filter(tissue!="large airway basal culture") %>%
  mutate(submission_date=ymd(submission_date),
         submission_year=year(submission_date),
         submission_month=month(submission_date)) %>%
  unite(year_month, c(submission_year, submission_month), remove=FALSE, sep="-")
ae_input_expr2 <- ae_input_expr[,ae_input_meta5$geo_accession]

# train model w/o covar
design <- model.matrix(~smoking+sex+sex*smoking, data=ae_input_meta5) # model
fit <- lmFit(ae_input_expr2, design)
fit <- eBayes(fit)
topTable(fit, coef="smokingS:sexmale", number=nrow(ae_only)) %>% head()


# train model w/ covar
#  gene ~ sex + smok + sex*smok + pkyrs + race_ethnicity + age
ae_input_meta6 <- ae_input_meta5 %>%
  mutate(race_ethnicity=as.character(race_ethnicity)) %>%
  mutate(race_ethnicity=ifelse(race_ethnicity=="black/hispanic", "black", 
                               race_ethnicity)) %>%
  mutate(age=scale(age)) %>%
  mutate(race_ethnicity=as.factor(race_ethnicity)) %>%
  dplyr::rename(race=race_ethnicity)

ae_input_meta6 %>%
  group_by(year_month, smoking, sex) %>%
  dplyr::count()
ggplot(ae_input_meta6 %>% unite("study_group", c(sex, smoking)), 
       aes(x=factor(submission_date), fill=study_group))+
  geom_bar(stat="count")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  ylab("number of samples")
ggsave("figures/ae_study_grp_breakdown_by_date.png")

complete_grps <- ae_input_meta6 %>%
  distinct(submission_date, sex, smoking)  %>%
  group_by(submission_date) %>%
  summarize(ngrps=n()) %>%
  filter(ngrps==4)


ggplot(ae_input_meta6, 
       aes(x=factor(submission_date), fill=race))+
  geom_bar(stat="count")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  ylab("number of samples")

ae_grp_complete <- ae_input_meta6 %>% 
  filter(submission_date %in% c(complete_grps$submission_date))
ae_input_grp <- ae_input_expr2[,ae_grp_complete$geo_accession]
dim(ae_input_grp)
design_c <- model.matrix(~smoking+sex+sex*smoking+race+age, 
                       data=ae_grp_complete) # model
fit_c <- lmFit(ae_input_grp, design_c)
fit_c <- eBayes(fit_c)
top_tabs <- lapply(colnames(fit_c$coefficients),
       function(x){
         tt <- data.frame(topTable(fit_c, coef=x));
       tt$gene <- rownames(tt);
       tt$coef<-x;
       return(tt)})
tt_comb <- do.call(rbind, top_tabs)
names(top_tabs) <- colnames(fit_c$coefficients)
top_tabs[["smokingS:sexmale"]]

tt_comb %>% 
  filter(coef!="(Intercept)") %>% arrange(adj.P.Val) %>%
  as_tibble() %>%
  View()


topTable(fit_c, coef="racehispanic", number=nrow(ae_only)) %>% head()
topTable(fit_c, coef="racewhite", number=nrow(ae_only)) %>% head()
topTable(fit_c, coef="raceblack", number=nrow(ae_only)) %>% head()
topTable(fit_c, coef="smokingS", number=nrow(ae_only)) %>% head()
topTable(fit_c, coef="pack_years", number=nrow(ae_only)) %>% head()
topTable(fit_c, coef="smokingS:sexmale", number=nrow(ae_only)) %>% head()

topTable(fit_c, coef="smokingS:sexmale", number=nrow(ae_only)) %>% head()

# --- how to acct for packyears?? b/c related to smoking status
# --- remove subjects in groups < 5?
# --- what terms to add (age, age^2, race_ethnicity)
# --- normalize age




# look at models

# compare models?
# -- model fit comparison for gene expression data?


# compare output
