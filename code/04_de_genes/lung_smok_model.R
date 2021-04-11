
# CHALLENGES:
# - initial MetaIntegrator didn't show anything

# POSSIBILITIES:
# - need to adjust by other covariates - extract these?
# - smoking labels are incorrect
# - sex labels are incorrect (how many conflict, visualize within dataset)
# - look at toptable data for each of the results -- what does it look like?
# - is there one dataset that is a huge outlier
#
# - use rb data because it is in the same space and try a glmnet model
# - use larger rb dataset and pull smok labels another way
# - add in adenocarcinoma


library(tidyverse)
library(MetaIntegrator)
source("code/00_utils.R")

lung_kept = read_csv("data/lung_sample_labels.csv")
head(lung_kept)
unique(lung_kept$study_acc) # 11
lung_kept %>%
  mutate(sex_lab=ifelse(is.na(sex_lab),'unlabeled',sex_lab)) %>%
  filter(!is.na(smok)  & smok != 'FS') %>%
  group_by(study_acc, smok, sex_lab) %>%
  count() %>%
  unite(grp,c(smok,sex_lab)) %>%
  pivot_wider(names_from=grp, values_from=n, values_fill=0) 
lung_kept1 <- lung_kept %>%
  mutate(sex_lab=ifelse(is.na(sex_lab),'unlabeled',sex_lab)) %>%
  filter(!is.na(smok)  & smok != 'FS') 
(lung_kept2 <- lung_kept1 %>%
  group_by(study_acc, smok) %>%
  count() %>%
  pivot_wider(names_from=smok, values_from=n, values_fill=0) %>%
  filter(NS >= 5  & S >= 5))
# six datasets

lung_kept2 %>% pull(study_acc)

# download with MetaIntegrator 
# what are the platforms?
# -- all affy so we might not be able to generalize to other non-affy but still a start!
# application ds: GPL10379 - Rosetta/Merck Human RSTA Custom Affymetrix 2.0 microarray [HuRSTA-2a520709]
# ... it is a bit of a problem we don't have this platform anywhere, but oh well
gse <- getGEOData("GSE103174") # GPL13667 -  [HG-U219] Affymetrix Human Genome U219 Array
gse2 <- getGEOData("GSE31210") # GPL570
gse3 <- getGEOData("GSE31908") # GPL96, 97, 570
gse4 <- getGEOData("GSE32539") # 6244, 8786(miR) -  	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version],
gse5 <- getGEOData("GSE37768") # GPL570
#gse6 <- getGEOData("GSE40364") # this failed download - mb validation?, GPL570 + gene
# this dataset is actually BLOOD


# add labels for smoking status
poss_samples <- lung_kept1 %>% semi_join(lung_kept2, by="study_acc")

add_smok_lab <- function(ds, my_acc) {
  my_samples = poss_samples %>% 
    filter(study_acc==my_acc) %>% 
    filter(sample_acc %in% ds$pheno$geo_accession) %>%
    arrange(sample_acc) 
  new_pheno = ds$pheno %>% filter(geo_accession %in% my_samples$sample_acc)
  new_expr = ds$expr[,my_samples$sample_acc]
  new_class = ifelse(my_samples$smok=="S", 1, 0)
  names(new_class) = my_samples$sample_acc
  ds$pheno = new_pheno
  ds$expr = new_expr
  ds$class = new_class
  stopifnot(checkDataObject(ds, "Dataset", "Pre-Analysis"))
  return(ds)
}


# there are three sub-datasets
gse3_gpl570 <- gse3$originalData$GSE31908_GPL570$pheno$geo_accession
gse3_gpl96 <- gse3$originalData$GSE31908_GPL96$pheno$geo_accession
gse3_gpl97 <- gse3$originalData$GSE31908_GPL97$pheno$geo_accession
poss_samples %>% 
  filter(sample_acc %in% gse3_gpl570) # all tumor

poss_samples %>% 
  filter(sample_acc %in% gse3_gpl96) # <-- let's just use this
poss_samples %>% 
  filter(sample_acc %in% gse3_gpl97)

gse1.1 <- add_smok_lab(gse$originalData$GSE103174, "GSE103174")
gse2.1 <- add_smok_lab(gse2$originalData$GSE31210, "GSE31210")
gse5.1 <- add_smok_lab(gse5$originalData$GSE37768, "GSE37768")

gse3.1_96 <- add_smok_lab(gse3$originalData$GSE31908_GPL96, "GSE31908")
gse3.1_97 <- add_smok_lab(gse3$originalData$GSE31908_GPL97, "GSE31908")

# TCGA lung data?

# divide into training and validation
# validation: GSE40364, GSE32539 (new platform + former)
gse4.1 <- add_smok_lab(gse4$originalData$GSE32539_GPL6244, "GSE32539")

# check expression
boxplot(gse1.1$expr)
boxplot(gse2.1$expr) # add to this one
min(gse2.1$expr)
gse2.1$expr <- gse2.1$expr-min(gse2.1$expr)+0.1
boxplot(gse2.1$expr) 

boxplot(gse5.1$expr)
boxplot(gse3.1_96$expr)
boxplot(gse3.1_97$expr)
boxplot(gse4.1$expr)
gse3.1_96$pheno$
gse5.1$class
# and create the metaObject
discovery_datasets <- list(gse1.1, gse2.1, gse5.1)
names(discovery_datasets) = c(gse1.1$formattedName, gse2.1$formattedName, 
                              gse5.1$formattedName)
metaObj=list() 
metaObj$originalData <- discovery_datasets
checkDataObject(metaObj, "Meta", "Pre-Analysis")


# save the metaObject

# run meta-integrator on the training data
metaObj <- runMetaAnalysis(metaObj)
metaObj <- filterGenes(metaObj, isLeaveOneOut = FALSE, FDRThresh = 0.05)
summarizeFilterResults(metaObj, getMostRecentFilter(metaObj))
str(metaObj, 2)
head(metaObj$metaAnalysis$pooledResults)


# - COVARIATES?
# - get last dataset
# Not relevant, actually airway epithelium
# gse_lung <- getGEO("GSE40364")
# pDat_lung=pData(gse_lung$`GSE40364-GPL570_series_matrix.txt.gz`)
# pDat_lung %>%
#   fct_summ()


# which of these data are included in refineio
s1 = read_csv("data/supp_tables/supp_table_1_annot.csv")
head(s1)

s1 %>%
  group_by( tissue) %>%
  filter(tissue %in% c('blood','lung'),
         included=='yes') %>%
  count()
s2=read_csv("data/supp_tables/s2_annot_counts.csv") 
head(s2)
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
head(smok_sl)

# get the counts for each of these?

# -- lung -- #
lung_studies = s1 %>% 
  filter(tissue=="lung") %>% 
  left_join(smok_sl, by='study_acc') %>%
  group_by(study_acc, expression) %>%
  count() %>%
  pivot_wider(names_from='expression', values_from='n') %>%
  filter(female>=5,male>=5)

# 02_phenotyping/lung_parse_rb_phe.R:lung_kept %>% write_csv("data/lung_sample_labels.csv")
# Which were included in PCs?
# 13 blood, 11 lung
# this reduced to five after filtering for sufficient samples per category

###   alternate set up  ###
# 1) try building a model with refine.bio data because it is normalized
load("data/compendia_pre_pc.RData")
dim(compendia_expr5)
list_studies = setdiff(lung_kept2 %>% pull(study_acc), "GSE40364")
lung_kept_filtered = lung_kept %>%
   filter(study_acc  %in% list_studies) %>%
  filter(smok != 'FS') #153
lung_kept_f <- lung_kept_filtered %>% 
  filter(sample_acc %in% colnames(compendia_expr5)) #102
(counts_study = lung_kept_f %>%
  group_by(study_acc, smok, sex_lab) %>%
  count() %>% 
  unite(category,c(smok,sex_lab))  %>%
  pivot_wider(names_from=category,values_from=n))

setdiff(lung_studies$study_acc, lung_kept$study_acc )

s1 %>% 
  filter(tissue=="lung") %>% 
  left_join(smok_sl, by='study_acc')  %>% 
  filter(sample_acc %in% colnames(compendia_expr5)) %>%
  group_by(study_acc, expression) %>%
  count()
# --> 26
# no smoker labels... hmm


# GSE31908 no NS_male
# many of the rest do not have sufficient samples anymore
# with the exception of GSE32539

# -- alternate idea: use all samples in refine.bio regardless study breakdown -- #
# A) try ones where we have smok annot
lung_kept2 <- lung_kept %>%
  filter(sample_acc %in% colnames(compendia_expr5)) %>%
  filter(!is.na(smok), smok != 'FS')

(smoke_counts=lung_kept2 %>%
  group_by(study_acc, smok) %>%
  count() %>%
  pivot_wider(names_from=smok,values_from=n))
lung_kept3 <- lung_kept2 %>% inner_join(
  smoke_counts %>%
  filter(!is.na(NS), !is.na(S)), by="study_acc")
# --> 7 studies

# possible validation
study_valid=smoke_counts %>%
  filter(NS>5,S>5)

# remove valid
valid_ds = study_valid$study_acc[[1]]
lung_expr <- compendia_expr5[,lung_kept2$sample_acc] # 182
dim(lung_expr)
lung_kept_train = lung_kept2 %>%
  filter(study_acc != valid_ds)
lung_kept_valid = lung_kept3 %>%
  filter(study_acc == valid_ds)

# look at top DE
design_t <- model.matrix(~smok + study_acc, data=lung_kept_train) # model
fit_t <- lmFit(as.matrix(lung_expr[,lung_kept_train$sample_acc]), 
               design_t)
fit_t <- eBayes(fit_t)
topTable(fit_t, coef="smokS") %>% add_gene_info("rb")

design_v <- model.matrix(~smok, data=lung_kept_valid) # model
fit_v <- lmFit(as.matrix(lung_expr[,lung_kept_valid$sample_acc]), 
               design_v)
fit_v <- eBayes(fit_v)
topTable(fit_v, coef="smokS") %>% add_gene_info("rb")

# train model
library(glmnet)

# partition by group
#fold_list <- partition_group_data(train_part %>% 
# select(-partition), nfolds=NFOLDS) %>% arrange(class, sample_acc)
#samp_to_fold <- fold_list$partition
#names(samp_to_fold) <- fold_list$sample_acc

x_train = t(lung_expr[,lung_kept_train$sample_acc])
y_train = ifelse(lung_kept_train$smok == 'NS',0,1)
x_valid = t(lung_expr[,lung_kept_valid$sample_acc])
y_valid = ifelse(lung_kept_valid$smok == 'NS',0,1)

fit = cv.glmnet(x_train,y_train, family="binomial")
plot(fit)
coef_tab = data.frame(as.matrix(coef(fit, s = 'lambda.1se')))
coef_tab %>%
  filter(X1 != 0) %>%
  add_gene_info('rb')

# test model - not good!
p_train = predict(fit, newx = x_train, s = 'lambda.1se', type= 'response')
table(round(p_train),y_train)
sum(round(p_train)  == y_train)/length(y_train) # 88.16 %
p_valid = predict(fit, newx = x_valid, s = 'lambda.1se')
table(round(p_valid),y_valid)
sum(round(p_valid)  == y_valid)/length(y_valid) # 68.75 %

# TODO: try adjusting for study with combat

# modcombat = model.matrix(~1, data=pheno)
# combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

study_adj = lm(expr ~ study_acc, data=)

# TODO: use SVA

# visualize adjusted data


# B) try ones where we don't - figure out how to grab this info
#  -- use MetaIntegrator

lung_kept3 %>% 
  #group_by(smok, sex_lab) %>% 
  #dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% 
  filter(metadata_sex!="unlabeled" &
           metadata_sex!=sex_lab) %>% nrow()

length(unique(lung_kept3$study_acc))


start_time = Sys.time(); 
lung_bp1 <- big.PCA(lung_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

lung_pcs2 <- data.frame(lung_bp1$PCs)
lung_pcs2$smoking <- lung_kept3$smok
lung_pcs2$sex <- lung_kept3$sex_lab
lung_pcs2$study <- lung_kept3$study_acc

ggplot(lung_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()


# 2) consider including adenocarcinoma data, including tcga

