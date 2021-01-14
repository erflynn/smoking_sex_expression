# Code for exploring the RB data
#
# Gets counts, plots PCs, and does prelim DE analysis + dedup of AE 
#
# TODO:
# - make PC analysis into functions
# - clean up AE DE analysis

library(vroom)
library(bigpca)
library(limma)


# --- PC PLOT --- #
compendia_expr <- vroom("data/compendia_smoking_expr.csv")
rownames(compendia_expr) <- compendia_expr$rid
compendia_expr2 <- compendia_expr %>% select(-rid) %>% as.matrix.data.frame()
rownames(compendia_expr2) <- compendia_expr$rid

head(kept_samples)
kept_s2 <- kept_samples %>% 
  inner_join(incl_studies2 %>% select(study_acc, tissue), by="study_acc") %>%
  filter(!is.na(expression)) %>%
  distinct(sample_acc, expression, tissue)

# ---- deduplicate and identify duplicates ---- #
compendia_expr3 <- compendia_expr2[,!duplicated(t(compendia_expr2))] # -->11599 to 11326

# remove missing
row_sums <- apply(compendia_expr3, 1, function(x) sum(is.na(x)))
col_sums <- apply(compendia_expr3, 2, function(x) sum(is.na(x)))
compendia_expr4 <- compendia_expr3[row_sums==0,]
# figure out which are duplicated...
dup_dat <- data.frame(t(compendia_expr3[1:8,order(compendia_expr3[1,])]))
dup_dat$sample <- rownames(dup_dat)
dup_split <- dup_dat %>% 
  as_tibble() %>% 
  group_by(ENSG00000000003, ENSG00000000005, ENSG00000000419, 
           ENSG00000000457, ENSG00000000460,
           ENSG00000000938, ENSG00000000971, ENSG00000001036) %>% 
  summarize(n=n(), sample=paste(sample, collapse=";")) %>%
  filter(n>=2) %>%
  arrange(desc(n), sample) %>%
  ungroup() %>%
  dplyr::select(n, sample) %>%
  mutate(grp=1:n()) %>%
  select(-n)


meta_sex <- sample_metadata_filt %>% 
  filter(label_type=="metadata") %>%
  dplyr::select(sample_acc, sex_lab) %>%
  dplyr::rename(metadata_sex=sex_lab)
kept_s3 <- kept_s2 %>% filter(sample_acc %in% colnames(compendia_expr4))
kept_s3 %>% separate_rows(tissue, sep=";") %>% group_by(tissue) %>% count() %>%
  arrange(desc(n))

compendia_expr5 <- compendia_expr4[,kept_s3$sample_acc]
save(compendia_expr5, file="data/compendia_pre_pc.RData")

# --- run PCA for individual tissue groups --- #

blood_kept3 <- blood_kept2 %>% filter(sample_acc %in% colnames(compendia_expr4))
blood_expr <- compendia_expr4[,blood_kept3$sample_acc]

# grab some counts
blood_kept3 %>% 
  left_join(meta_sex) %>%
  group_by(smok, sex_lab) %>%
  count()

blood_kept3 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex != sex_lab & metadata_sex != "unlabeled") %>%
  nrow() # 17 mismatch


blood_kept3 %>% 
  #group_by(smok, sex_lab) %>% 
  #dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% 
  filter(metadata_sex!="unlabeled" &
           metadata_sex!=sex_lab) %>% nrow()
length(unique(blood_kept3$study_acc))

start_time = Sys.time(); 
bp1 <- big.PCA(blood_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time


pcs2 <- data.frame(bp1$PCs)
pcs2$smoking <- blood_kept3$smok
pcs2$sex <- blood_kept3$sex_lab
pcs2$study <- blood_kept3$study_acc

ggplot(pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()


ggsave("figures/blood_pcs.pdf")
ggplot(pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=study))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/blood_pcs_study.pdf")






# ---- Lung --- #

lung_kept3 <- lung_kept %>% filter(sample_acc %in% colnames(compendia_expr4))
lung_expr <- compendia_expr4[,lung_kept3$sample_acc]
lung_kept3 %>% left_join(meta_sex)


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


ggsave("figures/lung_pcs.pdf")
ggplot(lung_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=study))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/lung_pcs_study.pdf")




# --- run PCA on AE processed by RB --- #

ae_kept <- ae_df %>% filter(sample_acc %in% colnames(compendia_expr4)) %>% distinct() %>%
  filter(!str_detect(smok, "copd"))
ae_kept2 <- ae_kept %>% 
  group_by(sample_acc) %>% 
  summarize(study_acc=paste(unique(study_acc), collapse=";"),
            sex_lab=paste(unique(sex_lab), collapse=";"),
            smok=paste(unique(smok), collapse=";")) 
ae_kept2 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex != sex_lab,
         metadata_sex != "unlabeled") %>%
  nrow() # 9

ae_kept2 %>% left_join(meta_sex) %>% 
  filter(metadata_sex=="unlabeled") %>% nrow()

ae_kept2 %>%
  left_join(meta_sex) %>%
  filter(metadata_sex == sex_lab,
         metadata_sex != "unlabeled") %>%
  nrow() # 237
# 649
ae_kept2 %>% 
  filter(smok!="FS", sex_lab!="unlabeled") %>%
  left_join(meta_sex) %>% nrow()
length(unique(ae_kept2$study_acc))
group_by(smok, sex_lab) %>% 
  dplyr::count() %>%
  filter(smok!="FS", sex_lab!="unlabeled")

ae_expr <- compendia_expr4[,ae_kept2$sample_acc]
start_time = Sys.time(); 
ae_bp1 <- big.PCA(ae_expr, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

ae_pcs2 <- data.frame(ae_bp1$PCs) 
ae_pcs2$smoking <- ae_kept2$smok
ae_pcs2$sex <- ae_kept2$sex_lab
ae_pcs2$study <- ae_kept2$study_acc

ggplot(ae_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/ae_pcs.pdf")

ggplot(ae_pcs2 %>%
         filter(sex!="unlabeled", smoking!="FS"), 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()+
  xlim(-0.03, 0.0)+
  ylim(-0.05, 0.05)
ggsave("figures/ae_pcs_zoom.pdf")



# ggplot(ae_pcs2 %>%
#          filter(sex!="unlabeled", smoking!="FS"), 
#        aes(x=PC1, y=PC2, col=study))+
#   geom_point(alpha=0.5)+
#   theme_bw()
# ggsave("figures/ae_pcs_study.pdf")
# 
# 
# ggplot(ae_pcs2 %>%
#          filter(sex!="unlabeled", smoking!="FS"), 
#        aes(x=PC1, y=PC2, col=study))+
#   geom_point(alpha=0.5)+
#   theme_bw()+
#   xlim(-0.02, 0.015)+
#   ylim(-0.02, 0.03)
# 
# ggsave("figures/ae_pcs_study_zoom.pdf")

# -- REMOVE SAMPLES THAT ARE IDENTICAL IN PC SPACE -- #
# todo - figure out PC space cutoff
ae_pcs2$sample_acc <- rownames(ae_pcs2)
ae_pcs3 <- ae_pcs2 %>% 
  #arrange(PC1, PC2, PC3) %>%
  mutate(across(c(PC1, PC2, PC3), ~round(., digits=4))) %>%
  group_by(PC1, PC2, PC3) %>%
  summarize(study=paste(unique(study), collapse=";"),
            sample_acc=paste(unique(sample_acc), collapse=";"),
            sex=paste(unique(sex), collapse=";"),
            smoking=paste(unique(smoking), collapse=";")) 
ae_pcs3 %>% filter(str_detect(sample_acc, ";")) 

ae_pcs3 %>% filter(str_detect(sample_acc, ";"), 
                   str_detect(sex, ";"), str_detect(smoking, ";")) 
# at a 4 digit PC cutoff rounding, we have no conflicting labels
# 54 pairs+ of samples are identical (!)

ae_pcs3.2 <- ae_pcs3 %>%
  filter(PC1 > -0.03 , 
         PC1 < -0.015,
         PC2 > -0.05, 
         PC2 < 0.05)
# --> 407

ae_pcs3.2$sample_acc <- sapply(ae_pcs3.2$sample_acc, function(x) strsplit(x, ";")[[1]][[1]])

# NOW TRY our model. 
ae_only_meta <- ae_pcs3.2 %>% ungroup() %>% select(sample_acc, sex, smoking) %>%
  filter(!str_detect(smoking, ";"))
ae_only <- ae_expr[,ae_only_meta$sample_acc]

smok <- factor(ae_only_meta$smoking) 
sex <- factor(ae_only_meta$sex_lag) 
design <- model.matrix(~smok+sex+sex*smok) # model
fit <- lmFit(ae_only, design)
fit <- eBayes(fit)
