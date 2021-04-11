# explore_tcga_data.R
# Code for exploring what data is available via TCGA
#  - also does a quick DE analysis
#
# TODO:
#  - save clinical data intermediately
#  - double check on other available expr data
#  - check concordance btw smoking labels

library(tidyverse)
library(TCGAbiolinks)
library(limma)

source("code/00_utils.R")

# get a list of projects
proj <-TCGAbiolinks:::getGDCprojects()
proj2 <- proj %>% filter(str_detect(id, "TCGA")) # filter for TCGA ones

# download all the clinical data
query <- GDCquery(project = proj2$id, 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.all <- GDCprepare(query)

pt_tabs <- names(clinical.all)[str_detect(names(clinical.all), "patient")]
cancer_types = pt_tabs
counts_by_type <- lapply(pt_tabs, function(pt_tab){
        df <- clinical.all[[pt_tab]];
        if ( "tobacco_smoking_history_indicator" %in% colnames(df)){
          df2 <- df %>% 
            group_by(tobacco_smoking_history_indicator, gender) %>% 
            count() %>%
            mutate(cancer_type=pt_tab);
          return(df2)
        } else {
          return(NA)
        }})

# remove ones that do not contain smoking data
count_df <- do.call(rbind, counts_by_type[!is.na(counts_by_type)])

# TODO: save the patient data so we don't have to re-extract

# make a cleaned up df for counts
# KEY:
# 1 - lifelong non-smoker (<100 cigarette)
# 2 - current smoker
# 3 - reformed smoker >15y
# 4 - reformed smoker <=15y
# 5 - reformed smoker 
count_df2 <- count_df %>% 
  rename(smok=tobacco_smoking_history_indicator,
         sex=gender) %>%
  filter(smok %in% c("1","2","3","4","5")) %>%
  mutate(sex=tolower(sex),
         cancer_type=str_replace_all(cancer_type, "clinical_patient_", ""),
         smok=case_when(
           smok == "1" ~ "never",
           smok == "2" ~ "current",
           TRUE ~ "former"
         )) %>%
  select(cancer_type, everything()) %>%
  group_by(cancer_type, smok, sex) %>%
  summarize(n=sum(n)) 
count_df2 %>% write_csv("data/tcga_counts.csv")
count_df2 <- read_csv("data/tcga_counts.csv")

# plot breakdown of sex/smok counts by TCGA study
ggplot(count_df2 %>%
         mutate(cancer_type=toupper(cancer_type)), 
       aes(x=smok, y=n, fill=sex))+
  geom_bar(stat="identity")+
  facet_grid(.~cancer_type)+
  theme_bw()+
  xlab("smoking status")+
  ylab("number of subjects")+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
  
ggsave("figures/tcga_smoking_counts.pdf") 

# ---- Look specifically at LUAD data ---- #
# TODO: figure out if there is other expression workflow data available

luad_query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(luad_query)
luad_data <- GDCprepare(luad_query)

# the clinical data has different labels than the data assoc w clinical
luad_query_c <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(luad_query_c)
clinical.luad <- GDCprepare(luad_query_c)


# TODO - this duplicates code from earlier
# just add existing data
pt_df <- clinical.luad$clinical_patient_luad %>%
  as_tibble() %>%
  dplyr::rename(smok=tobacco_smoking_history_indicator,
       sex=gender) %>%
  filter(smok %in% c("1","2","3","4","5")) %>%
  mutate(sex=tolower(sex),
         smok=case_when(
           smok == "1" ~ "never",
           smok == "2" ~ "current",
           TRUE ~ "former"
         )) %>%
  dplyr::select(bcr_patient_barcode, bcr_patient_uuid, smok, sex)

length(unique(clinical.luad$clinical_patient_luad$bcr_patient_barcode))
luad_data2 <- luad_data[,!is.na(luad_data$gender)]
expr_mat <- luad_data2@assays$data$`HTSeq - Counts`
rownames(expr_mat) <- luad_data2@rowRanges$ensembl_gene_id
colnames(expr_mat) <- luad_data2@colData@rownames

# grab the pheno table from the expression data
pheno_dat <- luad_data2@colData

# put together so that we have smoking status
# TODO: make sure that smoking status info isn't conflicting
pDat2 <- pheno_dat %>% 
  as_tibble() %>%
  dplyr::select(barcode, patient, gender, subtype_Smoking.Status) %>%
  as_tibble() %>%
  dplyr::select(barcode, patient, gender, subtype_Smoking.Status) %>% 
  left_join(pt_df, by=c("patient"="bcr_patient_barcode"))

pDat3 <- pDat2 %>% filter(smok %in% c("never", "former"))

# run model
smok <- factor(pDat3$smok) # smok
sex <- factor(pDat3$sex) # expr_sex
design_t <- model.matrix(~smok+sex+sex*smok) # model
colnames(exppr_mat) <- pheno_dat$barcode
v <- voom(as.matrix(expr_mat[,pDat3$barcode]), design=design_t, plot=TRUE) #, normalize="quantile"
fit_t <- lmFit(v, design_t)
fit_t <- eBayes(fit_t)
colnames(fit_t$coefficients)

# look at results
tt <- topTable(fit_t, coef="smoknever:sexmale", n=nrow(expr_mat)) %>%
  add_gene_info("tcga") 

volcano_plot_de(tt, pcut=0.1, num_display=10)+
  scale_color_manual(values=c("red", "gray"))

topTable(fit_t, coef="smoknever", n=nrow(expr_mat)) %>%
  add_gene_info("tcga") %>%
  volcano_plot_de()

topTable(fit_t, coef="sexmale", n=nrow(expr_mat)) %>%
  add_gene_info("tcga") 
  


