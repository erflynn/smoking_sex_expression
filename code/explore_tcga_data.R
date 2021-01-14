
library(tidyverse)
clinical <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/clinical.tsv")
exposure <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/exposure.tsv")
fh <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/family_history.tsv") # NO ROWS

# ... hmm there are only 40 clinical rows and 20 exposure rows? what is going on here
# ok it's only selectinng the first page...

exposure2 <- exposure %>% mutate(across(everything(), ~ifelse(.=="'--", NA, .)))
clinical2 <- clinical %>% mutate(across(everything(), ~ifelse(.=="'--", NA, .)))
count_missing <- apply(exposure, 2, function(x) sum(x=="'--"))
exposure2 %>% 
  select(age_at_onset, cigarettes_per_day, 
         contains("smok"), contains("tobacco")) %>%
  mutate(across(c(contains("age"), contains("year"), 
                  contains("time"), contains("day")), as.numeric)) %>%
  summary()

clinical %>% select(gender, ethnicity, race, age_at_diagnosis)


# ---- 
library(TCGAretriever)
all_studies <- get_cancer_studies() # weird that there are 300 exactly

# Find published TCGA datasets
keep <- grepl("tcga_pub$", all_studies[,1]) # ... only 19 studies
tcga_studies <- all_studies[keep, ]

fetch_all_tcgadata()


###
library(TCGAbiolinks)
proj <-TCGAbiolinks:::getGDCprojects()
proj2 <- proj %>% filter(str_detect(id, "TCGA"))
proj2$id

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

count_df <- do.call(rbind, counts_by_type[!is.na(counts_by_type)])
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
count_df2 %>% write_csv("~/tcga_counts.csv")
count_df2 <- read_csv("~/tcga_counts.csv")
ggplot(count_df2 %>%
         mutate(cancer_type=toupper(cancer_type)), 
       aes(x=smok, y=n, fill=sex))+
  geom_bar(stat="identity")+
  facet_grid(.~cancer_type)+
  theme_bw()+
  xlab("smoking status")+
  ylab("number of subjects")+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
  
ggsave("tcga_smoking_counts.pdf") 

# -- need to figure out if this is what is available
luad_query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(luad_query)
luad_data <- GDCprepare(luad_query)

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.luad <- GDCprepare(query)

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

#library("DESeq2")
length(unique(clinical.luad$clinical_patient_luad$bcr_patient_barcode))
luad_data2 <- luad_data[,!is.na(luad_data$gender)]
expr_mat <- luad_data2@assays$data$`HTSeq - Counts`
rownames(expr_mat) <- luad_data2@rowRanges$ensembl_gene_id
colnames(expr_mat) <- luad_data2@colData@rownames

ibrary('biomaRt')
# download the data if it doesn't exist, note this takes some time
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convert_genes_tcga <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id",
                                      "chromosome_name"), 
                       filters = "ensembl_gene_id", values = rownames(expr_mat),
                       mart=mart)


#convert_genes <- getBM(attributes = c("hgnc_symbol", "entrezgene_id",
#"chromosome_name"), 
#                       filters = "entrezgene_id", values = rownames(ae_only),
#                       mart=mart)
save(convert_genes_tcga, file="data/tcga_ensembl_to_hgnc.RData")

add_gene_info_tcga <- function(df){
  df <- data.frame(df)
  #df$entrezgene_id <- rownames(df)
  df$ensembl_gene_id <- rownames(df)
  df2 <- df %>% 
    #mutate(entrezgene_id=as.numeric(entrezgene_id)) %>%
    left_join(convert_genes_tcga, by="ensembl_gene_id") %>% 
    dplyr::rename(gene=hgnc_symbol) %>%
    #dplyr::select(-entrezgene_id) %>%
    dplyr::select(-ensembl_gene_id) %>%
    as_tibble()
  
  # todo - check this, should remove NAs, duplicates
  df2 %>% 
    filter(!is.na(gene)) %>%
    filter(!duplicated(gene))
}

pheno_dat <- luad_data2@colData
length(intersect(pheno_dat$patient, clinical.luad$clinical_patient_luad$bcr_patient_barcode))



pDat2 <- pheno_dat %>% 
  as_tibble() %>%
  dplyr::select(barcode, patient, gender, subtype_Smoking.Status) %>%
  as_tibble() %>%
  dplyr::select(barcode, patient, gender, subtype_Smoking.Status) %>% 
  left_join(pt_df, by=c("patient"="bcr_patient_barcode"))
  #mutate(smok=case_when(
  #  subtype_Smoking.Status=="Lifelong Non-smoker" ~ 0,
  #  str_detect(subtype_Smoking.Status, "reformed") ~ 1)) %>%
  #filter(!is.na(smok)) %>%
  #mutate(sex=ifelse(gender=="female", 0, 1)) 
table(pDat2$smok)
pDat3 <- pDat2 %>% filter(smok %in% c("never", "former"))
table(pheno_dat$subtype_Smoking.Status)
#design <-tibble("sex"= pDat2$sex, "smok"=pDat2$smok, "sex*smok"=pDat2$smok*pDat2$sex)
smok <- factor(pDat3$smok) # smok
sex <- factor(pDat3$sex) # expr_sex
design_t <- model.matrix(~smok+sex+sex*smok) # model
colnames(expr_mat) <- pheno_dat$barcode
v <- voom(as.matrix(expr_mat[,pDat3$barcode]), design=design, plot=TRUE) #, normalize="quantile"
fit_t <- lmFit(v, design_t)
fit_t <- eBayes(fit_t)
colnames(fit$coefficients)
tt <- topTable(fit_t, coef="smoknever:sexmale", n=nrow(expr_mat)) %>%
  add_gene_info_tcga()

tt %>% inner_join(adj_res_int, by="gene") %>%
  arrange(adj.P.Val.y) %>%
  ggplot(aes(x=logFC.x, y=logFC.y))+geom_point(alpha=0.5)

topTable(fit, coef="smoknever", n=nrow(expr_mat)) %>%
  filter(adj.P.Val < 0.05) %>%
  add_gene_info_tcga()

topTable(fit, coef="sexmale", n=nrow(expr_mat)) %>%
  filter(adj.P.Val < 0.05)%>%
  add_gene_info_tcga()

add_gene_info_tcga()

head(tt , 30) %>%
# ---- TRY LIMMA + VOOM ---- #

# ddsSE <- DESeqDataSet(luad_data2, design = ~ gender)
# 
# keep <- rowSums(counts(ddsSE)) >= 10
# ddsSE <- ddsSE[keep,]
# ddsSE <- DESeq(ddsSE, parallel=TRUE)
# resultsNames(ddsSE)
# res <- results(ddsSE, name = "gender")
# dea <- as.data.frame(res)
# summary(res)

# --- STOP ----  #

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.LUAD.all <- GDCprepare(query)

luad_cols <- colnames(clinical.LUAD.all$clinical_patient_luad)[str_detect(colnames(clinical.LUAD.all$clinical_patient_luad), 
                                                             "nicotine|tobacco|smok|substance|pack|pkyr|pyr")]

# 75 non-smokers, 122 current, 411 former
# 524 TOTAL

# 1 - lifelong non-smoker (<100 cigarette)
# 2 - current smoker
# 3 - reformed smoker >15y
# 4 - reformed smoker <=15y
# 5 - reformed smoker 

summary(clinical.LUAD.all$clinical_patient_luad %>% select(luad_cols) %>% 
          mutate(across(everything(), as.factor)))
query <- GDCquery(project = "TCGA-LUSC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.LUSC.all <- GDCprepare(query)

tobacco_cols <- colnames(clinical.LUSC.all$clinical_patient_lusc)[
  str_detect(colnames(clinical.LUSC.all$clinical_patient_lusc),  "nicotine|tobacco|smok|substance|pack|pkyr|pyr")]
"tobacco_smoking_history_indicator"
summary(clinical.LUSC.all$clinical_patient_lusc %>% select(tobacco_cols) %>% 
         mutate(across(everything(), as.factor)))
# 18 non-smokers, 134 current, 340 former
# 506 total
summary(as.numeric(clinical.LUSC.all$clinical_patient_lusc$tobacco_smoking_pack_years_smoked))


# head and neck "TCGA-HNSC"
# ... should consider HPV
query <- GDCquery(project = "TCGA-HNSC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.HNSC.all <- GDCprepare(query)
tob_str <- "nicotine|tobacco|smok|substance|pack|pkyr|pyr"

summary(clinical.HNSC.all$clinical_patient_hnsc %>% select(contains("tobacco")) %>% 
          mutate(across(everything(), as.factor)))
# 122 non, 178 current out of

# esophageal "TCGA-ESCA"
query <- GDCquery(project = "TCGA-ESCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.ESCA.all <- GDCprepare(query)
summary(clinical.ESCA.all$clinical_patient_esca %>% select(contains("tobacco")) %>% 
          mutate(across(everything(), as.factor)))
# 56 non, 37 current

# bladder: "TCGA-BLCA"
query <- GDCquery(project = "TCGA-BLCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BLCA.all <- GDCprepare(query)
summary(clinical.BLCA.all$clinical_patient_blca %>% select(contains("tobacco")) %>% 
          mutate(across(everything(), as.factor)))
# 11 ns, 90 current

# head/neck, kidney, lung, uterus: CPTAC-3

query <- GDCquery(project = "CPTAC-3", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab") # doesn't work..
#GDCdownload(query)
#clinical.cptac3.all <- GDCprepare(query)

#kidney, pancreas, liver, bladder, cervix, colon and rectum,
#and a type of leukemia (acute myeloid leukemia).


# general - not sure how to query
# "GENIE-MSK", "GENE-VICC", "GENIE-DFCI", "GENE-NKI", "VAREPOP-APOLLO"
# "GENIE-GRCC", "FM-AD", "GENIE-UHN", "GENIE-MDA"

