
library('tidyverse')
library('GEOquery')
library('limma')
library('bigpca')
library('jetset')
library('vroom')


source("code/00_utils.R")

# ----- GSE36868 - simvastatin on LCLs ----- #

gse36868 <- getGEO("GSE36868")
pData(gse36868[[1]]) %>% fct_summ()
pDat1 <- pData(gse36868[[1]]) %>%
  select(geo_accession, title, source_name_ch1, #contains("characteristics_ch1"),
      contains(":ch1")) 
colnames(pDat1) <- str_replace_all(colnames(pDat1), ":ch1", "")
pDat1 %>% 
  filter(source_name_ch1=="Control buffer exposed LCLs") %>%
  group_by(`smoking status`, gender) %>%
  count()
# not a ton of smokers vs NS - could just look at NS + simvastation

pDat1_ns <- pDat1 %>% filter(`smoking status` != "1")
eDat <- exprs(gse36868[[1]])[,pDat1_ns$geo_accession]

start_time = Sys.time(); 
pc1 <- big.PCA(eDat, pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time

pcs2 <- data.frame(pc1$PCs)
pcs2$geo_accession <- rownames(pcs2) 
pcs2.2 <- pcs2 %>% left_join(pDat1, by=c("geo_accession"))
colnames(pDat1)
plotPC3(pcs2.2, gender)
plotPC3(pcs2.2, source_name_ch1)
plotPC3(pcs2.2, as.factor(`rna labeling batch`))
cor.test(pcs2.2$PC1, as.numeric(pcs2.2$bmi))
cor.test(pcs2.2$PC1, as.numeric(pcs2.2$`cell growth rate`))
cor.test(pcs2.2$PC1, as.numeric(pcs2.2$`age (yrs)`))

cor.test(pcs2.2$PC1, as.numeric(pcs2.2$bmi))
cor.test(pcs2.2$PC2, as.numeric(pcs2.2$`cell growth rate`))

summary(aov(PC1 ~ as.factor(`hybridization batch`), data=pcs2.2))
summary(aov(PC1 ~ as.factor(`rna labeling batch`), data=pcs2.2))
summary(aov(PC1 ~ as.factor(`exposure`), data=pcs2.2))
summary(aov(PC1 ~ as.factor(`gender`), data=pcs2.2))
summary(aov(PC1 ~ as.factor(`source_name_ch1`), data=pcs2.2))

pDat1_ns2 <- pDat1_ns %>%
  mutate(across(c(bmi, `age (yrs)`, `cell growth rate`), as.numeric)) %>%
  mutate(source_name_ch1=ifelse(str_detect(source_name_ch1, "Simvastatin"), 1, 0)) %>%
  rename(sex=gender,
         trt=source_name_ch1,
         age=`age (yrs)`,
         exp_batch=exposure,
         rna_batch=`rna labeling batch`,
         hybrid_batch=`hybridization batch`,
         cell_growth=`cell growth rate`) %>%
  mutate(trt=as.factor(trt))
design <- model.matrix(~sex + trt + sex*trt + exp_batch + rna_batch + hybrid_batch +
                         bmi + age + cell_growth,
                         data=pDat1_ns2) # model
boxplot(eDat[,1:10]) # do we need to log transform? looks normalized but maybe overly so
fit <- lmFit(eDat, design)
fit <- eBayes(fit)
topTable(fit, coef="sexMale:trt1")
topTable(fit, coef="sexMale") # there should be some stuff here??

# ----- GSE10846 - rituximab ----- #
# (GSE31312 is also rituximab)
library('GEOquery')
gse10846 <- getGEO("GSE10846")
pData(gse10846[[1]]) %>% fct_summ()
pDat2 <- pData(gse10846[[1]]) %>%
  select(geo_accession, submission_date, contains("characteristics_ch1")) %>%
  as_tibble() %>%
  filter(submission_date != "Jan 16 2009") %>%
  select(-submission_date) %>%
  pivot_longer(-geo_accession, names_to="col_lab", values_to="value") %>%
  mutate(value=str_replace_all(value, "Clinical info: ", "")) %>%
  separate(value, into=c("tag", "value"), sep=": ", extra="merge") %>%
  select(-col_lab) %>%
  pivot_wider(names_from=tag, values_from=value) %>%
  mutate(Age=as.numeric(Age))

pDat2.1 <- pDat2 %>% filter(Gender!="NA") %>% mutate(Gender=as.factor(Gender)) %>%
  mutate(across(c(Age, `Number of extranodal sites`, `LDH ratio`, 
                  `ECOG performance status`, `Stage`), as.numeric)) %>%
  distinct() %>%
  filter(!is.na(Stage))
pDat2.1 %>% group_by(Gender, Chemotherapy) %>% count()
eDat2 <- exprs(gse10846[[1]])[,pDat2.1$geo_accession]
design2 <- model.matrix(~Gender+Stage+Chemotherapy+Gender:Chemotherapy+Age,
                       data=pDat2.1) # model
summ_df <- pDat2.1 %>% group_by(Gender) %>%
  summarize(Age=mean(Age,na.rm=T), 
            sd_Age=sd(Age, na.rm=T), 
            n=n(),
            num_ABC=sum(`Final microarray diagnosis`=="ABC DLBCL"),
            num_GBC=sum(`Final microarray diagnosis`=="GCB DLBCL"),
            num_Unc=sum(`Final microarray diagnosis`=="Unclassified DLBCL"),
            num_r_chemo=sum(Chemotherapy=="R-CHOP-Like Regimen"),
            num_extra = mean(`Number of extranodal sites`, na.rm=T),
            sd_extra = sd(`Number of extranodal sites`, na.rm=T),
            ldh = mean(`LDH ratio`, na.rm=T),
            sd_ldh = sd(`LDH ratio`, na.rm=T),    
            ecog = mean(`ECOG performance status`, na.rm=T),
            sd_ecog = sd(`ECOG performance status`, na.rm=T),
            num_alive=sum(`Follow up status`=="ALIVE"),
            stage1 = sum(Stage==1),
            stage2 = sum(Stage==2),
            stage3 = sum(Stage==3),
            stage4 = sum(Stage==4))

test_numeric <- function(df, col){
  t.test(df %>% filter(Gender=="female", !is.na({{col}})) %>% pull({{col}}),
         df %>% filter(Gender=="male", !is.na({{col}})) %>% pull({{col}}))
}
test_categorical <- function(df, col){
  cont.tib <- df %>% group_by(Gender, {{col}}) %>% 
    count() %>% 
    pivot_wider(names_from={{col}}, values_from=n)
  cont.tab <- data.frame(cont.tib)
  cont.tab$Gender <- NULL
  rownames(cont.tab) <- cont.tib$Gender
  chisq.test(cont.tab)
}

t.test(pDat2.1 %>% filter(Gender=="female") %>% pull(Age),
       pDat2.1 %>% filter(Gender=="male") %>% pull(Age))

test_numeric(pDat2.1, Age)
test_numeric(pDat2.1, `LDH ratio`)
test_numeric(pDat2.1, `ECOG performance status`)
test_numeric(pDat2.1, `Number of extranodal sites`)

test_categorical(pDat2.1, `Follow up status`)
test_categorical(pDat2.1, `Final microarray diagnosis`)
test_categorical(pDat2.1, Chemotherapy)
#test_categorical(pDat2.1, Stage) # TODO - ordinal test

# TODO - try from RAW data
boxplot(eDat2[,1:10]) # do we need to log transform? looks normalized but maybe overly so
fit <- lmFit(eDat2, design2)
fit <- eBayes(fit)
colnames(fit$coefficients)
tt_int = data.frame(topTable(fit, coef="Gendermale:ChemotherapyR-CHOP-Like Regimen", 
                  n=nrow(eDat2)))
tt_int$probes <- rownames(tt_int)
tt_int %>% filter(adj.P.Val < 0.10)  %>% left_join(probe_gene)

load("ref/gpl570_probe_gene.RData")

topTable(fit, coef="Gendermale:ChemotherapyR-CHOP-Like Regimen", 
         n=nrow(eDat2)) %>% 
  add_gene() %>% anti_join(tt_sex, by="probes") 
tt_sex = topTable(fit, coef="Gendermale",  n=nrow(eDat2))   %>% 
  filter(adj.P.Val < 0.05) %>% add_gene()


# map to genes and then to STRING for STAMS
# TODO - check this is correct mapping?
load("data/affy_entrez_to_hgnc.RData")
res <- jmap("hgu133plus2", eg=unique(probe_gene$gene))
gene_to_probe = tibble("entrezgene_id"=names(res), "probe"=unlist(res)) %>%
  left_join(convert_genes %>%mutate(entrezgene_id=as.character(entrezgene_id)) %>%
              filter(chromosome_name %in% c(1:22, "X", "Y")))
df_int <- data.frame(tt_int)
df_int$probe <- rownames(df_int)
df_int2 <- df_int %>% left_join(gene_to_probe, by="probe")
df_int3 <- df_int2 %>% filter(!is.na(hgnc_symbol) & hgnc_symbol!="")
length(unique(df_int3$probe))
df_int3$adjP <- p.adjust(df_int3$P.Value, method="fdr")
string_db <- STRINGdb$new() 
gene_values <- df_int3 %>% 
  dplyr::rename(Pvalue=P.Value, Gene=hgnc_symbol) %>%
  filter(!duplicated(Gene))

# TODO - check STAMS mapping
gene_values2=gene_values %>% 
  filter(Pvalue != 1) %>% 
  filter(!is.na(Gene)) 
gene_mapped = string_db$map(data.frame(gene_values2), "Gene", removeUnmappedRows = TRUE )
save(gene_mapped, file="data/gene_mapped_ritux.RData")



#gse100833 <- getGEO("GSE100833")


# ------ next study ------- #
gse_usket <- getGEO("GSE100833")
eDat <- exprs(gse_usket$GSE100833_series_matrix.txt.gz)
pDat <- pData(gse_usket$GSE100833_series_matrix.txt.gz)

# title, source_name_ch1, 
pDat2 <- pDat %>% select(geo_accession,contains("characteristics_ch1")) %>%
  pivot_longer(-geo_accession) %>%
  separate(value, into=c("column", "value"), sep=": ") %>%
  mutate(column=tolower(column)) %>%
  mutate(column=ifelse(str_detect(column, "age"), "age", column)) %>%
  filter(column != "") %>%
  select(-name) %>%
  pivot_wider(names_from=column, values_from=value)

fct_summ(pDat2)
# so many different trt grps
