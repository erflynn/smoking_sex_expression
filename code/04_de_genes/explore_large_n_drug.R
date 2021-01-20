
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
                  `ECOG performance status`, `Stage`), as.numeric))
pDat2.1 %>% group_by(Gender, Chemotherapy) %>% count()
eDat2 <- exprs(gse10846[[1]])[,pDat2.1$geo_accession]
design <- model.matrix(~Gender+Stage+Chemotherapy+Gender:Chemotherapy+Age,
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
fit <- lmFit(eDat2, design)
fit <- eBayes(fit)
colnames(fit$coefficients)
tt_int = topTable(fit, coef="Gendermale:ChemotherapyR-CHOP-Like Regimen", n=nrow(eDat2))
topTable(fit, coef="Gendermale:ChemotherapyR-CHOP-Like Regimen", n=nrow(eDat2)) %>% 
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

# ---- smoking E-TABM-305 ---- #
etabm305_phe <- read_tsv("data/etabm305/E-TABM-305.sdrf.txt")
fct_phe <- etabm305_phe  %>% mutate(across(everything(), as.factor))
fct_l <- sapply(colnames(etabm305_phe), function(x) length(levels(fct_phe[[x]])))
colnames(etabm305_phe)[fct_l==1]
var_cols <- colnames(etabm305_phe)[fct_l!=1]
etabm305_phe2 <- etabm305_phe[,var_cols] 
etabm305_phe2 %>% fct_summ()

colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "Characteristics \\[", "c_") 
colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "Factor Value \\[", "f_") 
colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "\\[|\\]", "")
colnames(etabm305_phe2) <-tolower(str_replace_all(colnames(etabm305_phe2), " ", "_"))

etabm305_phe2 %>% dplyr::select(-f_age, -f_sex, -contains("term_accession"), -hybridization_name)
etabm305_phe3 <- etabm305_phe2 %>% 
  dplyr::select(f_age, c_sex, f_individual, c_test, c_clinicalinformation) %>%
  filter(c_clinicalinformation %in% c("smoker", "non-smoker")) %>%
  separate(c_test, into=c("test", "value"), sep="=") %>%
  dplyr::rename(hdl=value, age=f_age, sex=c_sex, individual=f_individual, smoking=c_clinicalinformation) %>%
  dplyr::select(-test)

etabm305_expr <- vroom("data/etabm305/E-TABM-305-processed-data-1770937661.txt")
colnames(etabm305_expr) <- str_replace_all(colnames(etabm305_expr), "Hybridization_individual ", "")
etabm305_expr <- etabm305_expr[-1,]
scan_id <- etabm305_expr$`Scan REF`
etabm305_expr$`Scan REF` <- NULL
etabm305_expr2 <- as.data.frame(etabm305_expr)
rownames(etabm305_expr2) <- scan_id

# grep -n "^\[main]" A-MEXP-691.adf.txt
etabm305_adf <- vroom("data/etabm305/A-MEXP-691.adf.txt", skip=20)
map_to_gene <- etabm305_adf %>%
  dplyr::select("Composite Element Name", "Reporter Group[role]", "Composite Element Database Entry[hugo]") %>%
  dplyr::rename(probe=`Composite Element Name`,
                gene=`Composite Element Database Entry[hugo]`,
                grp=`Reporter Group[role]`)
maps2 <- map_to_gene %>% filter(!is.na(gene))
etabm305_expr3 <- apply(etabm305_expr2[intersect(maps2$probe,rownames(etabm305_expr2)),], c(1,2), as.numeric)
etabm305_phe4 <- etabm305_phe3 %>% mutate(hdl=as.numeric(hdl)) %>% filter(!is.na(hdl))
design_blood <- model.matrix(~age+sex+hdl+smoking+smoking*sex, data=etabm305_phe4)
fit_b <- lmFit(etabm305_expr3[,etabm305_phe4$individual], design_blood)
fit_b <- eBayes(fit_b)
colnames(fit_b$coefficients)
b_sex <- topTable(fit_b, coef="sexmale", n=nrow(etabm305_expr3))
b_sex <- data.frame(b_sex)
b_sex$probe <- rownames(b_sex)
b_sex2 <- b_sex %>% left_join(map_to_gene, by="probe")
b_tt <- topTable(fit_b, coef="sexmale:smokingsmoker", n=nrow(etabm305_expr3))
b_tt <- data.frame(b_tt)
b_tt$probe <- rownames(b_tt)
b_tt2 <- b_tt %>% left_join(map_to_gene, by="probe") 

b_tt2 %>% filter(gene %in% c("ACAS2", "ACE")) %>% arrange(gene)
# TODO: use beadarray to look at!

# TODO"
# E-MTAB-6559 - tobacco heating system
# GSE71220 - statins + COPD + smoking

# -- ---- lung data ----- #
# missing smoking informating :/ 
# GSE23546  - super series
# GSE23352 , GSE23545, GSE23529
lung_studies <- getGEO("GSE23546")
pDat_lung <- pData(lung_studies[[1]])
pDat_lung %>% fct_summ()
lung1 <- getGEO("GSE23352")
lung2 <- getGEO("GSE23529")
lung3 <- getGEO("GSE23545")
samples1 <- pData(lung1[[1]])$geo_accession
samples2 <- pData(lung2[[1]])$geo_accession
samples3 <- pData(lung3[[1]])$geo_accession

pDat_lung2 <- pDat_lung %>% 
  mutate(ds=case_when(
    geo_accession %in% samples1 ~ "ds1",
    geo_accession %in% samples2 ~ "ds2",
    geo_accession %in% samples3 ~ "ds3"
    ) )
pDat_lung1 %>%  dplyr::select(characteristics_ch1, characteristics_ch1.1)

exprs(lung1)
