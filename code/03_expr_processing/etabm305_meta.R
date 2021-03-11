library('tidyverse')
library('vroom')


source("code/00_utils.R")

# ---- smoking E-TABM-305 ---- #
etabm305_phe <- read_tsv("data/etabm305/E-TABM-305.sdrf.txt")
fct_phe <- etabm305_phe  %>% mutate(across(everything(), as.factor))
fct_l <- sapply(colnames(etabm305_phe), function(x) length(levels(fct_phe[[x]])))
colnames(etabm305_phe)[fct_l==1]
var_cols <- colnames(etabm305_phe)[fct_l!=1]
etabm305_phe2 <- etabm305_phe[,var_cols] 
etabm305_phe2 %>% fct_summ()

head(etabm305_phe2) %>%
  select(`Array Data File`, contains('FTP'))

colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "Characteristics \\[", "c_") 
colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "Factor Value \\[", "f_") 
colnames(etabm305_phe2) <- str_replace_all(colnames(etabm305_phe2), "\\[|\\]", "")
colnames(etabm305_phe2) <- tolower(str_replace_all(colnames(etabm305_phe2), " ", "_"))
colnames(etabm305_phe2)
etabm305_phe2 %>% dplyr::select(-f_age, -f_sex, -contains("term_accession"), -hybridization_name)
etabm305_phe3 <- etabm305_phe2 %>% 
  dplyr::select(f_age, c_sex, f_individual, c_test, c_clinicalinformation, contains('ftp')) %>%
  filter(c_clinicalinformation %in% c("smoker", "non-smoker")) %>%
  separate(c_test, into=c("test", "value"), sep="=") %>%
  dplyr::rename(hdl=value, age=f_age, sex=c_sex, individual=f_individual, 
                smoking=c_clinicalinformation, ftp=comment_arrayexpress_ftp_file) %>%
  dplyr::select(-test)
etabm305_phe3 %>% write_csv("data/etabm305_phe_processed.csv")

etabm305_expr <- vroom("data/etabm305/E-TABM-305-processed-data-1770937661.txt")
colnames(etabm305_expr) <- str_replace_all(colnames(etabm305_expr), "Hybridization_individual ", "")
etabm305_expr <- etabm305_expr[-1,]
scan_id <- etabm305_expr$`Scan REF`
etabm305_expr$`Scan REF` <- NULL
etabm305_expr2 <- as.data.frame(etabm305_expr)
rownames(etabm305_expr2) <- scan_id
# 19648 x 1240

# this is an illumina array Illumina 'Sentrix' Human Whole Genome Series 1 BeadChip
# (47k transcripts)

# to get the starting line
# grep -n "^\[main]" A-MEXP-691.adf.txt
etabm305_adf <- vroom("data/etabm305/A-MEXP-691.adf.txt", skip=20)
setdiff(scan_id, etabm305_adf$`Composite Element Name`) #  all of the probes are in this set!

map_to_gene <- etabm305_adf %>%
  dplyr::select("Composite Element Name", "Reporter Group[role]", "Composite Element Database Entry[hugo]") %>%
  dplyr::rename(probe=`Composite Element Name`,
                gene=`Composite Element Database Entry[hugo]`,
                grp=`Reporter Group[role]`)
maps2 <- map_to_gene %>% filter(!is.na(gene))
etabm305_expr3 <- apply(etabm305_expr2[intersect(maps2$probe,rownames(etabm305_expr2)),], 
                        c(1,2), as.numeric)
etabm305_phe4 <- etabm305_phe3 %>% mutate(hdl=as.numeric(hdl)) %>% filter(!is.na(hdl)) # 1219 --> 1192


# map genes to chromosomes
if (!file.exists("data/etabm305_adf_gene_metadata.csv")){
  list_genes <- map_to_gene %>% distinct(gene) %>% pull(gene) # 32651
  library('biomaRt')
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  convert_genes <- getBM(attributes = c("hgnc_symbol",
                                        "chromosome_name"), 
                         filters = "hgnc_symbol", values = list_genes,
                         mart=mart)
  
  adf_genes <- map_to_gene %>% left_join(convert_genes, by=c("gene"="hgnc_symbol"))
  adf_genes %>% write_csv("data/etabm305_adf_gene_metadata.csv")
}
adf_genes <- read_csv("data/etabm305_adf_gene_metadata.csv")
adf_genes2 <- adf_genes %>% mutate(present=(probe %in% scan_id))

# it appears X + Y genes are removed during processing --> redo 
adf_genes2 %>% 
  group_by(chromosome_name, present) %>% 
  count() %>% 
  pivot_wider(names_from = "present", values_from="n") %>%
  filter(chromosome_name %in% c("1","X", "Y"))

