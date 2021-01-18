
library(tidyverse)
library(miceadds)
sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")


fct_summ <- function(df) {
  summary(df %>% 
            mutate(across(everything(), as.factor)))
}
add_sl <- function(df) {
  df %>% left_join(sl %>% 
                     dplyr::rename(sex_lab=expression) %>% 
                     select(sample_acc, sex_lab))
}

# TODO: gene convert should ALSO update adjusted pvals
# - GENE CONVERT NEEDS TO DEAL W DUPLICATES
#   select Min?? or meta-analyze?
# - double check this all works
#  

# --- code for converting to hgnc --- #

load_genes <- function(list_genes, my_f, id_type){
  stopifnot(id_type %in% c("entrezgene_id", "ensembl_gene_id"))
  # warning - this takes ~30min
  if (!file.exists(my_f)){
    library('biomaRt')
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    convert_genes <- getBM(attributes = c("hgnc_symbol", id_type,
                                          "chromosome_name"), 
                           filters = id_type, values = list_genes,
                           mart=mart)
    # TODO: double check that var substitutions work
    save(convert_genes, file=my_f)
  }
  load.Rdata(my_f, "convert_genes")
  return(convert_genes)
}

load_gene_convert <- function(dataset, list_genes){
  stopifnot(dataset %in% c("affy", "rb", "tcga"))
  if (dataset=="affy"){
    my_f <-"data/affy_entrez_to_hgnc.RData"
    return(load_genes(list_genes, my_f, "entrezgene_id"))
  }
  if (dataset=="rb"){
    my_f <- "data/rb_ensembl_to_hgnc.RData"
    return(load_genes(list_genes, my_f, "ensembl_gene_id"))
  } 
  if (dataset=="tcga"){
    my_f <- "data/tcga_ensembl_to_hgnc.RData"
    return(load_genes(list_genes, my_f, "ensembl_gene_id"))
  }
}

# TODO double check that join works!
add_gene_info <- function(df, dataset="affy"){
  stopifnot(dataset %in% c("affy", "rb", "tcga"))
  df <- data.frame(df)
  convert_genes <- load_gene_convert(dataset, rownames(df))
  df$orig_gene <- rownames(df)
  gene_type <- ifelse(dataset=="affy", "entrezgene_id", "ensembl_gene_id")
  if (dataset=="affy") {
    df <- df %>% mutate(entrezgene_id=as.numeric(entrezgene_id)) 
  }
  df2 <- df %>% 
    left_join(convert_genes, by=c("orig_gene"=gene_type))%>% 
    dplyr::rename(gene=hgnc_symbol) %>%
    dplyr::select(-orig_gene) %>%
    as_tibble()
  
  df2 %>% 
    filter(!is.na(gene)) %>%
    filter(!duplicated(gene))
  
  # -- TODO: adjust FDR! -- #
}


