
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


# GPL570 PROBE TO GENE
library(hgu133plus2.db)
x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])

probe_gene <- tibble(
  "probes"=names(xx),
  "gene" = unlist(xx)
)
save(probe_gene, file="ref/probe_gene.RData")
# TODO save this object, slow


convertToGenes <- function(expData, key.df, gene.list){
  # TODO
  #  optimize, this is slow!! <-- possibly switch to MetaIntegrator function
  #  make sure the object contains expression data, keys
  #  change so that this is INPUT expr matrix + key mapping
  #  should work with both GEOQuery -AND- MetaIntegrator Objects
  
  list.keys <- key.df[key.df$gene %in% gene.list,]
  
  gene.to.probe <- split(key.df$probes,  key.df$gene) # this is slow... mb store for each platform
  expData2 <- do.call(cbind, lapply(1:length(gene.to.probe), function(x) {
    # get the gene and the probe
    g <- names(gene.to.probe)[x]
    p <- unlist(gene.to.probe[g])
    if (length(p)>1){
      expD <- expData[p,]
      df <- (data.frame(colMeans(expD, na.rm=TRUE)))
      return(df)
    }
    else {
      df <- data.frame(expData[p,])
      return(df)
    }})) ### ALSO SLOW...
  
  colnames(expData2) <- names(gene.to.probe)
  expData2.2 <- data.frame(t(expData2)) # columns are samples, rows are genes
  
  # create a data fram of NAs for missing genes
  missing.genes <- setdiff(gene.list, list.keys)
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)
  
  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expData2.3 <- expDataPlusMiss[gene.list,] # REORDER so it matches other data
  
  return(expData2.3)
}


