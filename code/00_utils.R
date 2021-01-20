
library(tidyverse)
library(miceadds)


# ---- general summary  ---- #

fct_summ <- function(df) {
  summary(df %>% 
            mutate(across(everything(), as.factor)))
}

sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
add_sl <- function(df) {
  df %>% left_join(sl %>% 
                     dplyr::rename(sex_lab=expression) %>% 
                     select(sample_acc, sex_lab))
}

# --- code for DE viz --- #

# plot chromosome distribution
plotChrDist <- function(df,  pcut=0.01){
  df %>% filter(adj.P.Val < pcut & 
                  chromosome_name %in% c(1:22, "X", "Y")) %>%
    mutate(chromosome_name=factor(chromosome_name, 
                                  levels= c(1:22, "X", "Y"))) %>%
    ggplot(aes(x=chromosome_name))+
    geom_histogram(stat="count")+
    theme_bw()+
    ylab("number of genes")
}

# volcano plot
# TODO - fix if missing category
volcano_plot_de <- function(df, pcut=0.01, num_display=20){
  # labels the top "num_display" genes
  df2 <- df %>% 
    mutate(gene_grp=case_when(
      adj.P.Val > pcut ~ "not sig",
      logFC < 0 ~ "down",
      logFC > 0 ~ "up")) %>%
    mutate(gene_grp=factor(gene_grp, levels=c("down", "up", "not sig")))
  ggplot(df2, aes(x=logFC, y=-log10(adj.P.Val)))+
    geom_point(aes(col=gene_grp), alpha=0.5)+
    geom_label_repel(data=df2 %>% head(num_display), 
                     aes(label=gene), size=3)+
    scale_color_manual(values=c("red", "blue", "gray"))+
    theme_bw()
}

# --- code for visualizing PCs --- #
plotPC3 <- function(df, my_col){
  df2 <- df %>%
    mutate(PC2_copy=PC2) %>%
    pivot_longer(c(PC1,PC2), names_to="PCA", values_to="value1") %>%
    dplyr::rename(PC2=PC2_copy) %>%
    pivot_longer(c(PC2,PC3), names_to="PCB", values_to="value2") %>%
    filter(PCA!=PCB) %>%
    unite(PC_pair, c(PCA, PCB), sep=" vs ")
  
  ggplot(df2, 
         aes(x=value1, y=value2, col={{my_col}}))+
    geom_point(alpha=0.5)+
    theme_bw()+
    ylab("")+
    xlab("")+
    facet_wrap(.~PC_pair, nrow=3, scales="free")
}


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
add_gene <- function(df){
  if (file.exists("ref/probe_gene.RData")){
    df <- data.frame(df)
    df$probes <- rownames(df)
    left_join(df, probe_gene, by="probes")
  } else {
    library(hgu133plus2.db)
    x <- hgu133plus2SYMBOL
    mapped_probes <- mappedkeys(x)
    xx <- as.list(x[mapped_probes])
    probe_gene <- tibble(
      "probes"=names(xx),
      "gene" = unlist(xx)
    )
    save(probe_gene, file="ref/probe_gene.RData")
  }

}



