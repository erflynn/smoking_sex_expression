
library(tidyverse)
library(miceadds)
library(meta)
library(ggrepel)
library(variancePartition)

# ---- general summary  ---- #

fct_summ <- function(df) {
  summary(df %>% 
            mutate(across(everything(), as.factor)))
}

sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
add_sl <- function(df) {
  df %>% left_join(sl %>% 
                     dplyr::rename(sex_lab=expression) %>% 
                     dplyr::select(sample_acc, sex_lab))
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
  ggplot(df2, aes(x=logFC, y=-log10(P.Value)))+
    geom_point(aes(col=gene_grp), alpha=0.5)+
    geom_label_repel(data=df2 %>% head(num_display), 
                     aes(label=gene), size=3)+
    scale_color_manual(values=c("red", "blue", "gray"))+
    theme_bw()
}

# --- summarize probes to genes --- #
# TODO - rename ID, geneSymbol
prep_data_ma <- function(df){
  df2 <- df %>%
    mutate(SD=(logFC-CI.L)/1.96) %>% # SE for effect size - we have 95%  CI = 1.96*SD
    dplyr::select(ID, geneSymbol, chromosome, logFC, SD, P.Value) %>%
    filter(geneSymbol!="") %>%
    group_by(chromosome,geneSymbol) %>%
    mutate(n=n()) 
  return(df2)
}

ma_probes_genes <- function(df){
  # generic inverse variance methods
  ma <- metagen(df %>% pull(logFC), # treatment estimate
                df %>% pull(SD), # standard error
                studlab=df %>% pull(ID),
                comb.fixed = TRUE,
                comb.random = FALSE,
                method.tau = "DL", # method for between study variance
                hakn = FALSE,
                prediction = FALSE,
                sm = "MD") # underlying summary measure (RR, OR, ASD, ROM, HR, MD)
  return(list("chromosome"=unique(unlist(df$chromosome)),
              "gene"=unique(df$geneSymbol),
              "logFC"=ma$TE.fixed,
              "logFC.l"=ma$lower.fixed,
              "logFC.u"=ma$upper.fixed,
              "p"=ma$pval.fixed,
              "n"=unique(df$n)))
  
}
# TODO update src to n
run_clean_ma <- function(df){
  df2 <- prep_data_ma(df)
  multi_probe_gene <- df2 %>% filter(n>1) %>% group_split(geneSymbol)
  if (length(multi_probe_gene) <= 1){
    return(df2)
  }
  ma_vals <- lapply(multi_probe_gene, ma_probes_genes)
  mult_gene_df <- data.frame(apply(data.table::rbindlist(ma_vals), 
                                   c(1,2), unlist)) %>%
    mutate(across(c(contains("logFC"), p), ~as.numeric(as.character(.)) )) %>%
    mutate(n=as.numeric(as.character(n)))
  comb_ma_dat <- df2 %>% filter(n==1) %>% 
    dplyr::rename(gene=geneSymbol, p=P.Value) %>%
    mutate(logFC.l=logFC-1.96*SD,
           logFC.u=logFC+1.96*SD) %>%
    dplyr::select(colnames(mult_gene_df)) %>%
    mutate(src="single") %>%
    bind_rows(mult_gene_df %>% mutate(src="mult")) %>% 
    arrange(p)
  comb_ma_dat$adj.p <- p.adjust(comb_ma_dat$p, method="fdr")
  return(comb_ma_dat)
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

grab_chr <- function(df){
   
      library('biomaRt')
      mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      chr_map <- getBM(attributes = c("hgnc_symbol",
                                            "chromosome_name"),
                             filters = "hgnc_symbol", values =unique(df$gene) ,
                             mart=mart)
    return(chr_map)
}


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

  
load_gene_convert <- function(list_genes, dataset){
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
    if (!"probes" %in% colnames(df)){
      df$probes <- rownames(df)
    }
    load("ref/probe_gene.RData")
    
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
  return(left_join(df, probe_gene, by="probes"))
}

probe_gene_entrez <- function(df){
  if (file.exists("ref/probe_gene_entrez.RData")){
    df <- data.frame(df)
    df$probes <- rownames(df)
    load("ref/probe_gene_entrez.RData")
    
  } else {
    library(hgu133plus2.db)
    x <- hgu133plus2ENTREZID
    mapped_probes <- mappedkeys(x)
    xx <- as.list(x[mapped_probes])
    probe_gene <- tibble(
      "probes"=names(xx),
      "gene" = unlist(xx)
    )
    save(probe_gene, file="ref/probe_gene_entrez.RData")
  }
  return(left_join(df, probe_gene, by="probes"))
}


## ----- code for mapping to STAMS ----- ##
map_to_stams <- function(df, fname){
  library(STRINGdb)
  string_db <- STRINGdb$new() 
  gene_values <- df %>% 
    dplyr::rename(Pvalue=p, Gene=gene) 
  gene_values2=gene_values %>% 
    filter(Pvalue != 1) %>% 
    filter(!is.na(Gene)) 
  gene_mapped = string_db$map(data.frame(gene_values2), "Gene", 
                              removeUnmappedRows = TRUE )
  plot(density(gene_mapped$Pvalue))
  # plot p dist
  min_p = 0
  if (any(gene_mapped$Pvalue==0)){
    min_p <- min(gene_mapped %>% filter(Pvalue!=0) %>% pull(Pvalue))
  } 
  gene_mapped2 <- gene_mapped %>% mutate(ifelse(Pvalue==0, min_p, Pvalue))
  save(gene_mapped2, file=fname)
  return(gene_mapped2)
}


####### variance partition ##### 
varPartPC <- function(expr, cutoff=0.8){
  expr_pcs <- prcomp(t(expr))
  last_kept <- which(summary(expr_pcs)$importance["Cumulative Proportion",] < cutoff)
  if (length(last_kept)==0){
    last_pc <- 2
  } else {
    last_pc <- max(last_kept)
  }
  expr_pcs2 <- expr_pcs$x[,1:last_pc]
  props <- summary(expr_pcs)$importance["Proportion of Variance", 1:last_pc]
  return(list("pcs"=expr_pcs2, "props"=props))
}

get_pvca <- function(expr_pcs_t, phe, form, props, cutoff=0.8){
  varPart <- fitExtractVarPartModel(expr_pcs_t, form, phe)
  return( colSums(varPart*props)/cutoff)
}

is.error <- function(x) inherits(x, "try-error")


rand_pvca1 <- function(expr_pcs_t, phe, form, props, cutoff=0.8){
  phe2 <- phe
  phe2$sex <- sample(phe$sex, nrow(phe2), replace=F)
  phe2$smoking <- sample(phe$smoking, nrow(phe2), replace=F)
  return(try(get_pvca(expr_pcs_t, phe2, form, props, cutoff)))
}

rand_pvca <- function(n, expr_pcs_t, phe, form, props, cutoff=0.8){
  res <- lapply(1:n,  function(i) 
    rand_pvca1(expr_pcs_t, phe, form, props, my_var, cutoff))
  succeeded <- !vapply(res, is.error, logical(1))
  
  df <- do.call(cbind, lapply(res[succeeded], data.frame)) 
  data.frame(t(df) )  %>% 
    as_tibble() %>%
    mutate(n=1:n()) %>%
    pivot_longer(-n, names_to="covariate", values_to="variance") %>%
    select(-n) %>%
    mutate(src="random") %>%
    mutate(covariate=str_replace_all(covariate, "\\.", ":"))
}

plot_pvca_rand <- function(est_var, rand_var){
  tibble("covariate"=names(est_var),
         "variance"=unlist(est_var)) %>%
    mutate(src="estimated") %>%
    bind_rows(rand_var) %>%
    filter(covariate!="Residuals") %>%
    ggplot(aes(x=covariate, y=variance, col=src))+
    geom_boxplot(position=position_dodge(0), alpha=0.5)+
    theme_bw()+
    ylab("fraction of variance")+
    xlab("")+
    scale_color_manual(values=c("red", "black"))+
    theme(legend.position = "None")
}
add_annot_mi <- function(df, platform){
  df2 <- data.frame(df  )
  df2$probes <- rownames(df2)
  load(sprintf("ref/%s_probe_gene.RData", tolower(platform))) #probe_gene
  df3 <-df2 %>% left_join(probe_gene, by="probes")
  unmapped <- df3 %>% filter(is.na(gene)) %>% nrow()
  nmulti <- df3 %>% filter(str_detect(gene, ",")) %>% nrow()
  print(sprintf("%s probes did not map out of %s; %s multi-mapped.", 
                unmapped, nrow(df3), nmulti))
  return(df3)
}
