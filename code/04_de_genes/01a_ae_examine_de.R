# 01a_ae_examine_de.R
# 11/22/2020
# Code for examining differentially expressed genes in AE data
#
# TODO
#  - write the list of genes to a file

library('tidyverse')
library('limma')
library('ggrepel')


load("data/ae_eset_pheno.RData")

# ---  code/functions for adding in HGNC symbol --- #
if(!file.exisits("data/affy_entrez_to_hgnc.RData")){
  library('biomaRt')
  # download the data if it doesn't exist, note this takes some time
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  convert_genes <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "chromosome_name"), 
                         filters = "entrezgene_id", values = rownames(ae_only),
                         mart=mart)
  save(convert_genes, file="data/affy_entrez_to_hgnc.RData")
}
load("data/affy_entrez_to_hgnc.RData")

add_gene_info <- function(df){
  df <- data.frame(df)
  df$entrezgene_id <- rownames(df)
  df2 <- df %>% 
    mutate(entrezgene_id=as.numeric(entrezgene_id)) %>%
    left_join(convert_genes) %>% 
    dplyr::rename(gene=hgnc_symbol) %>%
    dplyr::select(-entrezgene_id) %>%
    as_tibble()
  
  # todo - check this, should remove NAs, duplicates
  df2 %>% 
    filter(!is.na(hgnc_symbol)) %>%
    filter(!duplicated(hgnc_symbol))
}



# separate out: bronchial epithelium, alveolar macrophages
ae_only_meta <- comb_metadata3 %>% 
  filter(!tissue %in% 
           c("bronchial epithelium","alveolar macrophages")) 
ae_only <- eset_comb3[, ae_only_meta$geo_accession] # 539 vs 671
save(ae_only, ae_only_meta, file="data/ae_only_eset.RData")

# -- run differential expression analysis -- #
smok <- factor(ae_only_meta$smok)
sex <- factor(ae_only_meta$expr_sex)
design <- model.matrix(~smok+sex+sex*smok) # model
fit <- lmFit(ae_only, design)
fit <- eBayes(fit)

# get the lists of genes for each of the terms:
#  smoking, sex, sex*smoking interaction
# note: to get these terms, you'll have to list the fit$coef
smok_ci <- topTable(fit, coef="smokS", number=nrow(ae_only), confint = TRUE) %>%
  add_gene_info()
ggplot(smok_ci %>% head(20), aes(x=gene, y=logFC))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R))+
  coord_flip()+
  theme_bw()
  

# TODO: update to confint=TRUE
res_smok <- topTable(fit, coef="smokS", number=nrow(ae_only)) %>% 
  add_gene_info()
res_sex <- topTable(fit, coef="sexmale", number=nrow(ae_only)) %>% 
  add_gene_info()

res_int <- topTable(fit, coef="smokS:sexmale", number=nrow(ae_only)) %>% 
  add_gene_info()


# filter by p-value cutoff

plotChrDist <- function(df,  pcut=0.01){
  df %>% filter(adj.P.Val < pcut & 
                  chromosome_name %in% c(1:22, "X", "Y")) %>%
    # clean up the chromosome names
    mutate(chromosome_name=factor(chromosome_name, 
                                  levels= c(1:22, "X", "Y"))) %>%
    ggplot(aes(x=chromosome_name))+
    geom_histogram(stat="count")+
    theme_bw()+
    ylab("number of genes")
}


# --- plot the chromosome distribution of DE genes --- #
plotChrDist(res_smok)+ 
  ggtitle("smoking-related genes (P<0.01)")
ggsave("figures/ae_chr_dist_smok.png")

plotChrDist(res_sex)+
  ggtitle("sex-related genes (P<0.01)")
ggsave("figures/ae_chr_dist_sex.png")

plotChrDist(res_int)+
  ggtitle("smoking*sex genes (P<0.01)")
ggsave("figures/ae_chr_dist_smok_sex.png")

# --- plot the DE genes in a volcano plot --- #
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


volcano_plot_de(res_smok)+ggtitle("smoking x sex interaction effects")
ggsave("figures/volcano_ae_sex_smoking.pdf")

volcano_plot_de(res_sex)+ggtitle("sex effects")
ggsave("figures/volcano_ae_sex.pdf")

volcano_plot_de(res_int)+ggtitle("smoking effects")
ggsave("figures/volcano_ae_smoking.pdf")


# TODO: add code to write the lists of genes to a file
