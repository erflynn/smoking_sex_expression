# 01a_ae_examine_de.R
# 11/22/2020
# Code for examining differentially expressed genes in AE data
#
# TODO:
#  - change so GPL570 only

library('tidyverse')
library('limma')
library('ggrepel')
source("code/00_utils.R") # for converting


load("data/ae_eset_pheno.RData")

# separate out: bronchial epithelium, alveolar macrophages
ae_only_meta <- comb_metadata3 %>%
 filter(!tissue %in%
          c("bronchial epithelium","alveolar macrophages"))
ae_only <- eset_comb3[, ae_only_meta$geo_accession] # 539 vs 671
save(ae_only, ae_only_meta, file="data/ae_only_eset.RData")

# add submision date
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
dates <- dbGetQuery(con, sprintf("SELECT gsm, submission_date FROM gsm WHERE gsm IN ('%s')", 
                       paste(ae_only_meta$geo_accession, collapse="','")))
dbDisconnect(con)
ae_only_meta2 <- ae_only_meta %>% left_join(dates, by=c("geo_accession"="gsm"))

# -- run differential expression analysis -- #
smok <- factor(ae_only_meta2$smok) # smok
sex <- factor(ae_only_meta2$expr_sex) # expr_sex
submission_date <- factor(ae_only_meta2$submission_date)
design <- model.matrix(~smok+sex+sex*smok+submission_date) # model
fit <- lmFit(as.matrix(ae_only[,ae_only_meta2$geo_accession]), design)
fit <- eBayes(fit)

# get the lists of genes for each of the terms:
#  smoking, sex, sex*smoking interaction
# note: to get these terms, you'll have to list the fit$coef
res_smok <- topTable(fit, coef="smokS", number=nrow(ae_only)) %>% 
  add_gene_info("affy")
res_sex <- topTable(fit, coef="sexmale", number=nrow(ae_only)) %>% 
  add_gene_info("affy")

res_int <- data.frame(topTable(fit, coef="smokS:sexmale", n=nrow(ae_only)))%>% 
  add_gene_info("affy")

res_int$entrezgene_id <- rownames(res_int)
res_int2 <- res_int %>% 
  left_join(convert_genes %>% mutate(entrezgene_id=as.character(entrezgene_id))) 

# write out for STAMS
res_int2$entrezgene_id <- rownames(res_int2)
res_int2 %>% dplyr::select(entrezgene_id, P.Value) %>% write_csv("data/ae_interaction_pvals.csv")
res_int2 %>% dplyr::select(gene, P.Value) %>% write_csv("data/ae_interaction_pvals_hgnc.csv")

# plot some of the CI and individual genes
smok_ci_int <- topTable(fit, coef="smokS:sexmale", number=nrow(ae_only), 
                    confint = TRUE) %>%
  add_gene_info() %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(chromosome=case_when(
    chromosome_name %in% c("X", "Y") ~ chromosome_name,
    TRUE ~ "autosomal"))

ggplot(smok_ci_int %>% head(40) %>%
         mutate(gene=factor(gene, levels=rev(c(smok_ci_int %>% head(40))$gene))),
       aes(x=gene, y=logFC, fill=chromosome))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R))+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values=c("gray", "blue", "yellow"))
ggsave("figures/smok_int_ci_top.png")

# make barplots of expr level in groups to help with
# intuitive understanding
ent_ids <- load_gene_convert(ae_only, "affy") %>%
  filter(hgnc_symbol %in% smok_ci_int$gene) %>%
  filter(hgnc_symbol!="",
         chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(entrezgene_id=as.character(entrezgene_id))

ae_only$entrezgene_id <- rownames(ae_only)
ae_long <- ae_only %>% 
  pivot_longer(-entrezgene_id,
               names_to="geo_accession", values_to="expr") %>%
  inner_join(ent_ids, by="entrezgene_id") %>%
  left_join(ae_only_meta, by=c("geo_accession"))

# TODO - figure out how to add lines!
ae_long %>%
  filter(hgnc_symbol %in% head(smok_ci_int$gene, 12)) %>%
  ggplot(aes(y=expr, x=smok, fill=expr_sex))+
  geom_boxplot()+
  facet_wrap(~hgnc_symbol, scales="free")+
  theme_bw()+
  ylab("")+
  xlab("")
ggsave("figures/smok_int_ci_top_expr_levels.png")


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
volcano_plot_de(res_smok)+
  ggtitle("smoking effects")
ggsave("figures/volcano_ae_smoking.pdf")

volcano_plot_de(res_sex)+
  ggtitle("sex effects")
ggsave("figures/volcano_ae_sex.pdf")

volcano_plot_de(res_int)+
  ggtitle("smoking x sex interaction effects")
ggsave("figures/volcano_ae_sex_smoking.pdf")


