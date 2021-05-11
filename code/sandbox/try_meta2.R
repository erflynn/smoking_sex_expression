library('MetaIntegrator')
library('limma')
library('tidyverse')

gse <- getGEOData("GSE103174") # GPL13667 -  [HG-U219] Affymetrix Human Genome U219 Array
gse2 <- getGEOData("GSE31210") # GPL570

gse4 <- getGEOData("GSE32539") # 6244, 8786(miR) -  	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version],
phe2 <- read_csv("data/pdata_filt/GSE31210.csv")
head(gse2$originalData$GSE31210$pheno)
phe <- read_csv("data/pdata_filt/GSE103174.csv")


add_smok_lab <- function(ds, pheno) {
  
  new_pheno = ds$pheno %>% dplyr::select(geo_accession) %>% 
    inner_join(pheno, by=c("geo_accession"="sample_acc")) %>%
    dplyr::filter(smok %in% c("S", "NS"), !is.na(smok))
  new_expr = ds$expr[,new_pheno$geo_accession]
  new_class = ifelse(new_pheno$smok=="S", 1, 0)
  names(new_class) = new_pheno$geo_accession
  ds$pheno = new_pheno
  rownames(ds$pheno) <- new_pheno$geo_accession
  ds$expr = new_expr
  ds$class = new_class
  stopifnot(checkDataObject(ds, "Dataset", "Pre-Analysis"))
  return(ds)
}
meta <- add_smok_lab(gse$originalData$GSE103174, phe)
design <- model.matrix(~ smok, data=meta$pheno) 
fit <- lmFit(meta$expr, design)
fit <- eBayes(fit)
smok_out <- topTable(fit, coef="smokS", number=nrow(fit), 
                     confint = TRUE) %>% add_gene()
smok_out$gene <- gse$originalData$GSE103174$keys[smok_out$probes]
smok_ma0 <- run_clean_ma(smok_out %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome=""))


meta2 <- add_smok_lab(gse2$originalData$GSE31210, phe2)
design <- model.matrix(~ smok, data=meta2$pheno) 
fit <- lmFit(meta2$expr, design)
fit <- eBayes(fit)
smok_out <- topTable(fit, coef="smokS", number=nrow(fit), 
                          confint = TRUE) %>% add_gene()
head(smok_out)

smok_out_na <- smok_out %>% filter(is.na(gene)) %>% 
  select(probes, logFC, CI.L, CI.R, P.Value) %>%
  rename(gene=probes, logFC.l=CI.L, logFC.u=CI.R, p=P.Value) %>%
  mutate(src="single")
  
smok_out %>% filter(is.na(gene)) %>% nrow() # 12743 probes removed that don't map to genes


dim(smok_ma) # 20183
table(smok_ma$src) # 9426 genes have only one probe
top_ma <- smok_out %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome="") %>% 
  filter(geneSymbol %in% head(smok_out$gene))
top_ma2 <- run_clean_ma(top_ma)
smok_ma <- run_clean_ma(smok_out %>% rename(ID=probes, geneSymbol=gene) %>% mutate(chromosome=""))

# add in the NAs
smok_ma_w_na <- smok_ma %>% 
  ungroup() %>% select(-chromosome, -adj.p) %>%
  bind_rows(smok_out_na) %>%
  arrange(p) %>%
  mutate(adj.p=p.adjust(p, "fdr")) 


smok_ma %>% filter(adj.p < 0.05) %>% nrow() # 43
smok_ma_w_na  %>% filter(adj.p < 0.05) %>% nrow() # 32

# --- what abt top p, or top expr? --- #
top_p <- smok_out %>% 
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  arrange(P.Value) %>%
  mutate(n=n()) %>%
  top_n(1, wt=-P.Value) %>%
  ungroup() %>%
  mutate(adj.p=p.adjust(P.Value, "fdr")) 
head(top_p)

top_expr <- smok_out %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  arrange(desc(AveExpr)) %>%
  mutate(n=n()) %>%
  top_n(1, wt=AveExpr) %>%
  ungroup() %>%
  mutate(adj.p=p.adjust(P.Value, "fdr")) 

# jetset

smok_ma %>% 
  ungroup() %>% 
  filter(adj.p < 0.05 & abs(logFC)>=0.4) %>% 
  select(-chromosome) %>%
  View()

# volcano plot - does it look ok?

smok_ma %>% filter(gene %in% c("AHRR", "SERPIND1", "FUCA1", "DNASE2B"))

plot(meta2$expr[c("205576_at"),], meta2$expr[c("229354_at"),])


volcano_plot_de(smok_ma %>% rename(adj.P.Val=adj.p), pcut=0.05)+
  ylab("-log10 P value")
