# example analysis for a small dataset

library('GEOquery')
library('tidyverse')
library('limma')


gse = getGEO("GSE42057") # get expression data here or from rb?
phe = read_csv("data/pdata_filt/gse42057.csv")
exp = exprs(gse$GSE42057_series_matrix.txt.gz)

phe2 = phe %>% rename(sex=sex_lab, smoking=smok, age=age_enroll, pack_years=ats_packyears) %>%
  filter(sex %in% c("male", "female"))

design = model.matrix(~sex + smoking + sex:smoking + age + pack_years + bmi, data=phe2)
fit = lmFit(exp[,phe2$sample_acc], design)
fit = eBayes(fit)
topTable(fit, coef="sexmale:smokingS") %>% add_gene()
smok_all <- topTable(fit, coef="smokingS", number=nrow(exp), confint=T) 
# 12531 transcripts, 7605 genes
head(smok_all, 30)


ma_smok <- data.frame(smok_all) %>% 
  #probe_gene_entrez() %>%
  add_gene() %>%
  mutate(chromosome="") %>%
  dplyr::rename(geneSymbol=gene, ID=probes) %>%
  run_clean_ma()
(ma_smok)
ma_smok %>% filter(gene=="LRRN3")
