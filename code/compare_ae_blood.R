# compare AE and blood results

# ----- load the data for comparison ------ #
load("data/results/blood_int.RData") # smok_ci_int_p, comb_ma_dat
load("data/results/blood_smok.RData") # smok_ma, sex_ma, smok_ci_smok, smok_ci_sex
smok_int_blood <- smok_ci_int_p
meta_int_blood <- comb_ma_dat
smok_blood <- smok_ci_smok 
meta_blood <- smok_ma
rm("smok_ci_int_p", "comb_ma_dat", "smok_ma", "sex_ma", "smok_ci_smok", "smok_ci_sex")

load("data/results/ae_int.RData") 
load("data/results/ae_smok.RData")
smok_int_ae <- smok_ci_int_p
meta_int_ae <- comb_ma_dat
smok_ae <- smok_ci_smok 
meta_ae <- ma_smok # this is diff
rm("smok_ci_int_p", "comb_ma_dat", "ma_smok", "ma_sex", "smok_ci_smok", "smok_ci_sex")

# TODO write the results out to the results tables?

# ---- replication of blood datasets? ---- #


# ----- Smoking overlap ---- #
# how many DE genes found in AE are also blood and vv? w/ corr ES, compare P

ae_genes <- meta_ae %>% 
  filter(adj.p < 0.05 & abs(logFC) >=0.3) # 932 genes
blood_genes <- meta_blood %>% 
  filter(adj.p < 0.05 & abs(logFC) >=0.3) # 13 genes
#smok_blood %>% filter(adj.P.Val < 0.05 & abs(logFC) >= 0.3) %>% nrow()

comb_both <- meta_ae %>% 
  filter(gene %in% c(ae_genes$gene, blood_genes$gene)) %>%
  inner_join(meta_blood, by="gene")

combined_dat <- ae_genes %>% ungroup() %>% select(-chromosome) %>% 
  inner_join(ae_in_blood_g %>% ungroup() %>% select(-chromosome), by="gene") # 482
combined_dat %>% 
  filter(logFC.x*logFC.y > 0) # 173
combined_dat %>% 
  filter(logFC.x*logFC.y > 0, adj.p.y < 0.05) # only 1!

comb_both %>% 
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(alpha=0.7)+
  #geom_label_repel(size=3)+
  theme_bw()+
  xlab("log fold change in airway epithelium")+
  ylab("log fold change in lymphocytes")

# compare to small ds?
load("data/results/GSE20189_de.RData") # df_smok1.1, df_int1, df_smok2.1, smok_ma1, smok_ma2, smok_int
ds1_genes <- smok_ma2 %>% filter(adj.p < 0.05, abs(logFC) >= 0.3) # 76
comb_blood <- smok_ma2 %>% 
  filter(gene %in% c(ds1_genes$gene, blood_genes$gene)) %>%
  inner_join(meta_blood, by="gene")

comb_blood %>% 
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(alpha=0.7)+
  #geom_label_repel(size=3)+
  theme_bw()+
  xlab("log fold change in whole blood rep")+
  ylab("log fold change in lymphocytes")
