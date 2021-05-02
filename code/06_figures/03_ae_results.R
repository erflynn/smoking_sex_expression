
library('GEOquery')

# examine output for AE
load( file="data/results/ae_int.RData") # --> smok_ci_int_p, comb_ma_dat,
volcano_plot_de(comb_ma_dat %>% rename(adj.P.Val=adj.p, P.Value=p), pcut=0.05)+
  ylab("-log10 P value")
ggsave("figures/paper_figs/ae_int_volcano.png")
load(file="data/results/ae_smok.RData") # --> ma_smok, ma_sex, smok_ci_smok, smok_ci_sex,
volcano_plot_de(ma_smok %>% rename(adj.P.Val=adj.p, P.Value=p), pcut=0.05)+
  ylab("-log10 P value")
ggsave("figures/paper_figs/ae_smok_volcano.png")

# add chromosome and write out the results
probes_per_gene <- smok_ci_int_p %>% group_by(gene) %>% count()
load("data/affy_entrez_to_hgnc.RData") # --> convert_genes
clean_genes <- function(x){
  annot_dat <- x %>% ungroup() %>%
    select(-chromosome, -src)  %>% 
    left_join(probes_per_gene) %>%
    rename(`n_probes`=n) %>%
    left_join(convert_genes %>% select(hgnc_symbol, chromosome_name) %>% filter(chromosome_name %in% c(1:22, "X", "Y")), 
              by=c("gene"="hgnc_symbol")) %>%
    rename(chromosome=chromosome_name) %>%
    select(chromosome, gene, logFC, logFC.l, logFC.u, n_probes, p, adj.p)
  annot_dat %>% 
    filter(adj.p < 0.05 & abs(logFC)>=0.3) %>%
    mutate(across(where(is.numeric), ~signif(., 3)))
}
comb_ma_dat %>% clean_genes() %>%
  write_csv("data/results/ae_sig_genes_int.csv")

ma_smok %>% clean_genes() %>%
  write_csv("data/results/ae_sig_genes_smok.csv")

# check replication
gse_ae <- getGEO("GSE4302")
str(gse_ae$GSE4302_series_matrix.txt.gz,2)
expr_ae <- exprs(gse_ae$GSE4302_series_matrix.txt.gz)
phe_ae <- read_csv("data/pdata_filt/gse4302.csv")
phe_ae2 <- phe_ae %>% filter(smok %in% c("S", "NS"), sex_lab %in% c("male", "female")) %>%
  rename(smoking=smok, sex=sex_lab)
design <- model.matrix(~ sex + smoking + smoking*sex , data=phe_ae2) 
fit <- lmFit(expr_ae[,phe_ae2$sample_acc], design)
fit <- eBayes(fit)
smok_int_ae <- topTable(fit, coef="sexmale:smokingS", number=nrow(fit), 
                          confint = TRUE) %>% add_gene()
smok_ae <- topTable(fit, coef="smokingS", number=nrow(fit), 
                        confint = TRUE) %>% add_gene()

dim(smok_int_ae) # 54675
smok_ma_ae <- run_clean_ma(smok_int_ae %>% 
                             rename(ID=probes, geneSymbol=gene) %>% 
                             mutate(chromosome=""))

smok_ma_ae_smok <- run_clean_ma(smok_ae %>% 
                             rename(ID=probes, geneSymbol=gene) %>% 
                             mutate(chromosome=""))
# save the results
save(smok_int_ae, smok_ae, smok_ma_ae, smok_ma_ae_smok, file="data/results/gse4302_out.RData")
disc_genes <- comb_ma_dat %>% 
  filter(adj.p < 0.05 & abs(logFC) >=0.3)

rep_genes <- smok_ma_ae %>% filter(gene %in% disc_genes$gene)
rep_genes %>% mutate(p.adj=p.adjust(p, "fdr"))
combined_dat <- disc_genes %>% ungroup() %>% select(-chromosome) %>% 
  left_join(rep_genes %>% ungroup() %>% select(-chromosome), by="gene")
combined_dat %>% filter(logFC.x*logFC.y > 0) # 21 out of 30

combined_dat %>% 
  filter(adj.p.x < 0.05 & abs(logFC.x) >= 0.3) %>%
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(alpha=0.7)+
  geom_label_repel(size=3)+
  theme_bw()+
  ylab("log fold change in replication")+
  xlab("log fold change in discovery")
ggsave("figures/paper_figs/ae_discovery_rep.png")

# ----- compare smoking ------ #
# ma_smok, smok_ci_smok, 


disc_smok <- ma_smok %>% filter(adj.p < 0.05, abs(logFC) >= 0.3) %>% ungroup() %>% select(-chromosome)
valid_smok <- smok_ma_ae_smok %>% filter(gene %in% disc_smok$gene) %>% 
  ungroup() %>% select(-chromosome) %>%
  mutate(p.adj=p.adjust(p, "fdr"))
valid_smok %>% filter(p.adj < 0.05) # 59/932

both_smok <- disc_smok %>% left_join(valid_smok, by="gene") 
both_smok %>% filter(logFC.x*logFC.y > 0) %>% nrow() # 754/932

both_smok %>% 
  filter(adj.p.x < 0.05 & abs(logFC.x) >= 0.3) %>%
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(alpha=0.5)+
 # geom_label_repel(size=3)+
  theme_bw()+
  ylab("log fold change in replication")+
  xlab("log fold change in discovery")
ggsave("figures/paper_figs/ae_discovery_rep_smoking.png")


# other AE/BE gses
gse_ae2 <- getGEO("GSE7895")
phe_gse7895 <- read_csv("data/pdata_filt/gse7895.csv")

my_expr <- exprs(gse_ae2$GSE7895_series_matrix.txt.gz)
  my_phe2 <- phe_gse7895 %>% filter(smok %in% c("S", "NS"), 
                            sex_lab %in% c("male", "female")) %>%
    mutate(pack_years=ifelse(smok=="NS", 0, pack_years)) %>%
    rename(smoking=smok, sex=sex_lab)
my_expr_sm <- my_expr[,my_phe2$sample_acc]
  my_design <- model.matrix(~ sex + smoking + smoking*sex + age + pack_years, data=my_phe2) 
  my_fit <- lmFit(my_expr_sm, my_design)
  my_fit <- eBayes(my_fit)

my_fit <- get_fit(, phe_gse7895)

smok_int_gse7895 <- topTable(my_fit, coef="sexmale:smokingS", number=nrow(my_fit), 
                        confint = TRUE) %>% add_gene()
smok_gse7895 <- topTable(my_fit, coef="smokingS", number=nrow(fit), 
                    confint = TRUE) %>% add_gene()

dim(smok_int_gse7895) 
ma_int_gse7895 <- run_clean_ma(smok_int_gse7895 %>% 
                             rename(ID=probes, geneSymbol=gene) %>% 
                             mutate(chromosome=""))

ma_smok_gse7895 <- run_clean_ma(smok_gse7895 %>% 
                                  rename(ID=probes, geneSymbol=gene) %>% 
                                  mutate(chromosome=""))

# save the results
save(smok_int_gse7895, smok_gse7895, ma_int_gse7895, ma_smok_gse7895, 
     file="data/results/gse7895_out.RData")

load("data/results/gse7895_out.RData")
ae_rep <- ma_smok_gse7895 %>% filter(gene %in% ae_de_gene$gene) %>%
  inner_join(ae_de_gene, by="gene") 
ae_rep %>%
  filter(logFC.x*logFC.y > 0, p.x < 0.05)

cor.test(ae_rep$logFC.x, ae_rep$logFC.y) # 63.0
ae_rep_int <- ma_int_gse7895 %>% filter(gene %in% ae_de_gene_i$gene) %>% 
  inner_join(ae_de_gene_i, by="gene")
ae_rep_int %>%
  filter(logFC.x*logFC.y > 0, p.x < 0.05)

cor.test(ae_rep_int$logFC.x, ae_rep_int$logFC.y) # -0.04, p=0.855

# compare smok
valid_smok <- ma_smok_gse7895 %>% filter(gene %in% disc_smok$gene) %>% 
  ungroup() %>% select(-chromosome) %>%
  mutate(p.adj=p.adjust(p, "fdr"))
valid_smok %>% filter(p.adj < 0.05) # 52/645

both_smok <- disc_smok %>% inner_join(valid_smok, by="gene") %>%
  mutate(replicated=(logFC.x*logFC.y > 0 & p.adj < 0.05))
both_smok %>% filter(logFC.x*logFC.y > 0) %>% nrow() # 475/645

both_smok %>%  
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(aes(col=replicated), alpha=0.5)+
  geom_label_repel(data=head(both_smok %>% 
                               arrange(desc(abs(logFC.x))) %>%
                               filter(replicated), 20), size=3)+
  theme_bw()+
  ylab("log fold change in replication (GSE9875)")+
  xlab("log fold change in discovery")+
  scale_color_manual(values=c("black", "red"))+
  theme(legend.position ="None")
ggsave("figures/paper_figs/ae_rep_smok_gse7895.png")

# compare smok*sex


rep_genes <- ma_int_gse7895 %>% filter(gene %in% disc_genes$gene) %>% 
  mutate(p.adj=p.adjust(p, "fdr"))

combined_dat <- disc_genes %>% ungroup() %>% select(-chromosome) %>% 
  inner_join(rep_genes %>% ungroup() %>% select(-chromosome), by="gene")  %>%
  mutate(replicated=(logFC.x*logFC.y > 0 & p.adj < 0.05))
combined_dat %>% filter(p.adj < 0.1) # 8

combined_dat %>% 
  filter(adj.p.x < 0.05 & abs(logFC.x) >= 0.3) %>%
  ggplot(aes(x=logFC.x, y=logFC.y, label=gene)) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=0, lty=2)+
  geom_point(aes(col=replicated), alpha=0.7)+
  geom_label_repel(size=3)+
  theme_bw()+
  ylab("log fold change in replication (GSE9875)")+
  xlab("log fold change in discovery")+
  scale_color_manual(values=c("black", "red"))+
  theme(legend.position ="None")
ggsave("figures/paper_figs/ae_rep_int_gse7895.png")
save(both_smok, combined_dat, file="data/results/ae_rep.RData")

# probability by chance?
 # 475
both_smok2 <- both_smok
set.seed(416)
sim_counts_s <- sapply(1:1000, function(x){
  both_smok2$logFC.y <- sample(both_smok$logFC.y, nrow(both_smok), replace=T)
  both_smok2 %>% filter(logFC.x*logFC.y>0) %>% nrow()
})

ggplot(tibble("num_probes"=sim_counts_s), aes(x=num_probes))+
  geom_histogram()+theme_bw()+
  geom_vline(xintercept=both_smok %>% filter(logFC.x*logFC.y>0) %>% nrow(), col="red")+
  xlab("Number of same direction genes")+
  ylab("Number of randomized runs")
ggsave("figures/ae_rep_random_smok.png")


combined_dat %>% filter(logFC.x*logFC.y>0) %>% nrow() 
both_sex <- combined_dat

set.seed(417)
sim_counts_s2 <- sapply(1:1000, function(x){
  both_sex$logFC.y <- sample(combined_dat$logFC.y, nrow(combined_dat), replace=T)
  both_sex %>% filter(logFC.x*logFC.y>0) %>% nrow()
})

ggplot(tibble("num_probes"=sim_counts_s2), aes(x=num_probes))+
  geom_histogram(binwidth=1)+
  theme_bw()+
  geom_vline(xintercept=combined_dat %>% filter(logFC.x*logFC.y>0) %>% nrow(), col="red")+
  xlab("Number of same direction genes")+
  ylab("Number of randomized runs")
ggsave("figures/ae_rep_random_int.png")


