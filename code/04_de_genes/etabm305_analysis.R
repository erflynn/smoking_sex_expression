# Code for performing differential expression analysis
# prior to running picture to run the meta, load, and process files
#
# Datasets:
# - E-TABM-305 (Illumina, n=1240), Charlesworth et al., lymphocytes



# To do: 
# Charlesworth identified many more DE genes in smokers (greater than 300), while we're getting about 60 
# - try doing analysis without the interaction term
# - try adding an age by sex term 
# - perform model comparison
# - try going through their preprocessing steps (in etabm305_process.R)
# - compare to their genes

library('tidyverse')
library('limma')
library('bigpca')
library('STRINGdb')


# --- 0. view in PC space --- #
load("data/blood/blood_proc.RData")
pDat2 <- pData(filt_ds)
start_time = Sys.time(); 
blood_bp1 <- big.PCA(exprs(filt_ds), pcs.to.keep=3); 
end_time=Sys.time()
end_time-start_time # 4 MINUTES

blood_pcs2 <- data.frame(blood_bp1$PCs)
blood_pcs2$smoking <- pDat2$smoking
blood_pcs2$sex <- pDat2$sex
blood_pcs2$age <- pDat2$age
blood_pcs2$hdl <- pDat2$hdl

ggplot(blood_pcs2, 
       aes(x=PC1, y=PC2, col=sex, shape=smoking))+
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/blood_pc.png")
ggplot(blood_pcs2, 
       aes(x=PC1, y=PC2, col=hdl))+
  geom_point(alpha=0.5)+
  theme_bw()

ggplot(blood_pcs2, 
       aes(x=PC1, y=PC2, col=age))+
  geom_point(alpha=0.5)+
  theme_bw()


# --- 1. remove undetectable probes  --- #
# compare to ctls... hmm?
exp_d <- exprs(filt_ds)
# plot(density(exp_d["trpF",]), xlim=c(0, 10))
# lines(density(exp_d["lysA",]))
# lines(density(exp_d["pheA",]))
# lines(density(exp_d["thrB",]))
# 
# lines(density(exp_d["GI_10047089-S",]), col="red")
# lines(density(exp_d["GI_10047091-S",]), col="red")
# lines(density(exp_d["GI_10047093-S",]), col="red")
# lines(density(exp_d["GI_10047099-S",]), col="red")
# lines(density(exp_d["GI_10047103-S",]), col="red")


# remove probes in the bottom 1% of variances ?
#row_vars <- apply(exp_d, 1, sd)
#quantile(row_vars, 0.01)

# include genes that are present
presentCount <- detectionCall(filt_ds)
exp_d2 <- exp_d[presentCount > 0,]
probeList <- rownames(exp_d2)

# double check if it's disproportionately X+Y probes
# ... nope, not bad - yay
gene_metadata <- read_csv("data/etabm305_adf_gene_metadata.csv")
gene_metadata %>%
  mutate(is_retained=(probe %in% probeList)) %>% 
  group_by(chromosome_name, is_retained) %>%
  count() %>% 
  filter(chromosome_name %in% c("1", "5", "20","22", "X", "Y")) %>%
  pivot_wider(names_from=is_retained, values_from=n, 
              names_prefix="kept") %>%
  mutate(frac_retained=keptTRUE/(keptTRUE+keptFALSE))

my_collapse <- function(x){
  y <- x[!is.na(x)];
  z <- sort(unique(y));
  return(paste(z, collapse=";"))
}

gene_metadata2 <- gene_metadata %>% 
  filter(probe %in% probeList) %>%
  mutate(chromosome_name=ifelse(chromosome_name %in% c(1:22, "X","Y"), 
                                chromosome_name, NA)) %>%
  group_by(probe, grp) %>%
  summarize(gene=my_collapse(gene),
            chromosome_name=my_collapse(chromosome_name)) 
stopifnot(nrow(gene_metadata2)==length(probeList))
gene_metadata_sort <- data.frame(gene_metadata2)
rownames(gene_metadata_sort) <- gene_metadata_sort$probe


# ----- 2. covariate analysis  ----- #
# - count of smoker/non-smoker by sex
pDat2 %>% 
  group_by(sex, smoking) %>% 
  count() %>%
  pivot_wider(names_from="smoking", values_from="n") %>%
  mutate(prop_smoker=smoker/(smoker+`non-smoker`))

# two proportions Z-test
# group A
# pA = fraction of female smokers, 0.166, n=721
# pB = fraction of male smokers, 0.355 n=498
# p = proportion of smokers, 0.244
# q = proportion of non-smokers, 0.756

# null hypothesis: pA=pB (two-tailed)

prop.test(x=c(120, 177), n=c(721, 498))

# - violin plots for age, HDL
ggplot(pDat2 %>% 
         unite( "grp", c(sex, smoking)), 
       aes(x=grp, y=age))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  geom_point(position=position_jitter(0.15), alpha=0.2)+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab("")
ggsave("figures/blood_age_dist.png")

ggplot(pDat2 %>% 
         unite( "grp", c(sex, smoking)), 
       aes(x=grp, y=hdl))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  geom_point(position=position_jitter(0.15), alpha=0.2)+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab("")+
  ylab("HDL")
ggsave("figures/blood_hdl_dist.png")


aov_in <- pDat %>%
  dplyr::select(smok, expr_sex, age) %>%
  filter(!is.na(age))

two.way <- aov(age ~ smoking+sex+smoking*sex, data = pDat2 %>% 
                 select(smoking, sex, age))
summary(two.way) 

two.way.hdl <- aov(hdl ~ smoking+sex+smoking*sex, data = pDat2 %>% 
                 select(smoking, sex, hdl))
summary(two.way.hdl) 

# ---- 3. run DE analysis at the probe level ---- #
design <- model.matrix(~ sex + smoking + smoking*sex + age, data=pDat2) 
fit <- lmFit(exp_d2, design)
fit <- eBayes(fit)
fit$genes <- data.frame(ID= probeList, geneSymbol=gene_metadata_sort$gene, 
                        chromosome=gene_metadata_sort$chromosome_name, 
                        stringsAsFactors=FALSE)

df_smok <- data.frame(topTable(fit, coef="smokingsmoker",  number=nrow(fit))) %>% 
  rename(gene=geneSymbol) 
df_smok %>% filter(adj.P.Val < 0.05) %>% nrow() 
# 65 probes, this is many fewer than was found before
# should we compare with a model that does not include the interaction term

# --- 4. plot DE results --- #
volcano_plot_de(data.frame(topTable(fit, coef="sexmale",  number=nrow(fit)) %>% 
                             rename(gene=geneSymbol) ))
ggsave("figures/blood_volcano_sex.png")

volcano_plot_de(data.frame(topTable(fit, coef="smokingsmoker",  number=nrow(fit)) %>% 
                             rename(gene=geneSymbol) ))
ggsave("figures/blood_volcano_smok.png")


volcano_plot_de(data.frame(topTable(fit, coef="sexmale:smokingsmoker",  number=nrow(fit)) %>% 
                             rename(gene=geneSymbol) ), num_display=10)+
  scale_color_manual(values=c("dark gray"))
ggsave("figures/blood_volcano_interaction.png")

plotChrDist(data.frame(topTable(fit, coef="sexmale",  number=nrow(fit))) %>%
              rename(gene=geneSymbol, chromosome_name=chromosome), pcut=0.01)
ggsave("figures/blood_chr_sex.png")

plotChrDist(data.frame(topTable(fit, coef="smokingsmoker",  number=nrow(fit))) %>%
              rename(gene=geneSymbol, chromosome_name=chromosome))
ggsave("figures/blood_chr_smok.png")

# --- 5. run meta-analysis to convert probes to genes --- #
smok_ci_int_p <- topTable(fit, coef="sexmale:smokingsmoker", number=nrow(fit), 
                          confint = TRUE)
smok_ci_int_p <- data.frame(smok_ci_int_p)
comb_ma_dat <- run_clean_ma(smok_ci_int_p)
comb_ma_dat %>% select(-src) %>% head(10)

smok_ci_sex <- topTable(fit, coef="sexmale", number=nrow(fit), 
                          confint = TRUE)
sex_ma <- run_clean_ma(data.frame(smok_ci_sex))
sex_ma %>% select(-src) %>% head(10)

smok_ci_smok <- topTable(fit, coef="smokingsmoker", number=nrow(fit), 
                        confint = TRUE)
smok_ma <- run_clean_ma(data.frame(smok_ci_smok))
smok_ma %>% select(-src) %>% head(10)
smok_ma %>% filter(adj.p < 0.05) # 45 (their FDR cutoff)
smok_ma %>% filter(p < 0.001) # 68 (their nominal P cutoff)


# ---- 6. visualize some interaction effects --- #
expDat4 <- data.frame(exp_d2)
expDat4$probe <- rownames(exp_d2)
exp_long <- expDat4 %>% 
  left_join(gene_metadata_sort, by="probe") %>%
  select(-grp, -chromosome_name, -probe) %>%
  pivot_longer( -gene, names_to="sampleID", values_to="expr") %>%
  left_join(pDat2, by=c("sampleID"))

# TODO - figure out how to add lines!
df_exp <- exp_long %>%
  filter(!is.na(gene), !is.na(sex), !is.na(smoking)) %>%
  filter(gene %in% head(smok_ci_int_p$geneSymbol, )) %>%
  unite(grp, c("smoking", "sex"), remove=FALSE) 
df_exp %>%
  ggplot(aes(y=expr, x=smoking, fill=sex, group=grp))+
  geom_boxplot()+
  facet_wrap(~gene, scales="free")+
  theme_bw()+
  ylab("")+
  xlab("")

df_exp %>%
  ggplot(aes(y=expr, x=smoking))+
  geom_boxplot()+
  facet_wrap(~gene, scales="free")+
  theme_bw()+
  ylab("")+
  xlab("")

smok_ci_int_p %>% filter(geneSymbol=="SLC1A5")


design <- model.matrix(~ smoking, data=pDat2) 
fit <- lmFit(exp_d2, design)
fit <- eBayes(fit)
fit$gene <- data.frame(ID= probeList, geneSymbol=gene_metadata_sort$gene, 
                        chromosome=gene_metadata_sort$chromosome_name, stringsAsFactors=FALSE)

df_smok <- data.frame(topTable(fit, coef="smokingsmoker",  number=nrow(fit), confint=T)) %>% 
  rename(gene=geneSymbol)  
df_smok %>% filter(adj.P.Val < 0.05) %>% nrow() 
df_smok %>% filter(gene=="SLC1A5")
df_exp %>% 
  filter(gene=="SLC1A5") %>%
  group_by(smoking, sex) %>%
  summarize(mean_expr=mean(expr), sd_expr=sd(expr))
# mean_f(s-ns) -0.15, mean_m(s-ns): 0.06
# diff means: -0.15-0.06=0.21
df_exp %>% 
  filter(gene=="PFAS") %>%
  group_by(smoking, sex) %>%
  summarize(mean_expr=mean(expr), sd_expr=sd(expr))
smok_ci_int_p %>% filter(geneSymbol=="PFAS")

# mean_f(s-ns) 0.07, mean_m(s-ns): -0.06
# diff means: 0.07--0.06=0.13

df_exp %>% 
  filter(gene=="OLR1") %>%
  group_by(smoking, sex) %>%
  summarize(mean_expr=mean(expr), sd_expr=sd(expr))
smok_ci_int_p %>% filter(geneSymbol=="OLR1")
#(4.95-5.09)-(5.22-4.97)
#-0.14-0.25=0.39


# DATA is *already* log-transformed
sd(exp_long$expr)

#  ---- 7. prepare data for STAMS ---- #
# TODO - check STAMS mapping - though everything mapped?

map_to_stams <- function(df, fname){
  string_db <- STRINGdb$new() 
  gene_values <- df %>% 
    dplyr::rename(Pvalue=p, Gene=gene) 
  gene_values2=gene_values %>% 
    filter(Pvalue != 1) %>% 
    filter(!is.na(Gene)) 
  gene_mapped = string_db$map(data.frame(gene_values2), "Gene", removeUnmappedRows = TRUE )
  save(gene_mapped, file=fname)
}
map_to_stams(comb_ma_dat, "data/gene_mapped_blood_int.RData")
map_to_stams(sex_ma, "data/gene_mapped_blood_sex.RData")
map_to_stams(smok_ma, "data/gene_mapped_blood_smok.RData")


# ---- 8. run immune cell deconvolution ---- #

# TODO - move this to a separate file
runIS <- function(gseMetaObj){
  dt <- MetaIntegrator:::.extractDataForGenesDT(gseMetaObj, promiscProbes = F)
  outDT <- as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix, dt), keep.rownames = T)
  outDT[, `:=`(natural_killer_cell, CD56bright_natural_killer_cell + 
                 CD56dim_natural_killer_cell)]
  outDT[, `:=`(monocyte, CD14_positive_monocyte + CD16_positive_monocyte)]
  outDT[, `:=`(B_cell, naive_B_cell + memory_B_cell + 
                 plasma_cell)]
  outDT[, `:=`(T_cell, CD8_positive_alpha_beta_T_cell + 
                 CD4_positive_alpha_beta_T_cell + gamma_delta_T_cell)]
  outDT[, `:=`(granulocyte, neutrophil + eosinophil + 
                 basophil)]
  
  return(outDT %>% as_tibble())
}

prepDeIS <- function(is_out, phe_data){
  is_out %>%
    left_join(phe_data, by=c("rn"="geo_accession")) %>%
    select(-`P-value`, -Correlation, -RMSE, -rn) %>%
    pivot_longer(-time_pt, names_to="cell_type", values_to="proportion") %>%
    group_by(cell_type) %>%
    mutate(sd=sd(proportion)) %>%
    filter(sd > 0)
}

runDeIS <- function(dat_mat){
  split_dat <- dat_mat %>% group_split(cell_type)
  cell_types <-unique(dat_mat$cell_type)
  names(split_dat) <- cell_types
  mod_res <- lapply(split_dat, function(x){
    lmod <- lm(proportion ~ factor(time_pt), data=x)
    return(summary(lmod)$coefficients[2,])
  })
  mod_res2 <- do.call(rbind, mod_res) %>% as_tibble()
  mod_res2$cell_type <- cell_types
  mod_res2 %>% 
    arrange(`Pr(>|t|)`)
}

talbi_is_out <- runIS(talbi_v2$originalData$GSE4888)
dat_mat  <- prepDeIS(talbi_is_out, data.frame(talbi_v2$originalData$GSE4888$pheno))
mod_res2 <- runDeIS(dat_mat)
mod_res2 %>% filter(`Pr(>|t|)`< 0.05/(nrow(mod_res2)))



