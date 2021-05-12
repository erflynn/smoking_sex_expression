
prep_meta_study <- function(ds, study){
  ds %>%
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(SD=(logFC-logFC.l)/1.96) %>%
    mutate(ID=study) %>%
    dplyr::select(gene, logFC, SD, ID)
}

run_meta <- function(df){
  df<- df %>% filter(!is.na(logFC))
  
  return(meta::metagen(df %>% pull(logFC), # treatment estimate
                       df %>% pull(SD), # standard error
                       studlab=df %>% pull(ID),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       method.tau = "DL", # method for between study variance
                       hakn = FALSE,
                       prediction = FALSE,
                       sm = "MD") )
}
meta_study <- function(df){
  ma <- run_meta(df)
  # ma <- meta::metagen(df %>% pull(logFC), # treatment estimate
  #                           df %>% pull(SD), # standard error
  #                           studlab=df %>% pull(ID),
  #                           comb.fixed = FALSE,
  #                           comb.random = TRUE,
  #                           method.tau = "DL", # method for between study variance
  #                           hakn = FALSE,
  #                           prediction = FALSE,
  #                           sm = "MD") 
  return(list("gene"=unique(df$gene),
              "logFC"=ma$TE.random,
              "logFC.l"=ma$lower.random,
              "logFC.u"=ma$upper.random,
              "p"=ma$pval.random,
              "n"=unique(df$n)))
}
valid_studies <- c("GSE7895", "GSE27002", "GSE21862", "E-MTAB-5279")

meta_in <- data.table::rbindlist(lapply(all_studies2, function(x) 
  prep_meta_study(x$gene_smok, x$study)))

meta_in_disc <- meta_in %>% filter(!ID %in% valid_studies)
meta_in_valid <- meta_in %>% filter(ID %in% valid_studies)

meta_in2 <- meta_in_disc %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  arrange(gene)
meta_in3 <- meta_in2 %>% filter(n >=15) %>% group_split(gene)
names(meta_in3) <- unique(meta_in2 %>% filter(n>=15) %>% pull(gene)) # 14367

ma_vals <- lapply(meta_in3, meta_study)

# 5min
mult_study <- data.frame(apply(data.table::rbindlist(ma_vals) ,
                               c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study$adj.p <- p.adjust(mult_study$p, method="fdr")
save(mult_study, file="data/results/meta_across_sm2.RData")

mult_study %>% filter(adj.p < 0.05) %>% nrow() # 48
sig_genes <- mult_study %>% filter(adj.p < 0.05 & abs(logFC)>=0.3) # 7
valid_g <- do.call(rbind, 
                   lapply(valid_studies,
                          function(x) all_studies2[[x]]$gene_smok %>% 
                            filter(gene %in% as.character(sig_genes$gene)) %>%
                            mutate(study=x) %>%
                            ungroup() %>%
                            dplyr::select(gene, logFC, p, adj.p, study)))
valid_g %>% filter(p < 0.05/nrow(valid_g))

#gene   logFC        p    adj.p study      
# NQO1   1.99  4.83e-20 6.00e-16 GSE7895    
# CYP1B1 2.61  7.71e- 5 8.47e- 2 GSE27002   
# LRRN3  1.58  1.88e- 6 1.16e- 2 E-MTAB-5279
# AHRR   0.878 7.84e- 4 1.26e- 1 E-MTAB-5279

valid_g %>% group_by(gene) %>% arrange(gene, p) %>%
  filter(p < 0.05) %>%
  dplyr::select(-adj.p) %>%
  write_csv("data/supp_tables/meta_smok_valid.csv")
sig_genes %>%
  write_csv("data/supp_tables/meta_smok.csv")




# sex
meta_in_sex <- data.table::rbindlist(lapply(all_studies2[! names(all_studies2) %in% c("GSE16149", "GSE18723", "GSE87072")], function(x) 
  prep_meta_study(x$gene_sex, x$study)))

meta_in_disc_s <- meta_in_sex %>% filter(!ID %in% valid_studies)
meta_in_valid_s <- meta_in_sex %>% filter(ID %in% valid_studies)

meta_in2_s <- meta_in_disc_s %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  arrange(gene)
summary(meta_in2_s$n)
meta_in3_s <- meta_in2_s %>% filter(n >=10) %>% group_split(gene)
names(meta_in3_s) <- unique(meta_in2_s %>% filter(n>=10) %>% pull(gene)) # 14367

ma_vals_s <- lapply(meta_in3_s, meta_study)

# 5min
mult_study_s <- data.frame(apply(data.table::rbindlist(ma_vals_s) ,
                                 c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study_s$adj.p <- p.adjust(mult_study_s$p, method="fdr")
save(mult_study_s, file="data/results/meta_across_sm2_sex.RData")

mult_study_s %>% filter(adj.p < 0.05) %>% nrow() # 60
mult_study_s %>% filter(adj.p < 0.05 & abs(logFC)>=0.3) #22
mult_study_s %>% 
  filter(adj.p < 0.05) %>%
  left_join(xy_chr_map, by=c("gene"="hgnc_symbol")) %>%
  filter(!chromosome_name %in% c("X", "Y"))
sig_genes_sex <- mult_study_s %>% filter(adj.p < 0.05 & abs(logFC)>=0.3)





plot_gene(meta_in2_s, "DRAM2")


# forestPlot
plot_gene <- function(my_meta, my_gene){
  df <- my_meta %>% filter(gene==my_gene)
  ma <- run_meta(df)
  meta::forest.meta(ma)
}

plot_gene(meta_in2, "GZMA")

plot_gene(meta_in2, "LRRN3")


# 5. meta grouped by tissue
get_meta_out <- function(list_g){
  ma_vals <- lapply(list_g, meta_study)
  
  # 5min
  mult_study <- data.frame(apply(data.table::rbindlist(ma_vals) ,
                                 c(1,2), unlist)) %>%
    mutate(across(c(contains("logFC"), p), 
                  ~as.numeric(as.character(.)) )) %>%
    mutate(n=as.numeric(as.character(n))) %>% 
    arrange(p)
  mult_study$adj.p <- p.adjust(mult_study$p, method="fdr")
  return(mult_study)
}

# AE/BE
# TODO - remove validations studies
ae_studies <-list_studies %>% 
  filter(!study %in% valid_studies) %>%
  filter(str_detect(tissue,"airway") | str_detect(tissue, "trachea")) %>%
  pull(study)
ae_meta <- meta_in %>% filter(ID %in% ae_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=4) %>%
  group_split(gene)
ae_meta_out <- get_meta_out(ae_meta)
ae_sig <- ae_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
ae_sig # 18

# blood or blood component
blood_studies <-list_studies %>% 
  filter(!study %in% valid_studies) %>%
  filter(str_detect(tissue, "blood"), !str_detect(tissue, "b cell") ) %>%
  pull(study) # PBMC, whole blood
blood_meta <- meta_in %>% 
  filter(ID %in% blood_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=5) %>%
  group_split(gene)
blood_meta_out <- get_meta_out(blood_meta)
blood_meta_out %>% filter(adj.p < 0.05) # 15
blood_sig <- blood_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
blood_sig
plot_gene(meta_in %>% filter(ID %in% blood_studies), "FAM81B")


# compare - 314 genes
meta2 <- ae_meta_out %>% 
  left_join(blood_meta_out, by="gene") %>%
  filter(adj.p.x < 0.1 | adj.p.y < 0.1) 
meta2 %>% 
  filter(logFC.x*logFC.y > 0 & p.x < 0.05 & 
           p.y < 0.05) %>%
  dplyr::select(gene, logFC.x, logFC.y, p.x, p.y) # 22
meta2 %>% filter(logFC.x*logFC.y < 0) %>% nrow()

# meta-analyze by sex
ae_meta_sex <- meta_in_sex %>% 
  filter(ID %in% ae_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=4) %>%
  group_split(gene)
ae_meta_out_s <- get_meta_out(ae_meta_sex)
ae_sig_sex <- ae_meta_out_s %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
ae_sig_sex

blood_meta_sex <- meta_in_sex %>% 
  filter(ID %in% blood_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=5) %>%
  group_split(gene)
blood_meta_out_s <- get_meta_out(blood_meta_sex)
blood_meta_out_s %>% filter(adj.p < 0.05) # 15
blood_sig_sex <- blood_meta_out_s %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)

# save
save(ae_meta_out, blood_meta_out, ae_meta_out_s, blood_meta_out_s, file="data/tissue_spec_ma.RData")


# add chromosomes to all the sig genes
sig_df <- do.call(rbind, 
                  list(sig_genes %>% mutate(term="smoking", ds="all"),
                       sig_genes_sex  %>% mutate(term="sex", ds="all"),
                       blood_sig  %>% mutate(term="smoking", ds="blood"),
                       blood_sig_sex  %>% mutate(term="sex", ds="blood"),
                       ae_sig %>% mutate(term="smok", ds="ae"),
                       ae_sig_sex %>% mutate(term="sex", ds="ae")))
chr_map <- grab_chr(sig_df) 
sig_df2 <- sig_df %>% 
  left_join(chr_map %>% filter(chromosome_name %in% c(1:22, "X", "Y")), by=c("gene"="hgnc_symbol")) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene == "TTTY15" ~ "Y",
    gene=="MARCH2" ~ "19"
  ))

# - write out output

sig_meta_comb <- sig_df2 %>% 
  as_tibble() %>%
  dplyr::select(ds, term, everything()) %>%
  mutate(across(c(contains("logFC"), contains("p")), ~signif(., 3)))
sig_meta_comb %>% write_csv("data/supp_tables/sig_meta_comb.csv")


ae_de_gene_s %>% 
  filter(adj.p < 0.05 & abs(logFC) >= 0.3) %>% 
  left_join(chr_ae %>% filter(chromosome_name %in% c(1:22, "X", "Y")),
            by=c("gene"="hgnc_symbol")) %>% 
  write_csv("data/supp_tables/ae_genes_sex.csv")

# look at validation

valid_g <- do.call(rbind, 
                   lapply(valid_studies,
                          function(x) all_studies2[[x]]$gene_smok %>% 
                            filter(gene %in% as.character(sig_genes$gene)) %>%
                            mutate(study=x) %>%
                            ungroup() %>%
                            dplyr::select(gene, logFC, p, adj.p, study)))
valid_g2 <- valid_g %>% 
  group_by(gene) %>% 
  arrange(gene, p) %>%
  filter(p < 0.05/nrow(valid_g)) %>%
  dplyr::select(-adj.p)
sig_genes %>% inner_join(valid_g2, by="gene"  ) %>%
  dplyr::select(gene, logFC.x, logFC.y) %>%
  filter(logFC.x*logFC.y > 0)
valid_g2 %>% write_csv("data/supp_tables/meta_smok_valid.csv")

valid_g_sex <- do.call(rbind, 
                       lapply(valid_studies,
                              function(x) all_studies2[[x]]$gene_sex %>% 
                                filter(gene %in% as.character(sig_genes_sex$gene)) %>%
                                mutate(study=x) %>%
                                ungroup() %>%
                                dplyr::select(gene, logFC, p, adj.p, study)))
valid_g_sex2 <- valid_g_sex %>% 
  filter(p < 0.05/nrow(valid_g_sex)) %>%
  arrange(gene, p) %>%
  dplyr::select(-adj.p)

valid_g_sex2 %>% 
  left_join(sig_genes_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

ae_valid <- all_studies2[["GSE7895"]]$gene_smok %>% 
  filter(gene %in% as.character(ae_sig$gene)) %>%
  ungroup() %>%
  dplyr::select(gene, logFC, p) 
ae_valid2 <- ae_valid %>% 
  filter(p < 0.05/nrow(ae_valid)) %>% 
  arrange(gene, p) %>%
  mutate(study="GSE7895")# 21
ae_valid2 %>% 
  left_join(ae_sig, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

ae_valid_sex <- all_studies2[["GSE7895"]]$gene_sex %>% 
  filter(gene %in% as.character(ae_sig_sex$gene)) %>%
  ungroup() %>%
  dplyr::select(gene, logFC, p)
ae_valid_sex2 <- ae_valid_sex %>%
  filter(p < 0.05/nrow(ae_valid_sex)) %>% 
  arrange(gene, p) %>%
  mutate(study="GSE7895")

ae_valid_sex2 %>%
  left_join(ae_sig_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

blood_valid <- do.call(rbind, 
                       lapply(c("GSE21862", "E-MTAB-5279"),
                              function(x) all_studies2[[x]]$gene_smok %>% 
                                filter(gene %in% as.character(blood_sig$gene)) %>%
                                mutate(study=x) %>%
                                ungroup() %>%
                                dplyr::select(gene, logFC, p, study)))
blood_valid2 <- blood_valid %>% 
  filter(p < 0.05/nrow(blood_valid)) %>%
  arrange(gene, p) # 3

blood_valid2 %>%
  left_join(blood_sig, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none


blood_valid_sex <- do.call(rbind, 
                           lapply(c("GSE21862", "E-MTAB-5279"),
                                  function(x) all_studies2[[x]]$gene_sex %>% 
                                    filter(gene %in% as.character(blood_sig_sex$gene)) %>%
                                    mutate(study=x) %>%
                                    ungroup() %>%
                                    dplyr::select(gene, logFC, p, study)))
blood_valid_sex2 <- blood_valid_sex %>% 
  filter(p < 0.05/nrow(blood_valid_sex)) %>%
  arrange(gene, p) 

blood_valid_sex2 %>%
  left_join(blood_sig_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

all_smok <- sig_genes %>% 
  left_join(valid_g2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

all_sex <- sig_genes_sex %>% 
  left_join(valid_g_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")



blood_smok <- blood_sig %>% 
  left_join(blood_valid2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

blood_sex <- blood_sig_sex %>% 
  left_join(blood_valid_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")

ae_smok <- ae_sig %>% 
  left_join(ae_valid2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

ae_sex <- ae_sig_sex %>% 
  left_join(ae_valid_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")



sig_df_v <- do.call(rbind, 
                    list(all_smok %>% mutate(term="smoking", ds="all"),
                         all_sex  %>% mutate(term="sex", ds="all"),
                         blood_smok  %>% mutate(term="smoking", ds="blood"),
                         blood_sex  %>% mutate(term="sex", ds="blood"),
                         ae_smok %>% mutate(term="smok", ds="ae"),
                         ae_sex %>% mutate(term="sex", ds="ae")))
sig_df2_v <- sig_df_v %>% 
  left_join(chr_map %>% filter(chromosome_name %in% c(1:22, "X", "Y")), by=c("gene"="hgnc_symbol")) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene == "TTTY15" ~ "Y",
    gene=="MARCH2" ~ "19"
  ))

# - write out output

sig_meta_comb_v <- sig_df2_v %>% 
  as_tibble() %>%
  rename(valid_study=study) %>%
  dplyr::select(ds, term, gene, chromosome_name, everything()) %>%
  mutate(across(c(contains("logFC"), contains("p")), ~signif(., 3)))
sig_meta_comb_v %>% write_csv("data/supp_tables/sig_meta_comb_valid.csv")
sig_meta_comb_v %>% 
  group_by(ds, term) %>%
  dplyr::summarize(n=n(), n_valid=sum(!is.na(logFC.valid)))
# ds    term        n n_valid
# <chr> <chr>   <int>   <int>
#   1 ae    sex         6       4
# 2 ae    smok       66      21
# 3 all   sex        22      20
# 4 all   smoking     7       4
# 5 blood sex        32      13
# 6 blood smoking    19       3

bs_auto <- sig_meta_comb_v %>% 
  filter(ds=="blood", term=="sex") %>%
  mutate(auto=(!chromosome_name %in% c("X", "Y"))) %>%
  dplyr::select(-term, -ds, -logFC.l, -logFC.u, -chromosome_name, -n, -p)  
table(bs_auto$auto)
# overlap?
intersect(blood_sig$gene, sig_genes$gene) # GZMH
intersect(ae_sig$gene, sig_genes$gene) # NQO1
intersect(blood_sig$gene, ae_sig$gene) # AKR1C3 -- opp dir
intersect(ae_sig$gene, sig_genes$gene) # NQO1
intersect(blood_sig_sex$gene, ae_sig_sex$gene) # 4



# heatmap for gene data  #

blood_plot_dat <- meta_in %>% filter(gene %in% blood_sig$gene, 
                                     ID %in% blood_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(blood_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 

blood_plot_dat_w <- blood_plot_dat %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
blood_plot_dat2 <- blood_plot_dat_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=blood_plot_dat_w$gene[order(blood_plot_dat_w$pooled)]))


ggplot(blood_plot_dat2, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC", 
                       limit=c(-1.75, 0.75))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/blood_meta_heatmap2.png")


ae_plot_dat <- meta_in %>% 
  filter(gene %in% ae_sig$gene, ID %in% ae_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(ae_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
ae_plot_dat_w <- ae_plot_dat %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
ae_plot_dat2 <- ae_plot_dat_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=ae_plot_dat_w$gene[order(ae_plot_dat_w$pooled)]))

ggplot(ae_plot_dat2, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/ae_meta_heatmap2.png")


# by sex
blood_plot_dat_s <- meta_in_sex %>% filter(gene %in% blood_sig_sex$gene, 
                                           ID %in% blood_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(blood_sig_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 

blood_plot_dat_s_w <- blood_plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
blood_plot_dat_s_2 <- blood_plot_dat_s_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=blood_plot_dat_s_w$gene[order(blood_plot_dat_s_w$pooled)]))


ggplot(blood_plot_dat_s_2 %>% 
         mutate(logFC=ifelse(logFC < -6, logFC+4, logFC)), aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/blood_meta_heatmap2_sex.png")


ae_plot_dat_s <- meta_in_sex %>% 
  filter(gene %in% ae_sig_sex$gene, ID %in% ae_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(ae_sig_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
ae_plot_dat_w_s <- ae_plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
ae_plot_dat2_s <- ae_plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=ae_plot_dat_w_s$gene[order(ae_plot_dat_w_s$pooled)]))

ggplot(ae_plot_dat2_s, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/ae_meta_heatmap2_sex.png")


plot_dat_s <- meta_in %>% 
  filter(gene %in% sig_genes$gene ) %>%
  dplyr::select(-SD) %>%
  bind_rows(sig_genes %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
plot_dat_w_s <- plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
plot_dat_w_s2 <- plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, 
                     levels=plot_dat_w_s$gene[order(plot_dat_w_s$pooled)])) 
plot_dat <- plot_dat_w_s2 %>% 
  rename(study=ID) %>%
  left_join(list_studies %>% 
              dplyr::select(study, tissue, samples) ) %>%
  mutate(study=as.character(study), tissue=as.character(tissue)) %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study)
plot_dat$study <- factor(plot_dat$study, levels=unique(plot_dat$study)  )

ggplot(plot_dat, aes(y=study, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/heatmap_all.png")
color_bar(plot_dat)
ggsave("figures/paper_figs/heatmap_all_colorbar.png")


plot_dat_s <- meta_in_sex %>% 
  filter(gene %in% sig_genes_sex$gene ) %>%
  dplyr::select(-SD) %>%
  bind_rows(sig_genes_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
plot_dat_w_s <- plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
plot_dat_w_s2 <- plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, 
                     levels=plot_dat_w_s$gene[order(plot_dat_w_s$pooled)])) 
plot_dat <- plot_dat_w_s2 %>% 
  rename(study=ID) %>%
  left_join(list_studies %>% 
              dplyr::select(study, tissue, samples) ) %>%
  mutate(study=as.character(study), tissue=as.character(tissue)) %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study)
plot_dat$study <- factor(plot_dat$study, levels=unique(plot_dat$study)  )

ggplot(plot_dat, aes(y=study, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")

ggsave("figures/paper_figs/heatmap_all_sex.png")
color_bar(plot_dat)
ggsave("figures/paper_figs/heatmap_all_sex_colorbar.png")


### FOREST PLOTS
my_gene <- "AHRR"
plot_gene2 <- function(my_gene){
  gene_df <- meta_in %>% filter(gene %in% my_gene ) %>%
    dplyr::rename(study=ID) %>%
    bind_rows(mult_study %>% 
                mutate(SD=(logFC-logFC.l)/1.96) %>%
                dplyr::select(gene, logFC, SD) %>% 
                filter(gene %in% my_gene) %>%
                mutate(study="pooled"))  %>%
    as_tibble() %>% left_join(list_studies %>% 
                                dplyr::select(study, tissue, samples) ) %>%
    mutate(study=as.character(study), tissue=as.character(tissue)) %>%
    filter(study!="GSE44456") %>%
    mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
    mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                          "nasal epithelium", "oral cavity",
                                          "buccal mucosa", "sputum", "alveolar macrophages", 
                                          "lung", 
                                          "blood - b cells", "blood - pbmcs", "blood - whole",
                                          "liver", "kidney", 
                                          "brain - prefrontal cortex"))) %>%
    arrange(tissue, study)
  gene_df$study <- factor(gene_df$study, levels=unique(gene_df$study)  )
  ggplot(gene_df, aes(x=logFC, y=study))+
    geom_vline(xintercept=0, lty=2, col="darkgray") +
    geom_point(aes(size=samples), alpha=0.6)+
    geom_errorbarh(aes(xmin=logFC-1.96*SD, xmax=logFC+1.96*SD))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ylab("")+
    facet_grid(.~gene, scales="free")
}
plot_gene2(sig_genes$gene)
ggsave("figures/all_forest.png")

ggplot(gene_df)+
  geom_bar(mapping = aes(x = study, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=c("black", paired3))
ggsave("figures/smok_bar.png")

plot_gene2("AHRR")
ggsave("figures/ahrr_forest_plot2.png")  

plot_gene2("LRRN3")
ggsave("figures/lrrn3_forest_plot2.png")  

plot_gene2("CYP1B1")
ggsave("figures/cy1b1_forest_plot2.png")  

plot_gene2("NQO1")
ggsave("figures/nqo1_forest_plot2.png") 

plot_gene2("CCL4")
ggsave("figures/nccl4_forest_plot2.png")

plot_gene2("GZMH")
ggsave("figures/gzmh_forest_plot2.png") 

#facet_grid(tissue~., scales="free")
plot_gene(meta_in, "AHRR")
