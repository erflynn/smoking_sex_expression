library(tidyverse)
load("data/all_studies.RData") # all_studies


fix_colnames_ds <- function(ds){
  if (!"logFC.l" %in% colnames(ds$gene_smok)){
    
    ds$gene_smok <- ds$gene_smok %>%
      dplyr::rename(gene=geneSymbol, p=P.Value,
                    probes=ID) %>%
      mutate(adj.p=p.adjust(p, "fdr")) %>%
      mutate(logFC.l=logFC-SD*1.96)
    if (!is.na(ds$gene_int)){
      ds$gene_int <- ds$gene_int %>%
        dplyr::rename(gene=geneSymbol, p=P.Value,
                      probes=ID) %>%
        mutate(adj.p=p.adjust(p, "fdr")) %>%
        mutate(logFC.l=logFC-SD*1.96) 
    }
  
  }
    return(ds)
}
all_studies <- lapply(all_studies, fix_colnames_ds)
names(all_studies) <-  lapply(all_studies, function(x) x$study)
# table info:
# study, prop_var_smok, prop_var_sex, prop_var_int, # de probes, # de genes, # de probes int, # de genes int

table_ds <- function(ds){
  num_tot_probes <- ds$probe_smok %>% nrow()
  num_tot_genes <- ds$gene_smok %>% nrow()
  num_de_probes <- ds$probe_smok %>%
    filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) %>%
    nrow()

  num_de_genes <- ds$gene_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    nrow()
  
  if (!is.na(ds$probe_int)){
    num_de_probes_int <- ds$probe_int %>%
      filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) %>%
      nrow()

    num_de_genes_int <- ds$gene_int %>%
      filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
      nrow()
  } else {
    num_de_probes_int <- NA
    num_de_genes_int <- NA
  }
  
  if (!is.na(ds$var)){
    var_smok =  ds$var[["smoking"]]
    if ("sex" %in% names(ds$var)){
      var_sex = ds$var[["sex"]]
      var_int= ds$var[["sex:smoking"]]
    } else {
      var_sex =NA
      var_int= NA
    }
  } else {
    var_smok = NA
    var_sex = NA
    var_int = NA
  }
  return(list(ds$study, var_smok, var_sex, var_int, num_tot_probes, num_tot_genes, num_de_probes, 
              num_de_genes, num_de_probes_int, num_de_genes_int))
}

my_tab <- do.call(rbind, lapply(all_studies, table_ds))
data_df <- data.frame(my_tab)
colnames(data_df) <- c("study",  "var_smok", "var_sex", "var_int", 
                       "nprobes", "ngenes", "de_probes", 
                       "de_genes", "de_probes_int", "de_genes_int")
data_df2 <- data_df %>% 
  mutate(across(-study, ~as.numeric(as.character(.)))) %>%
  as_tibble() %>%
  mutate(study=unlist(study))
list_studies <- read_csv("data/supp_tables/s4_list_smok_studies_formatted.csv")

list_studies <- list_studies %>% mutate(
  tissue=case_when(
    tissue == "pbmc" ~ "blood - pbmc",
    tissue=="whole blood" ~ "blood - whole",
    tissue=="b cell" ~ "blood - b cell",
    str_detect(tissue, "bronchial") ~ "airway epithelium",
    str_detect(tissue, "trachea") ~ "trachea epithelium",
    tissue=="oral" ~ "oral cavity",
    TRUE ~ tissue
    )
)

list_studies <- list_studies %>% filter(study != "GSE55962")
list_studies$tissue <- factor(list_studies$tissue, levels=c("airway epithelium",
                                           "trachea epithelium",
                                           "nasal epithelium",
                                           "oral cavity",
                                           "buccal mucosa",
                                           "alveolar macrophages",
                                           "lung",
                                           "blood - b cell", 
                                           "blood - pbmc",
                                           "blood - whole",
                                           "kidney", 
                                           "hippocampus"))
list_studies$tissue
data_df3 <- data_df2 %>%
  mutate(across(contains("var"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, `follow up analysis by sex`)) %>%
  dplyr::select(study, tissue, `follow up analysis by sex`, everything()) %>%
  arrange(tissue, study) 
data_df3 %>% write_csv("data/supp_tables/small_study_summary.csv")

df <- ae_l2 %>% left_join(list_studies %>% 
                            dplyr::select(study, samples, females, smokers)) %>%
  mutate(frac_smokers=smokers/samples) %>%
  filter(!is.na(samples))
cor.test(df$samples, df$de_genes)
cor.test(df$frac_smokers, df$de_genes)

ae_l2 %>% left_join(list_studies %>% 
                      dplyr::select(study, samples, females, smokers)) %>%
  mutate(frac_smokers=smokers/samples) %>%
  mutate(`fraction female`=females/samples) %>%
  ggplot(aes(x=frac_smokers, y=de_genes, col=`fraction female`))+
  geom_point(aes(size=samples), alpha=0.9)+
  theme_bw()+
  scale_fill_continuous(low = "purple", high = "green", limit = c(0,1))+
  xlab("fraction smokers")+
  ylab("number of DE genes")
ggsave("figures/paper_figs/sample_frac.png")  


ae_l2 %>% left_join(list_studies %>% 
                      dplyr::select(study, samples, females)) %>%
  mutate(frac_females=females/samples) %>%
  ggplot(aes(x=frac_females, y=de_genes, size=samples))+
  geom_point()+
  theme_bw()+
  xlab("fraction female")+
  ylab("number of DE genes")

l_b <- read_csv("data/supp_tables/s4b_breakdown_smok_sex.csv")
ggplot(l_b %>% pivot_longer(contains("male"), names_to="group", values_to="count") %>%
         mutate(group=str_replace_all(group,"_NS", " non-smokers"),
                group=str_replace_all(group, "_S", " smokers")), 
       aes(x=study_acc, y=count,fill=group))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("figures/paper_figs/sex_breakdown_sm.png")
# fraction male vs fraction female


ae_l <-  list("Grouped AE", "airway epithelium", "Y", 0.0215, 0.0211, 0.00666, 54675,20183,2625,932,1,30)
names(ae_l) <- colnames(data_df3)
ae_l2 <- rbind(data_df3, ae_l) %>% arrange(tissue, study)
# ADD AE HERE!

ae_l2$study <- factor(ae_l2$study, levels=ae_l2$study)
ae_l2 %>% 
  dplyr::select(study, tissue, contains("var")) %>%
  pivot_longer(contains("var"), names_to="term", values_to="variance") %>%
  mutate(term=case_when(
    term=="var_int" ~ "smoking*sex",
    term=="var_sex" ~ "sex",
    term=="var_smok" ~ "smoking"
  )) %>%
  ggplot(aes(x=1, y=variance, fill=term)) +
  geom_bar(stat="identity", position=position_dodge(0.8))+
  xlab("")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  coord_flip()+
  facet_grid(study~., switch="y")+
  theme(strip.text.y.left = element_text(angle = 0))
ggsave("figures/paper_figs/var_frac.png")

ggplot(ae_l2, aes(x=study, y=de_genes))+
  geom_bar(stat="identity")+
  theme_bw()

ggplot(ae_l2, aes(x=study, y=de_genes_int))+
  geom_bar(stat="identity")+
  theme_bw()

# add AE to the list?
data_df3 %>% filter(de_probes==0 | de_genes==0) # 4 have both, 9 have one of them


# --- ae data --- #
load("data/results/ae_int.RData") 
load("data/results/ae_smok.RData")
probe_int_ae <- smok_ci_int_p
gene_int_ae <- comb_ma_dat %>%
  ungroup() %>%
  dplyr::select(-chromosome)
probe_smok_ae <- smok_ci_smok 
gene_smok_ae <- ma_smok %>%
  ungroup() %>%
  dplyr::select(-chromosome)


# overlap + correlation of ES?
ae_de_probe <- probe_smok_ae %>% 
  filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) # 2625
# 932

ae_de_probe_i <- probe_int_ae %>% 
  filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 )  # 1
ae_de_gene_i <- gene_int_ae %>%
  filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) # 30

# // TODO cleanup and write these out again?

get_pair_overlap <- function(ds1, study1, ds2, study2){
  # look at gene overlap
  ds_smok <- ds1 %>% 
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(id=1:nrow(ds1)) %>% 
    separate_rows(gene, sep=",") 
  ds_de_gene <- ds_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    pull(gene)
  ds2_de_gene <- ds2 %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 )  
  
  overlap <- ds_smok %>%
    filter(gene %in% ds_de_gene | gene %in% ds2_de_gene$gene) %>%
    inner_join(ds2, by="gene") # 702

  ggplot(overlap, aes(x=logFC.x, y=logFC.y))+
    geom_point(alpha=0.5)+
    theme_bw()
  if (nrow(overlap) < 10){
   cor <- NA
   corp <- NA
  } else {
    pcor <- cor.test(overlap$logFC.x, overlap$logFC.y, method="pearson") # 0.00042, p=0.99
    cor <- pcor$estimate[["cor"]]
    corp <- pcor$p.value
  }

  same_dir <- overlap %>% filter(logFC.x*logFC.y > 0) # 345/702
  rep_sig <- same_dir %>% filter(p.x < 0.05/nrow(overlap) & p.y < 0.05/nrow(overlap))
  return(list("study1"=study1, "study2"=study2, "tot"=nrow(overlap), "cor"=cor, "cor.p"=corp,
              "same_dir"=nrow(same_dir), "rep_sig"=nrow(rep_sig)))
}



res <- lapply(all_studies, function(x) get_pair_overlap(x$gene_smok,
                                                      x$study, 
                                                      gene_smok_ae,
                                                      "AE"))
ae_overlap <- data.table::rbindlist(res)
ae_overlap_df <- ae_overlap %>%
  dplyr::select(-study2) %>%
  mutate(across(-study1, ~as.numeric(as.character(.)))) %>%
  mutate(study=unlist(study1)) %>%
  mutate(across(contains("cor"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, platform)) %>%
  dplyr::select(study, tissue, platform, everything()) %>%
  arrange(tissue, study)
ae_overlap_df %>%
  #mutate(same_dir=round(same_dir/tot, digits=3)) %>%
  #dplyr::select(-cor.p) %>%
  write_csv("data/supp_tables/ae_overlap_results_smok.csv")

# make a correlation plot

all_pairs <- combn(names(all_studies),2)
all_pair_overlap <- apply(all_pairs, 2, function(x) {
  study1 <- x[1]
  study2 <- x[2]
  print(sprintf("%s,%s", study1, study2))
  get_pair_overlap(all_studies[[study1]]$gene_smok, study1, 
                   all_studies[[study2]]$gene_smok, study2)
})

all_pair_df <- data.table::rbindlist(all_pair_overlap)
full_pair <- ae_overlap %>% 
  dplyr::select(-study2) %>%
  dplyr::rename(study2=study1)  %>%
  mutate(study1="Grouped AE") %>% 
  dplyr::select(study1, everything()) %>%
  bind_rows(all_pair_df)
# TODO:
# - add AE
# - reorder by tissue

self_cor <- data.table::rbindlist(lapply(c(names(all_studies),"Grouped AE"), function(x)
  list("study1"=x, "study2"=x, "cor"=1)))


full_pair2 <- full_pair %>%
  bind_rows(full_pair %>% 
              mutate(study3=study1, study1=study2, study2=study3) %>%
              dplyr::select(-study3) %>%
              dplyr::select(study1, study2, everything())) %>%
  mutate(cor=ifelse(tot < 30, NA, cor )) %>%
  bind_rows(self_cor) 

counts2 <- full_pair2 %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) 
  
counts2$tissue =factor(counts2$tissue, levels=unique(list_studies %>% arrange(tissue) %>% pull(tissue)))
counts2 <- counts2 %>% arrange(tissue,study)
counts2$study1=factor(counts2$study1, levels=unique(counts2$study1))
counts2$study2=factor(counts2$study2, levels=unique(counts2$study1))

library('RColorBrewer')
set3 <- brewer.pal(12, "Paired")
ggplot(counts2)+
  geom_bar(mapping = aes(x = study1, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=set3)+
  theme_void()
 ggplot(counts2,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=rep_sig)) +
  geom_text(aes(label=rep_sig), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,35), 
                       space = "Lab", 
                       name="Number of replicated genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/gene_overlap.png")



write_csv(all_pair_df, "data/results/all_pair_cors.csv")

#  filter(study1 < study2)
 # Create a ggheatmap

long_reord2 <- full_pair2 %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", 1, tissue)) %>%
  arrange(tissue) 
long_reord2$study1=factor(long_reord2$study1, levels=unique(long_reord2$study1))
long_reord2$study2=factor(long_reord2$study2, levels=unique(long_reord2$study1))
ggplot(long_reord2, aes(study1, study2, fill = cor))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/correlation_plot.png")
# add airway epithelium

# plot fraction same direction? or n????

# meta-analyze?

prep_meta_study <- function(ds, study){
  ds %>%
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(SD=(logFC-logFC.l)/1.96) %>%
    mutate(ID=study) %>%
    dplyr::select(gene, logFC, SD, ID)
}

meta_study <- function(df){
  ma <- meta::metagen(df %>% pull(logFC), # treatment estimate
                            df %>% pull(SD), # standard error
                            studlab=df %>% pull(ID),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "DL", # method for between study variance
                            hakn = FALSE,
                            prediction = FALSE,
                            sm = "MD") 
  return(list("gene"=unique(df$gene),
              "logFC"=ma$TE.fixed,
              "logFC.l"=ma$lower.fixed,
              "logFC.u"=ma$upper.fixed,
              "p"=ma$pval.fixed,
              "n"=unique(df$n)))
}
meta_in <- data.table::rbindlist(lapply(all_studies, function(x) 
  prep_meta_study(x$gene_smok, x$study)))
meta_in2 <- meta_in %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup()
meta_in3 <- meta_in2 %>% filter(n >=10) %>% group_split(gene)
# 20183

ma_vals <- lapply(meta_in3, meta_study)

# 5min
mult_study <- data.frame(apply(data.table::rbindlist(ma_vals) ,
                                 c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study$adj.p <- p.adjust(mult_study$p, method="fdr")
save(mult_study, file="data/results/meta_across_sm.RData")

mult_study %>% filter(adj.p < 0.05) %>% nrow() # 787
mult_study %>% filter(adj.p < 0.05 & abs(logFC)>=0.3) # 2

mult_study %>% filter(adj.p < 0.05 & abs(logFC)>=0.2) 

# forestPlot
plot_gene <- function(my_meta, my_gene){
  df <- my_meta %>% filter(gene==my_gene)
  ma <- meta::metagen(df %>% pull(logFC), # treatment estimate
                      df %>% pull(SD), # standard error
                      studlab=df %>% pull(ID),
                      comb.fixed = FALSE,
                      comb.random = TRUE,
                      method.tau = "DL", # method for between study variance
                      hakn = FALSE,
                      prediction = FALSE,
                      sm = "MD") 
  meta::forest.meta(ma)
}
plot_gene(meta_in2, "AHRR")

# ....... #
# [x] 0. fix IDs
# 1. results for sex-related studies
sex_r_s <- setdiff(list_studies %>% 
                     filter(`follow up analysis by sex`=="Y") %>% 
                     pull(study),"GSE55962")
res_int <- lapply(all_studies[sex_r_s], function(x) get_pair_overlap(x$gene_int,
                                                        x$study, 
                                                        gene_int_ae,
                                                        "AE"))
int_overlap <- data.frame(do.call(rbind, res_int))
int_overlap_df <- int_overlap %>%
  dplyr::select(-study2) %>%
  mutate(across(-study1, ~as.numeric(as.character(.)))) %>%
  mutate(study=unlist(study1)) %>%
  mutate(across(contains("cor"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, platform)) %>%
  dplyr::select(study, tissue, platform, everything()) %>%
  arrange(tissue, study)
int_overlap_df %>%
  mutate(same_dir=round(same_dir/tot, digits=3)) %>%
  dplyr::select(-cor.p, -study1) %>%
  write_csv("data/supp_tables/ae_int_overlap_results_smok.csv")

# [x] pairwise correlation?
# - is correlation associated with platform
# 3. divide into disc/valid for meta --> re-run
# 4. meta sex-related

meta_int <- data.table::rbindlist(lapply(all_studies[sex_r_s], function(x) 
  prep_meta_study(x$gene_int, x$study)))
meta_int2 <- meta_int %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup()
table(meta_int2$n)
meta_int3 <- meta_int2 %>% filter(n >=4) %>% group_split(gene)
# 20183

# 5min
ma_vals_int <- lapply(meta_int3, meta_study)

mult_study_int <- data.frame(apply(data.table::rbindlist(ma_vals_int) ,
                               c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study_int$adj.p <- p.adjust(mult_study_int$p, method="fdr")
plot_gene(meta_int2, "SNX3")


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

save(mult_study, file="data/results/meta_across_sm.RData")

# AE/BE
ae_studies <-list_studies %>% 
  filter(str_detect(tissue,"airway") | str_detect(tissue, "trachea")) %>%
  pull(study)
ae_meta <- meta_in %>% filter(ID %in% ae_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n==4) %>%
  group_split(gene)
ae_meta_out <- get_meta_out(ae_meta)
ae_sig <- ae_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)

# blood or blood component
blood_studies <-list_studies %>% 
  filter(str_detect(tissue, "blood"), !str_detect(tissue, "b cell") ) %>%
  pull(study) # PBMC, whole blood
blood_meta <- meta_in %>% filter(ID %in% blood_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n==5) %>%
  group_split(gene)
blood_meta_out <- get_meta_out(blood_meta)
blood_sig <- blood_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)



intersect(ae_sig$gene, blood_sig$gene) # AKR1C3
plot_gene( meta_in %>% filter(ID %in% ae_studies) %>%
             group_by(gene) %>% 
             mutate(n=n())  %>% 
             ungroup() %>%
             filter(n==4), "AKR1C3")

plot_gene( meta_in %>% filter(ID %in% blood_studies) %>%
             group_by(gene) %>% 
             mutate(n=n())  %>% 
             ungroup() %>%
             filter(n==5), "LRRN3")


plot_gene( meta_in %>% filter(ID %in% ae_studies) %>%
             group_by(gene) %>% 
             mutate(n=n())  %>% 
             ungroup() %>%
             filter(n==4), "AKR1C3") # up in smoking in oral cancer cells


# buccal mucosa
bc_studies <- list_studies %>% filter(tissue=="buccal mucosa") %>% pull(study)
bc_meta <- meta_in %>% 
  filter(ID %in% bc_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n==2) %>%
  group_split(gene)
bc_meta_out <- get_meta_out(bc_meta)
bc_sig <- bc_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)

# alveolar macrophages
am_studies <- list_studies %>% filter(tissue=="alveolar macrophages") %>% pull(study)
am_meta <- meta_in %>% filter(ID %in% am_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n==2) %>%
  group_split(gene)
am_meta_out <- get_meta_out(am_meta)
am_sig <- am_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)



# lung
lung_studies <- list_studies %>% filter(tissue=="lung") %>% pull(study)
lung_meta <- meta_in %>% filter(ID %in% lung_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n==3) %>%
  group_split(gene)
lung_meta_out <- get_meta_out(lung_meta)
lung_sig <- lung_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
ae_sig %>% filter(gene %in% lung_sig$gene)
blood_sig %>% filter(gene %in% lung_sig$gene)
plot_gene( meta_in %>% filter(ID %in% lung_studies) %>%
             group_by(gene) %>% 
             mutate(n=n())  %>% 
             ungroup(), "MS4A6A") # these two don't look good
save(ae_meta_out, blood_meta_out, lung_meta_out, file="data/meta_tiss_grp_out.RData")

# overlap btw studies where cor > 1


# heatmap for gene data  #

blood_plot_dat <- meta_in %>% filter(gene %in% blood_sig$gene, 
                   ID %in% blood_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(blood_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 
blood_plot_dat$gene <- factor(blood_plot_dat$gene, levels=blood_sig %>% arrange(logFC) %>% pull(gene))
blood_plot_dat$study <- factor(blood_plot_dat$study, 
                               levels=c("GSE20189", "GSE56768", "GSE21862", "GSE42057",
                                        "GSE87072", "pooled"))  
ggplot(blood_plot_dat, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-1.5,1.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/blood_meta_heatmap.png")


ae_plot_dat <- meta_in %>% filter(gene %in% ae_sig$gene, 
                                     ID %in% ae_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(ae_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 
ae_plot_dat$gene <- factor(ae_plot_dat$gene, levels=ae_sig %>% arrange(logFC) %>% pull(gene))

ggplot(ae_plot_dat, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-1.5,1.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/ae_meta_heatmap.png")



am_plot_dat <- meta_in %>% filter(gene %in% am_sig$gene, 
                                  ID %in% am_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(am_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 
am_plot_dat$gene <- factor(am_plot_dat$gene, 
                           levels=am_sig %>% arrange(logFC) %>% pull(gene))

ggplot(am_plot_dat, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-1.5,1.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/am_meta_heatmap.png")

