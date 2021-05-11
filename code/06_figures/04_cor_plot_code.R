
library(tidyverse)
library(psych) # for cor.wt

load("data/fig_color_scale.RData")

tiss_levels =c("airway epithelium", "trachea epithelium",
   "nasal epithelium", "oral cavity",
   "buccal mucosa", "sputum", "alveolar macrophages", 
   "lung", 
   "blood - b cells", "blood - pbmcs", "blood - whole",
   "liver", "kidney", 
   "brain - prefrontal cortex")
tiss_to_col <- c(paired3, "black")
names(tiss_to_col) <- c(tiss_levels, "all")


xy_genes <-  c(xy_chr_map$hgnc_symbol, "TTTY15", "KAL1", "ARSE", "ARSL","LOC389906")


get_cor <- function(df1, study1, df2, study2, cutoff, auto=FALSE){
  both_df <- df1 %>% mutate(SD=(logFC-logFC.l)/1.96)  %>%
    inner_join(df2 %>% mutate(SD=(logFC-logFC.l)/1.96), by="gene") %>% 
    dplyr::select(gene, logFC.x, logFC.y, SD.x, SD.y, adj.p.x, adj.p.y) %>%
    filter(!is.na(logFC.x), !is.na(logFC.y))
  both_df2 <- both_df %>% filter(adj.p.x < cutoff | adj.p.y < cutoff)
  if (auto){
    both_df2 <- both_df2 %>% filter(!gene %in% xy_genes)
  }
  
  if (nrow(both_df2) > 10){
    vals <- both_df2 %>% ungroup() %>% dplyr::select(logFC.x, logFC.y)
    wts <- both_df2 %>% ungroup() %>% dplyr::select(SD.x, SD.y)
    cor_wt <- cor.wt(vals, w=wts)
    my_cor <- cor_wt$r[1,2]
  } else {
    my_cor <- NA
  }

  return(list=c("study1"=study1, "study2"=study2, 
                "tot"=nrow(both_df), "filt"=nrow(both_df2), "cor"=my_cor))
}


create_cor_mat <- function(my_dfs, cutoff, auto=FALSE){
  all_pairs <- combn(names(my_dfs),2)
  all_cors <- apply(all_pairs, 2, function(x) {
    study1 <- x[1]
    study2 <- x[2]
    print(sprintf("%s,%s", study1, study2));
    get_cor(my_dfs[[study1]], study1, 
                     my_dfs[[study2]], study2, cutoff, auto)
  })
  all_pair_dat<- data.frame(t(all_cors)) %>%
    mutate(across(c(cor, tot), ~as.numeric(as.character(.))))
  
  self_cor <- data.table::rbindlist(lapply(names(my_dfs), 
                          function(x) list("study1"=x, "study2"=x, "cor"=1)))
  
  full_pair <- all_pair_dat %>%
    bind_rows(all_pair_dat %>% 
                mutate(study3=study1, study1=study2, study2=study3) %>%
                dplyr::select(-study3) %>%
                dplyr::select(study1, study2, everything())) %>%
    bind_rows(self_cor) 
  return(full_pair)
}



reform_for_plot <- function(full_pair, list_studies){
  list_studies2 <- list_studies %>% 
    dplyr::select(study, tissue) %>%
    mutate(across(c(tissue, study), as.character)) %>%
    bind_rows(tibble("study"="Grouped AE", "tissue"="airway epithelium"))
  
  ds_w_tiss <- full_pair %>% 
    left_join(list_studies2, 
              by=c("study1"="study"))  
  long_reord2 <- ds_w_tiss %>%
    mutate(tissue=factor(tissue, levels=intersect(tiss_levels,
                                                  unique(ds_w_tiss$tissue)))) %>%
    arrange(tissue, study1)
  long_reord2$study1=factor(long_reord2$study1, levels=unique(long_reord2$study1))
  long_reord2$study2=factor(long_reord2$study2, levels=unique(long_reord2$study1))
  return(long_reord2)
}

cor_plot <- function(long_reord2){
  ggplot(long_reord2, aes(study1, study2, fill = cor))+
    geom_tile()+
    scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                         midpoint = 0,  space = "Lab", 
                         name="Pearson\nCorrelation",
                         na.value="#DCDCDC") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust=1))+
    coord_fixed()+
    ylab("")+
    xlab("")
}




# // TODO fix this so they always match
color_bar <- function(ds){
  if ("study1" %in% colnames(ds)){
    ds <- ds %>% dplyr::rename(study=study1)
  }
  ggplot(ds)+
    geom_bar(mapping = aes(x = study, y = 1, fill = tissue), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values=
                        tiss_to_col[ds %>% 
                                      distinct(tissue, study) %>% 
                                      mutate(tissue=as.character(tissue)) %>%
                                      pull(tissue)])
}

# --- SMOKING --- #
all_dfs <- lapply(all_studies2 , function(x) x$gene_smok)
names(all_dfs) <- names(all_studies2)
all_dfs2 <- c(all_dfs[names(all_dfs)!="GSE44456"], 
              "Grouped AE"=list(gene_smok_ae))
cm <- create_cor_mat(all_dfs2, 0.1)
cm2 <- reform_for_plot(cm, list_studies)
cor_plot(cm2)
ggsave("figures/paper_figs/wcor_plot_smoking.pdf")
color_bar(cm2)
ggsave("figures/paper_figs/wcor_plot_smoking_colorbar.png")

# --- SEX --- #
all_dfs_s <- lapply(all_studies2 , function(x) x$gene_sex)
names(all_dfs_s) <- names(all_studies2)
missing <- names(which(is.na(all_dfs2_s)))
all_dfs2_s <- c(all_dfs_s[!names(all_dfs_s) %in% c("GSE44456", missing)], 
              "Grouped AE"=list(gene_sex_ae))
cm_sex<- create_cor_mat(all_dfs2_s, 0.1)
cm2_sex <- reform_for_plot(cm_sex, list_studies)
cor_plot(cm2_sex)
ggsave("figures/paper_figs/wcor_plot_sex.pdf")
color_bar(cm2_sex)
ggsave("figures/paper_figs/wcor_plot_sex_colorbar.png")

save(cm2, cm2_sex, file="data/data_cor_mat.RData")

# --- AUTOSOMAL SEX --- #

cm_sex_a<- create_cor_mat(all_dfs2_s, 0.1, auto=T)
cm2_sex_a <- reform_for_plot(cm_sex_a, list_studies)
cor_plot(cm2_sex_a)
ggsave("figures/paper_figs/wcor_plot_sex_auto.pdf")

color_bar(cm2_sex_a)
ggsave("figures/paper_figs/wcor_plot_sex_auto_colorbar.png")

save(cm2, cm2_sex, cm2_sex_a, file="data/data_cor_mat.RData")

# --- INTERACTION EFFECTS
all_dfs_int <- lapply(all_studies2 , function(x) x$gene_int)
all_dfs_int2 <- c(all_dfs_int[!names(all_dfs_int) %in%
                                c("GSE4302", "GSE19027", "GSE8987",
                                  "GSE42743", "GSE16149",
                                     "GSE2125", "GSE103174",
                                  "GSE56768", "GSE21862", "GSE87072",
                                  "GSE18723",
                                  "GSE44456")], "Grouped AE"=list(gene_int_ae))
cm_int <- create_cor_mat(all_dfs_int2, 0.1)
cm2_int <- reform_for_plot(cm_int, list_studies)
cor_plot(cm2_int)


cm_int2 <- create_cor_mat(all_dfs_int2, 1)
cm2_int2 <- reform_for_plot(cm_int2, list_studies)
cor_plot(cm2_int2)

