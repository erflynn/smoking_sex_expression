# compare final results to Yang et al
library(tidyverse)
library(readxl)
library(psych)
library(dmetar) # for grabbing SDs, download from github


# load our data
load("data/results/ae_smok.RData")
load("data/results/ae_int.RData")
probe_int_ae <- smok_ci_int_p
gene_int_ae <- comb_ma_dat %>% ungroup() %>%dplyr::select(-chromosome)
probe_smok_ae <- smok_ci_smok 
gene_smok_ae <- ma_smok %>%ungroup() %>%dplyr::select(-chromosome)
probe_sex_ae <- smok_ci_sex 
gene_sex_ae <- ma_sex %>% ungroup() %>% dplyr::select(-chromosome)


# load Yang et al data
fill_empty_cells <- function(df){
  # fill empty cells w previous values
  df2 <- df %>% 
    mutate(gene=case_when(
      is.na(gene) ~ lag(gene),
      TRUE ~ gene),
      chromosome=case_when(
        is.na(chromosome) ~ lag(chromosome),
        TRUE ~ chromosome
      )) 
  num_nas <- length(which(is.na(df2$gene)))
  if(num_nas==0){
    return(df2)
  }
  return(fill_empty_cells(df2))
}

format_sig_genes <- function(df){
  disc <- df[,1:6]
  rep <- df[,8:13] 
  colnames(disc) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
  colnames(rep) <- c("gene", "probe", "logFC", "pval", "FDR", "chromosome")
  disc1 <- fill_empty_cells(disc %>% filter(!is.na(pval)))
  rep1 <- fill_empty_cells(rep %>% filter(!is.na(pval)))
  return(list("disc"=disc1, "rep"=rep1))
}

disc_int <- read_csv("ref/yang_int_disc.csv")
rep_int <- read_csv("ref/yang_int_rep.csv")

supp_files2 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=2, skip=4, col_names=TRUE)
supp_files3 <- read_xlsx("ref/41598_2019_54051_MOESM2_ESM.xlsx", sheet=3, skip=4, col_names=TRUE)
#

smok_dr <- format_sig_genes(supp_files2)
sex_dr <- format_sig_genes(supp_files3)

disc_sex1 <- sex_dr$disc
disc_smok1 <- smok_dr$disc
rep_sex1 <- sex_dr$rep
rep_smok1 <- smok_dr$rep




setup_df <- function(rep_out, disc_out){
  rep_out %>% 
    select(probes, CI.L,logFC, P.Value, adj.P.Val) %>%
    rename(FDR=adj.P.Val, pval=P.Value) %>%
    inner_join(disc_out, by=c("probes"="probe")) %>%
    select(probes, gene, chromosome, everything()) %>%
    mutate(SD.x=(logFC.x - CI.L)/1.96) %>%
    mutate(SD.y=se.from.p(logFC.y, pval.y, 211, 
                          effect.size.type = 'difference',
                          calculate.g = FALSE)$StandardError) %>%
    as_tibble()
}

count_cor_rep <- function(probe_comb, fdr){
  num_r <- probe_comb %>% filter(FDR.x < fdr) %>% nrow()
  num_d <- probe_comb %>% filter(FDR.y < fdr) %>% nrow()
  num_int <- probe_comb %>% filter(FDR.x < fdr, FDR.y < fdr) %>% nrow()
  
  
  ct <- cor.test(probe_comb %>% 
                   filter(FDR.x < fdr | FDR.y < fdr) %>% 
                   pull(logFC.x),
                 probe_comb2 %>% 
                   filter(FDR.x < fdr | FDR.y < fdr) %>% 
                   pull(logFC.y))
  
  return(list("fdr"=fdr, "num_r"=num_r, "num_d"=num_d, "num_int"=num_int,
        "cor"=ct$estimate[["cor"]]))
}

calc_cor_wt <- function(probe_comb){
  vals <- probe_comb %>% ungroup() %>% dplyr::select(logFC.x, logFC.y)
  wts <- probe_comb %>% ungroup() %>% dplyr::select(SD.x, SD.y)
  cor_wt <- cor.wt(vals, w=wts)
  return(cor_wt$r[1,2])
}



# Smoking
probe_comb = setup_df(probe_smok_ae,disc_smok1)
smok_comparison <- do.call(rbind, 
        lapply(c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0001),
       function(fdr)
       count_cor_rep(probe_comb, fdr))) %>%
  apply(c(1,2), unlist) %>%
  as_tibble()
calc_cor_wt(probe_comb)

plot_rep <- function(probe_comb, fdr){
  ggplot(probe_comb %>% 
           mutate(`FDR < 0.05 in rep`=FDR.x < 0.05) %>%
           filter(FDR.y < fdr), aes(x=logFC.y, y=logFC.x, 
                                   col=`FDR < 0.05 in rep`))+
    geom_point(alpha=0.3)+
    theme_bw()+
    scale_color_manual(values=c("black", "red"))+
    xlab("logFC in Yang et al.")+
    ylab("logFC in expanded analysis")+
    ggtitle(sprintf("FDR < %s in Yang et al.", fdr))
}
plot_rep(probe_comb, 0.05)
plot_rep(probe_comb, 0.0001)


probe_comb %>% filter(FDR.y > 0.05 & FDR.x < 0.01)


# COMPARE SMOKING RESULTS
# FDR < 0.1
nrow(disc_smok1) # 8000
de_smok <- probe_smok_ae %>% filter(adj.P.Val < 0.1) %>% as_tibble() # 5681
length(intersect(de_smok$probes, disc_smok1$probe)) # 3112 probes overlap

# FDR < 0.05
disc_smok_filt <- disc_smok1 %>% filter(FDR < 0.05) # 5604
de_smok <- probe_smok_ae %>% filter(adj.P.Val < 0.05) %>% as_tibble() # 3610
length(intersect(de_smok$probes, disc_smok_filt$probe)) # 2018

