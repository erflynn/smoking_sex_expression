
library(tidyverse)
library(limma)
library(meta)
source("code/00_utils.R")

summarize_max_p <- function(ds){
  ds %>% 
    group_by(gene) %>%
    mutate(n=n()) %>%
    slice_min(P.Value) %>%
    arrange(P.Value) %>% 
    mutate(adj.P.Val=p.adjust(P.Value, "fdr"))
}
summarize_max_expr <- function(ds){
  # find the probe that has the max expression
  ds %>% 
    group_by(gene) %>%
    mutate(n=n()) %>%
    slice_max(AveExpr) %>%
    arrange(P.Value) %>% 
    mutate(adj.P.Val=p.adjust(P.Value, "fdr"))
  # return it
}
summarize_average <- function(my_expr, design){
  probe_df <- add_annot_mi(my_expr, "GPL570")
  expr_summ <- probe_df %>% 
    select(-probes) %>%
    pivot_longer( -gene, names_to="sample", values_to="expr") %>%
    group_by(gene, sample) %>%
    mutate(n=n()) %>%
    summarize(expr=mean(expr)) %>%
    pivot_wider(names_from=sample, values_from="expr")
  expr_summ2 <- expr_summ %>% ungroup() %>% 
    filter(!is.na(gene)) 
  expr_summ_mat <- as.matrix(expr_summ2 %>% 
                               select(-gene))
  rownames(expr_summ_mat) <- expr_summ2$gene
  fit_a <- lmFit(expr_summ_mat, design)
  fit_a <- eBayes(fit_a)
  df_a <- data.frame(topTable(fit_a, coef="smokingS",  
                                  number=nrow(fit_a), confint=T)) 
  df_a$gene <- rownames(df_a)
  return(df_a)
}

prep_data_ma <- function(df){
  df2 <- df %>%
    mutate(SD=(logFC-CI.L)/1.96) %>% # SE for effect size - we have 95%  CI = 1.96*SD
    dplyr::select(ID, geneSymbol, chromosome, logFC, SD, P.Value) %>%
    filter(geneSymbol!="") %>%
    group_by(chromosome,geneSymbol) %>%
    mutate(n=n()) 
  return(df2)
}

ma_probes_genes <- function(df){
  # generic inverse variance methods
  ma <- metagen(df %>% pull(logFC), # treatment estimate
                df %>% pull(SD), # standard error
                studlab=df %>% pull(ID),
                comb.fixed = TRUE,
                comb.random = FALSE,
                method.tau = "DL", # method for between study variance
                hakn = FALSE,
                prediction = FALSE,
                sm = "MD") # underlying summary measure (RR, OR, ASD, ROM, HR, MD)
  return(list("chromosome"=unique(unlist(df$chromosome)),
              "gene"=unique(df$geneSymbol),
              "logFC"=ma$TE.fixed,
              "logFC.l"=ma$lower.fixed,
              "logFC.u"=ma$upper.fixed,
              "p"=ma$pval.fixed,
              "n"=unique(df$n)))
  
}



get_top_genes <- function(eDat, pDat, frac){
  # select xx fraction
  sel_idx <- sample(nrow(pDat), floor(frac*nrow(pDat)))
  pDat2 <- pDat[sel_idx,]
  eDat2 <- eDat[,pDat2$geo_accession]

  
  
  # get probe results
  design2 <- model.matrix(~ sex + smoking + age + sex*smoking + 
                           race_ethnicity + submission_date, 
                         data=pDat2) 
  
  #df_a <- summarize_average(eDat2, design2)
  
  fit <- lmFit(eDat2, design2)
  fit <- eBayes(fit)
  df_smok1 <- data.frame(topTable(fit, coef="smokingS",  
                                  number=nrow(fit), confint=T)) 
  
  # add gene annotations
  probe_res <- add_annot_mi(df_smok1, "GPL570") %>%
    filter(!is.na(gene))
  
  smok_ma1 <- run_clean_ma(probe_res %>% 
                             rename(ID=probes, geneSymbol=gene) %>% 
                             mutate(chromosome=""))
  smok_ma <- smok_ma1 %>% 
    ungroup() %>% 
    select(logFC, gene, n, p, adj.p) %>%
    rename(P.Value=p, adj.P.Val=adj.p)
  # summarize to genes
  max_e <- summarize_max_expr(probe_res)
  min_p <- summarize_max_p(probe_res)

  # compare to "true"
  # - overlap in true lists of genes
  
  
  # output results
  calc_out <- function(df, tt_true, method1, method2){
    do.call(rbind, list(compare_out(df, tt_true, method1, method2, 0.05),
                        compare_out(df, tt_true, method1, method2, 0.01)))}
  compare_out <- function(df, tt_true, method1, method2, pcut){
    true_cut <- tt_true %>% filter(adj.P.Val < pcut)
    df_cut <- df %>% filter(adj.P.Val < pcut)
    
    # counts
    num_true = nrow(true_cut)
    num_df = nrow(df_cut)
    true_df_overlap = length(intersect(true_cut$gene, df_cut$gene))
    
    # calculate stats
    true_df = true_cut %>% select(gene, logFC, P.Value) %>% 
      inner_join(df_cut  %>% select(gene, logFC, P.Value), by="gene") 
    
    est_cor = cor.test(true_df$logFC.x, true_df$logFC.y)
    p_corS =cor.test(true_df$P.Value.x, true_df$P.Value.y, 
                         method="spearman")
    p_corP =cor.test(log10(true_df$P.Value.x), log10(true_df$P.Value.y))

    return(list(
      "frac_miss"= frac,
      "method1"=method1,
      "method2"=method2,
      "cutoff" = pcut,
      "num_true" = num_true,
      "num_df" = num_df,
      "true_imp_overlap" = true_df_overlap,
      "est_cor_imp"= est_cor$estimate[[1]],
      "est_cor_p" =est_cor$p.value,
      "p_corS"= p_corS$estimate[[1]],
      "p_corS_p" =p_corS$p.value,
      "p_corP"= p_corP$estimate[[1]],
      "p_corP_p" =p_corP$p.value

    ))
  }
  out1 <- calc_out(smok_ma, tt_true, "ma", "true")
  out2 <- calc_out(max_e, tt_true, "max_e", "true")
  out3 <- calc_out(min_p, tt_true, "min_p", "true")
  #out4 <- calc_out(df_a, tt_true, "average", "true")
  out5 <- calc_out(smok_ma, max_e, "ma", "max_e")
  out6 <- calc_out(smok_ma, min_p, "ma", "min_p")
  out7 <- calc_out(max_e, min_p, "max_e", "min_p")
  return(do.call(rbind, list(out1, out2, out3, out5, out6, out7)))
}


# setup
load("data/ae_full_exp.RData") # expDat5, pDat5.1
complete_pDat <- pDat5.1 %>%
  filter(!is.na(age), !is.na(race_ethnicity), !is.na(pack_years)) # 244
complete_eDat <- expDat5[,complete_pDat$geo_accession]
eDat <- complete_eDat
pDat <- complete_pDat

# create tt_true
design <- model.matrix(~ sex + smoking + age + sex*smoking + 
                          race_ethnicity + submission_date, 
                        data=pDat) 

fit <- lmFit(eDat, design)
fit <- eBayes(fit)
df_smok <- data.frame(topTable(fit, coef="smokingS",  
                                number=nrow(fit), confint=T)) 

# add gene annotations
tt_true <- add_annot_mi(df_smok, "GPL570") %>%
  filter(!is.na(gene))

start = Sys.time()
res <- get_top_genes(eDat,pDat, 0.5) # 3min
end = Sys.time()
print(end-start)

# res2 <- data.frame(apply(all_df, c(1,2), unlist))
# 
# res2 %>% filter(method2=="true") %>%
#   mutate(across(c(num_true, num_df, true_imp_overlap),
#                 ~ as.numeric(as.character(.)))) %>%
#   mutate(tpr=true_imp_overlap/num_true,
#          fpr=(num_df-true_imp_overlap)/num_df) %>%
#   ggplot(aes(x=tpr, y=fpr, shape=cutoff, col=method1)) +
#   geom_point()+
#   facet_grid(. ~frac_miss)
# 
# res2 %>% filter(method2=="true") %>%
#   mutate(across(c(est_cor_imp),
#                 ~ as.numeric(as.character(.)))) %>%
#   ggplot(aes(x=frac_miss, y=est_cor_imp, shape=cutoff, col=method1)) +
#   geom_point()+
#   facet_grid(.~cutoff)


ds_fraction <- c( 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
all_df <- do.call(rbind, lapply(ds_fraction, function(frac){
  do.call(rbind, lapply(1:10, function(i) get_top_genes(eDat, pDat, frac)))
}))
save(all_df, file="data/probe_gene_comp.RData")


# test spurious genes
# - how many replicate have multiple direction probes?
#  - vs. if we use the meta-analyzed 


#####
load("data/probe_gene_comp.RData")
all_df <- as_tibble(apply(all_df, c(1,2), unlist))
all_df2 <- all_df %>% mutate(across(c(-method1, -method2), as.numeric)) 
  

all_df2 %>% filter(method2=="true") %>%
         mutate(recall=true_imp_overlap/num_true,
                #fpr=(num_df-true_imp_overlap)/num_df,
                precision=true_imp_overlap/num_df) %>%
  mutate(method=case_when(
    method1=="ma" ~ "meta-analysis",
    method1=="max_e" ~ "max expression",
    method1=="min_p" ~ "min p-value"),
    cutoff=sprintf("FDR < %s", cutoff)) %>%
         ggplot(aes(x=recall, y=precision, col=method)) +
         geom_point(alpha=0.5)+
          #geom_label(aes(label=num_df))+
         facet_grid(frac_miss~cutoff)+
  xlim(0,1)+
  ylim(0.2,1)+
  theme_bw()
ggsave("figures/probe_summary_comparison.png")

