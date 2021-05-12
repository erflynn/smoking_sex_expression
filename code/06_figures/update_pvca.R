# updates:
# - separate into X + Y and autosomal genes
# - bootstrap CI for estimate
# - permutation model for interaction
source("code/00_utils.R")
load("data/ae_full_exp.RData")  
# pDat5.1, expDat5  


load("ref/probe_gene.RData")
load("ref/chr_map.RData")
eDat_r <- tibble(probes=rownames(expDat5)) %>%
  left_join(probe_gene, by=c("probes")) %>%
  left_join(chr_map2, by=c("gene"="hgnc_symbol"))
num_no_gene <- eDat_r %>% filter(is.na(gene)) %>% nrow()
num_no_chr <- eDat_r %>% 
  filter(!is.na(gene) & is.na(chromosome_name)) %>% nrow()
print(sprintf("%s did not map to gene and %s did not map to chromosome out of %s.",
              num_no_gene, num_no_chr, nrow(eDat_r)))
auto_probes <-eDat_r %>% filter(!chromosome_name %in% c("X", "Y")) %>% pull(probes)
xy_probes <- eDat_r %>% filter(chromosome_name %in% c("X", "Y")) %>% pull(probes)

set.seed(1014)

calc_var1 <- function(eDat, pDat){
  rs <- rowSums(is.na(eDat))
  eDat_complete <- eDat[rs==0,]
  print(sprintf("Removed %s rows out of %s because of NAs", 
                nrow(eDat)-nrow(eDat_complete),
                nrow(eDat)))
  pcs_l <- varPartPC(eDat_complete)
  expr_pcs_t <- t(pcs_l$pcs)
  props <- pcs_l$props
  my_model <- as.formula("~ sex + smoking + sex:smoking")
  
  (est_var <- get_pvca(expr_pcs_t, pDat, my_model, props))
  return(est_var)
}

bootstrap_ci <- function()


auto_var <- calc_var1(expDat5[auto_probes,], pDat5.1)
#xy_var <- calc_var1(expDat5[xy_probes,], pDat5.1)
all_var <- calc_var1(expDat5, pDat5.1)
all_var_xy = all_var - auto_var
# the better way to do this would be to go directly


save(est_var,
     file=sprintf("data/results/%s_pvca.RData", ds))

rand_var <- rand_pvca(100, expr_pcs_t, pDat, my_model, props)
save(rand_var, 
     file=sprintf("data/results/%s_rand.RData", ds))