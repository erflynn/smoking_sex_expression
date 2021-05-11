
# TODO:
# - add in control probes?
# - try vst


# Charlesworth pre-processing:
# 0. chi-square test to identify significantly expressed transcripts (20k)
# 1. the score within individuals using decile percentage bins grouped by average log transform signal across individuals
# 2. linear regression against average log transform signal and it's squared value
# 3. normalized residuals with an inverse gaussian across individuals
#  eventually tested allowing for residual genetic effects because of family relationships
#  used an LRT test to see if smoking was associated with gene expression levels
#  all FDR cut off 0.05



# -- use LUMI -- # 
library('lumi')

load("data/blood/processed/proc1.RData") 
full_L  <- myL
for (i in 2:13){
  print(i)
  load(sprintf("data/blood/processed/proc%s.RData", i) )
  full_L <- combine(full_L, myL)
}

full_L.T <- lumiT(full_L, method="log2")
full_L.NT <- lumiN(full_L.T, method="quantile")

save(full_L.NT, file="data/blood/blood_full.RData")

# --- add pheno data ---- #
phe_dat <- read_csv("data/etabm305_phe_processed.csv")
phe_dat2 <- phe_dat %>% 
  dplyr::select(-ftp) %>% 
  dplyr::rename(sampleID=individual) %>%
  mutate(hdl=as.numeric(hdl)) # 27 have missing HDL
empty_pDat <- pData(full_L.NT)
pDat2 <- empty_pDat %>% 
  inner_join(phe_dat2, by="sampleID")


# remove samples
filt_ds <- full_L.NT[,pDat2$sampleID]
# TODO - double check this is the RIGHT order after filt

pData(filt_ds) <- pDat2
save(filt_ds, file="data/blood/blood_proc.RData")