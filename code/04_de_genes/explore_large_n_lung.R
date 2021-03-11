# -- ---- lung data ----- #
# missing smoking informating :/ 
# GSE23546  - super series
# GSE23352 , GSE23545, GSE23529
lung_studies <- getGEO("GSE23546")
pDat_lung <- pData(lung_studies[[1]])
pDat_lung %>% fct_summ()
lung1 <- getGEO("GSE23352")
lung2 <- getGEO("GSE23529")
lung3 <- getGEO("GSE23545")
# what are the platforms?

samples1 <- pData(lung1[[1]])$geo_accession
samples2 <- pData(lung2[[1]])$geo_accession
samples3 <- pData(lung3[[1]])$geo_accession

pDat_lung2 <- pDat_lung %>% 
  mutate(ds=case_when(
    geo_accession %in% samples1 ~ "ds1",
    geo_accession %in% samples2 ~ "ds2",
    geo_accession %in% samples3 ~ "ds3"
  ) )
pDat_lung1 %>%  dplyr::select(characteristics_ch1, characteristics_ch1.1)

exprs(lung1)