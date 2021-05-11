# WHICH ARE THE BLOOD STUDIES?
# Parker
# Armili
# Su
# Zeller
# Na
# Vink
# Morrow
# Georgiadis - Northern Sweden Health + Disease, Epic
# 


# 0. identify a blood signature 

# Huan et al?

# Beneike, how tested

# Martin
gse <- getGEOData("GSE20680")

pDat <- gse$originalData$GSE20680$pheno
eDat <- gse$originalData$GSE20680$expr
gse$originalData$GSE20680$platform

pDat %>% 
  select(geo_accession, characteristics_ch1) %>%
  filter(str_detect(characteristics_ch1, "Control")) %>%
  nrow() #52 samples but no smoking data 


gse2 <- getGEOData("GSE20681")

pDat2 <- gse2$originalData$GSE20681$pheno
eDat2 <- gse$originalData$GSE20681$expr
gse2$originalData$GSE20681$platform


pDat2.2 <- pDat2 %>% 
  select(geo_accession, source_name_ch1, contains("characteristics_ch1")) %>%
  pivot_longer(c(-geo_accession,-source_name_ch1)) %>%
  separate(value, into=c("column", "value"), sep=": ") %>%
  mutate(column=tolower(column)) %>%
  mutate(column=ifelse(str_detect(column, "age"), "age", column)) %>%
  filter(column != "") %>%
  select(-name) %>%
  pivot_wider(names_from=column, values_from=value)

pDat2.2 %>% 
  filter(`disease state`=="Control (0)", `smoking status` %in% c('Current', 'Never')) %>%
  group_by(`sex (gender)`, `smoking status`) %>%
  count() # mostly male

# 1. test w/ small datasets
# 2. do we see sex differences in these genes?
# test Chatz dataset (which dataset was it found on)


# 3. apply to 

