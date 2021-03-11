library(tidyverse)
library(MetaIntegrator)


lung_kept = read_csv("data/lung_sample_labels.csv")
head(lung_kept)
lung_kept %>%
  mutate(sex_lab=ifelse(is.na(sex_lab),'unlabeled',sex_lab)) %>%
  filter(!is.na(smok)  & smok != 'FS') %>%
  group_by(study_acc, smok, sex_lab) %>%
  count() %>%
  unite(grp,c(smok,sex_lab)) %>%
  pivot_wider(names_from=grp, values_from=n, values_fill=0) 
lung_kept1 <- lung_kept %>%
  mutate(sex_lab=ifelse(is.na(sex_lab),'unlabeled',sex_lab)) %>%
  filter(!is.na(smok)  & smok != 'FS') 
(lung_kept2 <- lung_kept1 %>%
  group_by(study_acc, smok) %>%
  count() %>%
  pivot_wider(names_from=smok, values_from=n, values_fill=0) %>%
  filter(NS > 5  & S > 5))
# six datasets

lung_kept2 %>% pull(study_acc)

# download with MetaIntegrator 
# what are the platforms?
# -- all affy so we might not be able to generalize to other non-affy but still a start!
# application ds: GPL10379 - Rosetta/Merck Human RSTA Custom Affymetrix 2.0 microarray [HuRSTA-2a520709]
# ... it is a bit of a problem we don't have this platform anywhere, but oh well
gse <- getGEOData("GSE103174") # GPL13667 -  [HG-U219] Affymetrix Human Genome U219 Array
gse2 <- getGEOData("GSE31210") # GPL570
gse3 <- getGEOData("GSE31908") # GPL96, 97, 570
gse4 <- getGEOData("GSE32539") # 6244, 8786(miR) -  	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version],
gse5 <- getGEOData("GSE37768") # GPL570
gse6 <- getGEOData("GSE40364") # this failed download - mb validation?, GPL570 + gene
# figure out how to download this


# add labels for smoking status
poss_samples <- lung_kept1 %>% semi_join(lung_kept2, by="study_acc")

add_smok_lab <- function(ds, my_acc) {
  my_samples = poss_samples %>% 
    filter(study_acc==my_acc) %>% 
    filter(sample_acc %in% ds$pheno$geo_accession) %>%
    arrange(sample_acc) 
  new_pheno = ds$pheno %>% filter(geo_accession %in% my_samples$sample_acc)
  new_expr = ds$expr[,my_samples$sample_acc]
  new_class = ifelse(my_samples$smok=="S", 1, 0)
  names(new_class) = my_samples$sample_acc
  ds$pheno = new_pheno
  ds$expr = new_expr
  ds$class = new_class
  stopifnot(checkDataObject(ds, "Dataset", "Pre-Analysis"))
  return(ds)
}


# there are three sub-datasets
gse3_gpl570 <- gse3$originalData$GSE31908_GPL570$pheno$geo_accession
gse3_gpl96 <- gse3$originalData$GSE31908_GPL96$pheno$geo_accession
gse3_gpl97 <- gse3$originalData$GSE31908_GPL97$pheno$geo_accession
poss_samples %>% 
  filter(sample_acc %in% gse3_gpl570) # all tumor

poss_samples %>% 
  filter(sample_acc %in% gse3_gpl96) # <-- let's just use this
poss_samples %>% 
  filter(sample_acc %in% gse3_gpl97)

gse1.1 <- add_smok_lab(gse$originalData$GSE103174, "GSE103174")
gse2.1 <- add_smok_lab(gse2$originalData$GSE31210, "GSE31210")
gse5.1 <- add_smok_lab(gse5$originalData$GSE37768, "GSE37768")

gse3.1_96 <- add_smok_lab(gse3$originalData$GSE31908_GPL96, "GSE31908")
gse3.1_97 <- add_smok_lab(gse3$originalData$GSE31908_GPL97, "GSE31908")

# TCGA lung data?

# divide into training and validation
# validation: GSE40364, GSE32539 (new platform + former)
gse4.1 <- add_smok_lab(gse4$originalData$GSE32539_GPL6244, "GSE32539")

# check expression
boxplot(gse1.1$expr)
boxplot(gse2.1$expr) # add to this one
min(gse2.1$expr)
gse2.1$expr <- gse2.1$expr-min(gse2.1$expr)+0.1
boxplot(gse2.1$expr) 

boxplot(gse5.1$expr)
boxplot(gse3.1_96$expr)
boxplot(gse3.1_97$expr)
boxplot(gse4.1$expr)
gse3.1_96$pheno$
gse5.1$class
# and create the metaObject
discovery_datasets <- list(gse1.1, gse2.1, gse5.1)
names(discovery_datasets) = c(gse1.1$formattedName, gse2.1$formattedName, 
                              gse5.1$formattedName)
metaObj=list() 
metaObj$originalData <- discovery_datasets
checkDataObject(metaObj, "Meta", "Pre-Analysis")

# run meta-integrator on the training data
metaObj <- runMetaAnalysis(metaObj)
metaObj <- filterGenes(metaObj, isLeaveOneOut = FALSE, FDRThresh = 0.05)
summarizeFilterResults(metaObj, getMostRecentFilter(metaObj))
str(metaObj, 2)
head(metaObj$metaAnalysis$pooledResults)


# - COVARIATES?
# - get last dataset