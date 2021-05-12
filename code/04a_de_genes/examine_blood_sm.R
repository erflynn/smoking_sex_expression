# smaller?
# - GSE20686 (n=209, Agilent), Beineke et al., whole blood
# - E-TABM-5278 (Affy), Martin et al., whole blood
# - Huan meta-analysis GSE33828, GSE48152, E-MTAB-1708, GSE48348, GSE36382
# GSE42057 PBMC

# Most missing labels
# - which of these are in?


#This code we attempt to examine the additional blood datasets to see if any 
#have sufficient number or easily downloadable to add to our analysis
# steps:
#-Look at phenotype categories, Identify samples
#-Run de analysis

library('tidyverse')
library('GEOquery')
gse1 <- getGEO("GSE20686")
gse2 <- getGEO("GSE33828")
gse3 <- getGEO("GSE48152")
gse4 <- getGEO("GSE48348")
gse5 <- getGEO("GSE36382")

pDat1 = pData(gse1$GSE20686_series_matrix.txt.gz)

head(pDat1)
nrow(pDat1)


# Beineke - 41 S, 100 NS
pDat1 %>%
  fct_summ()
#filter out disease state
#Look at pair information - 110 case controlled by age and sex
#age is in dob
#look at source name
pDat1 %>%
  rename(smoking=`smoking status:ch1`, 
         sex=`sex (gender):ch1`,
       disease=`disease state:ch1`,
       age_date=`date of birth:ch1`,
       pair=`pair:ch1`,
       timing=growth_protocol_ch1,
       date=last_update_date) %>%
  filter(str_detect(disease,'Control')) %>%
  group_by(smoking,sex) %>%
  count()
"E-TABM-5278"
"E-MTAB-1708"

pDat1 %>%
  rename(smoking=`smoking status:ch1`, 
         sex=`sex (gender):ch1`,
         disease=`disease state:ch1`,
         age_date=`date of birth:ch1`,
         pair=`pair:ch1`) %>%
  filter(is.na(smoking)) %>%
  select(title,source_name_ch1,sex,disease) %>%
  fct_summ()
#It appears that the blood was collected either 
#prior or during catheterization

pDat2=pData(gse2$GSE33828_series_matrix.txt.gz)
pDat2 %>%
  fct_summ()
#age,gender
pDat2 %>%
  select(title,source_name_ch1) %>%
  View()
#No smoking information provided
#Maas, rotterdam study https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5662692/

pDat3=pData(gse3$GSE48152_series_matrix.txt.gz)
pDat3 %>%
  fct_summ()
# No metadata, InChianti cohort.

pDat4=pData(gse4$GSE48348_series_matrix.txt.gz)
pDat4 %>%
  fct_summ()
# No metadata, Estonian Biobank

pDat5=pData(gse5$GSE36382_series_matrix.txt.gz)
pDat5%>%
  fct_summ()
# No metadata, SHIP-TREND

# E-TABM-5278 -does not appear to exist

# E-MTAB-1708 
# These data exist but have no metadata labels related to smoking


# GSE42057 PBMC
gse6=getGEO("GSE42057")
