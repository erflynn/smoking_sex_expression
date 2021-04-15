# 03_covariate_tables.R
# 11/24/2020
#  
# Code for creating covariate tables and performing significance tests

require('tidyverse')
require('googlesheets4')
require('googledrive')
require('lubridate')


covar_table <- function(ds, my_var){
  ds2 <- ds %>%
    unite(grp, c("smoking", "sex"), sep=" - ") %>%
    group_by(grp, get(my_var)) %>%
    dplyr::count() %>%
    ungroup() %>%
    pivot_wider(names_from=grp, 
                values_from = n, 
                values_fill = list(n = 0)) 
  colnames(ds2)[1] <- my_var
  return(ds2)
}

#ae_only_meta <- read_csv("data/ae_only_meta.csv") %>% rename(smoking=smok, sex=expr_sex)
ae_only_meta <- read_csv("data/ae_pdat_full.csv") 

chisq.test(table(ae_only_meta$sex, ae_only_meta$smoking))

# --- race/ethnicity breakdown --- #
# counts per race/ethnicity group
(counts_per_re <- ae_only_meta %>%
    dplyr::select(smoking, sex, race_ethnicity) %>%
    #mutate(race_ethnicity=
    #         ifelse(str_detect(as.character(race_ethnicity), "black"),
    #                "black", as.character(race_ethnicity))) %>%
    covar_table("race_ethnicity"))

re_filt <- ae_only_meta %>%
  dplyr::select(smoking, sex, race_ethnicity) %>%
  mutate(race_ethnicity=ifelse(str_detect(as.character(race_ethnicity), "black"), "black", as.character(race_ethnicity)),
         expr_sex=as.character(sex)) %>%
  filter(!race_ethnicity %in% c("asian"), !is.na(race_ethnicity))


# race vs smok status
(chisq.race_smok <- chisq.test(table(re_filt$race_ethnicity, re_filt$smoking))) # p=0.090
# p=0.01
# race vs sex
(chisq.race_sex <- chisq.test(table(re_filt$race_ethnicity, re_filt$sex))) # p=0.835


# -- is the missingness signficiant?
missing_re <- ae_only_meta %>%
  dplyr::select(smoking, sex, race_ethnicity) %>%
  mutate(sex=as.character(sex)) %>%
  mutate(race_missing=ifelse(race_ethnicity=="(Missing)" | 
                               is.na(race_ethnicity), "y", "n"))
missing_re %>%
  covar_table("race_missing")

(chisq.race_m_smok <- chisq.test(table(missing_re$smoking, missing_re$race_missing))) # p = 0.217
# more smokers with missing data than expected

(chisq.race_m_sex <- chisq.test(table(missing_re$sex, missing_re$race_missing))) # p = 0.597



# --- look at age distributions --- #
ae_only_meta %>%
  dplyr::select(smoking, sex, age) %>%
  unite(grp,c("smoking", "sex"), sep=" - ") %>%
  ggplot(aes(x=factor(grp), y=age))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  geom_point(position=position_jitter(0.15), alpha=0.5)+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab("")

# get a list of means
ae_only_meta %>%
  dplyr::select(smoking, sex, age) %>%
  filter(!is.na(age)) %>%
  unite(grp,c("smoking", "sex"), sep=" - ") %>%
  group_by(grp) %>%
  summarize(mean_age=mean(age), median_age=median(age), sd_age=sd(age))

# perform a statistical test
aov_in <- ae_only_meta %>%
  dplyr::select(smoking, sex, age) %>%
  filter(!is.na(age))

two.way <- aov(age ~ smoking+sex+smoking*sex, data = aov_in)
summary(two.way) 
# sex: 0.01, smok: 0.01, sex*smok: 0.336



# -- missing age -- #
missing_age <- ae_only_meta %>%
  dplyr::select(smoking, sex, age) %>%
  mutate(sex=as.character(sex)) %>%
  mutate(age_missing=ifelse(is.na(age), "y", "n"))


missing_age %>% covar_table("age_missing")


(chisq.age_m_sex <- chisq.test(table(missing_age$sex, missing_age$age_missing))) 
# p=0.92

(chisq.age_m_smok <- chisq.test(table(missing_age$smoking, missing_age$age_missing))) 
# p=0.009


# -- pack_years -- #
smokers <- ae_only_meta %>%
  filter(smoking=="S")
missing_pkyrs <- ae_only_meta %>%
  dplyr::select(race_ethnicity, sex, pack_years) %>%
  mutate(miss_py=ifelse(is.na(pack_years), "y", "n"))

(chisq.py_m_sex <- chisq.test(table(missing_pkyrs$sex, missing_pkyrs$miss_py))) 
# 0.29

missing_pkyrs2 <- missing_pkyrs %>%
  mutate(race_ethnicity=ifelse(str_detect(as.character(race_ethnicity), "black"), "black", as.character(race_ethnicity)),
                         expr_sex=as.character(sex)) %>%
  filter(!race_ethnicity %in% c("asian"), !is.na(race_ethnicity))
(chisq.py_m_re <- chisq.test(table(missing_pkyrs2$race_ethnicity, missing_pkyrs2$miss_py))) 
# 0.009

pkyrs_dat <- smokers %>% 
  mutate(pack_years=as.numeric(pack_years)) %>% 
  filter(!is.na(pack_years)) %>%
  select(sex, pack_years)

t.test(pkyrs_dat %>% filter(sex=="female") %>% pull(pack_years),
       pkyrs_dat %>% filter(sex=="male") %>% pull(pack_years))

# create the table

age_py_data <- ae_only_meta %>%
  group_by(smoking, sex) %>%
  summarize(n_tot=n(),
         mean_age=mean(age, na.rm=T),
         age_sd=sd(age, na.rm=T), 
         missing_age=sum(is.na(age)),
         m_pack_years=mean(pack_years, na.rm=T),
         sd_pack_years=sd(pack_years, na.rm=T),
         missing_py=sum(is.na(pack_years))) %>% 
  unite(grp, c(smoking, sex)) %>%
  pivot_longer(-grp) %>%
  pivot_wider(names_from=grp, values_from=value, values_fill=0)


re_data <- ae_only_meta %>%
  group_by(smoking, sex, race_ethnicity) %>%
  summarize(n=n()) %>%
  unite(grp, c(smoking, sex)) %>%
  pivot_wider(names_from=grp, values_from=n, values_fill=0)

table_ae <- age_py_data %>% 
  bind_rows(re_data %>% 
              mutate(race_ethnicity=ifelse(is.na(race_ethnicity), 
                                           "missing_race", race_ethnicity)) %>%
              rename(name=race_ethnicity)) %>%
  mutate(across(everything(), ~ifelse(is.na(.) | is.nan(.), 0, .)))



age_py_data_s <- ae_only_meta %>%
  group_by(smoking) %>%
  summarize(n_tot=n(),
            mean_age=mean(age, na.rm=T),
            age_sd=sd(age, na.rm=T), 
            missing_age=sum(is.na(age)),
            m_pack_years=mean(pack_years, na.rm=T),
            sd_pack_years=sd(pack_years, na.rm=T),
            missing_py=sum(is.na(pack_years))) %>% 
  pivot_longer(-smoking) %>%
  pivot_wider(names_from=smoking, values_from=value, values_fill=0)


re_data_s <- ae_only_meta %>%
  group_by(smoking, race_ethnicity) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from=smoking, values_from=n, values_fill=0)

table_ae_s <- age_py_data_s %>% 
  bind_rows(re_data_s %>% 
              mutate(race_ethnicity=ifelse(is.na(race_ethnicity), 
                                           "missing_race", race_ethnicity)) %>%
              rename(name=race_ethnicity)) %>%
  mutate(across(everything(), ~ifelse(is.na(.) | is.nan(.), 0, .)))

table_ae2 <-table_ae_s %>% left_join(table_ae, by="name")
table_ae2 %>%
  select(name, S, "S_female", "S_male", "NS", "NS_female", "NS_male") %>%
  write_csv("data/supp_tables//ae_table_info.csv")
