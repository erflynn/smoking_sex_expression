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
    unite(grp, c("smok", "expr_sex"), sep=" - ") %>%
    group_by(grp, get(my_var)) %>%
    dplyr::count() %>%
    ungroup() %>%
    pivot_wider(names_from=grp, 
                values_from = n, 
                values_fill = list(n = 0)) 
  colnames(ds2)[1] <- my_var
  return(ds2)
}

ae_only_meta <- read_csv("data/ae_only_meta.csv")


# --- race/ethnicity breakdown --- #
# counts per race/ethnicity group
(counts_per_re <- ae_only_meta %>%
    dplyr::select(smok, expr_sex, race_ethnicity) %>%
    mutate(race_ethnicity=
             ifelse(str_detect(as.character(race_ethnicity), "black"),
                    "black", as.character(race_ethnicity))) %>%
    covar_table("race_ethnicity"))

re_filt <- ae_only_meta %>%
  dplyr::select(smok, expr_sex, race_ethnicity) %>%
  mutate(race_ethnicity=ifelse(str_detect(as.character(race_ethnicity), "black"), "black", as.character(race_ethnicity)),
         expr_sex=as.character(expr_sex)) %>%
  filter(!race_ethnicity %in% c("asian", "hispanic", "(Missing)"))


# race vs smok status
(chisq.race_smok <- chisq.test(table(re_filt$race_ethnicity, re_filt$smok))) # p=0.090

# race vs sex
(chisq.race_sex <- chisq.test(table(re_filt$race_ethnicity, re_filt$expr_sex))) # p=0.835


# -- is the missingness signficiant?
missing_re <- ae_only_meta %>%
  dplyr::select(smok, expr_sex, race_ethnicity) %>%
  mutate(expr_sex=as.character(expr_sex)) %>%
  mutate(race_missing=ifelse(race_ethnicity=="(Missing)" | 
                               is.na(race_ethnicity), "y", "n"))
missing_re %>%
  covar_table("race_missing")

(chisq.race_m_smok <- chisq.test(table(missing_re$smok, missing_re$race_missing))) # p = 0.217
# more smokers with missing data than expected

(chisq.race_m_sex <- chisq.test(table(missing_re$expr_sex, missing_re$race_missing))) # p = 0.597



# --- look at age distributions --- #
ae_only_meta %>%
  dplyr::select(smok, expr_sex, age) %>%
  unite(grp,c("smok", "expr_sex"), sep=" - ") %>%
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
  dplyr::select(smok, expr_sex, age) %>%
  filter(!is.na(age)) %>%
  unite(grp,c("smok", "expr_sex"), sep=" - ") %>%
  group_by(grp) %>%
  summarize(mean_age=mean(age), median_age=median(age), sd_age=sd(age))

# perform a statistical test
aov_in <- ae_only_meta %>%
  dplyr::select(smok, expr_sex, age) %>%
  filter(!is.na(age))

two.way <- aov(age ~ smok+expr_sex+smok*expr_sex, data = aov_in)
summary(two.way) 
# sex: 0.01, smok: 0.09, sex*smok: 0.67



# -- missing age -- #
missing_age <- ae_only_meta %>%
  dplyr::select(smok, expr_sex, age) %>%
  mutate(expr_sex=as.character(expr_sex)) %>%
  mutate(age_missing=ifelse(is.na(age), "y", "n"))



missing_age %>% covar_table("age_missing")


(chisq.age_m_sex <- chisq.test(table(missing_age$expr_sex, missing_age$age_missing))) # p=0.3797

(chisq.age_m_smok <- chisq.test(table(missing_age$smok, missing_age$age_missing))) # p=0.1894