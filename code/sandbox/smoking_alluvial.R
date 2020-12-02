# smoking_alluvial.R
# E Flynn
# 12/17/2019
#
# load the results of smoking sex labeling 
# get a bunch of counts (which studies failed?) 
# and make an alluvial diagram


require('tidyverse')
require('ggalluvial')

# read in the text metadata-based labels
dedup_smok2 <- read_csv("data/smok_dat/smokDat_dedup2.csv")
dedup_smok2$X <- NULL
ale_dat <- read_csv("data/ale_combined_data.csv")
smok_ale <- left_join(dedup_smok2 %>% select(-organism, -title), 
                      ale_dat %>% select(gse, gsm, text_sex), 
                      by=c("gse", "gsm")) # 15455

# get some platform counts, get some failed counts

# ---- reading in the exprsex labels ---- #

smok_sex_labels <- read_tsv("data/smok_dat/smok_sex_labels.txt")
smok_sex_labels2 <- smok_sex_labels %>% 
  rename(sex_label=res, study.sample=V2) 

# create a new column with the study/sample labels separated by an underscore
# so that we can separate these out
# (we could not before b/c sometimes study/gpl were also sep by ".")
smok_sex_labels2$new_col <- sapply(smok_sex_labels2$study.sample, 
                               function(x) 
                                 {y <- strsplit(x, "\\.")[[1]]; 
                                  len_y <- length(y);
                                  z <- paste(y[[1]], 
                                             y[[len_y]], sep="_");
                                  return(z)})

smok_sex_labels3 <- smok_sex_labels2 %>% 
  separate(new_col, into=c("gse", "gsm"), sep="_") %>%
  select(-study.sample) %>%
  filter(!is.na(sex_label))

# we have sex labels for 12,548 samples, 168 studies -- fairly good??

labeled <- left_join(smok_ale, smok_sex_labels3, by=c("gse", "gsm"))

labeled2 <- labeled %>% 
  rename(expr_sex = sex_label) %>%
  mutate(text_sex = 
           case_when(text_sex=="M" ~ "male",
                     text_sex=="F" ~"female"),
         expr_sex =
           case_when(expr_sex==1 ~ "male",
                     expr_sex==0 ~ "female")) 

write_csv(labeled2, "data/smok_dat/smoking_labels_reform.csv")

flow_freq_counts <- labeled2 %>% 
  ungroup() %>% 
  group_by(text_sex, expr_sex) %>% 
  mutate(Freq=n()) %>% 
  select(-gse, -gpl, -gsm) %>% 
  unique() %>%
  ungroup() %>%
  mutate(text_sex=ifelse(is.na(text_sex), "unlabeled", text_sex),
         expr_sex=ifelse(is.na(expr_sex), "unlabeled", expr_sex)) %>% 
  mutate(row_id=1:n()) %>%
  gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("text_sex", "expr_sex")),
         sex=as.factor(sex)) %>%
  unique() 

ggplot(flow_freq_counts,
       aes(x = labeling_method, 
           stratum = sex, alluvium = row_id,
           y = Freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of samples")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 



### --- GROUP --- ####

comb_labels_long <- labeled2 %>% 
  select(-gpl) %>%
  rename(metadata=text_sex, exprsex=expr_sex) %>%
  gather(key="labeling_method", value="sex", -gse, -gsm) %>%
  group_by(gse, labeling_method) %>%
  summarize(num_samples=n(),
            num_f=sum(sex=="female"),
            num_m=sum(sex=="male")) %>%
  mutate(study_type= case_when(
    (is.na(num_f) & is.na(num_m)) ~ "unlabeled",
    (!is.na(num_f) & !is.na(num_m) & num_f > 0 & num_m > 0 ) ~ "mixed",
    (!is.na(num_f) & num_f > 0 ) ~ "female-only",
    (!is.na(num_m) & num_m > 0 )~ "male-only")) %>%
  mutate(freq=1) %>% 
  ungroup(gse) %>%  
  mutate(study_type=as.factor(study_type), 
         labeling_method=factor(labeling_method, 
                                levels=c("metadata", "exprsex")), 
         gse=as.factor(gse)) %>% 
  rename(study=gse) %>%
  rename(sex=study_type)


ggplot(comb_labels_long,
       aes(x = labeling_method, 
           stratum = sex, 
           alluvium = study,
           y = freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of studies")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 


# 80/20
large_study <- labeled2 %>% 
  select(-gpl) %>%
  rename(metadata=text_sex, exprsex=expr_sex) %>%
  gather(key="labeling_method", value="sex", -gse, -gsm) %>%
  group_by(gse, labeling_method) %>%
  summarize(num_samples=n(),
            num_f=sum(sex=="female"),
            num_m=sum(sex=="male")) %>%
  mutate(study_type= case_when(
    (is.na(num_f) & is.na(num_m)) ~ "unlabeled",
    (!is.na(num_f) & !is.na(num_m) & num_f/num_samples > 0.8 & num_m > 0 ) ~ "mostly-female",
    (!is.na(num_f) & !is.na(num_m) & num_m/num_samples > 0.8 & num_f > 0 ) ~ "mostly-male",
    (!is.na(num_f) & !is.na(num_m) & num_f > 0 & num_m > 0 ) ~ "mixed",
    (!is.na(num_f) & num_f > 0 ) ~ "female-only",
    (!is.na(num_m) & num_m > 0 ) ~ "male-only")) %>%
  mutate(freq=1) %>% 
  ungroup(gse) %>%  
  mutate(study_type=factor(study_type, 
                           levels=c("female-only", "mostly-female", "mixed", "mostly-male", "male-only", "unlabeled")), 
         labeling_method=factor(labeling_method, 
                                levels=c("metadata", "exprsex")), 
         gse=as.factor(gse)) %>% 
  rename(study=gse) %>%
  rename(sex=study_type)


ggplot(large_study,
       aes(x = labeling_method, 
           stratum = sex, 
           alluvium = study,
           y = freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of studies")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 


### - remove the cell line studies - focus on the bin ###
