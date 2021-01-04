
require('tidyverse')
# read in the old data
smoking_labels <- read_csv("data/smok_dat/smoking_labels_reform.csv")
gse_gsm <- smoking_labels %>% select(gse, gsm) %>% unique()


smoking_labels2 <- smoking_labels %>%
  select(gsm, text_sex, expr_sex) %>% unique() %>%
  group_by(gsm) %>%
  mutate(text_sex=ifelse(all(is.na(text_sex)), NA, text_sex[!is.na(text_sex)]),
         expr_sex=ifelse(all(is.na(expr_sex)), NA, text_sex[!is.na(expr_sex)])) %>%
  unique()

# update labels
human_all_sl <- read_csv("../../drug_trt/data/02_labeled_data/human_all_sl.csv")
length(intersect(smoking_labels2$gsm, human_all_sl$id))

length(intersect(smoking_labels$gsm, human_all_sl$id))
# 8984
length(smoking_labels2$gsm)
# 15275
length(setdiff(smoking_labels$gsm, human_all_sl$id)) # 6291

options(stringsAsFactors=FALSE)

human_exp_to_samp <- read_csv("../../drug_trt/data/01_metadata/human_exp_to_sample.csv")
metadata <- read.csv("../../drug_trt/data/01_metadata/human_metadata.csv")

# missingness of gses
smok_gses_all <- read_csv("data/smok_dat/smoking_data.csv") %>% select(gse) # 220
smok_gsms <- smok_gses_all %>% left_join(human_exp_to_samp, by=c("gse"="study_acc"))
missing_gses <- setdiff(smok_gses_all$gse, human_exp_to_samp$study_acc) # 53 :(
## it appears that a large number of these *ARE* missing :(

# combine the metadata!


# possibly due to gene missingness, etc
labels2 <- smoking_labels2 %>% 
  left_join(human_all_sl, by=c("gsm"="id")) %>% 
  left_join(metadata %>% select(acc, sex) %>% 
              rename(text_sex2=sex), by=c("gsm"="acc"))
labels3 <- labels2 %>% 
  mutate(rb_sex=ifelse(is.na(pred), NA,
                       ifelse(pred > 0.5, "male", "female"))) %>%
  mutate(text_sex2=as.character(text_sex2)) %>%
  mutate(text_sex2=
           case_when(text_sex2 %in% c("na", "") ~ "NA",
                     text_sex2==1 ~ "male",
                     text_sex2==0 ~ "female",
                     TRUE ~ text_sex2)) %>%
  mutate(text_sex2=ifelse(text_sex2=="NA", NA, text_sex2))

labels3 %>% filter(is.na(pred)) %>% group_by(gse) %>% 
  summarize(missing=sum(is.na(pred)), present=sum(!is.na(pred)))%>% arrange(desc(missing))

# add metadata sex labels too

# how many studies change?
table(labels3$text_sex==labels3$text_sex2)
labels3 %>% filter(text_sex != text_sex2) 
# for these we trust text_sex2...
labels3 %>% filter(is.na(text_sex2) & !is.na(text_sex))

sum(labels3$expr_sex==labels3$rb_sex, na.rm=TRUE) # 1535  5989 (79.6% match)
labels3 %>% filter(!is.na(rb_sex) & expr_sex!=rb_sex)
table(labels3$expr_sex==labels3$rb_sex) # 84.5% 
labels3 %>% filter(pred > 0.9 | pred < 0.1) %>% 
  group_by(expr_sex, rb_sex) %>% count()
table(labels3$text_sex==labels3$rb_sex) # 84.5% 
table(labels3$text_sex==labels3$expr_sex) # 79.5% 
table(labels3$text_sex2==labels3$rb_sex) # 88.6%
table(labels3$text_sex2==labels3$expr_sex) # 80.7% 

# complicated consensus game
labels4 <- labels3 %>%
  mutate(consensus_sex=
           case_when(
             is.na(rb_sex) & !is.na(expr_sex) ~ expr_sex,
             is.na(expr_sex) & !is.na(rb_sex) ~ rb_sex,
             expr_sex==rb_sex ~ expr_sex,
             !is.na(text_sex2) & rb_sex==text_sex2 ~ rb_sex,
             !is.na(text_sex2) & expr_sex==text_sex2 ~ expr_sex,
             !is.na(text_sex) & rb_sex==text_sex ~ rb_sex,
             !is.na(text_sex) & expr_sex==text_sex ~ expr_sex,
             is.na(rb_sex) & is.na(expr_sex) ~ "unknown",
             (pred <= 0.2 | pred >=0.8) ~ rb_sex,
             TRUE ~ "mismatch"
           ))
labels4 %>% write_csv("data/smok_dat/sex_labels_w_rb.csv")

# make new alluvial diagram
require('ggalluvial')


label4.1 <- labels4 %>% 
  select(gsm, text_sex2, consensus_sex) %>%
  rename(text_sex=text_sex2, expr_sex=consensus_sex)
  
alluvPlotSample <- function(ds){
  flow_freq_counts <- ds %>% 
    ungroup() %>% 
    rename(metadata=text_sex, expression=expr_sex) %>%
    group_by(metadata, expression) %>% 
    mutate(Freq=n()) %>% 
    select(-gsm) %>% 
    unique() %>%
    ungroup() %>%
    mutate(row_id=1:n()) %>%
    gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
    mutate(row_id=as.factor(row_id), 
           labeling_method=factor(labeling_method, levels=c("metadata", "expression")),
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
}

alluvPlotSample(label4.1 %>% 
                  mutate(text_sex=ifelse(is.na(text_sex), "unknown", text_sex),
                         expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)))
label4.2 <- labels4 %>% 
  select(gsm, text_sex2, expr_sex) %>%
  rename(text_sex=text_sex2)
alluvPlotSample(label4.2 %>% 
                  mutate(text_sex=ifelse(is.na(text_sex), "unknown", text_sex),
                         expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)))
label4.3 <- labels4 %>% 
  select(gsm, text_sex2, rb_sex) %>%
  rename(text_sex=text_sex2, expr_sex=rb_sex)
alluvPlotSample(label4.3 %>% 
                  mutate(text_sex=ifelse(is.na(text_sex), "unknown", text_sex),
                         expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)))

###
alluvPlotStudy <- function(ds){
  large_study <- ds %>% 
    mutate(text_sex=ifelse(text_sex %in% c("unknown", "mixed"), NA, text_sex)) %>%
    rename(metadata=text_sex, expression=expr_sex) %>%
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
                                  levels=c("metadata", "expression")), 
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
  
}
# are samples mislabeled?




ggplot(label4.1 %>% filter(expr_sex!="unknown"), aes(x=1))+geom_bar(aes(fill=expr_sex))+ylab("number of samples")


swap <- labels4 %>% filter(consensus_sex!=text_sex2)
counts_per_study <- swap %>% inner_join(gse_gsm) %>% group_by(gse) %>% count() %>% arrange(desc(n))
# across 12 studies, (517 in 1)
ggplot(swap, aes(x=pred))+geom_histogram()+ylab("number of samples")+xlab("P(male)")


swap_df <- counts_per_study %>% ungroup() %>% left_join(gse_gsm) %>% 
  left_join(labels4 %>% select(gsm, pred, text_sex2)) %>%
  rename(text_sex=text_sex2, study=gse) %>%
  filter(!is.na(text_sex))

swap_stat <- swap_df %>% group_by(study, text_sex) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred)) %>%
  mutate(mu_l=mu-3*sigma, mu_u=mu+3*sigma) %>%
  mutate(mu_2l=mu-2*sigma, mu_2u=mu+2*sigma) %>%
  ungroup()
swap_stat2 <- swap_stat %>%
  pivot_longer(c("mu_2l", "mu_2u"), names_to="stat", values_to="stat_val") %>%
  mutate(stat=paste(stat,text_sex, sep="_"))  %>%
  select(-text_sex) %>%
  pivot_wider(id_cols = "study", names_from="stat", values_from="stat_val") %>%
  mutate(sex_sep=ifelse(mu_2u_female >= mu_2l_male, FALSE, TRUE))

likely_mislabeled <- 
  swap_df %>%
  left_join(swap_stat, by=c("study", "text_sex")) %>% 
  filter(pred < mu_2l | pred > mu_2u) %>%
  left_join(swap_stat2 %>% select(study, sex_sep)) %>%
  filter(sex_sep) 


ggplot(swap_df,
       aes(x=pred, y=0))+
  geom_point(alpha=0.3, aes(col=text_sex))+
  facet_grid(rows=vars(study), scales="free") +xlab("P(male)")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  theme(strip.text.y = element_text(angle = 0))

