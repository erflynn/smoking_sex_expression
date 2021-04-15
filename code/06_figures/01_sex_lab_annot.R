# 01_sex_lab_annot.R
# annotation breakdown
# grab the sex labeling results for all the data we kept and plot
#
# figures generated:
#   - alluvial for breakdown
#   - tissue bar plot
#   - bar plots for breakdown (s1?)
# tables generated:
#   - s2a annot breakdown counts
#   - s2b trt counts
#   - s2c tissue counts
#   [- s2d trt cl types?]
#   - s3a+b sample + study breakdown
# intermed files generated:
#  - data/smok_study_sex_lab.csv  -- list of all study sex lab


library(tidyverse)
library(ggalluvial)
library(RColorBrewer)
dark2 <- brewer.pal(8, "Dark2")
set2 <- brewer.pal(8, "Set2")
alluv_col <- c(dark2[1], set2[1], set2[4], set2[3], dark2[3], dark2[8])


# --- 1. read in the list of studies/samples that we are working with --- #
annot_studies <- read_csv("data/supp_tables/supp_table_1_annot.csv")
smok_sl <- read_csv("data/smok_samples_w_sl.csv", col_types="clddccccc")
(annot_counts <- annot_studies   %>%
  group_by(included, study_type) %>% 
  summarise(num_studies=n(), 
            num_samples=sum(num_samples)) %>%
  arrange(desc(included)))
annot_counts %>%
  write_csv("data/supp_tables/s2_annot_counts.csv") # SUPPLEMENTARY TABLE 2

# make a plot showing the sex breakdown / labeling of the filtered
smok_tab3 <- read_csv("data/supp_tables/s4_list_smok_studies.csv")
smok_filt_sl <- smok_sl %>% 
  filter(study_acc %in% smok_tab3$study_acc) %>%
  select(study_acc, sample_acc, expression, metadata)

# --- 4. breakdown by tissue type --- #
kept_history <- annot_studies %>% 
  filter(study_type=="smoking history") 
kept_history %>% group_by(source) %>% count()
(tiss <- kept_history %>%
  separate_rows(tissue, sep=";") %>%
  group_by(tissue) %>% count() %>%
  arrange(desc(n)) %>%
  rename( num_studies=n))
tiss %>% write_csv("data/supp_tables/s2_tiss_breakdown.csv")

tiss %>%
  filter(num_studies>1) %>%
  mutate(tissue=factor(tissue, levels=c(tiss %>% filter(num_studies>1) %>% pull(tissue)))) %>%
  ggplot(aes(x=tissue, y=num_studies, fill=tissue))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(legend.position="none")+
  theme(axis.text.x=element_text(angle = 90, vjust=0.5, hjust=1) )+
  ylab("number of studies")+
  xlab("")
ggsave("figures/paper_figs/s1_tiss_breakdown.pdf")

# --- 5. breakdown by smoke exposure + cell line --- #
treated_cells <- annot_studies %>%
  filter(study_type=="treated cells" & included=="yes") 


(trt_cls_count <- treated_cells %>% 
  separate_rows(treatment, sep="; ") %>% 
  group_by(treatment) %>% 
  count() %>%
  arrange(desc(n)) %>%
  rename(num_studies=n))
trt_cls_count %>%
  write_csv("data/supp_tables/s2_trt_breakdown.csv")



# ---- 6. divide up and make alluvials + tables ---- #
sex_lab_kept <- smok_sl %>% 
  filter(keep=="yes" & type != "treated cells")
sex_lab_tc <- smok_sl %>% 
  filter(keep=="yes" & type == "treated cells")
sex_lab_all <- smok_sl %>% 
  filter(keep=="no" & type == "all smokers")
sex_lab_non <- smok_sl %>% 
  filter(keep=="no" & type == "all nonsmokers")
# CHECK PMALE + NUM_READS
# --- 4. make alluvial diagrams --- #

# pick a color scheme?
# sample data
plotAlluvialSample <- function(df){
  flow_freq_counts <- df %>% 
    distinct(sample_acc, expression, metadata) %>%
    group_by(metadata, expression) %>% 
    mutate(Freq=n()) %>% 
    select(-sample_acc) %>%
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
    scale_fill_manual(values=dark2[c(1,3,8)])+
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    xlab("Label source")+ylab("Number of samples")+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + 
    theme(legend.position = "none") 
}
plotAlluvialSample(smok_filt_sl)
ggsave("figures/paper_figs/figs5_de_sample_breakdown.pdf")

plotAlluvialSample(sex_lab_kept)
ggsave("figures/paper_figs/fig2_human_smok.pdf")

plotAlluvialSample(sex_lab_tc)
ggsave("figures/paper_figs/fig2_tc_smok.pdf")

plotAlluvialSample(sex_lab_all)
ggsave("figures/paper_figs/s_fig2_all_smok.pdf")

plotAlluvialSample(sex_lab_non)
ggsave("figures/paper_figs/s_fig2_non_smok.pdf")

### --- GROUP --- ####
study_sex_lab <- function(df){
  df %>% 
    select(-present, -num_reads, -p_male, -keep, -type) %>%
    pivot_longer(c(expression, metadata), names_to="labeling_method", values_to="sex") %>%
    group_by(study_acc, labeling_method) %>%
    summarize(num_samples=n(),
              num_f=sum(sex=="female"),
              num_m=sum(sex=="male"),
              num_unlab=sum(sex=="unlabeled")) %>%
    mutate(study_type= case_when(
      num_unlab/num_samples > 0.5 & (num_samples-num_unlab) < 30 ~ "unlabeled",
      (!is.na(num_f) & !is.na(num_m) & num_f/num_samples > 0.8 & num_m > 0 ) ~ "mostly-female",
      (!is.na(num_f) & !is.na(num_m) & num_m/num_samples > 0.8 & num_f > 0 ) ~ "mostly-male",
      (!is.na(num_f) & !is.na(num_m) & num_f > 0 & num_m > 0 ) ~ "mixed sex",
      (!is.na(num_f) & num_f > 0 ) ~ "female-only",
      (!is.na(num_m) & num_m > 0 ) ~ "male-only"))
}

studyAlluvial <- function(df){
  comb_labels_long <- df  %>%
    study_sex_lab() %>%
    mutate(freq=1) %>% 
    ungroup(study_acc) %>%  
    mutate(study_type=as.factor(study_type), 
           labeling_method=factor(labeling_method, 
                                  levels=c("metadata", "expression")), 
           study_acc=as.factor(study_acc)) %>% 
    rename(study=study_acc) %>%
    rename(sex=study_type) %>%
    mutate(sex=factor(sex, levels=c("female-only", "mostly-female", "mixed sex", "mostly-male", "male-only", "unlabeled")))
  
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
    theme(legend.position = "none") +
    scale_fill_manual(values=alluv_col)
  
}
studyAlluvial(sex_lab_kept)
ggsave("figures/paper_figs/fig2_study_smok.pdf")
studyAlluvial(sex_lab_tc)+
  scale_fill_manual(values=alluv_col[c(1,3,5,6)])
ggsave("figures/paper_figs/fig2_tc_study_smok.pdf")
studyAlluvial(sex_lab_all)+
  scale_fill_manual(values=alluv_col[c(1,3,4,5,6)])
ggsave("figures/paper_figs/s_fig2_study_all_smok.pdf")


studyAlluvial(sex_lab_non)+
  scale_fill_manual(values=alluv_col[c(1,2,3,5,6)])
ggsave("figures/paper_figs/s_fig2_study_non_smok.pdf")


studyAlluvial(smok_filt_sl %>% mutate(present=NA, num_reads=NA, p_male=NA, keep=NA, type=NA))+
  scale_fill_manual(values=c(alluv_col[1], alluv_col[3], alluv_col[4], alluv_col[5], alluv_col[6]))
ggsave("figures/paper_figs/figs5_de_study_breakdown.pdf")

# --- TABLES --- #
get_sample_counts <- function(df){
  df %>% 
    distinct(sample_acc, expression, metadata) %>%
    pivot_longer(c(expression, metadata), names_to="label_type", values_to="sex") %>%
    group_by(label_type, sex) %>% 
    count()
}
sample_sl <- get_sample_counts(sex_lab_kept) %>% 
  mutate(dataset="smoking history") %>%
  bind_rows(get_sample_counts(sex_lab_tc) %>%
              mutate(dataset="treated cells")) %>%
  bind_rows(get_sample_counts(sex_lab_all) %>%
              mutate(dataset="all smokers")) %>%
  bind_rows(get_sample_counts(sex_lab_non) %>%
              mutate(dataset="all nonsmokers")) %>%
  ungroup() %>% 
  select(dataset, everything()) %>%
  pivot_wider(names_from="sex", values_from="n") %>%
  arrange(dataset, desc(label_type)) %>%
  mutate(total_samples=female+male+unlabeled) %>%
  mutate(across(c(female,male,unlabeled), ~./total_samples, .names="frac_{col}") ) %>%
  rename(num_female=female, num_male=male, num_unlabeled=unlabeled)

sample_sl %>%
  mutate(across(contains("frac"), ~round(., digits=3))) %>%
  select(dataset, label_type, total_samples, everything()) %>%
  write_csv("data/supp_tables/s3a_sample_sex_breakdown.csv")
  

tab_study_sex_lab <- study_sex_lab(sex_lab_kept) %>%
  mutate(dataset="smoking history") %>%
  bind_rows(study_sex_lab(sex_lab_tc) %>%
              mutate(dataset="treated cells")) %>%
  bind_rows(study_sex_lab(sex_lab_all) %>%
              mutate(dataset="all smokers")) %>%
  bind_rows(study_sex_lab(sex_lab_non) %>%
              mutate(dataset="all nonsmokers"))

tab_study_sex_lab %>% 
  select(dataset, everything()) %>%
  write_csv("data/smok_study_sex_lab.csv")
study_sex_lab_counts <- tab_study_sex_lab %>%
  ungroup() %>%
  rename(label_type=labeling_method) %>%
  group_by(dataset, label_type, study_type) %>%
  count() %>%
  group_by(dataset, label_type) %>%
  mutate(total_samples=sum(n) ) %>%
  pivot_wider(names_from="study_type", values_from="n", values_fill=0) %>%
  arrange(dataset, desc(label_type)) %>%
  mutate(across(`female-only`:`mostly-male`, ~./total_samples, .names="frac_{col}") ) 

study_sex_lab_counts %>%
  mutate(across(contains("frac"), ~round(., digits=3))) %>%
  write_csv("data/supp_tables/s3b_study_sex_breakdown.csv")

plot_sex_breakdown <- function(df) {
  df %>% 
    mutate(dataset=factor(dataset, levels=c("smoking history", "treated cells", "all smokers", "all nonsmokers"))) %>%
    mutate(study_type=factor(study_type, levels=c("female-only", "mostly-female", "mixed sex", "mostly-male", "male-only", "unlabeled"))) %>%
    rename(study_sex=study_type) %>%
    group_by(dataset, study_sex) %>%
    count() %>%
    ggplot(aes(x=dataset, y=n, fill=study_sex))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab("number of studies")+
    xlab("")+
    scale_fill_manual(values=alluv_col)+
    theme(axis.text.x=element_text(angle = 90, vjust=0.5, hjust=1))
}
tab_study_sex_lab %>% 
  filter(labeling_method=="expression") %>% 
  plot_sex_breakdown()
ggsave("figures/paper_figs/s1_annot_sex_breakdown.pdf")

tab_study_sex_lab %>% 
  filter(labeling_method=="metadata") %>% 
  plot_sex_breakdown()
ggsave("figures/paper_figs/s1_annot_sex_breakdown_metadata.pdf")


