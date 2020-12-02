

manual_annot <- readxl::read_excel("data/smok_dat/smok_data_manual_annot_0319.xlsx")

sex_lab <- read_csv("data/smok_dat/sex_labels_w_rb.csv") 
sex_lab2 <- sex_lab %>% select(gsm,consensus_sex) %>% inner_join(gse_gsm) %>%
  group_by(gse) %>%
  summarize(num_f=sum(consensus_sex=="female"),
            num_m=sum(consensus_sex=="male"),
            num_samples=n()) %>%
  mutate(study_sex=case_when(
    num_f==num_samples ~ "female-only",
    num_m==num_samples ~ "male-only",
    (num_m+num_m) < (0.5*num_samples) ~ "unknown",
    num_m/num_samples >= 0.8  ~ "mostly-male",
    num_f/num_samples >= 0.8  ~ "mostly-female",
    num_f/num_samples <= 0.2 ~ "mostly-male",
    num_m/num_samples <= 0.2 ~ "mostly-female",
    num_f==0 & num_m!=0 ~"male-only",
    num_f!=0 & num_m==0 ~"female-only",
    TRUE ~ "mixed"
  )) %>%
  mutate(study_sex=factor(study_sex, levels=c("female-only", 
                                              "mostly-female", 
                                              "mixed", "mostly-male", 
                                              "male-only", "unknown")))

# combine the treatments into fewer

adj_trt <- trt_cls %>% 
  select(treatment) %>% unique() %>%
  mutate(adj_trt=case_when(
    treatment %in% c("cigarette component","nicotine","tobacco carcinogens") ~ "cigarette component",
    treatment %in% c("cigarette smoke","whole cigarette smoke") ~ "whole cigarette smoke",
    treatment %in% c("cigarette smoke condensate", 
                       "cigarette smoke extract", 
                       "cigarette smoke extract; other components",
                       "tobacco smoke extract") ~ "cigarette smoke extract",
    treatment %in% c("tobacco smoke exposure") ~ "whole tobacco smoke",
    TRUE ~ treatment
  )) 

new_breakdown <- manual_annot %>% filter(keep=="yes") %>% left_join(sex_lab2)
ggplot(new_breakdown %>% filter(study_sex!="unknown"), aes(x=1))+geom_bar(aes(fill=study_sex))+ylab("num studies")

trt_cls <- manual_annot %>% filter(keep=="yes" & type=="treated cells") %>%
  left_join(adj_trt)
tissues <- manual_annot %>% filter(keep=="yes" & type=="smoking")

tissues2 <- tissues %>% select(gse, tissue2) %>% left_join(sex_lab2) %>% separate_rows(tissue2, sep=";")%>% group_by(tissue2) %>% mutate(num_studies=n())
tissues2 %>% filter(tissue2=="buccal mucosa" & study_sex=="mixed") %>%
  ungroup() %>%
  select(gse, num_f, num_m, num_samples)

be <- tissues2 %>% filter(tissue2=="bronchial epithelium" & study_sex=="mixed") %>%
  ungroup() %>%
  select(gse, num_f, num_m, num_samples)

full%>% select(gse,title, design, summary) %>% right_join(be, by="gse")  %>% write_csv("be_tmp.csv")



full%>% select(gse,title, design, summary) %>% 
  right_join(tissues2 %>% 
               filter(tissue2=="lung" & study_sex=="mixed") %>% 
               ungroup() %>%
               select(gse, num_f, num_m, num_samples), by="gse")  %>% write_csv("lung_tmp.csv")

full%>% select(gse,title, design, summary) %>% 
  right_join(tissues2 %>% 
               filter(tissue2=="airway epithelium" & study_sex=="mixed") %>% 
               ungroup() %>%
               select(gse, num_f, num_m, num_samples), by="gse")  %>% write_csv("ae_tmp.csv")


full%>% select(gse,title, design, summary) %>% 
  right_join(tissues2 %>% 
               filter(tissue2=="small airway epithelium" & study_sex=="mixed") %>% 
               ungroup() %>%
               select(gse, num_f, num_m, num_samples), by="gse")  %>% write_csv("sae_tmp.csv")

full%>% select(gse,title, design, summary) %>% 
  right_join(tissues2 %>% 
               filter(tissue2=="bronchial alveolar lavage" & study_sex=="mixed") %>% 
               ungroup() %>%
               select(gse, num_f, num_m, num_samples), by="gse")  %>% write_csv("bae_tmp.csv")



ggplot(tissues2 %>% filter(num_studies >=2), aes(x=tissue2))+geom_bar(aes(fill=study_sex))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("number of studies")+xlab("tissue")


trt_cl_breakdown <- trt_cls %>%
ggplot(trt_cl_breakdown %>% separate_rows(adj_trt, sep=";"), aes(x=adj_trt))+geom_bar(aes(fill=study_sex))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("number of studies")+xlab("treatment")

ggplot(trt_cl_breakdown, aes(x=1))+geom_bar(aes(fill=study_sex))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("number of studies")+xlab("treatment")

mixed_sex <- trt_cl_breakdown %>% 
  separate_rows(adj_trt, sep=";") %>%
  filter(study_sex=="mixed")

trt_cl_breakdown %>% 
  separate_rows(adj_trt, sep=";") %>%
  filter(study_sex=="mixed" & adj_trt=="cigarette smoke extract") %>%
  select(-adj_trt, -tissue2, -study_sex) %>%
  left_join(full%>% select(gse,title, design, summary)) %>%
  write_csv("tmp_cse.csv")
# NOW

require('GEOquery')
gse14383 <- getGEO("GSE14383")
gse14385 <- getGEO("GSE14385")

# they're male... BUT loss of Y chromsomes is a thing
pData1 <- pData(gse14383$GSE14383_series_matrix.txt.gz)
pData1.2 <- pData1 %>% 
  select(geo_accession, characteristics_ch1.1, characteristics_ch1.2, characteristics_ch1.3, characteristics_ch1.4) %>% 
  rename(gsm=geo_accession, trt=characteristics_ch1.1, month=characteristics_ch1.2, replicate=characteristics_ch1.3, batch=characteristics_ch1.4) %>%
  left_join(sex_lab %>% select(gsm, pred, rb_sex)) %>%
  mutate(gse="GSE14383")

pData1.2 %>% View()
# why / how are these sex split?

pData2 <- pData(gse14385$GSE14385_series_matrix.txt.gz)

pData2.2 <- pData2 %>% 
  select(geo_accession, characteristics_ch1.1, characteristics_ch1.2, characteristics_ch1.3) %>%
  rename(gsm=geo_accession, trt=characteristics_ch1.1, demethyl=characteristics_ch1.2, replicate=characteristics_ch1.3) %>%
  left_join(sex_lab %>% select(gsm, pred, rb_sex)) %>%
  mutate(gse="GSE14385")
pData2.2  %>% View()

ggplot(rbind(pData1.2 %>% select(gse, pred, rb_sex), 
             pData2.2%>% select(gse, pred, rb_sex)), 
       aes(x=pred))+geom_histogram(aes(fill=rb_sex))+
  facet_wrap(.~gse, ncol=1)+
  xlab("P(male)")+ylab("number of samples")

cl2 <- trt_cl_breakdown %>% 
  separate_rows(adj_trt, sep=";") %>%
  filter(study_sex=="mixed" & adj_trt=="cigarette smoke extract") %>%
  select(-adj_trt, -tissue2, -study_sex) %>%
  inner_join(gse_gsm) %>%
  left_join(sex_lab %>% select(gsm, pred, rb_sex)) 

ggplot(cl2,aes(x=pred))+geom_histogram(aes(fill=rb_sex))+
  facet_wrap(.~gse, ncol=1)+
  xlab("P(male)")+ylab("number of samples")+

plot(density(pData2.2$pred, na.rm=TRUE))

# how do I represent this??

require('MetaIntegrator')
# do buccal mucosa fun!


