library(tidyverse)

library('LaCroixColoR')
paired2 <- lacroix_palette(type = "paired")
paired3 <- c(paired2[[1]], paired2[[2]], paired2[[11]], paired2[[12]], 
             paired2[[5]], paired2[[6]], paired2[[3]], paired2[[4]],
             paired2[[10]],paired2[[13]], paired2[[14]],
             paired2[[7]], paired2[[8]], paired2[[9]])
save(paired3, file="data/fig_color_scale.RData")


counts2 <- full_pair2 %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)


counts2$study1=factor(counts2$study1, levels=unique(counts2$study1))
counts2$study2=factor(counts2$study2, levels=unique(counts2$study1))



counts2_sex <- full_pair2_sex %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)


counts2_sex$study1=factor(counts2_sex$study1, levels=unique(counts2_sex$study1))
counts2_sex$study2=factor(counts2_sex$study2, levels=unique(counts2_sex$study1))


counts2_sex_auto <- autosomal_sex_cts %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)


counts2_sex_auto$study1=factor(counts2_sex_auto$study1, levels=unique(counts2_sex_auto$study1))
counts2_sex_auto$study2=factor(counts2_sex_auto$study2, levels=unique(counts2_sex_auto$study1))


#library('RColorBrewer')
#set3 <- brewer.pal(12, "Paired")
#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A"
# "#FFFF99" "#B15928"

ggplot(counts2)+
  geom_bar(mapping = aes(x = study1, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=paired3)
ggsave("figures/smok_bar.png")

counts2.2 <- counts2 %>% 
  filter(as.numeric(study1) < as.numeric(study2))
# make it only the top half

ggplot(counts2.2,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=rep_sig)) +
  geom_text(aes(label=rep_sig), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
  scale_x_discrete(position = "top") +
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/gene_overlap2.png")

counts2_sex <- counts2_sex %>% 
  select(-tissue) %>%
  
  filter(study1!="GSE44456", study2!="GSE44456",
         !study1 %in% c("GSE87072", "GSE18723", "GSE16149"),
         !study2 %in% c("GSE87072", "GSE18723", "GSE16149")
  ) %>%
  inner_join(list_studies %>% dplyr::select(study, tissue), 
             by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", "blood - b cell",
                                        "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study1)
counts2_sex$study1=factor(counts2_sex$study1, levels=unique(counts2_sex$study1))
counts2_sex$study2=factor(counts2_sex$study2, levels=unique(counts2_sex$study1))



counts2.2_sex <- counts2_sex %>% 
  filter(as.numeric(study1) < as.numeric(study2))
ggplot(counts2.2_sex)+
  geom_bar(mapping = aes(x = study1, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=paired3[c(1:8,10:13)])
ggsave("figures/sex_bar.png")


ggplot(counts2.2_sex,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=rep_sig)) +
  geom_text(aes(label=rep_sig), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
  scale_x_discrete(position = "top") +  
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/gene_overlap2_sex.png")
# ... how many are autosomal?

counts2.2_sex_auto <- counts2_sex_auto %>% 
  filter(as.numeric(study1) < as.numeric(study2)) 
ggplot(counts2.2_sex_auto,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=n)) +
  geom_text(aes(label=n), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
  scale_x_discrete(position = "top") +  
  coord_fixed()+
  ylab("")+
  xlab("")
# "GSE20681", "E-MTAB-5278" overlap on 4
# GSE21862    GSE20681 on 3

