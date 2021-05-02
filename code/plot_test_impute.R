
library(tidyverse)
summ_imp <- read_csv("data/summary_imp.csv")
summ_imp2 <- summ_imp %>% rename(true_sm_overlap=imp_sm_overlap,
                    imp_sm_overlap=imp_sm_overlap.1) %>%
  mutate(covars=str_extract(run_id, "[A-z]+"))

summ_imp2 %>%
  mutate(imputed=true_imp_overlap/num_true,
         small=true_sm_overlap/num_true) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="fraction") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, fraction, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
ggplot(aes(x=factor(frac_miss),
           y=fraction,
           group=grp,
           col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("fraction of true probes identified")
ggsave("figures/imputed_data_overlap.png")


summ_imp2 %>%
  mutate(imputed=(num_imp-true_imp_overlap)/num_imp,
         small=(num_sm-true_sm_overlap)/num_sm) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="fraction") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, fraction, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
  ggplot(aes(x=factor(frac_miss),
             y=fraction,
             group=grp,
             col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("fraction of false probes identified")
ggsave("figures/imputed_data_fp.png")


summ_imp2 %>%
  mutate(imputed=true_imp_overlap/num_imp,
         small=true_sm_overlap/num_sm) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="fraction") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, fraction, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
  ggplot(aes(x=factor(frac_miss),
             y=fraction,
             group=grp,
             col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("fraction of TP")
ggsave("figures/imputed_data_tp.png")



summ_imp2 %>%
  mutate(imputed=est_cor_imp,
         small=est_cor_sm) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="es_cor") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, es_cor, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
  ggplot(aes(x=factor(frac_miss),
             y=es_cor,
             group=grp,
             col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("effect size correlation")
ggsave("figures/imputed_data_es_cor.png")


summ_imp2 %>%
  mutate(imputed=p_cor_impP,
         small=p_cor_smP) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="p_cor") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, p_cor, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
  ggplot(aes(x=factor(frac_miss),
             y=p_cor,
             group=grp,
             col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("p-value correlation")
ggsave("figures/imputed_data_p_cor.png")


summ_imp2 %>%
  mutate(imputed=p_cor_impS,
         small=p_cor_smS) %>%
  pivot_longer(c(imputed, small), names_to="source", 
               values_to="p_corS") %>%
  unite("grp", c(source, frac_miss), remove=F) %>%
  select(grp, frac_miss, covars, cutoff, p_corS, source) %>%
  mutate(covars=str_replace_all(covars, "_", " "),
         covars=str_replace_all(covars, "re", "race"),
         covars=str_replace_all(covars, "pk", "pkyrs"),
         covars=str_trim(covars),
         covars=str_replace_all(covars, " ", ", "),
         cutoff=paste("FDR <", cutoff)) %>%
  ggplot(aes(x=factor(frac_miss),
             y=p_corS,
             group=grp,
             col=source))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.4, 
                                           dodge.width = 0.9), alpha=0.5)+
  facet_grid(cutoff~covars)+
  theme_bw()+
  xlab("fraction of missing data")+
  ylab("p-value correlation (spearman)")
ggsave("figures/imputed_data_p_corS.png")

# sens vs fp
summ_imp2 %>%
  select(run_id, covars, frac_miss, cutoff, contains("num"), contains("true")) %>%
  group_by(run_id) %>%
  mutate(imp_sens=true_imp_overlap/num_true,
         sm_sens=true_sm_overlap/num_true,
         imp_fp=(num_imp-true_imp_overlap)/num_imp,
         sm_fp=(num_sm-true_sm_overlap)/num_sm) %>%
  pivot_longer(c(imp_fp, sm_fp), names_to="source", values_to="fpr") %>%
  mutate(source=ifelse(source=="imp_fp", "imputed", "small")) %>%
  pivot_longer(c(imp_sens, sm_sens), names_to="source1", values_to="sensitivity") %>%
  mutate(source1=ifelse(source1=="imp_sens", "imputed", "small")) %>%
  filter(source1==source) %>%
  ungroup() %>%
  filter(covars=="_age", cutoff==0.05) %>%
  ggplot(aes(x=sensitivity, y=fpr, col=source)) +
  geom_point(alpha=0.5)+
  theme_bw()
ggsave("figures/imp_sens_fpr.png")


