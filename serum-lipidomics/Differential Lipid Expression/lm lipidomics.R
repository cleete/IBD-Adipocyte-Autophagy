library(reshape2)
library(tibble)
library(tidyr)
library(emmeans)

view(lipid.data)

#Making a long table with two columns only
lipid.data.new.long = rbind(lipid.data.new,metadata) %>% 
  drop_na() %>%
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('id') %>% 
  melt(id.vars=c('id','genotype','sex','condition'),variable.name='lipid') %>%
  mutate(lipid=droplevels(lipid),value=as.numeric(value))

results = data.frame()
for (cond in unique(lipid.data.new.long$condition)) {
  for (lip in unique(lipid.data.new.long$lipid)) {
    data = lipid.data.new.long %>%
      filter(lipid==lip,condition==cond)
    fit=lm(log2(value)~genotype+sex,data=data) #maybe ignore interaction effect? (can probably be ignored since, interaction effects show no significant differences)
    results = rbind(
      results,
      broom::tidy(fit) %>%
        mutate(condition=cond,lipid=lip)
    )
  }
}

results %>% 
  filter(term=='genotypeWT', condition=='DSS') %>% view()


results2 = data.frame()
for (geno in unique(lipid.data.new.long$genotype)) {
  for (lip in unique(lipid.data.new.long$lipid)) {
    data = lipid.data.new.long %>%
      filter(lipid==lip,genotype==geno)
    fit=lm(log2(value)~condition*sex,data=data)
    results2 = rbind(
      results2,
      broom::tidy(fit) %>%
        mutate(genotype=geno,lipid=lip)
    )
  }
}



write.csv(results, "Results DEL.csv")
view(results)

#Vulcano plotting differences between 
ggplot(results %>% filter(term=='genotypeWT', condition=='DSS')) + #%>% mutate(p.value=p.adjust(p.value,method='fdr'))) +
  geom_point(aes(x=estimate,y=-log10(p.value))) +
  geom_hline(yintercept=-log10(0.05)) +
  xlab('Log2 Fold Change')+
  scale_x_continuous(limits = c(-2, 2))+
  scale_y_continuous(limits = c(0,5))+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  geom_vline(xintercept = c(-0.5,0.5))

#Beide p-value darstellen in einem plot (p.adj and p.)
ggplot(results %>% filter(term=='genotypeWT', condition=='DSS') %>% mutate(p.adj=p.adjust(p.value,method='fdr'))) +
  geom_point(aes(x=estimate,y=-log10(p.value))) +
  geom_point(aes(x=estimate,y=-log10(p.adj)), color='red') +
  geom_hline(yintercept=-log10(0.05)) +
  xlab('Log2 Fold Change')+
  scale_x_continuous(limits = c(-2, 2))+
  scale_y_continuous(limits = c(0,5))+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15), text = element_text(size=20),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  geom_vline(xintercept = c(-0.5,0.5))

ggplot(results %>% filter(term=='genotypeWT', condition=='NI') %>% mutate(p.adj=p.adjust(p.value,method='fdr'))) +
  geom_point(aes(x=estimate,y=-log10(p.value))) +
  geom_point(aes(x=estimate,y=-log10(p.adj)), color='red') +
  geom_hline(yintercept=-log10(0.05)) +
  xlab('Log2 Fold Change')+
  scale_x_continuous(limits = c(-2, 2))+
  scale_y_continuous(limits = c(0,5))+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15), text = element_text(size=20),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  geom_vline(xintercept = c(-0.5,0.5))

#Vulcano plotting differences between NI vs DSS
ggplot(results2 %>% filter(term=='conditionNI', genotype=='WT') %>% mutate(p.value=p.adjust(p.value,method='fdr'))) +
  geom_point(aes(x=estimate,y=-log10(p.value))) +
  geom_hline(yintercept=-log10(0.05)) +
  xlab('Log2 Fold Change')+
  scale_x_continuous(limits = c(-2, 2))+
  scale_y_continuous(limits = c(0,5))+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  geom_vline(xintercept = c(-0.5,0.5))



emmeans(fit,c('genotype','sex')) %>% #estimated marginal means (emmeans)
  contrast('pairwise') %>%
  summary(infer=T)


lipid.data.long = rbind(lipid.data.ni,metadata) %>% 
  drop_na() %>%
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('id') %>% 
  melt(id.vars=c('id','genotype','sex','condition'),variable.name='lipid') %>%
  mutate(lipid=droplevels(lipid),value=as.numeric(value))
