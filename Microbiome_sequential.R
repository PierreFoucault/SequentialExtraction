#### Import packages ####

library(BiocManager)
library(phyloseq) 
library(qiime2R) 
library(vegan) 
library(tidyverse)
library(tidytext)
library(ggh4x)
library(patchwork)
library(rstatix)
library(rcompanion)
library(reshape2)
library(foreign)
library(pairwiseAdonis)
library(ggvenn)

# load your session
load("SeqExt.RData")

# Save your session
save.image(file="SeqExt.RData")

#for qiime2R :
#install.packages("devtools") #if you don't have devtools installed
#devtools::install_github("jbisanz/qiime2R")

#pour rcompanion :
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/rcompanion/rcompanion_2.3.26.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#need ‘DescTools’, ‘multcompView’, ‘EMT’, ‘lmtest’, ‘nortest’

####|####

#### Palettes, strip & design ####

# To color the horizontal strips (=site name) according to the 3 levels of eutrophication 

palette_origin<-c("#DF7857","#E7AB9A","#CA955C","#1C82AD",
                  "#40513B","#82954B","#609966","#9DC08B")

palette_origin_blank<-c("#DF7857","#E7AB9A","#CA955C","#1C82AD","black",
                        "#40513B","#82954B","#609966","#9DC08B")

palette_origin_inv<-c("#9DC08B","#609966","#82954B","#40513B",
                      "#1C82AD","#CA955C","#E7AB9A","#DF7857")

palette_phylum<-c("#CCEBC5","#B3CDE3","#DECBE4","#FED9A6","#FDDAEC",
                  "#FFFFCC","#E5D8BD","Gray38","#609966","#DF7857",
                  "#CBD5E8","#F2F2F2","#E6F5C9","#FFF2AE","#F1E2CC")

palette_genus<-c("white","#CCEBC5","#F18B6C","#FDCDAC","#FF8C00","#E5D8BD","#FAC898","#EC5800","#F5761A","#E6812F",
                 "#D74826","#8A3324","#F5761A","Gray38","#DECBE4","#B3CDE3","#CBD5E8","#2E6D74","#B3E2CD","#88F393",
                 "lightgrey")

palette_red<-c("#FF8C00","#FC4C02","#FAC898","#EC5800","#F5761A","#E6812F","#D74826","#8A3324","#F5761A","#F18B6C",
               "#FDCDAC")

pie(rep(1,length(palette_phylum)), col=palette_phylum)

strip_color_origin<- strip_themed(
  background_x = elem_list_rect(fill = palette_origin,color = "black"),
  text_x = elem_list_text(colour = "white", face = "bold",size=9))

strip_color_origin_y<- strip_themed(
  background_y = elem_list_rect(fill = palette_origin,color = "black"),
  text_y = elem_list_text(colour = "white", face = "bold",size=11.5))

strip_color_white_y<- strip_themed(
  background_y = elem_list_rect(fill = "white",color = "white"),
  text_y = elem_list_text(colour = "white", face = "bold",size=11))

strip_color_origin_blank<- strip_themed(
  background_x = elem_list_rect(fill = palette_origin_blank,color = "black"),
  text_x = elem_list_text(colour = "white", face = "bold",size=10))

####|####

#### Import data ####

####__Qiime2 to phyloseq ####
SeqExt<-qza_to_phyloseq(features="Data_16S/filtered_table.qza",
                          taxonomy="Data_16S/taxonomy_DE.qza",
                          metadata="Data_16S/metadata.txt",
                          tree="Data_16S/rooted_tree.qza")
SeqExt@sam_data<-
  read_delim("Data_16S/metadata.txt","\t", escape_double = FALSE,
             trim_ws = TRUE,show_col_types = FALSE) %>%
  remove_rownames %>%
  column_to_rownames(var="sample-id") %>%
  dplyr::mutate(.,
                samples=rownames(.),
                sample_type=factor(sample_type,levels=c("LATIPES","FLUVIATILIS","BIOFILM","WATER","BLANK",
                                              "PMC_728","PMC_807","PMC_810","PMC_826")),
                extraction_method=factor(extraction_method,levels=c("DNA","SEQ","BLANK")),
                sample_names=factor(sample_names,levels=c("LATIPES","FLUVIATILIS","BIOFILM","WATER","BLANK",
                                                          "PMC-728.11","PMC-807.12","PMC-810.12","PMC-826.12"))) %>% sample_data()

comp_reads<-as.data.frame(colSums(SeqExt@otu_table))

####__Remove spurious ASV (lab contamination) ####

## ASVs from the extraction blank
ASV_blank <- SeqExt %>% subset_samples(.,sample_type == "BLANK")

to_remove_row<-NULL
for (i in 1:nrow(ASV_blank@otu_table)) {
  if (sum(ASV_blank@otu_table[i,])==0) {
    to_remove_row<-append(to_remove_row,i)
  }
}
ASV_blank@otu_table<-ASV_blank@otu_table[-(to_remove_row),]
# 116 ASVs

## cyanobacterial ASVs from the PER-E-GUT1 sample
ASV_PER_E_GUT1 <- SeqExt %>% subset_samples(.,samples == "PER-E-GUT1")

to_remove_row<-NULL
for (i in 1:nrow(ASV_PER_E_GUT1@otu_table)) {
  if (sum(ASV_PER_E_GUT1@otu_table[i,])==0) {
    to_remove_row<-append(to_remove_row,i)
  }
}
ASV_PER_E_GUT1@otu_table<-ASV_PER_E_GUT1@otu_table[-(to_remove_row),]
# 94 ASVs

## ASVs intersect from extraction blank and PER-E-GUT1
badTaxa = intersect(rownames(ASV_blank@otu_table),rownames(ASV_PER_E_GUT1@otu_table))

# 1 ASV from lab cultures and extraction blank
# b6fda6888c5dd373460d44bb5474fc2b #ASV_407

## plot its relative abundance across all samples
contamination_barplot<-
  SeqExt %>% transform_sample_counts(., function(x) 100 * x/sum(x)) %>%
  prune_taxa(badTaxa,.)  %>% psmelt(.) %>%
  ggplot(.,aes(x=method_replicate,y=Abundance,fill=Genus))+
  geom_col(width=0.8,color="black",size=0.2)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=15,face="bold"),axis.title.x = element_blank(),
        axis.text.x=element_text(size=10,angle=90),
        axis.ticks =element_blank(),
        legend.position="bottom",
        legend.title = element_text(face="bold",size=15),legend.text = element_text(size=10),
        legend.key.size = unit(0.8,"cm"),legend.direction = "horizontal")+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02),limits = c(NA,100))+
  scale_fill_manual(values = "#CCEBC5",labels="ASV_407 - Microcystis_PCC_7914")+
  guides(fill=guide_legend(ncol=1))+
  labs(fill="",y="Relative Abundance")+
  facet_wrap2(~ sample_names,scale="free_x",nrow = 2,strip = strip_color_origin_blank)
contamination_barplot

## remove the ASV from the dataset expect from cyanobacterial culture samples
list_culture_sample<-
  c("728-E-L1","728-E-L2","728-E-L3","728-ES-L1","728-ES-L2","728-ES-L3",
    "807-E-L1","807-E-L2","807-E-L3","807-ES-L1","807-ES-L2","807-ES-L3",
    "810-E-L1","810-E-L2","810-E-L3","810-ES-L1","810-ES-L2","810-ES-L3",
    "826-E-L1","826-E-L2","826-E-L3","826-ES-L1","826-ES-L2","826-ES-L3")

list_ASV<-rownames(SeqExt@otu_table)

final_otu_table<- SeqExt@otu_table %>% as.tibble(.) %>% t(.) %>%
  as.data.frame(.) %>% `colnames<-`(list_ASV) %>%
  dplyr::mutate(.,
                b6fda6888c5dd373460d44bb5474fc2b=
                  if_else(rownames(.) %in%
                  list_culture_sample,b6fda6888c5dd373460d44bb5474fc2b,0)) %>%
    t(.)

SeqExt_final<-SeqExt
SeqExt_final@otu_table<-final_otu_table %>% otu_table(.,taxa_are_rows = T)

SeqExt_final<-subset_samples(SeqExt_final,sample_type != "BLANK")

####|####

#### Rar Curve ####

# get rarefaction curve data as dataframe
rar_curve.df<-otu_table(SeqExt_final) %>% as.data.frame() %>%
  rarecurve(., step=500, cex=0.5,tidy = T) #View(rar_curve.df)

# Plot the rar. curve data using ggplot2
rar_curve_plot <- rar_curve.df %>%
  ggplot(.,aes(x=Sample,y=Species,group=Site))+
  geom_line()+theme_bw()+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size=15,face="bold",hjust=0.5,color="blue"),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.text=element_text(size=12))+
  #geom_vline(xintercept=raremax_B,color="blue")+
  #annotate('label', x=raremax_B, y=100,label=raremax_B,size=5)+
  scale_x_continuous(expand=c(0,0),limits=c(0,10000))+ 
  scale_y_continuous(expand=c(0,0),limits=c(0,max(rar_curve.df$Species+15)))+
  labs(y = "ASVs Richness\n", x = "\nNb of filtered reads")
rar_curve_plot

####|####

#### Alpha Diversity ####

## data frame
View(as.data.frame(SeqExt@sam_data) %>%
  cbind(.,
        estimate_richness(SeqExt, measures=c("Shannon")),
        estimate_richness(SeqExt, measures=c("Observed"))) %>%
  dplyr::rename("shannon"="Shannon","richness"="Observed") %>%
  dplyr::mutate(.,
                evenness=as.numeric(shannon/log(richness))))

alphaDiv_SeqExt<- as.data.frame(SeqExt_final@sam_data) %>%
  cbind(.,
        estimate_richness(SeqExt_final, measures=c("Shannon")),
        estimate_richness(SeqExt_final, measures=c("Observed"))) %>%
  dplyr::rename("shannon"="Shannon","richness"="Observed") %>%
  dplyr::mutate(.,
                evenness=as.numeric(shannon/log(richness)))

## into .csv
write.csv(alphaDiv_SeqExt,
          "alphaDiv_SeqExt.csv")

####__Shannon ####

# Stats
shanonn.stats<- alphaDiv_SeqExt %>% wilcox_test(shannon ~ sample_type_method,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "sample_type_method", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_shannon<-tibble(sample_type_method=unique(sort(alphaDiv_SeqExt$sample_type_method)),
                       letters=cldList(p.adj ~ as.character(comp),shanonn.stats,threshold  = 0.05)$Letter)
View(shanonn.stats)

View(alphaDiv_SeqExt %>% kruskal_test(shannon ~ extraction_method))
shanonn_type.stats<- alphaDiv_SeqExt %>% wilcox_test(shannon ~ sample_type,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "sample_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

View(evenness_type.stats)

letters_shannon_type<-tibble(sample_type=unique(sort(alphaDiv_SeqExt$sample_type)),
                        letters=cldList(p.adj ~ as.character(comp),shanonn_type.stats,threshold  = 0.05)$Letter)
## Plot

Fig_Shannon_diff<-alphaDiv_SeqExt %>%
  cbind(.,letter=letters_shannon_type$letters[match(.$sample_type,letters_shannon_type$sample_type)]) %>%
  ggplot(.,aes(x=extraction_method,y=shannon))+
  geom_boxplot(aes(fill=sample_type),color="black",show.legend = F)+theme_bw()+
  geom_text(mapping=aes(y=1.1*max(shannon),
                        label=letter),vjust=0,size=5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=12.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12))+
  scale_fill_manual(values =palette_origin)+
  labs(y="Shannon index")+
  facet_wrap2(~sample_names ,scale="free_y",ncol= 1,
              strip = strip_color_white_y,
              strip.position = "left")+coord_flip()
Fig_Shannon_diff

strip_color_white_y<- strip_themed(
  background_y = elem_list_rect(fill = "white",color = "white"),
  text_y = elem_list_text(colour = "white", face = "bold",size=11))

####__Richness ####

## Stats
# richness.stats<- alphaDiv_SeqExt %>% wilcox_test(richness ~ sample_type_method,p.adjust.method = "BH") %>%
#   add_significance("p.adj")  %>% add_xy_position(x = "sample_type_method", dodge = 0.8) %>%
#   mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
#   unite("comp",g1:g2,sep=" - ",remove = FALSE)
# 
# letters_richness<-tibble(sample_type_method=unique(sort(alphaDiv_SeqExt$sample_type_method)),
#                         letters=cldList(p.adj ~ as.character(comp),richness.stats,threshold  = 0.05)$Letter)

richness_type.stats<- alphaDiv_SeqExt %>% wilcox_test(richness ~ sample_type,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "sample_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_richness_type<-tibble(sample_type=unique(sort(alphaDiv_SeqExt$sample_type)),
                             letters=cldList(p.adj ~ as.character(comp),richness_type.stats,threshold  = 0.05)$Letter)

Fig_richness_diff<- alphaDiv_SeqExt %>%
  cbind(.,letter=letters_richness_type$letters[match(.$sample_type,letters_richness_type$sample_type)]) %>%
  ggplot(.,aes(x=extraction_method,y=richness))+
  geom_boxplot(aes(fill=sample_type),color="black",show.legend = F)+theme_bw()+
  geom_text(mapping=aes(y=1.1*max(richness),
                        label=letter),vjust=0,size=5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        axis.title.y=element_blank(),
        axis.text = element_text(size=12.5),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12))+
  scale_fill_manual(values =palette_origin)+
  scale_y_continuous(limits =c(NA,1.15*max(alphaDiv_SeqExt$richness)))+
  labs(y="ASV Richness")+
  facet_wrap2(~sample_names ,scale="free_y",ncol= 1,strip = strip_color_origin_y,strip.position = "left")+coord_flip()
Fig_richness_diff

####__Evenness ####

## Stats
# evenness.stats<- alphaDiv_SeqExt %>% wilcox_test(evenness ~ sample_type_method,p.adjust.method = "BH") %>%
#   add_significance("p.adj")  %>% add_xy_position(x = "sample_type_method", dodge = 0.8) %>%
#   mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
#   unite("comp",g1:g2,sep=" - ",remove = FALSE)
# 
# letters_evenness<-tibble(sample_type_method=unique(sort(alphaDiv_SeqExt$sample_type_method)),
#                          letters=cldList(p.adj ~ as.character(comp),evenness.stats,threshold  = 0.05)$Letter)

evenness_type.stats<- alphaDiv_SeqExt %>% wilcox_test(evenness ~ sample_type,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "sample_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_evenness_type<-tibble(sample_type=unique(sort(alphaDiv_SeqExt$sample_type)),
                              letters=cldList(p.adj ~ as.character(comp),evenness_type.stats,threshold  = 0.05)$Letter)

Fig_evenness_diff<- alphaDiv_SeqExt %>%
  cbind(.,letter=letters_evenness_type$letters[match(.$sample_type,letters_evenness_type$sample_type)]) %>%
  ggplot(.,aes(x=extraction_method,y=evenness))+
  geom_boxplot(aes(fill=sample_type),color="black",show.legend = F)+theme_bw()+
  geom_text(mapping=aes(y=1.1*max(evenness),
                        label=letter),vjust=0,size=5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=12.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12))+
  scale_fill_manual(values =palette_origin)+
  labs(y="Evenness")+
  facet_wrap2(~sample_names ,scale="free_y",ncol= 1,strip = strip_color_white_y,strip.position = "left")+coord_flip()
Fig_evenness_diff

####|####

# #### Venn diagram ####
# 
# ####__ SEQ ####
# # subset the full dataset
# SeqExt_final_SEQ <- SeqExt_final %>% subset_samples(.,method == "SEQ") 
# 
# ## generate 4 sets of ASVs
# SEQ_venn.list = list(GUT     = rownames(subset_samples(SeqExt_final_SEQ,origin %in% c("MEDAKA","PERCA"))@otu_table %>%
#                                         as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      BIOFILM = rownames(subset_samples(SeqExt_final_SEQ,origin %in% c("BIOFILM"))@otu_table %>%
#                                         as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      WATER   = rownames(subset_samples(SeqExt_final_MDK,method == "DNA")@otu_table %>%
#                                           as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      CULTURE = rownames(subset_samples(SeqExt_final_SEQ,origin %in% c("PMC_728","PMC_807",
#                                                                                       "PMC_810","PMC_826"))@otu_table %>%
#                                            as.data.frame(.) %>% .[which(rowSums(.) != 0),]))
#                      
# ## plot
# SEQ_venn.plot <-
#   ggvenn(SEQ_venn.list,
#          fill_color = c("#DF7857","#CA955C","#1C82AD","#609966"),stroke_size = 0.5,set_name_size = 4,show_percentage = F)
# SEQ_venn.plot
#
# ####__ DNA ####
# # subset the full dataset
# SeqExt_final_DNA <- SeqExt_final %>% subset_samples(.,method == "DNA") 
# 
# ## generate 4 sets of ASVs
# DNA_venn.list = list(GUT     = rownames(subset_samples(SeqExt_final_DNA,origin %in% c("MEDAKA","PERCA"))@otu_table %>%
#                                           as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      BIOFILM = rownames(subset_samples(SeqExt_final_DNA,origin %in% c("BIOFILM"))@otu_table %>%
#                                           as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      WATER   = rownames(subset_samples(SeqExt_final_DNA,method == "DNA")@otu_table %>%
#                                           as.data.frame(.) %>% .[which(rowSums(.) != 0),]),
#                      CULTURE = rownames(subset_samples(SeqExt_final_DNA,origin %in% c("PMC_728","PMC_807",
#                                                                                       "PMC_810","PMC_826"))@otu_table %>%
#                                           as.data.frame(.) %>% .[which(rowSums(.) != 0),]))
# 
# ## plot
# DNA_venn.plot <-
#   ggvenn(DNA_venn.list,
#          fill_color = c("#DF7857","#CA955C","#1C82AD","#609966"),stroke_size = 0.5,set_name_size = 4,show_percentage = F)
# DNA_venn.plot

####|####

#### Composition ####

####__ Phylum rank ####
barplot_replicate<- SeqExt_final %>% transform_sample_counts(., function(x) 100 * x/sum(x)) %>%
  tax_glom(.,taxrank ="Phylum") %>% psmelt(.) %>%
  rename(rel_abundance=Abundance) %>%
  dplyr::mutate(.,
                Phylum_legend=ifelse(rel_abundance<1,"Phylum < 1%",Phylum),
                method_replicate=ifelse(method_replicate =="DNA_1","D1",
                                        ifelse(method_replicate =="DNA_2","D2",
                                               ifelse(method_replicate =="DNA_3","D3",
                                                      ifelse(method_replicate =="SEQ_1","S1",
                                                             ifelse(method_replicate =="SEQ_2","S2","S3")))))) %>%
  ggplot(.,aes(x=method_replicate,y=rel_abundance,
               fill=fct_rev(fct_reorder(Phylum_legend,rel_abundance))))+
  geom_col(width=0.8,color="black",size=0.2)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=12),
        axis.ticks =element_blank(),
        legend.position="right",
        legend.title = element_text(face="bold",size=15),legend.text = element_text(size=12),
        legend.key.size = unit(0.8,"cm"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02))+
  scale_fill_manual(values=c("#CCEBC5","#B3CDE3","#DECBE4","#FED9A6","#FDDAEC",
                             "#E5D8BD","#FFFFCC","Gray38","#609966","#DF7857",
                             "#E6F5C9","#CBD5E8","#F2F2F2"))+
  guides(fill=guide_legend(ncol=1))+
  labs(fill="Microbial Phyla",y="Relative Abundance")+
  facet_wrap2(~ sample_names,scale="fixed",nrow = 2,strip = strip_color_origin)
barplot_replicate


strip_color_origin<- strip_themed(
  background_x = elem_list_rect(fill = palette_origin,color = "black"),
  text_x = elem_list_text(colour = "white", face = "bold",size=11))

####__ Genus rank ####
Top_Genera <- SeqExt_final %>% tax_glom(.,taxrank ="Genus") %>% transform_sample_counts(., function(x) 100 * x/sum(x))
Top_Genera <- names(sort(taxa_sums(Top_Genera), TRUE)[1:20])

barplot_Genus_replicate.df<- SeqExt_final %>% 
  tax_glom(.,taxrank ="Genus") %>%
  transform_sample_counts(., function(x) 100 * x/sum(x)) %>%
  psmelt(.) %>%
  rename(rel_abundance=Abundance) %>%
  dplyr::mutate(.,
                Genus_legend=ifelse(OTU %in% Top_Genera,paste0(Phylum,"_",Genus),"Other genera")) %>%
  dplyr::group_by(sample_names,method_replicate,Sample,Genus_legend) %>%
  dplyr::summarise(rel_abun_sum=sum(rel_abundance))

barplot_Genus_replicate <- barplot_Genus_replicate.df %>%
  ggplot(.,aes(x=method_replicate, y=rel_abun_sum,
               fill=fct_rev(fct_reorder(Genus_legend,rel_abun_sum))))+
  geom_col(width=0.8,color="black",size=0.2)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=15,face="bold"),axis.title.x = element_blank(),
        axis.text.x=element_text(size=10,angle=90),
        axis.ticks =element_blank(),
        legend.position="right",
        legend.title = element_text(face="bold",size=15),legend.text = element_text(size=10),
        legend.key.size = unit(0.8,"cm"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02))+
  scale_fill_manual(values=palette_genus)+
  guides(fill=guide_legend(ncol=1))+
  labs(fill="Microbial Genera",y="Relative Abundance")+
  facet_wrap2(~ sample_names,scale="free",nrow = 2,strip = strip_color_origin)
barplot_Genus_replicate

####|####

#### Beta Diversity ####

####__Bray-Curtis ####

## PCoA
SeqExt.BC=ordinate(SeqExt_final, "PCoA")

Fig_SeqExt.BC<-
  plot_ordination(SeqExt_final,
                  SeqExt.BC,shape="extraction_method",color="sample_names")+
  geom_point(size=3.5)+theme_bw()+
  theme(panel.grid = element_blank(),#aspect.ratio = 1,
        plot.title = element_text(hjust=0,size=15,face = "italic"),
        axis.title = element_text(size=15,face = "bold"),
        axis.text=element_text(size=12), axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_color_manual(values=palette_origin)+
  scale_shape_manual(labels=c('DNA', 'SEQ'),values = c(16,17))+
  annotate("text", label = "Bray-Curtis dissimilarity",x = -0.25, y = 0.5, size = 5, colour = "black",fontface="italic")+
  labs(shape="Extraction\nmethod",color="Sample Type",
       x=paste0("PCoA Axis 1 [",round(SeqExt.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("PCoA Axis 2 [",round(SeqExt.BC$values$Relative_eig[2],3)*100,"%]"))
Fig_SeqExt.BC

####____ Sample type ####
#merge the factor to the metabolites table
SeqExt_final.JC.dist=phyloseq::distance(SeqExt_final,method="jaccard",binary=T)
SeqExt_final.BC.dist=phyloseq::distance(SeqExt_final,method = "bray")
ASV_factor<-  data.frame(t(SeqExt_final@otu_table)) %>%
  dplyr::mutate(.,sample_type=Fig_SeqExt.BC$data$sample_type[match(rownames(.),Fig_SeqExt.BC$data$samples)])
#View(ASV_factor)
####______ permanova ####
permanova_ASV<-adonis2(SeqExt_final.JC.dist~ ASV_factor$sample_type,data=ASV_factor,permutations=999,method="bray")
permanova_ASV$`Pr(>F)`[[1]] #General P-value (if ANY sig. diff.)
permanova_ASV
####__regression DNA vs SEQ ####

Ext_SeqExt.df<- SeqExt_final %>% phyloseq::distance(.,method = "bray") %>% as.matrix() %>%
  as.data.frame() %>% rownames_to_column(.,var="sample1_BC") %>%
  pivot_longer(-sample1_BC,names_to = "sample2_BC",values_to = "BC_dist") %>%
  dplyr::mutate(.,
                origin1 = SeqExt_final@sam_data$sample_names[match(.$sample1_BC,SeqExt_final@sam_data$samples)],
                origin2 = SeqExt_final@sam_data$sample_names[match(.$sample2_BC,SeqExt_final@sam_data$samples)],
                method1 = SeqExt_final@sam_data$extraction_method[match(.$sample1_BC,SeqExt_final@sam_data$samples)],
                method2 = SeqExt_final@sam_data$extraction_method[match(.$sample2_BC,SeqExt_final@sam_data$samples)],
                origin_rep1 = SeqExt_final@sam_data$sample_type_replicate[match(.$sample1_BC,SeqExt_final@sam_data$samples)],
                origin_rep2 = SeqExt_final@sam_data$sample_type_replicate[match(.$sample2_BC,SeqExt_final@sam_data$samples)],
                comp_sample=paste0(sample1_BC,"vs",sample2_BC),
                comp_origin_rep = paste0(origin_rep1,"vs",origin_rep2)) %>%
  subset(sample1_BC != sample2_BC) %>%  subset(origin1 == origin2) %>% subset(method1 == method2) #View(Ext_SeqExt.df)

SeqExt.df <- Ext_SeqExt.df %>% subset(method1 == "SEQ")

Fig_comp_BC<-Ext_SeqExt.df %>% subset(method1 == "DNA") %>%
  dplyr::rename("BC_dist_DNA"="BC_dist","comp_DNA"="comp_sample") %>%
  cbind(
    comp_SEQ=SeqExt.df$comp_sample[match(.$comp_origin_rep,SeqExt.df$comp_origin_rep)],
    BC_dist_SEQ=SeqExt.df$BC_dist[match(.$comp_origin_rep,SeqExt.df$comp_origin_rep)]) %>%
  dplyr::filter(duplicated(BC_dist_DNA) == FALSE) %>%
  dplyr::mutate(.,diff_BC=abs(BC_dist_DNA-BC_dist_SEQ)) %>%
  ggplot(.,aes(BC_dist_DNA,BC_dist_SEQ,color=origin1,group=1))+theme_classic()+
  geom_abline(slope=1, intercept=0,color="black",linetype = "dashed",size=1.15)+
  geom_smooth(method = "lm",color="#C21807",show.legend = F,se = T,alpha=0.2)+
  geom_point(size=4,shape=16,show.legend = F)+
  theme(panel.grid = element_blank(),#aspect.ratio = 1,
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=12),
        axis.ticks = element_blank(),
        plot.caption = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12))+
  scale_color_manual(values =palette_origin)+
  scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  labs(x="DNA extraction dissimilarity",
       y="SEQ extraction dissimilarity",color="Sample Type")

spearman_DNA_SEQ<-cor(Fig_comp_BC$data$BC_dist_SEQ,Fig_comp_BC$data$BC_dist_DNA,method = "spearman")
lm_BC_model <- lm(BC_dist_SEQ~BC_dist_DNA,data=Fig_comp_BC$data)

var_BC=expression(paste("Spearman ",round(spearman_DNA_SEQ,2),R^2,
                        round(summary(lm_BC_model)[["adj.r.squared"]],2),italic("p"),"<0.01"))
View(summary(lm_BC_model))


Fig_comp_BC<-Fig_comp_BC+geom_label(label=expression(paste("Spearman 0.62, ",R^2," 0.76, ",italic("p"),"<0.01")),
                                    aes(x=0.3,y=0.8),fill="white",color="#C21807",size=4.2)
Fig_comp_BC

# Fig_diff_BC<- Fig_comp_BC$data %>%
#   ggplot(.,aes(fct_rev(origin1),diff_BC))+
#   geom_boxplot(aes(fill=origin1),color="black",show.legend = F)+theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title.x = element_text(size=15,face = "bold"),
#         axis.title.y=element_blank(),
#         axis.text.y = element_text(size=12,face="bold",color=palette_origin_inv),
#         axis.text.x = element_text(size=12),
#         axis.ticks = element_blank(),
#         legend.title = element_text(size=15,face = "bold"),
#         legend.text = element_text(size=12))+
#   scale_y_continuous(limits = c(NA,max(Fig_comp_BC$data$BC_dist_SEQ)))+
#   scale_fill_manual(values =palette_origin)+
#   labs(y="SeqExt vs DNA extraction\nBray-Curtis dissimilarities")+coord_flip()
# Fig_diff_BC

####__W.Unifrac ####
SeqExt.WU=ordinate(SeqExt_final, "PCoA","unifrac",weighted=T)

Fig_SeqExt.WU<-
  plot_ordination(SeqExt_final,
                  SeqExt.WU,shape="extraction_method",color="sample_names")+
  geom_point(size=3.5)+theme_bw()+
  theme(aspect.ratio = 1,panel.grid = element_blank(),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=12),axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12),
        plot.title = element_text(hjust=0,size=15,face = "italic"))+
  scale_color_manual(values=palette_origin)+
  scale_shape_manual(labels=c('DNA', 'SEQ'),values = c(16,17))+
  labs(shape="Extraction\nmethod",color="Sample Type",title="Weighted Unifrac dissimilarity",
       x=paste0("PCoA Axis 1 [",round(SeqExt.WU$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("PCoA Axis 2 [",round(SeqExt.WU$values$Relative_eig[2],3)*100,"%]"))
Fig_SeqExt.WU

####__Jaccard ####
SeqExt.jaccard=ordinate(SeqExt_final, "PCoA","jaccard",binary=T)

Fig_SeqExt.jaccard<-
  plot_ordination(SeqExt_final,
                  SeqExt.jaccard,shape="extraction_method",color="sample_names")+
  geom_point(size=3.5)+theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size=15,face="italic"),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=12),axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_color_manual(values=palette_origin)+
  scale_shape_manual(labels=c('DNA', 'SEQ'),values = c(16,17))+
  annotate("text", label = "Jaccard dissimilarity",x = -0.2, y = 0.27, size = 5, colour = "black",fontface="italic")+
  labs(shape="Extraction\nmethod",color="Sample Type",
       x=paste0("PCoA Axis 1 [",round(SeqExt.jaccard$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("PCoA Axis 2 [",round(SeqExt.jaccard$values$Relative_eig[2],3)*100,"%]"))
Fig_SeqExt.jaccard

####|####

#### Figures ####

####__ F1: Alpha-div ####

design_F1="
ABCDD
ABCEE"

Fig_SeqExt.jaccard_final<-
  Fig_SeqExt.jaccard+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(color=guide_legend(override.aes=list(size=6),nrow=4))

Fig_richness_diff+Fig_evenness_diff+Fig_Shannon_diff+Fig_SeqExt.jaccard_final+guide_area()+
  plot_layout(design = design_F1,guides = 'collect')+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold",size=17))

####__ F2: BC + Barplot ####

design_F2="
AACC
BBCC"

Fig_SeqExt.BC_final<-
  Fig_SeqExt.BC+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(color=guide_legend(override.aes=list(size=6),nrow=8))

Fig_comp_BC_final<-Fig_comp_BC+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(color=guide_legend(override.aes=list(size=6),nrow=4))

barplot_replicate_final<-
  barplot_replicate+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(fill=guide_legend(override.aes=list(size=6),ncol=1))

Fig_BC_barplot<-Fig_SeqExt.BC_final+Fig_comp_BC_final+barplot_replicate_final+
  plot_layout(design = design_F2,guides = 'collect')+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold",size=17),
        legend.position = 'right',legend.direction = "vertical")
Fig_BC_barplot


####______________________________________________________________________####