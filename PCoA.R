library(ade4)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(plyr)
library(ggalt)
library(pairwiseAdonis)
library(ggpubr)
library(patchwork)

### Read data
ab_table <- read.table("ab_table.txt", sep = "\t", header = T, row.names = 1)  # Taxa abundance table
group <- as.matrix(read.table("group.txt", sep = "\t", row.names = 1, header = T))  # Sample information

### PCoA analysis
pcoa_plot <- pcoa(ab_table, group)

### PCoA function
pcoa <- function(ab_table, group){
  ### Bray-Curtis distance calculation
  dist <- vegdist(ab_table,method='bray') # Other distance calculating methods are available
  data_pcoa <- cmdscale(dist,k=(nrow(ab_table)-1),eig=TRUE) 
  data_pcoa <- cmdscale(dist,k=(nrow(ab_table)-1),eig=TRUE,add=TRUE)
  data_pcoa_eig <- data_pcoa$eig
  data_pcoa_exp <- data_pcoa_eig/sum(data_pcoa_eig)
  pcoa1 <- paste(round(100*data_pcoa_exp[1],2),'%')
  pcoa2 <- paste(round(100*data_pcoa_exp[2],2),'%')
  
  ### Permanova test
  sample_site <- data.frame(data_pcoa$points)[1:2]
  sample_site$level<-factor(group$Group,levels = c(compare_a, compare_b))
  names(sample_site)[1:3] <- c('PCoA1', 'PCoA2','level')
  head(sample_site)
  summary(sample_site)
  
  level_order <- c('compared_variable_1', 'compared_variable_2')
  level_order <- factor(1:length(level_order),labels = level_order)
  sample_site$level <- factor(sample_site$level,levels = levels(level_order))
  find_hull <- function(sample_site) sample_site[chull(sample_site$PCoA1, sample_site$PCoA2),]
  hulls <- ddply(sample_site, "level", find_hull)
  hulls[,1]
  
  PCoA1_mean <- tapply(sample_site$PCoA1,sample_site$level, mean)
  PCoA2_mean <- tapply(sample_site$PCoA2,sample_site$level, mean)
  
  mean_point <- rbind(PCoA1_mean,PCoA2_mean)
  mean_point <- data.frame(mean_point)
  t1 <- t(mean_point)
  t2 <- as.data.frame(t1)
  t2
  pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2, color=level,shape=level)) +
    theme_classic()+
    geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
    geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
    geom_point(size = 4) +
    stat_ellipse(data=sample_site,geom = "polygon",level=0.9,linetype = 2,size=0.5,aes(fill=level),alpha=0.1,show.legend = T)+
    scale_color_manual(values = c('color for variable_1','color for variable_2')) +
    scale_fill_manual(values = c('color for variable_1','color for variable_2')) +
    scale_shape_manual(values = c(20,20)) +
    scale_x_continuous(limits = c(-0.5,0.6),breaks=round(seq(-0.5,0.6,0.1),2)) +
    scale_y_continuous(limits = c(-0.4,0.4),breaks=round(seq(-0.4,0.4,0.1),2)) +
    theme(panel.grid = element_line(color = 'black', linetype = 1, size = 0.1),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.title=element_blank())+
    labs(x = paste('PCoA1_value'), y = paste('PCoA2_value'), title=dune_adonis)
}