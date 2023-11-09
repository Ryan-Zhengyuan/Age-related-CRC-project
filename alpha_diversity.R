library(vegan)
library(ggplot2)
library(ggsignif)
library(patchwork)

### Read data
CRC_data <- t(read.table("CRC_table_species.txt",sep = "\t", row.names = 1, header = T))
CTR_data <- t(read.table("CTR_table_species.txt",sep = "\t", row.names = 1, header = T))


### Shannon index
alpha_diversity_shannon <- function(a_table, b_table, a, b){
  shannon_a <- data.frame(diversity(a_table, index = 'shannon'))
  shannon_b <- data.frame(diversity(b_table, index = 'shannon'))
  colnames(shannon_a)[1] <- "shannon"
  colnames(shannon_b)[1] <- "shannon"
  compare <- list(c(a, b))
  shannon_a$Group <- a
  shannon_b$Group <- b
  shannon <- rbind(shannon_a, shannon_b)
  diversity_plot <- ggplot(shannon,aes(Group,shannon,fill=Group))+
    geom_violin(width=0.5, trim = FALSE)+
    geom_boxplot(width=0.2, color="white")+
    geom_signif(comparisons = compare,map_signif_level = F,test =wilcox.test,
                test.args=list(alternative = "two.sided", var.equal = FALSE, paired=F))
  return(diversity_plot)
}

shannon <- alpha_diversity_shannon(CRC_data, CTR_data, 'CRC', 'CTR')

### Simpson index
alpha_diversity_simpson <- function(a_table, b_table, a, b){
  simpson_a <- data.frame(diversity(a_table, index = 'simpson'))
  simpson_b <- data.frame(diversity(b_table, index = 'simpson'))
  colnames(simpson_a)[1] <- "simpson"
  colnames(simpson_b)[1] <- "simpson"
  compare <- list(c(a, b))
  simpson_a$Group <- a
  simpson_b$Group <- b
  simpson <- rbind(simpson_a, simpson_b)
  diversity_plot <- ggplot(simpson,aes(Group,simpson,fill=Group))+
    geom_violin(width=0.5, trim = FALSE)+
    geom_boxplot(width=0.2, color="white")+
    geom_signif(comparisons = compare,map_signif_level = F,test =wilcox.test,
                test.args=list(alternative = "two.sided", var.equal = FALSE, paired=F))
  return(diversity_plot)
}

simpson <- alpha_diversity_simpson(CRC_data, CTR_data, 'CRC', 'CTR')

### Plot generation
shannon + simpson  +  plot_layout(design=c(area(1,1), area(1,2)))
