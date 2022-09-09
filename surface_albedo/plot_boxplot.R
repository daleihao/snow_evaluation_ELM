library(tidyverse)
library(R.matlab)
library(cowplot)
library(reshape2)
setwd('C:/Users/haod776/OneDrive - PNNL/Documents/work/E3SM/writting/Snow_evaluation/writing/figure/code/surface_albedo/')

# import data
all_datas <- readMat('Winter.mat')
all_data <- as.data.frame(all_datas['AllData'])
colnames(all_data) <- c('SnowLevel','SnowFrac', 'ELM', 'MCD43')
all_data$SnowLevel = as.factor(all_data$SnowLevel)
all_data <- melt(all_data,id = c("SnowLevel","SnowFrac"))

plot_winter <- ggplot(all_data, mapping = aes(x = SnowLevel, y = value, fill = variable)) + xlab(expression(italic(f)[sno])) + 
  ylab(expression(italic(alpha)[sno])) +  ggtitle('Winter') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.5, show.legend = FALSE) +
  geom_boxplot(notch=F, size = 0.5, show.legend = FALSE, width = 0.6, outlier.shape = NA) + 
  #geom_jitter( stroke = 0, alpha = .5, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F,position = position_dodge(width = 0.6)) +  
  #scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('0-0.2','0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0')) + 
  scale_fill_manual(values=c('#4DBEEE', '#EDB120')) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=10, angle=0,hjust = 0.5),
        axis.title.x=element_text( color="black", size=10, angle=0),
        axis.text.x = element_text(color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_text(face="bold", color="black",size=12),
        axis.line = element_line(size = 1)) +
  ylim(c(0,0.8))


# import data
all_datas <- readMat('Spring.mat')
all_data <- as.data.frame(all_datas['AllData'])
colnames(all_data) <- c('SnowLevel','SnowFrac', 'ELM', 'MCD43')
all_data$SnowLevel = as.factor(all_data$SnowLevel)
all_data <- melt(all_data,id = c("SnowLevel","SnowFrac"))

plot_spring <- ggplot(all_data, mapping = aes(x = SnowLevel, y = value, fill = variable)) + xlab(expression(italic(f)[sno])) +  ggtitle('Spring') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.5, show.legend = FALSE) +
  geom_boxplot(notch=F, size = 0.5, show.legend = TRUE, width = 0.6, outlier.shape = NA) + 
 # geom_jitter( stroke = 0, alpha = .5, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F,position = position_dodge(width = 0.6)) +  
  #scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('0-0.2','0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0')) + 
  scale_fill_manual(values=c('#4DBEEE', '#EDB120')) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=10, angle=0,hjust = 0.5),
        axis.title.x=element_text( color="black", size=10, angle=0),
        axis.text.x = element_text(color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        legend.title = element_blank()) +
  ylim(c(0,0.8))


plot_figure <- plot_grid(   plot_winter,plot_spring,
                            labels = c("(a)", "(b)"),
                            label_size = 12,
                            nrow = 1, ncol = 2, align = 'v')

ggsave("../../figure_all_tif/surfacealbedo_different_snowcover.tiff", plot = plot_figure, width = 20, height = 7, units = "cm", dpi = 300, limitsize = FALSE, compression = "lzw")