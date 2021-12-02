#!/usr/bin/env Rscript
library("ggplot2")
library("reshape2")
library(pROC)
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "data/Fatemeh_eval"))

experimental <- read.table('all-experimental.csv')
gtpredict <- read.table('all-gtpredict.csv')
experimental <- experimental[,c(17,29,26,5,13,31)]
colnames(experimental) <- colnames(gtpredict)
experimental <- experimental[ order(row.names(experimental)), ]
toplot <- melt(as.matrix(experimental))
colnames(toplot) <- c('enzyme', 'substrate', 'experimental')
melted <- melt(as.matrix(gtpredict))
toplot$gtpredict <- as.character(melted$value)
toplot$substrate <- as.character(toplot$substrate)

# pdf("GTPredict-evaluation.pdf")
p <- ggplot(toplot, aes(x=enzyme, y=experimental))+
  geom_point(aes(color=gtpredict,shape=substrate,size=5)) +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position='right') +
  theme(axis.text.x = element_text(size=14, angle=90), axis.text.y = element_text(size=14),
        axis.title = element_text(size=18, face="bold")) +
  ylim(-0.03,0.4) + 
  xlab("Enzyme") + ylab("Rate") + 
  guides(size = FALSE) + 
  scale_color_manual(values=c('#5ab4ac','#8c510a')) +
  scale_shape_manual(values=seq(0,5)) + 
  labs(color="GT-predict predictions", shape="Substrate")
print(p)
# dev.off()


toplot$pred = fifelse(toplot$gtpredict=="Yes", 1, 0)
cor(toplot$pred, toplot$experimental)
auc(roc(toplot$pred, toplot$experimental))

