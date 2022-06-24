#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("data.table")) # fifelse
library("ggplot2")
library("reshape2")
library("pROC")
suppressPackageStartupMessages(library("here"))
setwd(paste0(here(), "/data/Fatemeh_eval"))

experimental <- read.table('all-experimental.tsv')
gtpredict <- read.table('all-gtpredict.csv')
experimental <- experimental[,c(17,29,26,5,13,31)]
colnames(experimental) <- colnames(gtpredict)
experimental <- experimental[ order(row.names(experimental)), ]
toplot <- melt(as.matrix(experimental))
colnames(toplot) <- c('enzyme', 'substrate', 'experimental')
melted <- melt(as.matrix(gtpredict))
toplot$gtpredict <- as.character(melted$value)
toplot$substrate <- as.character(toplot$substrate)

# prettier naming for plot
toplot$substrate = sub('_', ' ', toplot$substrate)
toplot$gtpredict = fifelse(toplot$gtpredict=="Yes", "reactive", "not reactive")

p = ggplot(toplot, aes(x=enzyme, y=experimental)) +
  geom_point(aes(color=gtpredict, shape=substrate), size=5) +
  theme_linedraw() + 
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        ) +
  xlab("Enzyme") + ylab("Rate") + 
  # colors made with http://medialab.github.io/iwanthue/ selecting two colors with colorblind friendly and hard (force vector)
  scale_color_manual(values=c("#4842b4", "#ff9559")) + 
  # scale_stroke_manual(values=c(1, 1.5)) + 
  scale_shape_manual(values=0:5) + 
  labs(color="GT-predict prediction", shape="Substrate")

insetp = p +
    guides(color="none", shape="none") +
    scale_y_continuous(limits=c(min(toplot$experimental), 0.2), expand=c(0.01,0)) +
    theme(axis.title=element_blank(),
        panel.grid.major.y=element_line(color="gray"),
        panel.grid.minor.y=element_line(color="lightgray"))

p2 = p + annotation_custom(ggplotGrob(insetp), xmin=2.2, xmax=18.8, ymin=0.5, ymax=1.93)
ggsave("GTPredict-evaluation.pdf", p2)

toplot$pred = fifelse(toplot$gtpredict=="reactive", 1, 0)
cat("cor =", cor(toplot$pred, toplot$experimental), end='\n')
cat("auc =", auc(roc(toplot$pred, toplot$experimental)), end='\n')

