#!/usr/bin/env Rscript
# there has been openMP (multi-process) issues with data.table on Mac. 
install.packages("data.table", type="source", repos="https://Rdatatable.gitlab.io/data.table")
install.packages(c("webchem", "ggplot2", "Matrix", "tidyverse", "BiocManager", "optparse"), repos="https://cloud.r-project.org/")
install.packages("rentrez")
# R version of git root
install.packages("here")
# for data/Fatemeh_eval
install.packages(c("reshape2", "pROC"))
