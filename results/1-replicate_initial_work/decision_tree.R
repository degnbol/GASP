#!/usr/bin/env Rscript
library(data.table)
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(here))
# go to git root
setwd(here())

training_data = fread("data/acceptor_interaction_data.txt")
classes_data = training_data[,-(1:22)]
training_data = training_data[,3:22]
training_data[is.na(training_data)] = ""

classes = classes_data[[1]]
valid = classes < 2 & !is.na(classes)
# should actually be ntree=1 for a decision tree
forest = randomForest(training_data[valid,], factor(classes[valid]), ntree=100)

cor(classes[valid], predict(forest, training_data[valid,]) == 1)
cor(classes[valid], predict(forest, training_data[valid,], type="prob")[,"1"])



