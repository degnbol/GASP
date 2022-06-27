#!/usr/bin/env Rscript
library(data.table)
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(here))
# go to git root
setwd(here())

training_data = fread("data/GT-Predict/acceptor_interaction_data.txt")
classes_data = training_data[,-(1:22)]
training_data = training_data[,3:22]
training_data[is.na(training_data)] = ""

classes = classes_data[[1]]
valid = classes < 2 & !is.na(classes)
# ntree=1 for a decision tree
tree = randomForest(training_data[valid,], factor(classes[valid]), ntree=1)
forest = randomForest(training_data[valid,], factor(classes[valid]), ntree=100)

# This is a circular test, trainset == testset, but it shows that a single tree 
# is clearly worse, it can't even capture this artificially easy prediction well.
cor(classes[valid], predict(tree, training_data[valid,]) == 1)
cor(classes[valid], predict(tree, training_data[valid,], type="prob")[,"1"])
cor(classes[valid], predict(forest, training_data[valid,]) == 1)
cor(classes[valid], predict(forest, training_data[valid,], type="prob")[,"1"])


