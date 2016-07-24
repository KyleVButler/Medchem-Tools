require(rcdk)
require(Rcpi)
require(caret)
require(kernlab)

#setwd("/home/kb/R/Medchem")

source("QSARtools.R")

#This script will predict compound activity (delta T in thermal shift assay in this case)
#We will use the DSF heatmap data from http://www.nature.com/nbt/journal/v34/n1/abs/nbt.3374.html; doi:10.1038/nbt.3374
#Functions are found in QSARtools.R


#can start here from my csv file    
smilesheatmap <- read.csv("heatmapsmilescmpdid.csv", header = TRUE)
#Any target in the heatmapdata will work
target <- "PLK4"
smilesheatmap <- smilesheatmap[complete.cases(smilesheatmap[, target]), ]
#set data under 0.8 to 0
heatmapdata <- smilesheatmap[, 3:ncol(smilesheatmap)]
heatmapdata[heatmapdata < 0.8] <- 0


data.y <- heatmapdata[, target]
#start of input to function
set.seed(1)

data.x <- as.character(smilesheatmap$SMILES)
#getmodelfromfp gives a model and accepts two arguments, a vector of data and a vector of smiles
#Calculates models using various fingerprints ("circular", "standard", "extended") and various learning 
#methods ("gbm", "pls", "svmPoly", "ranger") and returns the model for the pair with the 
#lowest test MSE   -- MAY BE SLOW
bestmodel <- getmodelfromfp(data.y, data.x)
dev.new()
plot(bestmodel)

#some test compounds to predict - the first 4 should be PLK4 inhibitors, the last 5 should not
NewSmilesIN <- read.csv("PREDICTPLK4.csv", header = TRUE)

#this function will get a vector of predicted data from the original training SMILES, new SMILES for predictions, and the model
predictfromsmiles(data.x = data.x, NewSmilesIN = NewSmilesIN, bestmodel = bestmodel)

#Makes a model with random forest and circular fingerprints
bestmodel <- rfcircularmodel(data.y, data.x)
plot(bestmodel)
NewSmilesIN <- read.csv("PREDICTPLK4.csv", header = TRUE)
predictfromsmiles(data.x = data.x, NewSmilesIN = NewSmilesIN, bestmodel = bestmodel)


