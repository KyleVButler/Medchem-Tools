require(rcdk)
require(Rcpi)
require(caret)
require(kernlab)
require(dplyr)
require(kernlab)
require(reshape2)

setwd("/home/kb/R/Medchem")

source("QSARtools.R")

#This script will predict compound activity using classification 
#We will use the heatmap and nanosyn data from http://www.nature.com/nbt/journal/v34/n1/abs/nbt.3374.html; doi:10.1038/nbt.3374
#Functions are found in QSARtools.R
heatmap <- read.csv("heatmapsmilescmpdid.csv", header = TRUE)
nanosyn <- read.csv("nanosyn.csv", header = TRUE)
nanosyn <- nanosyn[!duplicated(nanosyn$SMILES), ]
#convert activity to A for active or N for not active
#trying to get similar ratios for heatmap and nanosyn, but many more actives are in heatmap
#right now, 5% are active in nanosyn and 11% are active in heatmap
nanosynclass <- select(nanosyn, -c(Compound.ID, Pubmed.ID, SMILES))
nanosynclass <- mutate_each(nanosynclass, funs(ifelse(. >=15, "A", "N")))
heatmapclass <- select(heatmap, -c(Compound.ID, SMILES))
heatmapclass <- mutate_each(heatmapclass, funs(ifelse(. >=2.5, "A", "N")))
table(as.matrix(nanosynclass))
table(as.matrix(heatmapclass))
heatmapclass <- data.frame(SMILES = heatmap$SMILES, heatmapclass)
nanosynclass <- data.frame(SMILES = nanosyn$SMILES, nanosynclass)

set.seed(1)

y <- as.factor(heatmapclass$PLK4)
x <- as.character(heatmapclass$SMILES)

classmodel <- gbmclassmodel(y,x)
plot(classmodel)
#what works? - gbm 2/4 fast, ada 1/4 medium, knn gets 2/4, gets the other chemotype, kknn 2/4 1fp
#extensively tuned gbm gets both chemotypes, but misses one active - best method 

NewSmilesIN <- read.csv("PREDICTPLK4.csv", header = TRUE)
classpredictfromsmiles(x, NewSmilesIN, classmodel)

#setting up the PCM model
#focus on heatmap, class imbalance in nanosyn
#want to remove colnames from nanosyn that appear in heatmap
#amino acids found at atp binding site are extracted from paper in Bioinformatics, 2009, in file protseq
heatmap <- read.csv("heatmapsmilescmpdid.csv", header = TRUE, stringsAsFactors = FALSE)
heatmapclass <- select(heatmap, -c(Compound.ID, SMILES))
heatmapclass <- mutate_each(heatmapclass, funs(ifelse(. >=2.0, "A", "N")))
table(as.matrix(heatmapclass))
heatmapclass <- data.frame(SMILES = heatmap$SMILES, heatmapclass)

SMILES <- as.character(heatmapclass$SMILES)
containers <- parse.smiles(SMILES)
circlefpmatrix <- circularfp(containers)

smileswithfp <- data.frame(SMILES, circlefpmatrix)
nrow(smileswithfp) == length(SMILES)
maindata <- melt(heatmapclass, id.vars = "SMILES", na.rm = TRUE)
maindata <- right_join(maindata, smileswithfp, by="SMILES")

prot_seq <- read.csv("protseq.csv", header = TRUE, stringsAsFactors = FALSE)
blosummatrix <- matrix(0, nrow = 60, ncol = 175)
for(i in 1:nrow(prot_seq)){
  blosummatrix[i, ] <- extractPCMBLOSUM(prot_seq[i, 2], k = 5, lag = 7, scale = TRUE)
}

blosummatrix <- data.frame(variable = prot_seq$Kinase, blosummatrix)
blosummatrix$variable <- as.character(blosummatrix$variable)
maindata$variable <- as.character(maindata$variable)
maindata <- right_join(maindata, blosummatrix, by = "variable")
maindata <- maindata[complete.cases(maindata), ]
maindata$value <- as.factor(maindata$value)
ctrl <- trainControl(method = 'cv', number = 10, repeats = 5, 
                     classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
#model saved
#gbm_model <- train(maindata[, 4:609], maindata$value, method = "gbm", trControl = ctrl, tuneLength = 10,
#                  metric = 'ROC')
roc0 <- roc(maindata$value, predict(gbm_model, maindata[, 4:609], type = "prob")[, 1], levels = rev(levels(maindata$value)))
plot(roc0)
#blosum data for s6k2, our target is in row 60 of blosummatrix, need to add it to all rows of fingerprint matrix - 367 rows
s6k2df <- as.data.frame(blosummatrix[60, -1])
s6k2df <- s6k2df[rep(1, 367), ]
predict_df <- data.frame(x = circlefpmatrix, y = s6k2df)
s6k2_classpredict <- predict(gbm_model, newdata = predict_df)
s6k2_classpredict_probs <- predict(gbm_model, newdata = predict_df, type = "prob")
actives <- data.frame(SMILES = SMILES, activity = s6k2_classpredict, s6k2_classpredict_probs)
actives <- actives[actives$activity == "A", ]
actives
actives_containers <- parse.smiles(as.character(actives$SMILES))
view.molecule.2d(actives_containers)
