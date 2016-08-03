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


#Building AKT activity in to BRAF inhibitors
akt_data <- read.csv("AKT Bioactivity.csv", header = TRUE)
akt_data <- select(akt_data, CMPD_CHEMBLID, CANONICAL_SMILES, STANDARD_VALUE, STANDARD_UNITS)
akt_data <- filter(akt_data, STANDARD_UNITS == "nM")
akt_data <- mutate(akt_data, pIC50 = -log10(STANDARD_VALUE * 10^(-9)))
akt_data$CANONICAL_SMILES <- as.character(akt_data$CANONICAL_SMILES)
akt_gbm_model <- gbmcircularregmodel(akt_data$pIC50, akt_data$CANONICAL_SMILES)


#predict AKT activity for potent BRAF inhibitors with <100nM IC50
braf_data <- read.csv("braf bioactivity.csv", header = TRUE)
braf_data <- select(braf_data, CMPD_CHEMBLID, CANONICAL_SMILES, STANDARD_VALUE, STANDARD_UNITS)
braf_data <- filter(braf_data, STANDARD_UNITS == "nM")
braf_data <- mutate(braf_data, pIC50 = -log10(STANDARD_VALUE * 10^(-9)))
braf_data$CANONICAL_SMILES <- as.character(braf_data$CANONICAL_SMILES)
braf_data <- filter(braf_data, pIC50 > 8)
braf_predict <- gbmcircularregpredict(data.x = akt_data$CANONICAL_SMILES, NewSmilesIN = braf_data$CANONICAL_SMILES, 
                                      bestmodel = akt_gbm_model)
braf_predict <- data.frame(SMILES = braf_data$CANONICAL_SMILES, CHEMBLID = braf_data$CMPD_CHEMBLID, 
                           predictedAKTpIC50 = braf_predict)
braf_predict
braf_akt_predicted_inhibitors <- braf_predict[braf_predict$predictedAKTpIC50 >= 7, ]
braf_akt_predicted_containers <- parse.smiles(as.character(braf_akt_predicted_inhibitors$SMILES))
view.molecule.2d(braf_akt_predicted_containers)
# predicts CHEMBL388978(staurosporine) and CHEMBL570873

#build BRAF model

#predict AKT activity for potent BRAF inhibitors with <100nM IC50
braf_data <- read.csv("braf bioactivity.csv", header = TRUE)
braf_data <- select(braf_data, CMPD_CHEMBLID, CANONICAL_SMILES, STANDARD_VALUE, STANDARD_UNITS)
braf_data <- filter(braf_data, STANDARD_UNITS == "nM")
braf_data <- mutate(braf_data, pIC50 = -log10(STANDARD_VALUE * 10^(-9)))
braf_data$CANONICAL_SMILES <- as.character(braf_data$CANONICAL_SMILES)
braf_gbm_model <- gbmcircularregmodel(braf_data$pIC50, braf_data$CANONICAL_SMILES)

#load new molecules for predictions
newmols <- load.molecules(molfile = "dualinhibs3.sdf", aromaticity = TRUE)
view.molecule.2d(newmols)
newmols_smiles <- rep("", length(newmols))
for(i in 1:length(newmols)){newmols_smiles[i] <- get.smiles(newmols[[i]])}
akt_predict <- gbmcircularregpredict(data.x = akt_data$CANONICAL_SMILES, NewSmilesIN = newmols_smiles, 
                                      bestmodel = akt_gbm_model)
braf_predict <- gbmcircularregpredict(data.x = braf_data$CANONICAL_SMILES, NewSmilesIN = newmols_smiles, 
                                      bestmodel = braf_gbm_model)
braf_akt_predict <- data.frame(SMILES = newmols_smiles, AKT_pIC50 = akt_predict, 
                           BRAF_pIC50 = braf_predict)
braf_akt_predict
braf_akt_predict[braf_akt_predict$AKT_pIC50 == max(braf_akt_predict$AKT_pIC50), ]
