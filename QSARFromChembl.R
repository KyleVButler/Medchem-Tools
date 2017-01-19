#function to build qsar model from chembl data for any protein
#requires chembl_22.db 
#This function uses circular fingerprints from rcdk package and 
#the xgbtree learning method.
#Returns a model with the columns of the fingerprint matrix
#removed because of near-zero variance in attribute "LostColumns".
#Two arguments are required - path_to_chembl (e.g. "/home/kb/Desktop/chembl/chembl_22.db"),
#and the target protein used to build the model in uniprot format, e.g. "P14416"

QSARFromChembl <- function(target, path_to_chembl){
library(dplyr)
library(caret)
library(xgboost)
require(rcdk)
require(Rcpi)


#link to chembl

chembl_db <- src_sqlite(path_to_chembl, create = TRUE)
#choose target


#establish chembl data frames
activities <- tbl(chembl_db, "activities")
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
assays <- tbl(chembl_db, "assays")
component_sequences <- tbl(chembl_db, "component_sequences")
target_components <- tbl(chembl_db, "target_components")
compound_structures <- tbl(chembl_db, "compound_structures")

#get bioactivities for target
component <- component_sequences %>% filter(accession == target)
component <- collect(component, n = Inf)
tid <- target_components %>% filter(component_id == component$component_id)
tid <- collect(tid, n = Inf)
assay_ids <- assays %>% filter(tid %in% tid$tid & confidence_score > 7)
assay_ids <- collect(assay_ids , n = Inf)
bioactivities <- activities %>% filter(assay_id %in% assay_ids$assay_id & standard_units == "nM" & standard_value > 0) %>%
  select(molregno, standard_value)
bioactivities <- collect(bioactivities, n = Inf)
bioactivities <- bioactivities %>% filter(standard_value > 0) %>% mutate(pIC50 = -log10(standard_value * 10^(-9)))
molecules <- compound_structures %>% select(molregno, canonical_smiles) %>% filter(molregno %in% bioactivities$molregno)
molecules <- collect(molecules, n = Inf)
bioactivities <- left_join(bioactivities, molecules, by = "molregno")

#only use the maximum bioactivity for each molecule, then remove those containers that are NA
bioactivities <- bioactivities %>% group_by(molregno) %>% arrange(desc(pIC50)) %>% slice(1)
containersTrain <- parse.smiles(bioactivities$canonical_smiles)
bioactivities <- bioactivities[!is.na(containersTrain),]
containersTrain <- containersTrain[!is.na(containersTrain)]


circularFP <- function(cmpdlist){
  fptest <- get.fingerprint(cmpdlist[[1]], type = "circular", fp.mode = 'bit')
  fpmatrix <- matrix(data = 0, nrow = length(cmpdlist), ncol = attr(fptest, "nbit"))
  for(i in 1:length(cmpdlist)){
    fp <- get.fingerprint(cmpdlist[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
    print(i)
  }
  lostColumns <- nearZeroVar(fpmatrix)
  fpmatrix <- fpmatrix[, -nearZeroVar(fpmatrix)]
  attr(fpmatrix, "LostColumns") <- lostColumns
  
  return(fpmatrix)
}

fingerprints <- circularFP(containersTrain)
nrow(fingerprints) == length(containersTrain)
str(fingerprints)
ctrl <- trainControl(method = 'cv', number = 5, repeats = 3, 
                     summaryFunction = defaultSummary)

xgboost.fit1 <- train(fingerprints, bioactivities$pIC50, method = "xgbTree", trControl = ctrl, tuneLength = 10,
                  metric = 'RMSE')
attr(xgboost.fit1, "fingerprintMethod") <- "circular"
attr(xgboost.fit1, "LostColumns") <- attr(fingerprints, "LostColumns")
return(xgboost.fit1)
}
