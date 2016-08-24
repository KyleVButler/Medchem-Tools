library(dplyr)
library(data.table)
chembl_db <- src_sqlite("/Users/kylebutler/Desktop/chembl_21_sqlite/chembl_21.db", create = TRUE)
#type .tables in sqlite3 to see tables
activities <- tbl(chembl_db, "activities")
#collect only compounds with some ic50 below 100 nM
activities <- activities %>% filter(standard_units == "nM" & standard_value < 100) %>% select(molregno)
cmpdstokeep <- collect(activities, n = Inf)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities <- tbl(chembl_db, "activities")
activities_collected <- activities %>% filter(molregno %in% cmpdstokeep) %>% select(activity_id, assay_id, molregno, standard_value,
                                                                                    standard_units, standard_type, data_validity_comment)
activities_collected <- collect(activities_collected, n = Inf)
length(unique(cmpdstokeep)); length(unique(activities_collected$molregno))
dim(activities_collected)

# so i don't have to repeat above
activities_collected_temp <- activities_collected

#remove bad data
activities_collected <- activities_collected[is.na(activities_collected$data_validity_comment), ]
dim(activities_collected)
#append compound ids and assay ids for checking
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
cmpdstokeep <- unique(activities_collected$molregno)
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "molregno")
dim(activities_collected)
filter(activities_collected, molregno == 1837801) 
assays <- tbl(chembl_db, "assays")
cmpdstokeep <- unique(activities_collected$assay_id)
join_vector <- assays %>% select(assay_id, assay_type, tid, bao_format, confidence_score, description) %>% filter(assay_id %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "assay_id")
dim(activities_collected)
filter(activities_collected, molregno == 1837801) 
target_dictionary <- tbl(chembl_db, "target_dictionary")
cmpdstokeep <- unique(activities_collected$tid)
join_vector <- target_dictionary %>% select(tid, target_type, pref_name, organism) %>%
  filter(tid %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "tid")
dim(activities_collected)
filter(activities_collected, molregno == 1837801) 
#all is good until here
activities_collected <- activities_collected %>% filter(organism == "Homo sapiens" | is.na(organism))

#get only compounds with sub 100nM assay for a single protein human targets only
cmpdstokeep <- activities_collected %>% filter(standard_value < 100 & standard_units == "nM" & confidence_score == 9 & 
                                                 organism == "Homo sapiens" & target_type == "SINGLE PROTEIN") %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected)
length(unique(activities_collected$molregno))
#make sure the compound has at least one activity below 1000 nM in a cell based activity
cmpdstokeep <- activities_collected %>% filter(standard_value < 1000 & standard_units == "nM" & bao_format == "BAO_0000219") %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected)
length(unique(activities_collected$molregno))
#add description
#62630 compounds now

#it appears that some functional data is labeled as B but under format: cell-based
#chembl can return something like this from assay report card: Format 	Cell-based (BAO_0000219) 
# benchmark compounds
# unc1999 CHEMBL3414619   
# fk506 CHEMBL269732  View(filter(activities_collected, chembl_id == "CHEMBL269732"))
# vx-745 CHEMBL119385

#end part 1
#maybe focus data on Homo sapiens and NA
#so how to find the compounds... there will be two groups, those that have selectivity data and those that do not...
#there are also two groups of probes, ones that are selective for a single protein, and those that are selective for 2-3 isoforms first group first
# remember to remove all value <= 1 when subsettimnmg selectivity values
cmpdstokeep <- activities_collected %>% filter(standard_type == "Selectivity ratio" | standard_type == "Ratio IC50" | 
                                                 standard_type == "Fold selectivity") %>% filter(standard_value > 30) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected_selectivity <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected_selectivity)
length(unique(activities_collected_selectivity$molregno))
#3507 compounds have shown >30 fold selectivity against any other protein
#find compounds that have no selectivity values between 1 and 30
cmpdstokeep <- activities_collected_selectivity %>% filter(standard_type == "Selectivity ratio" | standard_type == "Ratio IC50" | 
                                                             standard_type == "Fold selectivity") %>% 
  filter(standard_value < 30 & standard_value >= 1) %>% count(molregno) 
cmpdstokeep <- cmpdstokeep %>% filter(n >= 1) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected_selectivity <- filter(activities_collected_selectivity, !(molregno %in% cmpdstokeep))
#get average potency (nM) values for each molecule/target pair 
test <- activities_collected_selectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                      standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) 
test$substringtoremove <- substr(test$pref_name, 1, 7)
#remove cytochromes
test <- test %>% filter(!(substringtoremove == "Cytochr")) %>% select(-substringtoremove)
#keep only compounds with 1 activity below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% count(molregno) 
n_lower_100 <- n_lower_100[n_lower_100$n == 1, ]$molregno
test <- test %>% filter(molregno %in% n_lower_100) %>% arrange(molregno, min_value) 
#must have 30 fold selectivity for remaining targets
mols <- unique(test$molregno)
x <- rep(TRUE, length(mols))
for(i in 1:length(mols)){
  subset <- test[test$molregno == mols[i], ]
  if(nrow(subset) > 1){
    if(subset[2, 3]/subset[1, 3] < 30){
     x[i] <- FALSE
    }
  }
}
mols <- mols[x]
test <- test %>% filter(molregno %in% mols)
#match each molecule with main target
test_min <- data_frame(molregno = rep(NA, length(mols)), pref_name = rep(NA, length(mols)), min_value = rep(NA, length(mols)))
for(i in 1:length(mols)){
  subset <- test[test$molregno == mols[i], ]
  test_min[i, ] <- subset[which.min(subset$min_value),]
}
#for above, just keep the top 2 compounds for each target, get number of observations total from activities_collected
cmpdstokeep <- activities_collected %>% filter(molregno %in% test_min$molregno) %>% count(molregno)
test_min <- full_join(test_min, cmpdstokeep, by = "molregno")
test_min <- test_min %>% arrange(pref_name, min_value) %>% group_by(pref_name) %>% top_n(1, n) %>% group_by(pref_name) %>% top_n(-1, min_value)
potential_probes_group1 <- test_min
potential_probes_group1$group <- rep(1, nrow(potential_probes_group1))
#210 targets, 253 compounds
#



#now to get compounds active against 2 isoforms
cmpdstokeep <- activities_collected %>% filter(standard_type == "Selectivity ratio" | standard_type == "Ratio IC50" | 
                                                 standard_type == "Fold selectivity") %>% filter(standard_value > 30) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected_selectivity <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected_selectivity)
length(unique(activities_collected_selectivity$molregno))
#3507 compounds have shown >30 fold selectivity against any other protein
#find compounds that have no selectivity values between 1 and 30
cmpdstokeep <- activities_collected_selectivity %>% filter(standard_type == "Selectivity ratio" | standard_type == "Ratio IC50" | 
                                                             standard_type == "Fold selectivity") %>% 
  filter(standard_value < 30 & standard_value >= 1) %>% count(molregno) 
cmpdstokeep <- cmpdstokeep %>% filter(n > 1) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected_selectivity <- filter(activities_collected_selectivity, !(molregno %in% cmpdstokeep))
#get average potency (nM) values for each molecule/target pair 
test <- activities_collected_selectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                      standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) 
test$substringtoremove <- substr(test$pref_name, 1, 7)
#remove cytochromes
test <- test %>% filter(!(substringtoremove == "Cytochr")) %>% select(-substringtoremove)
#keep only compounds with 2 activities below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% count(molregno) 
n_lower_100 <- n_lower_100[n_lower_100$n == 2, ]$molregno
test <- test %>% filter(molregno %in% n_lower_100) %>% arrange(molregno, min_value) 
#must have 30 fold selectivity for any targets not in the two below 100nM
mols <- unique(test$molregno)
x <- rep(TRUE, length(mols))
for(i in 1:length(mols)){
  subset <- test[test$molregno == mols[i], ]
  if(nrow(subset) > 2){
    if(subset[3, 3]/subset[1, 3] < 30){
      x[i] <- FALSE
    }
  }
}
mols <- mols[x]
test <- test %>% filter(molregno %in% mols)
#now get the two targets below 100nM and see if they are isoforms (will have similar pref_name, probably same first or last 6 characters)
test <- test %>% filter(min_value < 100)
test$namecheck <- substr(test$pref_name, 1, 6)
library(stringr)
test$namecheck2 <- str_sub(test$pref_name,-6,-1)
mols <- unique(test$molregno)
x <- rep(FALSE, length(mols))
for(i in seq_along(mols)){
  subset <- test[test$molregno == mols[i], ]
  if(subset[1, 4] == subset[2, 4] | subset[1, 5] == subset[2, 5]){
    x[i] <- TRUE
  }
}
mols <- mols[x]
test <- test %>% filter(molregno %in% mols)
#match each molecule with main target
test_min <- data_frame(molregno = rep(NA, length(mols)), pref_name = rep(NA, length(mols)), min_value = rep(NA, length(mols)))
for(i in seq_along(mols)){
  subset <- test[test$molregno == mols[i], ]
  test_min[i, ] <- subset[which.min(subset$min_value), 1:3]
}
#for above, just keep the top compoundfor each target, get number of observations total from activities_collected
cmpdstokeep <- activities_collected %>% filter(molregno %in% test_min$molregno) %>% count(molregno)
test_min <- full_join(test_min, cmpdstokeep, by = "molregno")
#for same counts, tiebreak on lowest potency
test_min <- test_min %>% arrange(pref_name, min_value) %>% group_by(pref_name) %>% top_n(1, n) %>% group_by(pref_name) %>% top_n(-1, min_value)
#get both active targets for any compounds in above list
test_min <- filter(test, molregno %in% test_min$molregno)
test_min <- test_min[, 1:3]
cmpdstokeep <- activities_collected %>% filter(molregno %in% test_min$molregno) %>% count(molregno)
test_min <- full_join(test_min, cmpdstokeep, by = "molregno")
potential_probes_group2 <- test_min
potential_probes_group2$group <- rep(2, nrow(potential_probes_group2))
#remove probes from group 2 that are already in group 1
dim(potential_probes_group2)
cmpdstokeep <- potential_probes_group2 %>% filter(!(pref_name %in% potential_probes_group1$pref_name)) %>% select(molregno)
potential_probes_group2 <- filter(potential_probes_group2, molregno %in% cmpdstokeep$molregno)
dim(potential_probes_group2)
#make main list
potential_probes <- potential_probes_group1
potential_probes <- rbind(potential_probes, potential_probes_group2)
#pains filter - returns smiles in different format so you have to manually copy the offenders - can i use schrodinger canvas?
# connect molregno and smiles with compound id
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% potential_probes$molregno)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- full_join(potential_probes, join_vector, by = "molregno")
compound_structures <- tbl(chembl_db, "compound_structures")
join_vector <- compound_structures %>% select(molregno, canonical_smiles) %>% filter(molregno %in% potential_probes$molregno)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- full_join(potential_probes, join_vector, by = "molregno")
library(ChemmineR)
library(ChemmineOB)
smiset <- as(potential_probes$canonical_smiles, "SMIset") 
cid(smiset) <- potential_probes$molregno
write.SMI(smiset, file="sub.smi", cid=TRUE) 
# i copied all into the variable xx
length(unique(potential_probes$molregno))
potential_probes <- filter(potential_probes, !(molregno %in% xx))
length(unique(potential_probes$molregno))
length(unique(xx))
