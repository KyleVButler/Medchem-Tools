library(dplyr)
library(tibble)
library(stringr)
library(readr)
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
activities_collected <- filter(activities_collected, standard_value > 0)
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
#make sure the compound has at least one activity below 1000 nM in a cell based assay, not admet, and with high confidence on target
cmpdstokeep <- activities_collected %>% filter(standard_value < 1000 & standard_units == "nM" & bao_format == "BAO_0000219" &
                                                 assay_type != "A" & confidence_score == 9) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected)
length(unique(activities_collected$molregno))
#add component_id and accession
target_components <- tbl(chembl_db, "target_components")
join_vector <- target_components %>% select(tid, component_id) %>% filter(tid %in% activities_collected$tid)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "tid")
activities_collected
#join accession
component_sequences <- tbl(chembl_db, "component_sequences")
join_vector <- component_sequences %>% select(accession, component_id) %>% filter(component_id %in% activities_collected$component_id)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "component_id")
activities_collected
activities_collected_temp <- activities_collected
#add description
#62630 compounds now

#it appears that some functional data is labeled as B but under format: cell-based
#chembl can return something like this from assay report card: Format 	Cell-based (BAO_0000219) 
# benchmark compounds
# unc1999 CHEMBL3414619   
# fk506 CHEMBL269732  View(filter(activities_collected, chembl_id == "CHEMBL269732"))
# vx-745 CHEMBL119385

#end part 1
#so how to find the compounds... there will be two groups, those that have selectivity data and those that do not...
#there are also two groups of probes, ones that are selective for a single protein, and those that are selective for 2-3 isoforms first group first
# remember to remove all value <= 1 when subsettimnmg selectivity values
cmpdstokeep <- activities_collected %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                 assay_type == "B" & standard_value > 30) %>% select(molregno)
activities_collected_selectivity <- filter(activities_collected, molregno %in% unique(cmpdstokeep$molregno))
dim(activities_collected_selectivity)
length(unique(activities_collected_selectivity$molregno))
#find compounds that have no selectivity values between 1 and 30
cmpdstokeep <- activities_collected_selectivity %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                             assay_type == "B") %>% 
  filter(standard_value < 30 & standard_value > 1 | standard_value < 1) %>% select(molregno)
activities_collected_selectivity <- filter(activities_collected_selectivity, !(molregno %in% unique(cmpdstokeep$molregno)))
#get minimum potency (nM) values for each molecule/target pair and remove cytochromes  

##
## should also remove herg and pgp
##
test <- activities_collected_selectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                      standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) %>% filter(substr(pref_name, 1, 7) != "Cytochr")

#keep only compounds with 1 activity below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% group_by(molregno) %>% filter(n() == 1) %>% select(molregno)
test <- test %>% filter(molregno %in% n_lower_100$molregno) %>% arrange(molregno, min_value) 
#must have 30 fold selectivity for remaining targets - figure out how to do this with dplyr later
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
test_min <- test %>% group_by(molregno) %>% slice(which.min(min_value)) 
#see if compounds in test_min have a >1000nM cell-based assay for the main target
x <- rep(FALSE, nrow(test_min))
activities_collected_subset <- activities_collected %>% filter(molregno %in% unique(test_min$molregno))
for(i in 1:nrow(test_min)){
  if(nrow(activities_collected_subset %>% filter(molregno == test_min$molregno[i] & standard_units == "nM" & standard_value < 1000 & bao_format == 
                                              "BAO_0000219" & pref_name == test_min$pref_name[i])) >= 1){
    x[i] <- TRUE
  }
}
test_min <- test_min[x, ]
#for above, just keep the top 1 compounds for each target, get number of observations total from activities_collected
cmpdstokeep <- activities_collected %>% filter(molregno %in% test_min$molregno) %>% count(molregno)
test_min <- left_join(test_min, cmpdstokeep, by = "molregno")
test_min <- test_min %>% arrange(pref_name, desc(n), min_value) %>% group_by(pref_name) %>% slice(1)
#keep those with more than 3 observations
potential_probes_group1 <- test_min %>% filter(n > 5)
potential_probes_group1$group <- rep(1, nrow(potential_probes_group1))
#201 compounds


#now to get compounds active against 2 isoforms but have positive selectivity data
cmpdstokeep <- activities_collected %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                 assay_type == "B" & standard_value > 30) %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected_selectivity <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected_selectivity)
length(unique(activities_collected_selectivity$molregno))
#3507 compounds have shown >30 fold selectivity against any other protein
#find compounds that have no more than 1 selectivity values between 1 and 30
cmpdstokeep <- activities_collected_selectivity %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                             assay_type == "B") %>% 
  filter(standard_value < 30 & standard_value > 1 | standard_value < 1) %>% count(molregno) %>% filter(n > 1)
activities_collected_selectivity <- filter(activities_collected_selectivity, !(molregno %in% cmpdstokeep$molregno))
#get min potency (nM) values for each molecule/target pair and remove cytochromes
test <- activities_collected_selectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                      standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) %>% filter(substr(pref_name, 1, 7) != "Cytochr")
#keep only compounds with 2 activities below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% group_by(molregno) %>% filter(n() == 2) %>% select(molregno)
test <- test %>% filter(molregno %in% n_lower_100$molregno) %>% arrange(molregno, min_value) 
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
test <- test %>% filter(molregno %in% mols) %>% select(-namecheck, -namecheck2)
#remove activities that are in the first group of compounds
mols <- test %>% filter(!(pref_name %in% potential_probes_group1$pref_name)) 
test <- filter(test, molregno %in% mols$molregno) 
#see if there is a cell assay for the main targets
x <- rep(FALSE, nrow(test))
activities_collected_subset <- activities_collected %>% filter(molregno %in% unique(test$molregno))
for(i in 1:nrow(test_min)){
  if(nrow(activities_collected_subset %>% filter(molregno == test$molregno[i] & standard_units == "nM" & standard_value < 1000 & bao_format == 
                                          "BAO_0000219" & pref_name == test$pref_name[i])) >= 1){
    x[i] <- TRUE
  }
}
test_min <- test[x, ]
test <- test %>% filter(molregno %in% test_min$molregno)

#match each probe with number of observations
cmpdstokeep <- activities_collected %>% filter(molregno %in% test$molregno) %>% count(molregno)
test <- left_join(test, cmpdstokeep, by = "molregno")
#for each target, keep the probe with the highest number of observations, tiebreak by min_value
mols <- test %>% arrange(pref_name, desc(n), min_value) %>% group_by(pref_name) %>% slice(1)
#keep those probes that appear in mols
test <- test %>% filter(molregno %in% mols$molregno) %>% select(molregno, pref_name, min_value, n) %>% filter(n > 5)

potential_probes_group2 <- test
potential_probes_group2$group <- rep(2, nrow(potential_probes_group2))

#make main list
potential_probes <- potential_probes_group1
potential_probes <- rbind(potential_probes, potential_probes_group2)




#addition of group 3 below - compounds have no selectivity data
cmpdstokeep <- activities_collected %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                 assay_type == "B") %>% select(molregno)
activities_collected_noselectivity <- activities_collected %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))
#get average potency (nM) values for each molecule/target pair 
test <- activities_collected_noselectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                        standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) %>% filter(substr(pref_name, 1, 7) != "Cytochr")
#keep only compounds with 1 activity below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% group_by(molregno) %>% filter(n() == 1) %>% select(molregno)
test <- test %>% filter(molregno %in% n_lower_100$molregno) %>% arrange(molregno, min_value) 
#must have 30 fold selectivity against at least one other target - would two targets be better?
mols <- unique(test$molregno)
x <- rep(FALSE, length(mols))
for(i in 1:length(mols)){
  subset <- test[test$molregno == mols[i], ]
  if(nrow(subset) > 1){
    if(subset[2, 3]/subset[1, 3] > 30){
      x[i] <- TRUE
    }
  }
}
mols <- mols[x]
test <- test %>% filter(molregno %in% mols)
#match each molecule with main target
test_min <- test %>% group_by(molregno) %>% slice(which.min(min_value)) 
#check to see if a cell assay exists for main target
#see if compounds in test_min have a >1000nM cell-based assay for the main target -- this is really slow
x <- rep(FALSE, nrow(test_min))
activities_collected_subset <- activities_collected %>% filter(molregno %in% unique(test_min$molregno))
for(i in 1:nrow(test_min)){
  if(nrow(activities_collected_subset %>% filter(molregno == test_min$molregno[i] & standard_units == "nM" & standard_value < 1000 & bao_format == 
                                          "BAO_0000219" & pref_name == test_min$pref_name[i])) >= 1){
    x[i] <- TRUE
  }
}
test_min <- test_min[x, ]
#for above, just keep the top compound for each target, get number of observations total from activities_collected
cmpdstokeep <- activities_collected %>% filter(molregno %in% test_min$molregno) %>% count(molregno)
test_min <- left_join(test_min, cmpdstokeep, by = "molregno")
test_min <- test_min %>% arrange(pref_name, desc(n), min_value) %>% group_by(pref_name) %>% slice(1)
#how many observations constitutes a good probe? 
potential_probes_group3 <- test_min %>% filter(n > 5)
potential_probes_group3$group <- rep(3, nrow(potential_probes_group3))
potential_probes_group3 <- potential_probes_group3 %>% filter(!(pref_name %in% potential_probes$pref_name))
potential_probes <- rbind(potential_probes, potential_probes_group3)





#should i find more probes that may have two homologous targets?
#addition of group 4 below - compounds have no selectivity data
cmpdstokeep <- activities_collected %>% filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Fold selectivity") & 
                                                 assay_type == "B") %>% select(molregno)
activities_collected_noselectivity <- activities_collected %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#get average potency (nM) values for each molecule/target pair 
test <- activities_collected_noselectivity %>% filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 &
                                                        standard_units == "nM") %>% group_by(molregno, pref_name) %>% 
  summarize(min_value = min(standard_value)) %>% filter(substr(pref_name, 1, 7) != "Cytochr")
#keep only compounds with 2 activity below 100
n_lower_100 <- test %>% filter(min_value < 100) %>% group_by(molregno) %>% filter(n() == 2) %>% select(molregno)
test <- test %>% filter(molregno %in% n_lower_100$molregno) %>% arrange(molregno, min_value) 
#both targets must have 30 fold selectivity against at least two other targets
mols <- unique(test$molregno)
x <- rep(FALSE, length(mols))
for(i in 1:length(mols)){
  subset <- test[test$molregno == mols[i], ]
  if(nrow(subset) >= 4){
    if(subset[3, 3]/subset[1, 3] > 30 & subset[4, 3]/subset[1, 3] > 30 & subset[3, 3]/subset[2, 3] > 30 & 
       subset[4, 3]/subset[2, 3] > 30){
      x[i] <- TRUE
    }
  }
}
mols <- mols[x]
test <- test %>% filter(molregno %in% mols)
#now get the two targets below 100nM and see if they are isoforms (will have similar pref_name, probably same first or last 6 characters)
test <- test %>% filter(min_value < 100)
test$namecheck <- substr(test$pref_name, 1, 6)
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
test <- test %>% filter(molregno %in% mols) %>% select(-namecheck, -namecheck2)
#remove activities that are in the probe list
mols <- test %>% filter(!(pref_name %in% potential_probes$pref_name)) 
test <- filter(test, molregno %in% mols$molregno) 
#see if there is a cell assay for the main targets
x <- rep(FALSE, nrow(test))
activities_collected_subset <- activities_collected %>% filter(molregno %in% unique(test$molregno))
for(i in 1:nrow(test_min)){
  if(nrow(activities_collected_subset %>% filter(molregno == test$molregno[i] & standard_units == "nM" & standard_value < 1000 & bao_format == 
                                          "BAO_0000219" & pref_name == test$pref_name[i])) >= 1){
    x[i] <- TRUE
  }
}
test_min <- test[x, ]
test <- test %>% filter(molregno %in% test_min$molregno)

#match each probe with number of observations
cmpdstokeep <- activities_collected %>% filter(molregno %in% test$molregno) %>% count(molregno)
test <- left_join(test, cmpdstokeep, by = "molregno")
#for each target, keep the probe with the highest number of observations, tiebreak by min_value
mols <- test %>% arrange(pref_name, desc(n), min_value) %>% group_by(pref_name) %>% slice(1)
#keep those probes that appear in mols
potential_probes_group4 <- test %>% filter(molregno %in% mols$molregno) %>% select(molregno, pref_name, min_value, n) %>% filter(n > 5)
potential_probes_group4$group <- rep(4, nrow(potential_probes_group4))
potential_probes <- rbind(potential_probes, potential_probes_group4)


#add accession to potential_probes
#join protein information
target_dictionary <- tbl(chembl_db, "target_dictionary")
join_vector <- target_dictionary %>% select(tid, pref_name, organism) %>% filter(organism == "Homo sapiens") %>%
  filter(pref_name %in% potential_probes$pref_name)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- left_join(potential_probes, join_vector, by = "pref_name")
#join component id
target_components <- tbl(chembl_db, "target_components")
join_vector <- target_components %>% select(tid, component_id) %>% filter(tid %in% potential_probes$tid)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- left_join(potential_probes, join_vector, by = "tid")
#join accession
component_sequences <- tbl(chembl_db, "component_sequences")
join_vector <- component_sequences %>% select(accession, component_id) %>% filter(component_id %in% potential_probes$component_id)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- left_join(potential_probes, join_vector, by = "component_id")


# connect molregno and smiles with compound id
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% potential_probes$molregno)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- full_join(potential_probes, join_vector, by = "molregno")
compound_structures <- tbl(chembl_db, "compound_structures")
join_vector <- compound_structures %>% select(molregno, canonical_smiles) %>% filter(molregno %in% potential_probes$molregno)
join_vector <- collect(join_vector, n = Inf)
potential_probes <- full_join(potential_probes, join_vector, by = "molregno")

#could filter out MW>1000
#PAINS filter

# i copied all into the variable xx
length(unique(potential_probes$molregno))
potential_probes <- filter(potential_probes, !(molregno %in% xx))
length(unique(potential_probes$molregno))



#sgc probes
sgcprobes <- read_csv("accessions.csv")
#connect accessions to pref_name
#remove compounds already covered by sgc probes
cmpdstokeep <- potential_probes %>% filter(!(accession %in% sgcprobes$accession))
potential_probes <- potential_probes %>% filter(molregno %in% cmpdstokeep$molregno)

sgcprobes <- sgcprobes %>% mutate(group = 0)
probe_list <- potential_probes %>% mutate(CHEMICAL_ID = chembl_id) %>% dplyr::select(CHEMICAL_ID, pref_name, accession, group)
probe_list <- bind_rows(list(sgcprobes, probe_list))
probe_list <- dplyr::select(probe_list, -PROTEIN_NAME)
probe_list %>% group_by(group) %>% summarise(Unique_Elements = n_distinct(pref_name)) 
probe_list %>% group_by(group) %>% summarise(Unique_Elements = n_distinct(CHEMICAL_ID)) 

write_csv(probe_list, "PROBELIST.csv")
probe_list <- read_csv("PROBELIST.csv")


probe_list %>% group_by(group) %>% summarise(N_Targets = n_distinct(pref_name)) 
probe_list %>% group_by(group) %>% summarise(N_Probes = n_distinct(CHEMICAL_ID)) 

#different species information causes problems, check CHEMBL2204357
#check if this makes sense:   filter(standard_value < 30 & standard_value > 1 | standard_value < 1)

#make list for protein annotation
protein_accessions <- tibble(accessions = unique(probe_list$accession))

#add protein information to probe list
probe_list <- left_join(probe_list, protein_information, by = c("accession" = "UNIPROTKB"))
