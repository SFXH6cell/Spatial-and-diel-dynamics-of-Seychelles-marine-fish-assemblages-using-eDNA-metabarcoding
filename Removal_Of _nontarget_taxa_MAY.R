##################################################################################################
#Script for removal of non-target taxa, mock community species, and control samples
##################################################################################################

library(phyloseq)
install.packages("dplyr")
library(dplyr)
#####Remove Mammalia
###Manually counted 81 Mammalia taxa in taxonomy file that need to be removed

ps_final <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")

###To keep dolphins in the dataset
# Check what mammalian taxa you have in the original object
mammal_tax <- as.data.frame(tax_table(ps_final)) %>% filter(Class == "Mammalia")
table(mammal_tax$Order) ### 3 Artiodactyla and 35 Primates
table(mammal_tax$Family) ### 1 Delphinidae, 35 Hominidae, 2 Suidae


# Remove Mammalia EXCEPT dephinidae
ps_final <- subset_taxa(ps_final, Class != "Mammalia" | Family == "Delphinidae")

dolphin_check <- as.data.frame(tax_table(ps_final)) %>% filter(Family == "Delphinidae")
cat("Dolphin ASVs retained:", nrow(dolphin_check), "\n")
print(dolphin_check)

mammal_check <- as.data.frame(tax_table(ps_final)) %>% filter(Class == "Mammalia")
cat("Remaining mammalian ASVs:", nrow(mammal_check), "\n")
print(mammal_check)

cat("ASVs after removing Mammalia:", ntaxa(ps_final), "\n")

#####ASVs before removing Mammalia: 833
#ASVs after removing Mammalia: 776
#833 - 776 = 57 Mammalia supprimés ==> PAS possible, on est censé avoir 37 


###UPDATE : FIXED 
#####ASVs after removing Mammalia: 794 ==> 1 qui manque je sais pas pourquoi


###investiguer pourquoi : 

table(mammal_tax$Family)
table(mammal_tax$Order)


ps_original <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")
ntaxa(ps_original)


# Compter les ASVs avec Class = NA dans le ps_final original
ps_original <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")
tax_original <- as.data.frame(tax_table(ps_original))
sum(is.na(tax_original$Class))
###20 
### So, R deleted the NA as well, this is why we have 

tax_original[is.na(tax_original$Class), ]


###Encore 1 NA manquant : 
sum(is.na(as.data.frame(tax_table(ps_final))$Class))
table(is.na(as.data.frame(tax_table(ps_final))$Class))

unique(as.data.frame(tax_table(ps_final))$Class)[1:10]
table(as.data.frame(tax_table(ps_final))$Class)

tax_df <- as.data.frame(tax_table(ps_final))
tax_df[tax_df$Class == "unknown", ]


###On va voir la liste de ceux supprimé : 

ps_original <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")
tax_original <- as.data.frame(tax_table(ps_original))

asvs_supprimes <- setdiff(rownames(tax_original), taxa_names(ps_final))
tax_original[asvs_supprimes, ]

######

####C'est BOOOOON J'ai 795 ASV après avoir supprimé mammalia. Donc 38 supprimé (36 hom + 2 sudiae)

##on reprend : 

######## Update : étape pour corriger autres taxons présents dans le dataset 

cat("ASVs after removing Mammalia:", ntaxa(ps_final), "\n")

##### NOUVELLE ÉTAPE : Remove other non-fish contaminant classes
other_nontarget <- as.data.frame(tax_table(ps_final)) %>% 
  filter(Class %in% c("Aves", "Amphibia", "Lepidosauria"))
cat("Non-fish ASVs to remove (Aves/Amphibia/Lepidosauria):", nrow(other_nontarget), "\n")
print(table(other_nontarget$Class))

ps_final <- subset_taxa(ps_final, !Class %in% c("Aves", "Amphibia", "Lepidosauria"))
cat("ASVs after removing other non-fish taxa:", ntaxa(ps_final), "\n")



#List all known mock species in the dataset
known_mock_species <- c("Raja_clavata", 
                        "Mustelus_asterias", 
                        "Chirostoma_humboldtianum", 
                        "Xiphophorus_variatus", 
                        "Pseudoxiphophorus_bimaculatus", 
                        "Lepomis_macrochirus", 
                        "Lepomis_sp.", 
                        "Algansea_tincella", 
                        "Algansea_sp.", 
                        "Micralestes_acutidens", 
                        "Micralestes_sp.", 
                        "Micralestes_sp._NM-2010", 
                        "Oreochromis_amphimelas", 
                        "Oreochromis_niloticus", 
                        "Oreochromis_tanganicae", 
                        "Oreochromis_sp.", 
                        "Cichlidae_sp.", 
                        "Cichliformes_sp.", 
                        "Synodontis_multipunctatus", 
                        "Synodontis_sp.", 
                        "Mochokidae_sp.", 
                        "Chiloglanis_cf._cameronensis")

# Find ASVs matching known mock species anywhere in the dataset
tax_mock <- as.data.frame(tax_table(ps_final))
asvs_to_remove <- rownames(tax_mock)[tax_mock$Species %in% known_mock_species]

cat("Known mock species ASVs to remove:", length(asvs_to_remove), "\n")
#Known mock species ASVs to remove: 5

if (length(asvs_to_remove) > 0) {
  print(tax_mock[asvs_to_remove, ])
  ps_final <- prune_taxa(!taxa_names(ps_final) %in% asvs_to_remove, ps_final)
  cat("ASVs after removing mock species:", ntaxa(ps_final), "\n")
} else {
  cat("No mock species ASVs found — check species names match taxonomy table\n")
}
#ASVs before removing mock species: 795
#ASVs after removing mock species: 790
# 5 removed 


#####Remove negative controls and mocks
###There are 26 known negative controls and 3 known mock communities in the dataset ??

colnames(sample_data(ps_final))
## quasiment sure on peut adapter avec control mais à check
ps_final <- subset_samples(ps_final, !sample_type %in% c("Control", "Mock"))
cat("Samples after removing controls:", nsamples(ps_final), "\n")
#Samples before removing controls: 49
#Samples after removing controls: 42
#7 removed



#remove ASVs now absent after sample removal
ps_final <- prune_taxa(taxa_sums(ps_final) > 0, ps_final)
cat("ASVs after removing controls:", ntaxa(ps_final), "\n")

#ASVs before removing controls: 790
#ASVs after removing controls: 782
# 8 removed 
### ==> 23/06/2026 : Faire - 6 sur chaque ASV parce que on a enlevé 6 autres especes 
cat("Samples after removing controls:", nsamples(ps_final), "\n")
cat("ASVs after removing controls:", ntaxa(ps_final), "\n")




#####Final Summary

cat("\n── Clean dataset ready for analysis ─────────────────\n")
print(ps_final)
cat("\nSamples per treatment:\n")
print(table(sample_data(ps_final)$sample_type))


save(ps_final, file = "ps_seychelles_V1.RData")
saveRDS(ps_final, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final_clean_V2.rds")


##  Changer dans ps_final changer le unknown de la premiere colonne to Eukaryota (1st column)
## And change kingdom with domain 

table(as.data.frame(tax_table(ps_final))$Kingdom)

# Corriger les typos ET renommer en Eukaryota dans Kingdom
tax_mat <- as.matrix(tax_table(ps_final))

tax_mat[tax_mat[, "Kingdom"] %in% c("unknown", "nknown", "unkonwn"), "Kingdom"] <- "Eukaryota"

# Renommer la colonne Kingdom en Domain
colnames(tax_mat)[colnames(tax_mat) == "Kingdom"] <- "Domain"

# Remettre dans ps_final
tax_table(ps_final) <- tax_table(tax_mat)

# Vérifier
table(as.data.frame(tax_table(ps_final))$Domain)



#### and then anything else unknown "incertae_sedis"

colSums(as.data.frame(tax_table(ps_final)) == "unknown", na.rm = TRUE)

tax_mat <- as.matrix(tax_table(ps_final))

# Remplacer unknown par incertae_sedis dans toutes les colonnes sauf Domain
cols_to_fix <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
for (col in cols_to_fix) {
  tax_mat[tax_mat[, col] == "unknown", col] <- "incertae_sedis"
}

# Remettre dans ps_final
tax_table(ps_final) <- tax_table(tax_mat)

# Vérifier
colSums(as.data.frame(tax_table(ps_final)) == "unknown", na.rm = TRUE)


tax_postps <- as.data.frame(tax_table(ps_final))
View(tax_postps)


##### Taxonomic agglomeration at species level
ps_glom <- tax_glom(ps_final, taxrank = "Species")

cat("ASVs before tax_glom:", ntaxa(ps_final), "\n")
cat("Taxa after tax_glom:", ntaxa(ps_glom), "\n")

# Sauvegarder ps_final (nettoyé, avant glom) et ps_glom (après glom)
saveRDS(ps_final, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final_clean_V3.rds")
saveRDS(ps_glom, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_glom_clean_V3.rds")

ps_glom <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_glom_clean_V2.rds")





