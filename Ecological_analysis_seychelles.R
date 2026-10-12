####Ecological Analysis Seychelles###
head(sample_data(ps_final))


head(sample_data(ps_final)[, c("Island", "Time")])
table(sample_data(ps_final)$Island)
table(sample_data(ps_final)$Time)


# Nombre total d'ASVs
cat("Total ASVs:", ntaxa(ps_final), "\n")

####
cat("Total ASVs:", ntaxa(ps_glom), "\n")
###
##376, donc pb ca veut dire y a pas identité distcincte 

# Nombre total de samples
cat("Total samples:", nsamples(ps_final), "\n")

# Extraire la table taxonomique
tax_summary <- as.data.frame(tax_table(ps_glom))

# Nombre de familles uniques
cat("Number of families:", length(unique(tax_summary$Family)), "\n")

# Nombre de genres uniques
cat("Number of genera:", length(unique(tax_summary$Genus)), "\n")

# Nombre d'espèces uniques
cat("Number of species:", length(unique(tax_summary$Species)), "\n")


unique(tax_summary$Species)


### Séparer par île
ps_felicite <- subset_samples(ps_glom, Island == "Félicité")
ps_fregate <- subset_samples(ps_glom, Island == "Frégate")

# Retirer les ASVs absents dans chaque sous-dataset
ps_felicite <- prune_taxa(taxa_sums(ps_felicite) > 0, ps_felicite)
ps_fregate <- prune_taxa(taxa_sums(ps_fregate) > 0, ps_fregate)

# Résumé Félicité
tax_fe <- as.data.frame(tax_table(ps_felicite))
cat("=== Félicité ===\n")
cat("Species:", length(unique(tax_fe$Species)), "\n")
cat("Genera:", length(unique(tax_fe$Genus)), "\n")
cat("Families:", length(unique(tax_fe$Family)), "\n")
cat("ASVs:", ntaxa(ps_felicite), "\n")


# Résumé Frégate
tax_fr <- as.data.frame(tax_table(ps_fregate))
cat("\n=== Frégate ===\n")
cat("Species:", length(unique(tax_fr$Species)), "\n")
cat("Genera:", length(unique(tax_fr$Genus)), "\n")
cat("Families:", length(unique(tax_fr$Family)), "\n")
cat("ASVs:", ntaxa(ps_fregate), "\n")


########
# Détail par famille
cat("\nASVs per Family:\n")
print(sort(table(tax_summary$Family), decreasing = TRUE))

# Détail par genre
cat("\nASVs per Genus:\n")
print(sort(table(tax_summary$Genus), decreasing = TRUE))

# Détail par espèce
cat("\nASVs per Species:\n")
print(sort(table(tax_summary$Species), decreasing = TRUE))


####regarder combien d'espèce sûre on a : 
# Espèces incertaines (_sp.)
cat("Espèces incertaines (_sp.):", sum(grepl("_sp\\.", unique(tax_summary$Species))), "\n")

# Espèces certaines
cat("Espèces certaines:", sum(!grepl("_sp\\.", unique(tax_summary$Species))), "\n")

# Total
cat("Total espèces uniques:", length(unique(tax_summary$Species)), "\n")


# Liste espèces certaines
sort(unique(tax_summary$Species)[!grepl("_sp\\.", unique(tax_summary$Species))])

# Liste espèces incertaines
sort(unique(tax_summary$Species)[grepl("_sp\\.", unique(tax_summary$Species))])


# Genres incertains (_gen.)
cat("Genres incertains (_gen.):", sum(grepl("_gen\\.", unique(tax_summary$Genus))), "\n")

# Genres certains
cat("Genres certains:", sum(!grepl("_gen\\.", unique(tax_summary$Genus))), "\n")

# Familles incertaines (_fam.)
cat("Familles incertaines (_fam.):", sum(grepl("_fam\\.", unique(tax_summary$Family))), "\n")

# Familles certaines
cat("Familles certaines:", sum(!grepl("_fam\\.", unique(tax_summary$Family))), "\n")




###Avant d'utiliser Taxglom, on veut que chaque _sp. reste une identité distincte, on va donc les numeroter


tax_mat <- as.matrix(tax_table(ps_final))

# Pour chaque niveau taxonomique avec des _sp. ou _gen. ou _fam.
for (col in c("Species", "Genus", "Family")) {
  # Trouver les entrées avec _sp., _gen., _fam.
  uncertain <- grep("_sp\\.|_gen\\.|_fam\\.", tax_mat[, col])
  
  if (length(uncertain) > 0) {
    # Compter les occurrences de chaque valeur incertaine
    vals <- tax_mat[uncertain, col]
    # Numéroter les doublons
    tax_mat[uncertain, col] <- ave(vals, vals, FUN = function(x) {
      if (length(x) == 1) x
      else paste0(x, seq_along(x))
    })
  }
}

# Remettre dans ps_final 
tax_table(ps_final) <- tax_table(tax_mat)

# Vérifier
head(sort(unique(as.data.frame(tax_table(ps_final))$Species)), 20)

####Corriger le point
tax_mat <- gsub("_sp\\.\\.([0-9])", "_sp._\\1", tax_mat)
tax_mat <- gsub("_gen\\.\\.([0-9])", "_gen._\\1", tax_mat)
tax_mat <- gsub("_fam\\.\\.([0-9])", "_fam._\\1", tax_mat)

tax_table(ps_final) <- tax_table(tax_mat)

# Vérifier
sort(grep("_sp\\._", unique(as.data.frame(tax_table(ps_final))$Species), value = TRUE))

 # Voir les espèces avec _sp. numérotées
sort(grep("_sp\\.", unique(as.data.frame(tax_table(ps_final))$Species), value = TRUE))

###Maintenant fusionner ave tax_Glom

ps_glom <- tax_glom(ps_final, taxrank = "Species")

# Vérifier
cat("ASVs avant tax_glom:", ntaxa(ps_final), "\n")
cat("ASVs après tax_glom:", ntaxa(ps_glom), "\n")

####REPRISE 

library(phyloseq)
list.files("~/Desktop/Research_project/Research_data_seychelles/All_data", pattern = "\\.rds")

library(phyloseq)
library(dplyr)

ps_final <- readRDS("~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")

# Vérifier
ps_final
ntaxa(ps_final)
nsamples(ps_final)


##### POst tax glom Descirption dataset


tax_glom_df <- as.data.frame(tax_table(ps_glom))

cat("Total taxa:", ntaxa(ps_glom), "\n")
cat("Species:", length(unique(tax_glom_df$Species)), "\n")
cat("Genera:", length(unique(tax_glom_df$Genus)), "\n")
cat("Families:", length(unique(tax_glom_df$Family)), "\n")


r# Créer les sous-objets par île
ps_glom_felicite <- subset_samples(ps_glom, Island == "Félicité")
ps_glom_felicite <- prune_taxa(taxa_sums(ps_glom_felicite) > 0, ps_glom_felicite)

ps_glom_fregate <- subset_samples(ps_glom, Island == "Frégate")
ps_glom_fregate <- prune_taxa(taxa_sums(ps_glom_fregate) > 0, ps_glom_fregate)

tax_fe <- as.data.frame(tax_table(ps_glom_felicite))
tax_fr <- as.data.frame(tax_table(ps_glom_fregate))

cat("=== Félicité ===\n")
cat("Total taxa:", ntaxa(ps_glom_felicite), "\n")
cat("Species:", length(unique(tax_fe$Species)), "\n")
cat("Genera:", length(unique(tax_fe$Genus)), "\n")
cat("Families:", length(unique(tax_fe$Family)), "\n\n")

cat("=== Frégate ===\n")
cat("Total taxa:", ntaxa(ps_glom_fregate), "\n")
cat("Species:", length(unique(tax_fr$Species)), "\n")
cat("Genera:", length(unique(tax_fr$Genus)), "\n")
cat("Families:", length(unique(tax_fr$Family)), "\n")



tax_glom_df <- as.data.frame(tax_table(ps_glom))

cat("Species certain:", sum(!grepl("_sp\\.", unique(tax_glom_df$Species))), "\n")
cat("Species uncertain (_sp.):", sum(grepl("_sp\\.", unique(tax_glom_df$Species))), "\n\n")

cat("Genera certain:", sum(!grepl("_gen\\.", unique(tax_glom_df$Genus))), "\n")
cat("Genera uncertain (_gen.):", sum(grepl("_gen\\.", unique(tax_glom_df$Genus))), "\n\n")

cat("Families certain:", sum(!grepl("_fam\\.", unique(tax_glom_df$Family))), "\n")
cat("Families uncertain (_fam.):", sum(grepl("_fam\\.", unique(tax_glom_df$Family))), "\n")


##Per island : 

# Félicité
cat("=== Félicité ===\n")
cat("Species uncertain (_sp.):", sum(grepl("_sp\\.", unique(tax_fe$Species))), 
    "/", length(unique(tax_fe$Species)), 
    "(", round(100 * sum(grepl("_sp\\.", unique(tax_fe$Species))) / length(unique(tax_fe$Species)), 1), "%)\n")

cat("Genera uncertain (_gen.):", sum(grepl("_gen\\.", unique(tax_fe$Genus))), 
    "/", length(unique(tax_fe$Genus)), 
    "(", round(100 * sum(grepl("_gen\\.", unique(tax_fe$Genus))) / length(unique(tax_fe$Genus)), 1), "%)\n")

cat("Families uncertain (_fam.):", sum(grepl("_fam\\.", unique(tax_fe$Family))), 
    "/", length(unique(tax_fe$Family)), 
    "(", round(100 * sum(grepl("_fam\\.", unique(tax_fe$Family))) / length(unique(tax_fe$Family)), 1), "%)\n\n")

# Frégate
cat("=== Frégate ===\n")
cat("Species uncertain (_sp.):", sum(grepl("_sp\\.", unique(tax_fr$Species))), 
    "/", length(unique(tax_fr$Species)), 
    "(", round(100 * sum(grepl("_sp\\.", unique(tax_fr$Species))) / length(unique(tax_fr$Species)), 1), "%)\n")

cat("Genera uncertain (_gen.):", sum(grepl("_gen\\.", unique(tax_fr$Genus))), 
    "/", length(unique(tax_fr$Genus)), 
    "(", round(100 * sum(grepl("_gen\\.", unique(tax_fr$Genus))) / length(unique(tax_fr$Genus)), 1), "%)\n")

cat("Families uncertain (_fam.):", sum(grepl("_fam\\.", unique(tax_fr$Family))), 
    "/", length(unique(tax_fr$Family)), 
    "(", round(100 * sum(grepl("_fam\\.", unique(tax_fr$Family))) / length(unique(tax_fr$Family)), 1), "%)\n")



###Taxon unique

# Espèces dans Félicité mais pas dans Frégate
species_felicite_only <- setdiff(unique(tax_fe$Species), unique(tax_fr$Species))
cat("Species unique to Félicité:", length(species_felicite_only), "\n")
sort(species_felicite_only)

# Espèces dans Frégate mais pas dans Félicité
species_fregate_only <- setdiff(unique(tax_fr$Species), unique(tax_fe$Species))
cat("\nSpecies unique to Frégate:", length(species_fregate_only), "\n")
sort(species_fregate_only)










###requin et raies : 


tax_glom_df <- as.data.frame(tax_table(ps_glom))

# Chercher dans toutes les colonnes taxonomiques
shark_ray <- tax_glom_df[grepl("Chondrichthyes", tax_glom_df$Class, ignore.case = TRUE), ]

cat("Number of shark/ray taxa (Chondrichthyes):", nrow(shark_ray), "\n\n")
print(shark_ray)


unique(tax_glom_df$Class[tax_glom_df$Class != "Actinopteri"])
tax_glom_df <- as.data.frame(tax_table(ps_glom))

shark_ray <- tax_glom_df[grepl("Chondrichthyes", tax_glom_df$Class, ignore.case = TRUE), ]

unique(shark_ray$Species)


table(shark_ray$Order)


# Vérifier toutes les incohérences potentielles Class/Order dans tout le dataset
table(tax_glom_df$Class, tax_glom_df$Order)

##### le dataset a été nettoyé des incoherences 

###Maintenant pour le ray et shark, il y a toujours le probleme du poisson lanterne
tax_glom_df <- as.data.frame(tax_table(ps_glom))

tax_glom_df[tax_glom_df$Order == "Myctophiformes", ]

#### Remplacer la mauvaise classification par Actinopteri 
tax_mat <- as.matrix(tax_table(ps_glom))

tax_mat[tax_mat[, "Order"] == "Myctophiformes", "Class"] <- "Actinopteri"

tax_table(ps_glom) <- tax_table(tax_mat)

# Vérifier
as.data.frame(tax_table(ps_glom))[as.data.frame(tax_table(ps_glom))$Order == "Myctophiformes", ]

##poisson Raies et shark : 
  
tax_glom_df <- as.data.frame(tax_table(ps_glom))

cat("Actinopteri:", sum(tax_glom_df$Class == "Actinopteri"), "\n")


# Chondrichthyes : séparer requins vs raies par Order
chond <- tax_glom_df[tax_glom_df$Class == "Chondrichthyes", ]

# Requins (sharks) : Carcharhiniformes et autres ordres de requins
sharks <- chond[chond$Order %in% c("Carcharhiniformes"), ]
cat("Sharks:", nrow(sharks), "\n")

# Raies (rays) : Myliobatiformes, Rhinopristiformes
rays <- chond[chond$Order %in% c("Myliobatiformes", "Rhinopristiformes"), ]
cat("Rays:", nrow(rays), "\n")

dolphins <- tax_glom_df[tax_glom_df$Family == "Delphinidae", ]
cat("Dolphins:", nrow(dolphins), "\n")

# Identifier le taxon manquant
table(tax_glom_df$Class)

###Checker incoherence : 
tax_glom_df <- as.data.frame(tax_table(ps_glom))

# Pour chaque Order, vérifier s'il y a plus d'une Class associée
inconsistent_orders <- tax_glom_df %>%
  group_by(Order) %>%
  summarise(n_classes = n_distinct(Class), classes = paste(unique(Class), collapse = ", ")) %>%
  filter(n_classes > 1)


print(inconsistent_orders)



###Corriger incoherence

tax_glom_df[tax_glom_df$Order == "Carcharhiniformes", ]
tax_glom_df[tax_glom_df$Order == "Myliobatiformes", ] 

tax_mat <- as.matrix(tax_table(ps_glom))

tax_mat[tax_mat[, "Order"] %in% c("Carcharhiniformes", "Myliobatiformes"), "Class"] <- "Chondrichthyes"

tax_table(ps_glom) <- tax_table(tax_mat)

# Vérifier qu'il n'y a plus d'incohérence
tax_glom_df <- as.data.frame(tax_table(ps_glom))
table(tax_glom_df$Order[tax_glom_df$Order %in% c("Carcharhiniformes", "Myliobatiformes")], 
      tax_glom_df$Class[tax_glom_df$Order %in% c("Carcharhiniformes", "Myliobatiformes")])

table(tax_summary)


tax_mat <- as.matrix(tax_table(ps_glom))
tax_mat[tax_mat[, "Order"] %in% c("Carcharhiniformes", "Myliobatiformes"), "Class"] <- "Chondrichthyes"
tax_table(ps_glom) <- tax_table(tax_mat)

tax_glom_df <- as.data.frame(tax_table(ps_glom))

inconsistent_orders <- tax_glom_df %>%
  group_by(Order) %>%
  summarise(n_classes = n_distinct(Class), classes = paste(unique(Class), collapse = ", ")) %>%
  filter(n_classes > 1)

print(inconsistent_orders)

saveRDS(ps_glom, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_glom_final_V3.rds")

#########SUITE  27/07

### ENREGISTRER le ps_glom tout propre
# Voir toutes les espèces dans ps_glom
all_species <- unique(as.data.frame(tax_table(ps_glom))$Species)

# Chercher les typos connus
grep("Abudefduf_sexfasciatu$", all_species, value = TRUE)  # lettre manquante
grep("Lujtanus", all_species, value = TRUE)  # faute de frappe Lutjanus
grep("Myripristis_sp$", all_species, value = TRUE)  # sp sans point
grep("Murdjan", all_species, value = TRUE)  # majuscule

# Et voir tous les _sp sans point
grep("_sp$", all_species, value = TRUE)





saveRDS(ps_glom, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_glom_final.rds")


##Encore une correction (dernière inch)
tax_mat <- as.matrix(tax_table(ps_glom))
tax_mat[c("ASV_617", "ASV_2040", "ASV_2165"), ]

tax_mat <- as.matrix(tax_table(ps_glom))

asvs_fix <- c("ASV_617", "ASV_2040", "ASV_2165")

tax_mat[asvs_fix, "Order"] <- "incertae_sedis"
tax_mat[asvs_fix, "Family"] <- "Pomacentridae"
tax_mat[asvs_fix, "Genus"] <- "Pomacentridae_gen."
tax_mat[asvs_fix, "Species"] <- "Pomacentridae_sp."

tax_table(ps_glom) <- tax_table(tax_mat)

# Vérifier
as.data.frame(tax_table(ps_glom))[asvs_fix, ]

tax_mat <- as.matrix(tax_table(ps_glom))

# Numéroter les _sp., _gen., _fam.
for (col in c("Species", "Genus", "Family")) {
  uncertain <- grep("_sp\\.|_gen\\.|_fam\\.", tax_mat[, col])
  
  if (length(uncertain) > 0) {
    vals <- tax_mat[uncertain, col]
    tax_mat[uncertain, col] <- ave(vals, vals, FUN = function(x) {
      if (length(x) == 1) x
      else paste0(x, seq_along(x))
    })
  }
}

# Corriger le double point
tax_mat <- gsub("_sp\\.\\.([0-9])", "_sp._\\1", tax_mat)
tax_mat <- gsub("_gen\\.\\.([0-9])", "_gen._\\1", tax_mat)
tax_mat <- gsub("_fam\\.\\.([0-9])", "_fam._\\1", tax_mat)

# Remettre dans ps_glom
tax_table(ps_glom) <- tax_table(tax_mat)

# Vérifier
sort(grep("Pomacentridae", unique(as.data.frame(tax_table(ps_glom))$Species), value = TRUE))

##Enregistrer à nouveau 

saveRDS(ps_glom, "~/Desktop/Research_project/Research_data_seychelles/All_data/Phyloseq/ps_glom_final_V2.rds")



##DEBUT ANALYSE : 

###On veut la liste de toutes les espèces certaines

tax_df <- as.data.frame(tax_table(ps_glom))

# Garder uniquement les espèces certaines (sans _sp.)
tax_certain <- tax_df[!grepl("_sp\\.", tax_df$Species), ]

View(tax_certain)
cat("Total espèces certaines:", nrow(tax_certain), "\n")
head(tax_certain)

#Garder que les actinopteri
tax_certain_actino <- tax_certain[tax_certain$Class == "Actinopteri", ]
cat("Total Actinopteri certains:", nrow(tax_certain_actino), "\n")
View(tax_certain_actino)


## On crée un objet phyloseq avec espèces certaines de actinopteri
# Créer l'objet phyloseq Actinopteri certains
ps_actino_certain <- subset_taxa(ps_glom, 
                                 Class == "Actinopteri" &
                                   !grepl("_sp\\.", Species))
##On a bien une dimension de 244 

df_actino <- psmelt(ps_actino_certain)
# Tu auras des colonnes : OTU, Sample, Abundance, Island, Time, Species, etc.
head(df_actino)

# Nombre d'espèces détectées par time point
df_actino %>%
  filter(Abundance > 0) %>%
  group_by(Time) %>%
  summarise(n_species = n_distinct(Species))



##Création 2 objets séparé pour chaque ile
# Félicité
ps_actino_felicite <- subset_samples(ps_actino_certain, Island == "Félicité")
ps_actino_felicite <- prune_taxa(taxa_sums(ps_actino_felicite) > 0, ps_actino_felicite)

# Frégate
ps_actino_fregate <- subset_samples(ps_actino_certain, Island == "Frégate")
ps_actino_fregate <- prune_taxa(taxa_sums(ps_actino_fregate) > 0, ps_actino_fregate)

# Vérifier
ntaxa(ps_actino_felicite)
ntaxa(ps_actino_fregate)



################ TEST AND TREND 
####FELICITE#####

df_felicite <- psmelt(ps_actino_felicite)
df_felicite %>%
  filter(Abundance > 0) %>%
  group_by(Time) %>%
  summarise(n_species = n_distinct(Species))

# Vérifier la structure
head(df_felicite)
colnames(df_felicite)

species_per_time_fe <- df_felicite %>%
  filter(Abundance > 0) %>%
  group_by(Time) %>%
  summarise(n_species = n_distinct(Species),
            total_reads = sum(Abundance)) %>%
  arrange(Time)

##Création graphique 
library(ggplot2)

# Créer un vecteur pour renommer les time points en heures réelles
time_labels <- c("Time 1" = "06:00", "Time 2" = "08:00", "Time 3" = "10:00",
                 "Time 4" = "12:00", "Time 5" = "14:00", "Time 6" = "16:00",
                 "Time 7" = "18:00")

# Graphique
ggplot(species_per_time_fe, aes(x = Time, y = n_species, group = 1)) +
  geom_line(color = "#2196F3", linewidth = 1) +
  geom_point(color = "#2196F3", size = 3) +
  scale_x_discrete(labels = time_labels) +
  labs(title = "Number of species detected over time - Félicité",
       x = "Time of day",
       y = "Number of species") +
  theme_classic()

# Félicité - richesse par sample
richness_fe <- estimate_richness(ps_actino_felicite, measures = "Observed")
richness_fe$Time <- sample_data(ps_actino_felicite)$Time

###KRUSKAL, MAIS EN FAIT NON
kruskal.test(Observed ~ Time, data = richness_fe)
view(richness_fe)

install.packages("MASS")
install.packages("emmeans")
library(emmeans)
library(MASS)

# Modèle binomial négatif - Félicité
glm_nb_fe <- MASS::glm.nb(Observed ~ Time, data = richness_fe)

# Test effet global
summary(glm_nb_fe)

# Comparaisons post-hoc
emmeans(glm_nb_fe, pairwise ~ Time, adjust = "tukey")

# Modèle nul (sans Time)
glm_nb_fe_null <- MASS::glm.nb(Observed ~ 1, data = richness_fe)

# Test du rapport de vraisemblance - effet global de Time
anova(glm_nb_fe_null, glm_nb_fe, test = "Chisq")



###FREGATE
# Modèle binomial négatif - Frégate
glm_nb_fr <- MASS::glm.nb(Observed ~ Time, data = richness_fr)

# Test effet global
summary(glm_nb_fr)

# Comparaisons post-hoc
emmeans(glm_nb_fr, pairwise ~ Time, adjust = "tukey")

# Modèle nul (sans Time)
glm_nb_fr_null <- MASS::glm.nb(Observed ~ 1, data = richness_fr)

# Test du rapport de vraisemblance - effet global de Time
anova(glm_nb_fr_null, glm_nb_fr, test = "Chisq")




# Frégate
richness_fr <- estimate_richness(ps_actino_fregate, measures = "Observed")
richness_fr$Time <- sample_data(ps_actino_fregate)$Time

kruskal.test(Observed ~ Time, data = richness_fr)

##ON save les PS

saveRDS(ps_actino_certain, "~/Desktop/Research_project/Research_data_seychelles/All_data/phyloseq/ps_actino_certain.V2.rds")
saveRDS(ps_actino_felicite, "~/Desktop/Research_project/Research_data_seychelles/All_data/phyloseq/ps_actino_felicite.V2.rds")
saveRDS(ps_actino_fregate, "~/Desktop/Research_project/Research_data_seychelles/All_data/phyloseq/ps_actino_fregate.V2.rds")



####Nombre de reads par heure 

# Graphique 2 - Total reads
ggplot(species_per_time_fe, aes(x = Time, y = total_reads, group = 1)) +
  geom_line(color = "#E91E63", linewidth = 1) +
  geom_point(color = "#E91E63", size = 3) +
  scale_x_discrete(labels = time_labels) +
  labs(title = "Number of reads detected over time - Félicité",
       x = "Time of day",
       y = "Number of reads")  +
  theme_classic()

##Carte NMDS - BRAY 

install.packages("vegan")
library(vegan)
library(ggplot2)

# 1. Extraire la table OTU de ton objet phyloseq Félicité
otu_fe <- as.data.frame(t(otu_table(ps_actino_felicite)))

# 2. Extraire les métadonnées (Time, Island etc.)
meta_fe <- as.data.frame(sample_data(ps_actino_felicite))

# 3. Faire le NMDS
nmds_fe <- metaMDS(otu_fe, 
                   distance = "bray",  # Bray-Curtis, standard pour abondance
                   k = 2,              # 2 dimensions
                   trymax = 100)       # nombre d'essais pour trouver la meilleure solution

# 4. Vérifier le stress
nmds_fe$stress


# 1. Extraire les coordonnées NMDS
nmds_coords <- as.data.frame(scores(nmds_fe, display = "sites"))
nmds_coords$Time <- meta_fe$Time

# 2. Graphique
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Time, fill = Time)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", alpha = 0.2, geom = "polygon") +  # ellipses par time point
  labs(title = "NMDS - Félicité (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fe$stress, 3))) +
  theme_classic() +
  scale_color_brewer(palette = "RdYlBu") +
  scale_fill_brewer(palette = "RdYlBu")



# Stressplot pour voir la qualité de la représentation
stressplot(nmds_fe)



###PCOA - BRAY

# 1. Faire le PCoA directement depuis ton objet phyloseq
pcoa_fe <- ordinate(ps_actino_felicite, 
                    method = "PCoA", 
                    distance = "bray")

# 2. Voir combien de variance est expliquée par chaque axe
plot_scree(pcoa_fe)

# 3. Faire le graphique
plot_ordination(ps_actino_felicite, 
                pcoa_fe,  
                color = "Time") +
  theme_classic() +
  scale_color_manual(values = c("#FF6B00", "#FF8C00", "#FFB347", 
                                "#FF69B4", "#FF1493", "#C71585", "#8B0057")) +
  labs(title = "PCoA - Félicité (Actinopteri)") + geom_point(size = 2)


# PERMANOVA Bray-Curtis (abondance)

# PERMANOVA Bray-Curtis
permanova_bray <- adonis2(dist_bray ~ sample_data(ps_actino_felicite)$Time, 
                          permutations = 999)
print(permanova_bray)






####### Matrice species x sample matrice binaire 


# STEP 1 Transformer en présence/absence
ps_actino_felicite_bi <- transform_sample_counts(ps_actino_felicite, 
                                                 function(x) ifelse(x > 0, 1, 0))

#STEP 2 PCOA - Jaccard 

# PCoA avec Jaccard
pcoa_fe_jac <- ordinate(ps_actino_felicite_bi, 
                        method = "PCoA", 
                        distance = "bray")

# Graphique
plot_ordination(ps_actino_felicite_bi, 
                pcoa_fe_jac, 
                color = "Time") +
  theme_classic() +
  scale_color_manual(values = c("#0D47A1", "#1976D2", "#42A5F5",
                                "#26A69A", "#66BB6A", "#AED581", "#F9A825")) +
  labs(title = "PCoA Jaccard - Félicité (Actinopteri)") +
geom_point(size = 2.5)


###NMDS - jaccard 
# Extraire la table OTU binaire
otu_fe_pa <- as.data.frame(t(otu_table(ps_actino_felicite_bi)))

# NMDS
nmds_fe_jac <- metaMDS(otu_fe_pa, 
                       distance = "jaccard", 
                       k = 2, 
                       trymax = 100,
                       binary = TRUE)

nmds_fe_jac$stress
# Extraire les coordonnées
nmds_jac_coords <- as.data.frame(scores(nmds_fe_jac, display = "sites"))
nmds_jac_coords$Time <- meta_fe$Time

# Graphique
ggplot(nmds_jac_coords, aes(x = NMDS1, y = NMDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_viridis_d() +
  labs(title = "NMDS Jaccard - Félicité (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fe_jac$stress, 3))) +
  theme_classic()


# Stressplot pour voir la qualité de la représentation
stressplot(nmds_fe_jac)
##Permanova - jaccard 
dist_jac <- phyloseq::distance(ps_actino_felicite_bi, method = "jaccard")

permanova_jac <- adonis2(dist_jac ~ sample_data(ps_actino_felicite)$Time, 
                         permutations = 999)
print(permanova_jac)


####FREGATE#####
df_fregate <- psmelt(ps_actino_fregate)
df_fregate %>%
  filter(Abundance > 0) %>%
  group_by(Time) %>%
  summarise(n_species = n_distinct(Species))

# Vérifier la structure
head(df_fregate)
colnames(df_fregate)

species_per_time_fr <- df_fregate %>%
  filter(Abundance > 0) %>%
  group_by(Time) %>%
  summarise(n_species = n_distinct(Species),
            total_reads = sum(Abundance)) %>%
  arrange(Time)

##création graphique 
library(ggplot2)

ggplot(species_per_time_fr, aes(x = Time, y = n_species, group = 1)) +
  geom_line(color = "#219", linewidth = 1) +
  geom_point(color = "#219", size = 3) +
  scale_x_discrete(labels = time_labels) +
  labs(title = "Number of species detected over time - Fregate",
       x = "Time of day",
       y = "Number of species") +
  theme_classic()

####Nombre de reads par heure 

# Graphique 2 - Total reads
ggplot(species_per_time_fr, aes(x = Time, y = total_reads, group = 1)) +
  geom_line(color = "#E91", linewidth = 1) +
  geom_point(color = "#E91", size = 3) +
  scale_x_discrete(labels = time_labels) +
  labs(title = "Number of reads detected over time - Fregate",
       x = "Time of day",
       y = "Number of reads")  +
  theme_classic()



##Carte NMDS - BRAY 

install.packages("vegan")
library(vegan)
library(ggplot2)

# 1. Extraire la table OTU de ton objet phyloseq Fregate
otu_fr <- as.data.frame(t(otu_table(ps_actino_fregate)))

# 2. Extraire les métadonnées (Time, Island etc.)
meta_fr <- as.data.frame(sample_data(ps_actino_fregate))

# 3. Faire le NMDS
nmds_fr <- metaMDS(otu_fr, 
                   distance = "bray",  # Bray-Curtis, standard pour abondance
                   k = 2,              # 2 dimensions
                   trymax = 100)       # nombre d'essais pour trouver la meilleure solution

# 4. Vérifier le stress
nmds_fr$stress

# 5. Extraire les coordonnées
nmds_coords_fr <- as.data.frame(scores(nmds_fr, display = "sites"))
nmds_coords_fr$Time <- meta_fr$Time

# 6. Graphique
ggplot(nmds_coords_fr, aes(x = NMDS1, y = NMDS2, color = Time, fill = Time)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", alpha = 0.2, geom = "polygon") +
  labs(title = "NMDS - Frégate (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fr$stress, 3))) +
  theme_classic() +
  scale_color_manual(values = c("#D62728", "#FF7F0E", "#FFDB58", 
                                "#2CA02C", "#17BECF", "#1F77B4", "#9467BD")) +
  scale_fill_manual(values = c("#D62728", "#FF7F0E", "#FFDB58",
                               "#2CA02C", "#17BECF", "#1F77B4", "#9467BD"))


# Stressplot pour voir la qualité de la représentation
stressplot(nmds_fr)


###PCOA - BRAY

# 1. Faire le PCoA directement depuis ton objet phyloseq
pcoa_fr <- ordinate(ps_actino_fregate, 
                    method = "PCoA", 
                    distance = "bray")

# 2. Voir combien de variance est expliquée par chaque axe
plot_scree(pcoa_fr)

# 3. Faire le graphique
plot_ordination(ps_actino_fregate, 
                pcoa_fr,  
                color = "Time") +
  theme_classic() +
  scale_color_manual(values = c("#2f6a8a", "#4a9aba", "#17BECF",
                                "#45b899", "#8bc34f", "#c8e058", "#FFDB58")) +
  labs(title = "PCoA - Fregate (Actinopteri)") +
  geom_point(size = 2)  


# PERMANOVA Bray-Curtis (abondance)
# Créer la bonne matrice de distance pour Frégate
dist_bray_fr <- phyloseq::distance(ps_actino_fregate, method = "bray")

# PERMANOVA Bray-Curtis
permanova_bray_fr <- adonis2(dist_bray_fr ~ sample_data(ps_actino_fregate)$Time, 
                          permutations = 999)
print(permanova_bray_fr)



####### Matrice species x sample matrice binaire 


# STEP 1 Transformer en présence/absence
ps_actino_fregate_bi <- transform_sample_counts(ps_actino_fregate, 
                                                 function(x) ifelse(x > 0, 1, 0))

#STEP 2 PCOA - Jaccard 

# PCoA avec Jaccard
pcoa_fr_jac <- ordinate(ps_actino_fregate_bi, 
                        method = "PCoA", 
                        distance = "jaccard")

# Graphique
plot_ordination(ps_actino_fregate_bi, 
                pcoa_fr_jac, 
                color = "Time") +
  theme_classic() +
  scale_color_manual(values = c("#0D47A1", "#1976D2", "#42A5F5",
                                "#26A69A", "#66BB6A", "#AED581", "#F9A825")) +
  labs(title = "PCoA Jaccard - Fregate (Actinopteri)") + 
  geom_point(size = 2.5) 



###NMDS - jaccard 
# Extraire la table OTU binaire
otu_fr_pa <- as.data.frame(t(otu_table(ps_actino_fregate_bi)))

# NMDS
nmds_fr_jac <- metaMDS(otu_fr_pa, 
                       distance = "jaccard", 
                       k = 2, 
                       trymax = 100,
                       binary = TRUE)

nmds_fr_jac$stress

# Extraire les coordonnées
nmds_jac_coords_fr <- as.data.frame(scores(nmds_fr_jac, display = "sites"))
nmds_jac_coords_fr$Time <- meta_fr$Time

# Graphique
ggplot(nmds_jac_coords_fr, aes(x = NMDS1, y = NMDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_viridis_d() +
  labs(title = "NMDS Jaccard - Fregate (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fr_jac$stress, 3))) +
  theme_classic()


stressplot(nmds_fr_jac)

##Permanova - jaccard 
dist_jac_fr <- phyloseq::distance(ps_actino_fregate_bi, method = "jaccard")

permanova_jac_fr <- adonis2(dist_jac_fr ~ sample_data(ps_actino_fregate_bi)$Time, 
                         permutations = 999)
print(permanova_jac_fr)




############Sequencing DEPTH##############

# Félicité
depth_fe <- data.frame(
  Sample = sample_names(ps_actino_felicite),
  Depth = sample_sums(ps_actino_felicite),
  Time = sample_data(ps_actino_felicite)$Time
)

# Frégate
depth_fr <- data.frame(
  Sample = sample_names(ps_actino_fregate),
  Depth = sample_sums(ps_actino_fregate),
  Time = sample_data(ps_actino_fregate)$Time
)

# Voir les tableaux
print(depth_fe)
print(depth_fr)

# Stats résumées
cat("Félicité - Min:", min(depth_fe$Depth), "Max:", max(depth_fe$Depth), "Mean:", round(mean(depth_fe$Depth)), "\n")
cat("Frégate - Min:", min(depth_fr$Depth), "Max:", max(depth_fr$Depth), "Mean:", round(mean(depth_fr$Depth)), "\n")



# Bar chart Félicité
ggplot(depth_fe, aes(x = Sample, y = Depth, fill = Time)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sequencing depth per sample - Félicité",
       x = "Sample",
       y = "Number of reads") +
  scale_fill_viridis_d()

# Bar chart Frégate
ggplot(depth_fr, aes(x = Sample, y = Depth, fill = Time)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sequencing depth per sample - Frégate",
       x = "Sample",
       y = "Number of reads") +
  scale_fill_viridis_d()




###Tous les read depth combiné 


# Read depth de tous les samples combinés
depth_all <- data.frame(
  Sample = sample_names(ps_actino_certain),
  Depth = sample_sums(ps_actino_certain),
  Time = sample_data(ps_actino_certain)$Time,
  Island = sample_data(ps_actino_certain)$Island
)

# Stats résumées
cat("Min:", min(depth_all$Depth), "\n")
cat("Max:", max(depth_all$Depth), "\n")
cat("Mean:", round(mean(depth_all$Depth)), "\n")

# Bar chart
ggplot(depth_all, aes(x = Sample, y = Depth, fill = Island)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(title = "Sequencing depth per sample - All samples",
       x = "Sample",
       y = "Number of reads") +
  scale_fill_manual(values = c("#2196F3", "#FF5722"))




# === CORRECTIONS TAXONOMIQUES (already applied in ps_final_corrected.rds) ===
# Extraire toutes les espèces de ps_final
all_species_final <- unique(as.data.frame(tax_table(ps_final))$Species)

# Chercher chaque typo
grep("Abudefduf_sexfasciatu$", all_species_final, value = TRUE)
grep("Lujtanus", all_species_final, value = TRUE)
grep("Myripristis_Murdjan", all_species_final, value = TRUE)
grep("Myripristis_sp$", all_species_final, value = TRUE)


# Extraire la tax table
tax_mat <- as.matrix(tax_table(ps_final))

# Corriger les typos

###FELICITE
tax_mat[tax_mat == "Abudefduf_sexfasciatu"] <- "Abudefduf_sexfasciatus"
tax_mat[tax_mat == "Lujtanus_rufolineatus"] <- "Lutjanus_rufolineatus"
tax_mat[tax_mat == "Myripristis_Murdjan"] <- "Myripristis_murdjan"
tax_mat[tax_mat == "Myripristis_sp"] <- "Myripristis_sp."
tax_mat[tax_mat == "Lujtanus_sp."] <- "Lutjanus_sp."
tax_mat[tax_mat == "Abudefduf_septemfiatus"] <- "Abudefduf_septemfasciatus"
tax_mat[tax_mat == "Platax_teira"] <- "Platax_boersii" ##Fausse correction
tax_mat[tax_mat == "Carangoides_fulvoguttatus"] <- "Turrum_fulvoguttatum" ##Fause correction



##FREGATE
tax_mat[tax_mat == "Pterocaesi_tile"] <- "Pterocaesio_tile"  # lettre manquante
tax_mat[tax_mat == "Eviota_sp"] <- "Eviota_sp."  # point manquant



##Autre correction
tax_mat[tax_mat == "Callionymu_sp."] <- "Callionymus_sp."

View(tax_table(ps_final))

tax_mat[tax_mat == "Chhaetodon_sp."] <- "Chaetodon_sp."



# Remettre dans ps_final
tax_table(ps_final) <- tax_table(tax_mat)

# Vérifier que les corrections ont bien été appliquées
grep("Abudefduf_septemfiatus", as.data.frame(tax_table(ps_final))$Species, value = TRUE)  # doit retourner rien
grep("Lujtanus_rufolineatus", as.data.frame(tax_table(ps_final))$Species, value = TRUE)  # doit retourner le nom corrigé





###mettre sur le même graphique les courbes espèces per time

# Combiner les deux datasets
species_per_time_fe$Island <- "Félicité"
species_per_time_fr$Island <- "Frégate"
species_combined <- rbind(species_per_time_fe, species_per_time_fr)

# Graphique combiné

# Combiner les deux datasets
species_per_time_fe$Island <- "Félicité"
species_per_time_fr$Island <- "Frégate"
species_combined <- rbind(species_per_time_fe, species_per_time_fr)

ggplot(species_combined, aes(x = Time, y = n_species, group = Island, color = Island)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Félicité" = "#2196F3", "Frégate" = "#9C27B0")) +
  scale_x_discrete(labels = time_labels) +
  labs(title = "Number of species detected over time",
       x = "Time of day",
       y = "Number of species",
       color = "Island") +
  theme_classic()
####PENSER A RE ENGITRER LES OBJETS PHYLOSEQ
#####Check que c'est bien dans ps_glom


# Combiner les deux objets phyloseq
ps_actino_all <- merge_phyloseq(ps_actino_felicite, ps_actino_fregate)

# PCoA Bray-Curtis sur tous les samples
pcoa_all <- ordinate(ps_actino_all, method = "PCoA", distance = "bray")

# Extraire coordonnées
pcoa_all_coords <- plot_ordination(ps_actino_all, pcoa_all, justDF = TRUE)

# Graphique
ggplot(pcoa_all_coords, aes(x = Axis.1, y = Axis.2, 
                            color = Time, shape = Island)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c("Félicité" = 17, "Frégate" = 15)) +
  scale_color_manual(values = c("#D62728", "#FF7F0E", "#FFDB58",
                                "#2CA02C", "#17BECF", "#1F77B4", "#9467BD")) +
  theme_classic() +
  labs(title = "PCoA - Bray-Curtis - All samples",
       x = paste0("PCo1 [", round(pcoa_all$values$Relative_eig[1]*100, 1), "%]"),
       y = paste0("PCo2 [", round(pcoa_all$values$Relative_eig[2]*100, 1), "%]"))




