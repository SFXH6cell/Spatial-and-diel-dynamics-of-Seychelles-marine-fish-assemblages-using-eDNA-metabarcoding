####Read depth vs richness plot → rapide à faire
##iNEXT curves → identifier les samples problématiques
##Décider si on retire des samples ou non
###Ensuite continuer les analyses écologiques

# Read depth vs richness plot
depth_richness <- data.frame(
  Sample = sample_names(ps_actino_certain),
  Depth = sample_sums(ps_actino_certain),
  Richness = estimate_richness(ps_actino_certain, measures = "Observed")$Observed,
  Island = sample_data(ps_actino_certain)$Island,
  Time = sample_data(ps_actino_certain)$Time
)

# Graphique
ggplot(depth_richness, aes(x = Depth, y = Richness, color = Island, label = Sample)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "grey50") +
  theme_classic() +
  scale_color_manual(values = c("#2196F3", "#FF5722")) +
  labs(title = "Read depth vs Species richness",
       x = "Number of reads (sequencing depth)",
       y = "Species richness (observed)",
       color = "Island")



###INEXT

install.packages("iNEXT")
library(iNEXT)
max(sample_sums(ps_actino_fregate))
# Extraire la table OTU pour Frégate
otu_fr_inext <- as.data.frame(t(otu_table(ps_actino_fregate)))

# Convertir en liste pour iNEXT
inext_list_fr <- as.list(otu_fr_inext)

# Lancer iNEXT
# Extraire uniquement les 3 samples FrT3
otu_fr_t3 <- as.data.frame(t(otu_table(ps_actino_fregate)))
otu_fr_t3_subset <- otu_fr_t3[grep("FrT3", rownames(otu_fr_t3)), ]

# Convertir en liste
inext_list_t3 <- as.list(t(otu_fr_t3_subset))

# iNEXT juste sur ces 3 samples
inext_t3 <- iNEXT(inext_list_t3,
                  q = 0,
                  datatype = "abundance",
                  endpoint = 400000)

# Graphique
ggiNEXT(inext_t3, type = 1) +
  theme_classic() +
  labs(title = "iNEXT - FrT3 samples (low depth)",
       x = "Number of reads",
       y = "Species richness")


# Extraire les 3 samples FrT3
otu_fr_t3 <- as.data.frame(t(otu_table(ps_actino_fregate)))
otu_fr_t3_subset <- otu_fr_t3[grep("FrT3", rownames(otu_fr_t3)), ]

# Convertir en liste ET retirer les espèces avec 0 reads
inext_list_t3 <- lapply(as.list(t(otu_fr_t3_subset)), function(x) x[x > 0])

# iNEXT
inext_t3 <- iNEXT(inext_list_t3,
                  q = 0,
                  datatype = "abundance",
                  endpoint = 300000)

# Graphique
ggiNEXT(inext_t3, type = 1) +
  theme_classic() +
  labs(title = "iNEXT - FrT3 samples (low depth)",
       x = "Number of reads",
       y = "Species richness")


# Extraire la table OTU pour Frégate
otu_fr <- as.data.frame(otu_table(ps_actino_fregate))

# Vérifier l'orientation - les samples sont en colonnes ou en lignes?
dim(otu_fr)
colnames(otu_fr)[1:5]
rownames(otu_fr)[1:5]


# Extraire les colonnes FrT3
otu_fr_t3 <- otu_fr[, grep("FrT3", colnames(otu_fr))]

# Retirer les lignes avec 0 dans tous les samples FrT3
otu_fr_t3 <- otu_fr_t3[rowSums(otu_fr_t3) > 0, ]

# Convertir en liste par sample
inext_list_t3 <- as.list(otu_fr_t3)

# Vérifier qu'il n'y a pas de zéros
lapply(inext_list_t3, function(x) sum(x == 0))

# iNEXT
inext_t3 <- iNEXT(inext_list_t3,
                  q = 0,
                  datatype = "abundance",
                  endpoint = 300000)

ggiNEXT(inext_t3, type = 1) +
  theme_classic() +
  labs(title = "iNEXT - FrT3 samples (low depth)",
       x = "Number of reads",
       y = "Species richness")

###POUR TOUS
# Extraire la table OTU pour tous les samples (les deux îles)
otu_all <- as.data.frame(otu_table(ps_actino_certain))

# Retirer les lignes avec 0 dans tous les samples
otu_all <- otu_all[rowSums(otu_all) > 0, ]

# Convertir en liste par sample
inext_list_all <- as.list(otu_all)

# iNEXT
inext_all <- iNEXT(inext_list_all,
                   q = 0,
                   datatype = "abundance",
                   endpoint = 477688)

# Graphique séparé par île pour plus de lisibilité
# Félicité
ggiNEXT(inext_all, type = 1) +
  theme_classic() +
  facet_wrap(~ Assemblage, scales = "free") +
  labs(title = "iNEXT - All samples",
       x = "Number of reads",
       y = "Species richness") +
  theme(legend.position = "none",
        strip.text = element_text(size = 6))




# Félicité seulement
otu_fe_inext <- otu_all[, grep("FeT", colnames(otu_all))]
otu_fe_inext <- otu_fe_inext[rowSums(otu_fe_inext) > 0, ]
inext_list_fe <- as.list(otu_fe_inext)
inext_fe <- iNEXT(inext_list_fe, q = 0, datatype = "abundance", endpoint = 477688)







####Hellinger 

# Transformation Hellinger
# Extraire les tables OTU
otu_fe <- as.data.frame(t(otu_table(ps_actino_felicite)))
otu_fr <- as.data.frame(t(otu_table(ps_actino_fregate)))


otu_fe_hel <- decostand(otu_fe, method = "hellinger")
otu_fr_hel <- decostand(otu_fr, method = "hellinger")


head(otu_fe_hel[, 1:5])

dim(otu_fe)
dim(otu_fr)


# NMDS Hellinger - Félicité
nmds_fe_hel <- metaMDS(otu_fe_hel,
                       distance = "euclidean",
                       k = 2,
                       trymax = 100)

nmds_fe_hel$stress

# NMDS Hellinger - Frégate
nmds_fr_hel <- metaMDS(otu_fr_hel,
                       distance = "euclidean",
                       k = 2,
                       trymax = 100)

nmds_fr_hel$stress


# Félicité
nmds_coords_fe_hel <- as.data.frame(scores(nmds_fe_hel, display = "sites"))
nmds_coords_fe_hel$Time <- meta_fe$Time

ggplot(nmds_coords_fe_hel, aes(x = NMDS1, y = NMDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#67001F", "#990000", "#CC0000",
                                "#E34A33", "#FC8D59", "#FDD49E", "#FFFF33")) +
  labs(title = "NMDS Hellinger - Félicité (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fe_hel$stress, 3))) +
  theme_classic()

# Frégate
nmds_coords_fr_hel <- as.data.frame(scores(nmds_fr_hel, display = "sites"))
nmds_coords_fr_hel$Time <- meta_fr$Time

ggplot(nmds_coords_fr_hel, aes(x = NMDS1, y = NMDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#053061", "#2166AC", "#4393C3",
                                "#74ADD1", "#ABD9E9", "#80CDC1", "#35978F")) +
  labs(title = "NMDS Hellinger - Frégate (Actinopteri)",
       subtitle = paste("Stress =", round(nmds_fr_hel$stress, 3))) +
  theme_classic()





###PCOA 

library(vegan)

# PCoA Hellinger - Félicité
pcoa_fe_hel <- capscale(otu_fe_hel ~ 1, distance = "euclidean")

# Extraire les coordonnées
pcoa_fe_hel_coords <- as.data.frame(scores(pcoa_fe_hel, display = "sites"))
pcoa_fe_hel_coords$Time <- meta_fe$Time

# % variance expliquée
eigenvalues <- pcoa_fe_hel$CA$eig
var_explained <- eigenvalues / sum(eigenvalues) * 100
cat("Axe 1:", round(var_explained[1], 1), "%\n")
cat("Axe 2:", round(var_explained[2], 1), "%\n")

# Graphique
ggplot(pcoa_fe_hel_coords, aes(x = MDS1, y = MDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#FF6B00", "#FF8C00", "#FFB347", 
                                 "#FF69B4", "#FF1493", "#C71585", "#8B0057")) +
  labs(title = "PCoA Hellinger - Félicité (Actinopteri)",
       x = paste0("PCo1 [", round(var_explained[1], 1), "%]"),
       y = paste0("PCo2 [", round(var_explained[2], 1), "%]")) +
  theme_classic()


# PCoA Hellinger - Frégate
pcoa_fr_hel <- capscale(otu_fr_hel ~ 1, distance = "euclidean")

# Extraire les coordonnées
pcoa_fr_hel_coords <- as.data.frame(scores(pcoa_fr_hel, display = "sites"))
pcoa_fr_hel_coords$Time <- meta_fr$Time

# % variance expliquée
eigenvalues_fr <- pcoa_fr_hel$CA$eig
var_explained_fr <- eigenvalues_fr / sum(eigenvalues_fr) * 100
cat("Axe 1:", round(var_explained_fr[1], 1), "%\n")
cat("Axe 2:", round(var_explained_fr[2], 1), "%\n")

# Graphique
ggplot(pcoa_fr_hel_coords, aes(x = MDS1, y = MDS2, color = Time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#08306B", "#08519C", "#2171B5",
                                "#4292C6", "#006D2C", "#238B45", "#41AB5D")) +
  labs(title = "PCoA Hellinger - Frégate (Actinopteri)",
       x = paste0("PCo1 [", round(var_explained_fr[1], 1), "%]"),
       y = paste0("PCo2 [", round(var_explained_fr[2], 1), "%]")) +
  theme_classic()

# PERMANOVA Hellinger - Félicité
permanova_fe_hel <- adonis2(otu_fe_hel ~ sample_data(ps_actino_felicite)$Time,
                            method = "euclidean",
                            permutations = 999)
print(permanova_fe_hel)

# PERMANOVA Hellinger - Frégate
permanova_fr_hel <- adonis2(otu_fr_hel ~ sample_data(ps_actino_fregate)$Time,
                            method = "euclidean",
                            permutations = 999)
print(permanova_fr_hel)




##SIMPE
# SIMPER - Félicité
simper_fe <- simper(otu_fe_hel, 
                    group = sample_data(ps_actino_felicite)$Time,
                    permutations = 999)

# Résumé
summary(simper_fe)

# Renommer les colonnes de otu_fe_hel avec les noms d'espèces
tax_fe <- as.data.frame(tax_table(ps_actino_felicite))
colnames(otu_fe_hel) <- tax_fe$Species[match(colnames(otu_fe_hel), rownames(tax_fe))]

# Relancer SIMPER
simper_fe <- simper(otu_fe_hel,
                    group = sample_data(ps_actino_felicite)$Time,
                    permutations = 999)

# Top 5 par comparaison
for(i in names(simper_fe)) {
  cat("\n===", i, "===\n")
  print(head(summary(simper_fe)[[i]], 5))
}


# Extraire top 1 espèce par comparaison
top_species <- sapply(names(simper_fe), function(x) {
  rownames(summary(simper_fe)[[x]])[1]
})

# Voir quelles espèces reviennent le plus
sort(table(top_species), decreasing = TRUE)


# Pour chaque espèce top 5, voir dans quels time points elle est abondante
top5_species <- names(sort(table(top_species), decreasing = TRUE))[1:5]

df_felicite %>%
  filter(Species %in% top5_species) %>%
  group_by(Time, Species) %>%
  summarise(mean_reads = mean(Abundance)) %>%
  ggplot(aes(x = Time, y = mean_reads, color = Species, group = Species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "Diel pattern - Top 5 SIMPER species - Félicité",
       x = "Time point",
       y = "Mean reads") +
  scale_color_viridis_d()

# Top 5 espèces qui reviennent le plus souvent comme top contributeur
sort(table(top_species), decreasing = TRUE)[1:5]


###FREGATE 

# Renommer les colonnes de otu_fr_hel avec les noms d'espèces
tax_fr <- as.data.frame(tax_table(ps_actino_fregate))
colnames(otu_fr_hel) <- tax_fr$Species[match(colnames(otu_fr_hel), rownames(tax_fr))]

# SIMPER Frégate
simper_fr <- simper(otu_fr_hel,
                    group = sample_data(ps_actino_fregate)$Time,
                    permutations = 999)

# Top espèce par comparaison
top_species_fr <- sapply(names(simper_fr), function(x) {
  rownames(summary(simper_fr)[[x]])[1]
})

sort(table(top_species_fr), decreasing = TRUE)

# Top 5 espèces SIMPER Frégate
top5_species_fr <- names(sort(table(top_species_fr), decreasing = TRUE))[1:5]

df_fregate %>%
  filter(Species %in% top5_species_fr) %>%
  group_by(Time, Species) %>%
  summarise(mean_reads = mean(Abundance)) %>%
  ggplot(aes(x = Time, y = mean_reads, color = Species, group = Species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "Diel pattern - Top 5 SIMPER species - Frégate",
       x = "Time point",
       y = "Mean reads") +
  scale_color_viridis_d()



## SIMPER moyenne et autres
#Felicite

# Moyenne de dissimilarité globale
mean_diss <- sapply(simper_fe, function(x) mean(x$overall))
cat("Mean dissimilarity range:", round(min(mean_diss)*100, 1), 
    "% to", round(max(mean_diss)*100, 1), "%\n")
cat("Overall mean:", round(mean(mean_diss)*100, 1), "%\n")

#Fregate
mean_diss_fr <- sapply(simper_fr, function(x) mean(x$overall))
cat("Frégate - Mean dissimilarity range:", round(min(mean_diss_fr)*100, 1), 
    "% to", round(max(mean_diss_fr)*100, 1), "%\n")
cat("Frégate - Overall mean:", round(mean(mean_diss_fr)*100, 1), "%\n")


##Kruskall - Wallis - SIMPER
df_selar <- df_felicite %>%
  filter(Species == "Selar_crumenophthalmus") %>%
  group_by(Sample, Time) %>%
  summarise(total_reads = sum(Abundance), .groups = "drop")

kruskal.test(total_reads ~ Time, data = df_selar)




##SHANNON 

# Shannon diversity - Félicité
shannon_fe <- estimate_richness(ps_actino_felicite, measures = "Shannon")
shannon_fe$Time <- sample_data(ps_actino_felicite)$Time

# Graphique
ggplot(shannon_fe, aes(x = Time, y = Shannon, fill = Time)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#67001F", "#990000", "#CC0000",
                               "#E34A33", "#FC8D59", "#D94701", "#8C2D04")) +
  labs(title = "Shannon diversity - Félicité",
       x = "Time point",
       y = "Shannon diversity index")

# Kruskal-Wallis
kruskal.test(Shannon ~ Time, data = shannon_fe)


# Shannon diversity - Frégate
shannon_fr <- estimate_richness(ps_actino_fregate, measures = "Shannon")
shannon_fr$Time <- sample_data(ps_actino_fregate)$Time

kruskal.test(Shannon ~ Time, data = shannon_fr)

















  