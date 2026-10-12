####### Diel cycle - functionnal trait #####



####Diet 

### STEP 1: Prepare the data

trophic_guilds <- read.csv("~/Downloads/FishFecesCNP-main/data/extrapolation_trophic_guilds.csv")

guild_labels <- c("1" = "Sessile invertivore",
                  "2" = "Herbivore/microvore/detritivore",
                  "3" = "Corallivore",
                  "4" = "Piscivore",
                  "5" = "Microinvertivore",
                  "6" = "Macroinvertivore",
                  "7" = "Crustacivore",
                  "8" = "Planktivore")

trophic_guilds$guild_name <- guild_labels[as.character(trophic_guilds$trophic_guild_predicted)]


###FELICITE
## Merge onto your tax table (using the same cleaned name format as before)
df_felicite_troph <- df_felicite %>%
  left_join(trophic_guilds %>%
              select(species, guild_name, trophic_guild_predicted),
            by = c("Species" = "species"))

####### Structure de trophic_guilds
colnames(trophic_guilds)
head(trophic_guilds$species)  # ou le nom de la colonne espèces

@####### Format des espèces dans df_felicite
head(df_felicite$Species)


## Check coverage
cat("Species with a trophic guild assigned:",
    n_distinct(df_felicite_troph$Species[!is.na(df_felicite_troph$guild_name)]), 
    "out of",
    n_distinct(df_felicite_troph$Species), "\n")

cat("Species WITHOUT trophic guild:",
    n_distinct(df_felicite_troph$Species[is.na(df_felicite_troph$guild_name)]),
    "\n")


# Voir les 32 espèces manquantes
missing_guild_fe <- df_felicite_troph %>%
  filter(is.na(guild_name)) %>%
  distinct(Species) %>%
  pull(Species)

print(sort(missing_guild_fe))

###32 espèces manquantes, et 171 assignées. 



########Je dois corriger tout ca, mais là je veux juste voir globalement : 

####Number of species per diet per time point
df_felicite_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Time, guild_name) %>%
  summarise(n_species = n_distinct(Species)) %>%
  ggplot(aes(x = Time, y = n_species, color = guild_name, group = guild_name)) +
  geom_line() + geom_point() +
  theme_classic()


###Number of reads 
df_felicite_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Time, guild_name) %>%
  summarise(total_reads = sum(Abundance)) %>%
  ggplot(aes(x = Time, y = total_reads, color = guild_name, group = guild_name)) +
  geom_line() + geom_point() +
  theme_classic()



# Voir les espèces herbivores
df_felicite_troph %>%
  filter(guild_name == "Herbivore/microvore/detritivore",
         Abundance > 0) %>%
  distinct(Species) %>%
  pull(Species)

# Préparer les données
heatmap_herb_fe <- df_felicite_troph %>%
  filter(guild_name == "Herbivore/microvore/detritivore",
         Abundance > 0) %>%
  group_by(Time, Species) %>%
  summarise(total_reads = sum(Abundance)) %>%
  ungroup()

# Heatmap
ggplot(heatmap_herb_fe, aes(x = Time, y = Species, fill = total_reads)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 7)) +
  labs(title = "Diel cycle - Herbivores - Félicité",
       x = "Time point",
       y = "Species",
       fill = "Total reads")



###Clustering huerarchique (en vrai à effacer surement)

install.packages("pheatmap")
library(pheatmap)

install.packages("tibble")
library(tibble)

install.packages("tidyr")
library(tidyr)

install.packages("viridis")
library(viridis)



# Préparer la matrice
herb_matrix <- df_felicite_troph %>%
  filter(guild_name == "Herbivore/microvore/detritivore") %>%
  group_by(Time, Species) %>%
  summarise(total_reads = sum(Abundance)) %>%
  pivot_wider(names_from = Time, values_from = total_reads, values_fill = 0) %>%
  column_to_rownames("Species") %>%
  as.matrix()



df_felicite_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Time, guild_name) %>%
  summarise(total_reads = sum(Abundance)) %>%
  ggplot(aes(x = Time, y = total_reads, color = guild_name, group = guild_name)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_viridis_d() +
  labs(title = "Diel cycle by trophic guild - Félicité",
       x = "Time point",
       y = "Total reads",
       color = "Trophic guild")







####FREGATE

## Merge Frégate
df_fregate_troph <- df_fregate %>%
  left_join(trophic_guilds %>%
              select(species, guild_name, trophic_guild_predicted),
            by = c("Species" = "species"))

## Check coverage
cat("Species with a trophic guild assigned:",
    n_distinct(df_fregate_troph$Species[!is.na(df_fregate_troph$guild_name)]), 
    "out of",
    n_distinct(df_fregate_troph$Species), "\n")
cat("Species WITHOUT trophic guild:",
    n_distinct(df_fregate_troph$Species[is.na(df_fregate_troph$guild_name)]),
    "\n")

# Voir les espèces manquantes
missing_guild_fr <- df_fregate_troph %>%
  filter(is.na(guild_name)) %>%
  distinct(Species) %>%
  pull(Species)
print(sort(missing_guild_fr))

# Graphique - Nombre d'espèces par guilde par time point
df_fregate_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Time, guild_name) %>%
  summarise(n_species = n_distinct(Species)) %>%
  ggplot(aes(x = Time, y = n_species, color = guild_name, group = guild_name)) +
  geom_line() + geom_point() +
  theme_classic() +
  labs(title = "Number of species per trophic guild - Frégate",
       x = "Time point", y = "Number of species", color = "Trophic guild")

# Graphique - Total reads par guilde par time point
df_fregate_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Time, guild_name) %>%
  summarise(total_reads = sum(Abundance)) %>%
  ggplot(aes(x = Time, y = total_reads, color = guild_name, group = guild_name)) +
  geom_line() + geom_point() +
  theme_classic() +
  labs(title = "Total reads per trophic guild - Frégate",
       x = "Time point", y = "Total reads", color = "Trophic guild")



####Est-ce que le guilde predit le Régime alimentaire d'un poisson ? "# Créer une matrice OTU par guilde

guild_matrix_fe <- df_felicite_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Sample, guild_name) %>%
  summarise(total_reads = sum(Abundance)) %>%
  pivot_wider(names_from = guild_name, 
              values_from = total_reads, 
              values_fill = 0) %>%
  column_to_rownames("Sample")

# Métadonnées
meta_guild_fe <- data.frame(
  Time = sample_data(ps_actino_felicite)$Time
)

# PERMANOVA
adonis2(guild_matrix_fe ~ Time,
        data = meta_guild_fe,
        method = "bray",
        permutations = 999)



###Fregate 

guild_matrix_fr <- df_fregate_troph %>%
  filter(!is.na(guild_name), Abundance > 0) %>%
  group_by(Sample, guild_name) %>%
  summarise(total_reads = sum(Abundance)) %>%
  pivot_wider(names_from = guild_name,
              values_from = total_reads,
              values_fill = 0) %>%
  column_to_rownames("Sample")

meta_guild_fr <- data.frame(
  Time = sample_data(ps_actino_fregate)$Time
)

adonis2(guild_matrix_fr ~ Time,
        data = meta_guild_fr,
        method = "bray",
        permutations = 999)




# Kruskal-Wallis par guilde - Félicité
guilds <- unique(df_felicite_troph$guild_name[!is.na(df_felicite_troph$guild_name)])

for(guild in guilds) {
  df_guild <- df_felicite_troph %>%
    filter(guild_name == guild) %>%
    group_by(Sample, Time) %>%
    summarise(total_reads = sum(Abundance), .groups = "drop")
  
  kw <- kruskal.test(total_reads ~ Time, data = df_guild)
  cat("\n===", guild, "===\n")
  cat("Chi-squared =", round(kw$statistic, 3), 
      "| df =", kw$parameter,
      "| p =", round(kw$p.value, 3), "\n")
}

##Fregate
guilds_fr <- unique(df_fregate_troph$guild_name[!is.na(df_fregate_troph$guild_name)])

for(g in guilds_fr) {
  df_guild <- df_fregate_troph[df_fregate_troph$guild_name == g & 
                                 !is.na(df_fregate_troph$guild_name), ]
  
  df_guild_rich <- aggregate(Species ~ Sample + Time,
                             data = df_guild[df_guild$Abundance > 0, ],
                             FUN = function(x) length(unique(x)))
  
  kw <- kruskal.test(Species ~ Time, data = df_guild_rich)
  cat("\n===", g, "===\n")
  cat("Chi-squared =", round(kw$statistic, 3),
      "| df =", kw$parameter,
      "| p =", round(kw$p.value, 3), "\n")
}

# Exemple pour les herbivores
df_herb <- df_felicite_troph %>%
  filter(guild_name == "Herbivore/microvore/detritivore") %>%
  group_by(Sample, Time) %>%
  summarise(total_reads = sum(Abundance), .groups = "drop")

kruskal.test(total_reads ~ Time, data = df_herb)

####CORRECTION ESPECES MANQUANTE : 


####on fait le PCOA des guildes 
##FELICTE
# Vérifier les colonnes de tax_table
colnames(tax_table(ps_actino_felicite))
# Extraire la tax table
tax_fe <- as.data.frame(tax_table(ps_actino_felicite))

# Ajouter les guildes
tax_fe <- tax_fe %>%
  left_join(trophic_guilds %>% select(species, guild_name),
            by = c("Species" = "species"))

# Vérifier
head(tax_fe)

# Ajouter les guildes dans la tax_table de ps_actino_felicite
tax_mat_fe <- as.matrix(tax_fe)
rownames(tax_mat_fe) <- rownames(tax_table(ps_actino_felicite))
tax_table(ps_actino_felicite) <- tax_table(tax_mat_fe)

# Vérifier
colnames(tax_table(ps_actino_felicite))


# Filtrer ps_actino_felicite pour les herbivores
ps_herb_fe <- subset_taxa(ps_actino_felicite, 
                          guild_name == "Herbivore/microvore/detritivore")
ps_herb_fe <- prune_taxa(taxa_sums(ps_herb_fe) > 0, ps_herb_fe)

# Vérifier
ntaxa(ps_herb_fe)

# PCoA Bray-Curtis
pcoa_herb_fe <- ordinate(ps_herb_fe, method = "PCoA", distance = "bray")

# % variance expliquée
cat("PCo1:", round(pcoa_herb_fe$values$Relative_eig[1]*100, 1), "%\n")
cat("PCo2:", round(pcoa_herb_fe$values$Relative_eig[2]*100, 1), "%\n")

# Graphique
plot_ordination(ps_herb_fe, pcoa_herb_fe, color = "Time") +
  theme_classic() +
  scale_color_viridis_d() +
  labs(title = "PCoA - Herbivores - Félicité")

##PErmanova pour verifier 
dist_herb_fe <- phyloseq::distance(ps_herb_fe, method = "bray")
meta_fe <- data.frame(sample_data(ps_herb_fe))

adonis2(dist_herb_fe ~ Time,
        data = meta_fe,
        permutations = 999)

##Les autres guils
guilds_list <- unique(tax_fe$guild_name[!is.na(tax_fe$guild_name)])

for(g in guilds_list) {
  ps_guild <- subset_taxa(ps_actino_felicite, guild_name == g)
  ps_guild <- prune_taxa(taxa_sums(ps_guild) > 0, ps_guild)
  
  if(ntaxa(ps_guild) < 2) next  
  
  dist_guild <- phyloseq::distance(ps_guild, method = "bray")
  meta_g <- data.frame(sample_data(ps_guild))
  
  perm <- adonis2(dist_guild ~ Time, data = meta_g, permutations = 999)
  
  cat("\n===", g, "===\n")
  cat("R² =", round(perm$R2[1], 3),
      "| F =", round(perm$F[1], 3),
      "| p =", perm$`Pr(>F)`[1], "\n")
}


##Picivore 

# Piscivores
ps_pisc_fe <- subset_taxa(ps_actino_felicite, guild_name == "Piscivore")
ps_pisc_fe <- prune_taxa(taxa_sums(ps_pisc_fe) > 0, ps_pisc_fe)
ntaxa(ps_pisc_fe)

dist_pisc_fe <- phyloseq::distance(ps_pisc_fe, method = "bray")
meta_fe <- data.frame(sample_data(ps_pisc_fe))

adonis2(dist_pisc_fe ~ Time, data = meta_fe, permutations = 999)






