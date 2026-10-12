# FIRST STEP : Create a metadata file
### We take our ASV count Table and we want another table thats sayys 
### Those Samples are Control if they have “Ex8” or Ex7” or “FB2” or “FB10” or “FB11” or “Pc5” or “Pc6” in their name 

##asv_count_lulu

# STEP 2 — Create metadata table
# Les noms des samples sont les noms des colonnes
sample_names <- colnames(asv_count_lulu)
control_pattern <- "Ex8|Ex7|FB2|FB10|FB11|Pc5|Pc6"

metadata <- data.frame(
  SampleID = sample_names,
  is.neg = grepl(control_pattern, sample_names),
  sample_type = ifelse(
    grepl(control_pattern, sample_names),
    "Control",
    "Real_sample"
  ))

rownames(metadata) <- metadata$SampleID

metadata
table(metadata$sample_type)

real_samples <- metadata$SampleID[metadata$is.neg == FALSE]

real_samples[!grepl("FrT|FeT", real_samples)]

write.table(
  metadata,
  file = "metadata_decontam_lulu.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



##Second step 

seqtab2 <- t(asv_count_lulu)

df <- metadata
df$LibrarySize <- rowSums(seqtab2)
df$Index <- seq_len(nrow(df))


ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_type)) +
  geom_point() +
  xlab("Sample Index") +
  ylab("Library Size")



# 4. Run decontam — prevalence method


### First we run it with default setting 

ps <- phyloseq(
  otu_table(asv_count_lulu, taxa_are_rows = TRUE),
  sample_data(metadata)
)


ps

sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "Control"

contamdf.prev2 <- isContaminant(ps, method = "prevalence", neg = "is.neg")

table(contamdf.prev2$contaminant)


###NNOW avec 0,5 threshold
contamdf.prev05bis <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)

table(contamdf.prev05bis$contaminant)

ps.pa2 <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))

ps.pa.neg2 <- prune_samples(sample_data(ps.pa2)$sample_type == "Control", ps.pa2)
ps.pa.pos2 <- prune_samples(sample_data(ps.pa2)$sample_type == "Real_sample", ps.pa2)


df.pa2 <- data.frame(
  pa.pos2 = taxa_sums(ps.pa.pos2),
  pa.neg2 = taxa_sums(ps.pa.neg2),
  contaminant = contamdf.prev05bis$contaminant
)

ggplot(data = df.pa2, aes(x = pa.neg2, y = pa.pos2, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)")

##### PLOT PART #####

control_samples <- metadata$SampleID[metadata$is.neg == TRUE]

# Extraire ces colonnes
controls_only <- asv_count_lulu[, control_samples]

# Total de reads par contrôle
reads_per_control <- data.frame(
  Sample = colnames(controls_only),
  TotalReads = colSums(controls_only)
)

# Bar chart
ggplot(reads_per_control, aes(x = Sample, y = TotalReads, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = TotalReads), vjust = -0.5, size = 3.5) +
  xlab("Negative Controls") +
  ylab("Total Reads") +
  ggtitle("Total reads per negative control") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )




##12/10###Check which species are in there : 


asv_count_lulu <- read.table(
  "~/Desktop/Research_project/Research_data_seychelles/All_data/Dada2output/Output2_ASVcounts.tsv",
  header = TRUE, sep = "\t", row.names = 1
)

# Noms des blanks d'après ton graphique
control_samples <- c("GC.JJ.1228.PLATE3.Ex7", "GC.JJ.1228.PLATE3.Ex8",
                     "GC.JJ.1228.PLATE3.FB10", "GC.JJ.1228.PLATE3.FB11",
                     "GC.JJ.1228.PLATE3.FB2", "GC.JJ.1228.PLATE3.Pc5",
                     "GC.JJ.1228.PLATE3.Pc6")

controls_only <- asv_count_lulu[, control_samples]
asv_in_controls <- controls_only[rowSums(controls_only) > 0, ]
cat("ASVs dans les blanks:", nrow(asv_in_controls), "\n")
print(asv_in_controls)


# Charger la table taxonomique
taxo_lulu <- read.table(
  "~/Desktop/Research_project/Research_data_seychelles/All_data/basta_output/ASVlulu_taxo_copie.txt",
  header = TRUE, sep = "\t"
)

# Vérifier la structure
head(tax_df)
colnames(tax_df)

# Croiser avec la taxonomie
species_in_controls <- tax_df[rownames(asv_in_controls), c("Species", "Family", "Class")]
species_in_controls$TotalReads <- rowSums(asv_in_controls)

print(species_in_controls[order(species_in_controls$TotalReads, decreasing = TRUE), ])

# Tableau long : pour chaque blank, quelles espèces et combien de reads




# Matrice présence/absence avec espèces
asv_in_controls$Species <- tax_df[rownames(asv_in_controls), "Species"]

aggregate(Species ~ ., data = asv_in_controls)

for(s in colnames(asv_in_controls)[-ncol(asv_in_controls)]) {
  cat("\n===", s, "===\n")
  cat(asv_in_controls$Species[asv_in_controls[,s] > 0], sep="\n")
}


# Créer un dataframe propre
result <- data.frame(
  Sample = rep(colnames(asv_in_controls)[-ncol(asv_in_controls)], 
               sapply(colnames(asv_in_controls)[-ncol(asv_in_controls)], 
                      function(s) sum(asv_in_controls[,s] > 0))),
  Species = unlist(lapply(colnames(asv_in_controls)[-ncol(asv_in_controls)], 
                          function(s) asv_in_controls$Species[asv_in_controls[,s] > 0]))
)

print(result)
