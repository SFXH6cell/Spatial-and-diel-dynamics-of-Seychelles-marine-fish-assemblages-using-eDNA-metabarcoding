##First, lets create metadata on computer


write.table (metadata,"~/Desktop/Research_project/Research_data_seychelles/All_data/Decontam/metadata.tsv",
             sep="\t", quote=FALSE)


###We want an R object with Taxo, Asv count, metadata, and adding a treatment column to metadata : To know which island it was sampled from and the time step 

# Lire le fichier Excel
install.packages("readxl")
library(readxl)
metadata_field <- read_excel("~/Desktop/Research_project/Research_data_seychelles/All_data/eDNA_Sample_Metadata.xlsx")

metadata_field <- metadata_field[, c("Sample", "Island", "Time")]


# Créer une colonne de correspondance dans ton metadata R existant
# Extraire FeT1-1 depuis GC-JJ-1228-PLATE3-FeT1-1


# Créer une clé dans metadata_field (déjà bien formaté)
metadata_field$join_key <- paste0(
  ifelse(grepl("Fé", metadata_field$Island), "Fe",
         ifelse(grepl("Frégate", metadata_field$Island), "Fr", NA)),
  "T",
  gsub(".*_T(\\d+)\\.(\\d+)$", "\\1-\\2", metadata_field$Sample)
)

metadata$join_key <- gsub(".*PLATE.-", "", rownames(metadata))
metadata$join_key <- gsub(".*PLATE2-", "", metadata$join_key)


head(metadata_field$join_key)
head(metadata$join_key)

head(rownames(metadata))

metadata$join_key <- gsub(".*PLATE\\d\\.", "", rownames(metadata))
metadata$join_key <- gsub("\\.", "-", metadata$join_key)


head(metadata$join_key)


# Jointure
metadata_final <- merge(metadata, 
                        metadata_field[, c("join_key", "Island", "Time")],
                        by = "join_key", 
                        all.x = TRUE)
head(metadata_final)
dim(metadata_final)

####Ajout de taxo and ASV count 

head(tax_lulu_clean)


ASV_taxo_clean <- readLines("~/Desktop/Research_project/Research_data_seychelles/All_data/Decontam/ASVlulu_taxo_clean.txt")


head(ASV_taxo_clean , 6)
length(ASV_taxo_clean )


# Séparer le nom ASV et la taxonomie
tax_split <- strsplit(ASV_taxo_clean, "\t")

# Extraire les noms ASV
asv_names <- sapply(tax_split, function(x) x[1])

# Extraire la taxonomie et séparer par ";"
tax_levels <- sapply(tax_split, function(x) x[2])
tax_matrix <- strsplit(tax_levels, ";")

# Mettre en tableau
tax_df <- do.call(rbind, lapply(tax_matrix, function(x) {
  length(x) <- 7  # forcer 7 colonnes (Kingdom à Species)
  return(x)
}))

# Nommer les colonnes et les lignes
colnames(tax_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_df) <- asv_names

# Vérifier
head(tax_df)
dim(tax_df)  # doit donner 833 x 7

# Vérifier que les ASV sont les mêmes dans les 3 objets
all(rownames(asv_count_clean) == rownames(tax_df))  # doit donner TRUE
all(colnames(asv_count_clean) == rownames(metadata_final)) # doit donner TRUE

dim(asv_count_clean)


##Aligner ASV count sur le sASV de taxo

asv_count_taxo <- asv_count_clean[rownames(asv_count_clean) %in% rownames(tax_df), ]
dim(asv_count_clean)
dim(asv_count_taxo)
dim(tax_df)


# ASV dans tax_df mais pas dans asv_count_clean
manquants_dans_count <- setdiff(rownames(tax_df), rownames(asv_count_clean))
length(manquants_dans_count)
manquants_dans_count

# ASV dans asv_count_clean mais pas dans tax_df
manquants_dans_taxo <- setdiff(rownames(asv_count_clean), rownames(tax_df))
length(manquants_dans_taxo)
manquants_dans_taxo

# Regarder les rownames de tax_df
head(rownames(tax_df), 10)

# Checker si y'a des espaces
rownames(tax_df)[1:5]
nchar(rownames(tax_df)[1:5])  # compter les caractères

# Chercher ASV_197 dans tax_df
grep("ASV_197", rownames(tax_df), value=TRUE)

# Chercher ASV_197 dans asv_count_clean
grep("ASV_197", rownames(asv_count_clean), value=TRUE)


# Nettoyer les rownames de tax_df qui ont absorbé la taxonomie
rownames(tax_df) <- gsub("^(ASV_\\d+).*", "\\1", rownames(tax_df))

# Vérifier que ASV_197 est corrigé
grep("ASV_197", rownames(tax_df), value=TRUE)


# Revérifier les dimensions
asv_communs <- intersect(rownames(asv_count_clean), rownames(tax_df))
length(asv_communs)

asv_count_taxo <- asv_count_clean[asv_communs, ]
tax_df_final <- tax_df[asv_communs, ]

dim(asv_count_taxo)   # doit donner 833 x 49
dim(tax_df_final)     # doit donner 833 x 7
all(rownames(asv_count_taxo) == rownames(tax_df_final)) # doit donner TRUE


# Mettre les bons rownames dans metadata_final
rownames(metadata_final) <- metadata_final$SampleID

# Construire l'objet phyloseq avec TOUT (vrais samples + contrôles)
ps_final <- phyloseq(
  otu_table(asv_count_taxo, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_df_final)),
  sample_data(metadata_final)
)

# Vérifier
ps_final


saveRDS(ps_final, "~/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")


readRDS("/Users/maylissaguer/Desktop/Research_project/Research_data_seychelles/All_data/ps_final.rds")
 
