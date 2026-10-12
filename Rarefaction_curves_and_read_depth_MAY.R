# Rarefaction curves — do samples reach saturation?

###### CHECK READ DEPTHS

read_depths <- sort(sample_sums(ps_final))
print(read_depths)
summary(read_depths)

# Ratio of max to min — if >> 10 rarefying is strongly recommended
cat("Max:min read depth ratio:", round(max(read_depths) / min(read_depths), 1), "\n")

# Visualise read depth per sample, coloured by treatment
depth_df <- data.frame(sample    = names(read_depths), depth = read_depths,
                       treatment = sample_data(ps_final)$sample_type[match(names(read_depths), rownames(sample_data(ps_final)))]
)

ggplot(depth_df, aes(x = reorder(sample, depth), y = depth, fill = treatment)) +
  geom_col() + coord_flip() +
  labs(title = "Read depth per sample", x = NULL, y = "Read count", fill = "Treatment") +
  theme_bw()


# Rarefaction curves — do samples reach saturation?
# If curves plateau, your sequencing depth is sufficient
# Force to a plain base R matrix
otu_mat <- matrix(as.numeric(otu_table(ps_final)), nrow = nsamples(ps_final),
                  ncol = ntaxa(ps_final), dimnames = list(sample_names(ps_final), 
                                                        taxa_names(ps_final))
)

# Check dimensions look right — should be samples x taxa
dim(otu_mat)

install.packages("vegan")
library(vegan)

rarecurve(otu_mat, step  = 100, col   = "steelblue", label = FALSE,
          main  = "Rarefaction Curves", ylab  = "ASVs", xlab  = "Read Count",
          xaxt  = "n")

# Add custom x axis with regular numbers
axis(1, at = seq(0, max(rowSums(otu_mat)), by = 50000),
     labels = format(seq(0, max(rowSums(otu_mat)), by = 50000), 
                     scientific = FALSE, big.mark = ","))




