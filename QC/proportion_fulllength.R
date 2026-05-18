# libraries
library(dplyr)
library(tidyr)
library(ggplot2)


# Read qc_bambu output files
dRNA_fulllength_counts <- read.delim(
  "dRNA_R10_4M_fullLengthCounts_transcript.txt",
  sep = "\t", stringsAsFactors = FALSE
)

cDNA_fulllength_counts <- read.delim(
  "cDNA_R10_4M_fullLengthCounts_transcript.txt",
  sep = "\t", stringsAsFactors = FALSE
)

dRNA_totalcounts <- read.delim(
  "dRNA_4M_total_reads.txt",
  stringsAsFactors = FALSE
)

cDNA_totalcounts <- read.delim(
  "cDNA_4M_total_reads.txt",
  stringsAsFactors = FALSE
)

# Transpose total counts table 

t_cDNA_totalcounts <- as.data.frame(t(cDNA_totalcounts$Counts), stringsAsFactors = FALSE)
colnames(t_cDNA_totalcounts) <- cDNA_totalcounts$Sample

t_dRNA_totalcounts <- as.data.frame(t(dRNA_totalcounts$Counts), stringsAsFactors = FALSE)
colnames(t_dRNA_totalcounts) <- dRNA_totalcounts$Sample

# Remove annotation columns from full-length count tables

dRNA_fulllength_counts <- dRNA_fulllength_counts[, !names(dRNA_fulllength_counts) %in% c("TXNAME", "GENEID")]
cDNA_fulllength_counts <- cDNA_fulllength_counts[, !names(cDNA_fulllength_counts) %in% c("TXNAME", "GENEID")]


# Calculate column sums of full-length counts table
dRNAfl_colSums <- colSums(dRNA_fulllength_counts, na.rm = TRUE)

# Reorder to match the sample order in t_dRNA_totalcounts
dRNAfl_colSums <- dRNAfl_colSums[ colnames(t_dRNA_totalcounts) ]

cDNAfl_colSums <- colSums(cDNA_fulllength_counts, na.rm = TRUE)
cDNAfl_colSums <- cDNAfl_colSums[ colnames(t_cDNA_totalcounts) ]

# Extract totals (as numeric)
dRNAtotals <- as.numeric(t_dRNA_totalcounts[1, ])
cDNAtotals <- as.numeric(t_cDNA_totalcounts[1, ])

# Build data.frames (ensure lengths match)
if (length(dRNAtotals) != length(dRNAfl_colSums)) stop("Mismatch in dRNA sample counts length!")
if (length(cDNAtotals) != length(cDNAfl_colSums)) stop("Mismatch in cDNA sample counts length!")

dRNA <- data.frame(
  sample = colnames(t_dRNA_totalcounts),
  total_counts = dRNAtotals,
  fulllength_counts = as.numeric(dRNAfl_colSums),
  stringsAsFactors = FALSE
)

cDNA <- data.frame(
  sample = colnames(t_cDNA_totalcounts),
  total_counts = cDNAtotals,
  fulllength_counts = as.numeric(cDNAfl_colSums),
  stringsAsFactors = FALSE
)


# Add protocol info to table

dRNA2 <- dRNA %>%
  mutate(protocol = "dRNA",
         percent_fullength = fulllength_counts / total_counts * 100)

cDNA2 <- cDNA %>%
  mutate(protocol = "cDNA",
         percent_fullength = fulllength_counts / total_counts * 100)

# Combine df together
combined <- bind_rows(dRNA2, cDNA2)

# Create summary stats per protocol
summary_per_protocol <- combined %>%
  group_by(protocol) %>%
  summarise(
    mean_percent = mean(percent_fullength, na.rm = TRUE),
    sd_percent   = sd(percent_fullength, na.rm = TRUE),
    n = n(),
    pooled_percent = sum(fulllength_counts, na.rm = TRUE) / sum(total_counts, na.rm = TRUE) * 100
  ) %>%
  ungroup()

print(summary_per_protocol)



# Creating stacked barplot
pooled <- pooled %>%
  mutate(
    label_y = ifelse(pooled_percent >= 12, pooled_percent / 2, pooled_percent + 4),
    label_colour = ifelse(pooled_percent >= 12, "white", "black")
  )

# Order to be stacked
desired_levels <- c("cDNA_other", "cDNA_full-length", "dRNA_other", "dRNA_full-length")
desired_levels <- desired_levels[desired_levels %in% unique(pooled_long$fill_group)]
pooled_long$fill_group <- factor(pooled_long$fill_group, levels = desired_levels)


# Plot
p1_polished <- ggplot() +
  geom_col(data = pooled_long,
           aes(x = protocol, y = percent, fill = fill_group),
           width = 0.5, color = NA, position = "stack") +
  geom_jitter(data = combined,
              aes(x = protocol, y = percent_fullength),
              width = 0.12, height = 0, size = 2.8,
              shape = 21, fill = "black", color = "black", stroke = 0.35, alpha = 0.95) +
  geom_text(data = pooled,
            aes(x = protocol, y = label_y,
                label = paste0(round(pooled_percent, 1), "%"),
                color = label_colour),
            inherit.aes = FALSE,
            fontface = "bold",
            size = 3.2) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"),
                     limits = c(0, 100),
                     expand = c(0, 0.01)) +
  scale_fill_manual(values = c(
    "cDNA_full-length" = "#2166ac",   # dark blue
    "cDNA_other"       = "#a6c8ff",   # light blue
    "dRNA_full-length" = "#b2182b",   # dark red
    "dRNA_other"       = "#f4a6a6"   # light red
  ),
  labels = c(
    "cDNA_full-length" = "cDNA full-length transcripts",
    "cDNA_other"       = "cDNA other transcripts",
    "dRNA_full-length" = "dRNA full-length transcripts",
    "dRNA_other"       = "dRNA other transcripts"
  ),
  name = NULL) +
  scale_color_identity() +
  labs(x = NULL, y = "Proportion of full-length transcripts (%)") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.45, "lines"),
    legend.spacing.x = unit(0.2, "lines"),
    legend.margin = margin(b = -6, t = 0, unit = "pt"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(face = "bold", size = 10),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 8)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))

# Print
print(p1_polished)

pdf("Proportion_fulllength_transcripts.pdf", width = 85/25.4, height = 60/24)
print(p1_polished)
dev.off()


ggsave(
  "Figure_full_length.pdf",
  plot = p1_polished,
  width = 10.5,
  height = 9.6,
  units = "in",
  device = cairo_pdf
)

