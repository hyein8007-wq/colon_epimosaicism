# FIGURE 4 - SOURCE-PRESERVING EXTRACT (NOT YET STANDALONE)
# Source: original_XCI_methylation.R; original line ranges: 286-375
# This extract preserves the original calculations and figure commands.
# It depends on objects from earlier sections and original server inputs.
# DO NOT publish as independently reproducible until dependencies are resolved.


# ===== ORIGINAL LINES 286-375 =====
###############################################################################

# Figure hypomethylated region enrichment

###############################################################################

# 1. Hypomethylated region 

# Sankey plot
library(ggalluvial)
df_transition <- cpg_x %>%
  dplyr::filter(par != "PAR" & ((male_WT == 0 & female_WT == 0) | (male_WT == 1 & female_WT == 0.5)) ) %>%
  mutate(type = paste0(male_WT, "_", female_WT)) %>%
  dplyr::select(type, genomic_region, cpg_context) %>%
  mutate(
    type = factor(type, levels = c("1_0.5", "0_0")),
    genomic_region = factor(genomic_region, levels = c("promoter", "gene_body", "intergenic")),
    cpg_context = factor(cpg_context, levels = c("Island", "Shore", "Shelf", "OpenSea"))
  ) %>%
  group_by(type, genomic_region, cpg_context) %>%
  summarise(Freq = n(), .groups = "drop")

g1 <- ggplot(df_transition %>% dplyr::filter(type == "1_0.5") %>% mutate(total = sum(Freq), prop = Freq / total, label = sprintf("%.1f%%", prop * 100)),
             aes(axis1 = genomic_region, axis2 = cpg_context, y = Freq)) +
  geom_alluvium(aes(fill = cpg_context), width = 1/6) +
  geom_stratum(width = 1/5, fill = "grey85", color = "black") +
  geom_text(stat = "alluvium", aes(label = label), size = 3, color = "black") +
  scale_fill_manual(values = c("Island" = "#1c284c", "Shore" = "#6f91bb", "Shelf"= "#bed2ec", "OpenSea" = "#e9e8e7")) +
  scale_x_discrete(
    limits = c("genomic_region", "cpg_context"),
    expand = c(.1, .1)
  ) +
  labs(y = "Count", fill = "cpg_context") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) + coord_flip()

g2 <- ggplot(df_transition %>% dplyr::filter(type == "0_0") %>% mutate(total = sum(Freq), prop = Freq / total, label = sprintf("%.1f%%", prop * 100)),
             aes(axis1 = genomic_region, axis2 = cpg_context, y = Freq)) +
  geom_alluvium(aes(fill = cpg_context), width = 1/6) +
  geom_stratum(width = 1/5, fill = "grey85", color = "black") +
  geom_text(stat = "alluvium", aes(label = label), size = 3, color = "black") +
  scale_fill_manual(values = c("Island" = "#1c284c", "Shore" = "#6f91bb", "Shelf"= "#bed2ec", "OpenSea" = "#e9e8e7")) +
  scale_x_discrete(
    limits = c("genomic_region", "cpg_context"),
    expand = c(.1, .1)
  ) +
  labs(y = "Count", fill = "cpg_context") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) + coord_flip()


tmp1 <- cpg_x %>% dplyr::filter(par != "PAR" & male_WT == 0 & female_WT == 0) %>% .$cpg_context %>% table %>% as.data.frame() %>% dplyr::rename("cpg_context" = ".") 
tmp2 <- cpg_x %>% dplyr::filter(par != "PAR") %>% .$cpg_context %>% table %>% as.data.frame() %>% dplyr::rename("cpg_context" = ".", "total" = "Freq") 
tmp3 <- left_join(tmp1, tmp2) %>% mutate(proportion = Freq/total)
tmp3$obs_frac <- tmp3$Freq / sum(tmp3$Freq)
tmp3$bg_frac <- tmp3$total / sum(tmp3$total)
tmp3$enrichment <- tmp3$obs_frac / tmp3$bg_frac
tmp3$log2_enrichment <- log2(tmp3$enrichment)
tmp3$se <- sqrt(tmp3$proportion * (1 - tmp3$proportion) / tmp3$total)
tmp3$cpg_context <- factor(tmp3$cpg_context,
                           levels = c("OpenSea", "Shelf", "Shore", "Island"))

gg <- ggplot(tmp3, aes(x = cpg_context, y = log2_enrichment, fill = cpg_context)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    y = "log2 enrichment (Observed / Expected)",
    x = "CpG context"
  ) +
  geom_errorbar(aes(
    ymin = log2((obs_frac - se) / bg_frac),
    ymax = log2((obs_frac + se) / bg_frac)
  )) + coord_flip()

pdf("output/demethylation_cpg_context_enrichment.pdf")
print(gg)
dev.off()

################################################################################
