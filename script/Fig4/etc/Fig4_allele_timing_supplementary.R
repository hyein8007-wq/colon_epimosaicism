# FIGURE 4 - SOURCE-PRESERVING EXTRACT (NOT YET STANDALONE)
# Source: original_XCI_methylation.R; original line ranges: 2181-2350
# This extract preserves the original calculations and figure commands.
# It depends on objects from earlier sections and original server inputs.
# DO NOT publish as independently reproducible until dependencies are resolved.


# ===== ORIGINAL LINES 2181-2350 =====
# 3. Supplementary table for timing ============================================

# 3.1 Sheet 1

library(tidyverse)

get_timing_stats <- function(df, xmin = -7, xmax = 30, adjust = 1.0) {
  
  d <- density(df$mutation_point, weights = df$weight / sum(df$weight), adjust = adjust, from = xmin, to = xmax, n = 4096)
  dx <- d$x[2] - d$x[1]
  kde_df <- tibble(
    x = d$x,
    density = d$y,
    mass = d$y * dx
  )
  # Peak
  peak <- d$x[which.max(d$y)]
  
  # 50% highest-density interval
  tmp <- kde_df %>%
    dplyr::arrange(desc(density)) %>%
    dplyr::mutate(cum_mass = cumsum(mass))
  
  cutoff_50 <- tmp %>%
    dplyr::filter(cum_mass >= 0.50) %>%
    dplyr::slice(1) %>%
    dplyr::pull(density)
  
  hdi50 <- kde_df %>%
    dplyr::filter(density >= cutoff_50)
  
  # 95% highest-density interval
  cutoff_95 <- tmp %>%
    dplyr::filter(cum_mass >= 0.95) %>%
    dplyr::slice(1) %>%
    dplyr::pull(density)
  
  hdi95 <- kde_df %>%
    dplyr::filter(density >= cutoff_95)
  
  tibble(
    peak_timing = peak,
    density50_lower = min(hdi50$x),
    density50_upper = max(hdi50$x),
    density95_lower = min(hdi95$x),
    density95_upper = max(hdi95$x)
  )
}

sheet1_donor <- mutation_density_complete %>%
  group_by(sample_id) %>%
  group_modify(~ get_timing_stats(.x)) %>%
  ungroup()

sheet1_overall <- get_timing_stats(
  mutation_density_total
) %>%
  mutate(sample_id = "Overall") %>%
  dplyr::select(sample_id, everything())

sheet1 <- bind_rows(sheet1_donor, sheet1_overall) %>%
  mutate(analysis_type = case_when(sample_id == "HC50" ~ "100 random subsamplings", sample_id == "Overall" ~ "Donor-weighted overall", TRUE ~ "Full phylogeny" )) %>%
  mutate(across( c(peak_timing, density50_lower, density50_upper, density95_lower, density95_upper ), ~ round(.x, 2)))
sheet1

# 3.2 sheet2
sheet2 <- df_shared_length_allele_final %>%
  dplyr::transmute(
    donor = sample,
    clone_1 = clone1,
    clone_2 = clone2,
    shared_somatic_mutations = shared_length,
    XCI_state_clone_1 = allele1,
    XCI_state_clone_2 = allele2,
    XCI_concordance = shared_allele
  ) %>%
  dplyr::arrange(
    donor,
    dplyr::desc(shared_somatic_mutations),
    clone_1,
    clone_2
  )

sheet2

# 3.3 sheet3
sheet3 <- hc50_subtree_complete %>%
  transmute(
    subsampling_id = n_subtree,
    mutation_point = mutation_point,
    probability_weight = weight_norm
  ) %>%
  arrange(
    as.numeric(subsampling_id),
    mutation_point
  )

sheet3

sheet3 %>%
  group_by(subsampling_id) %>%
  summarise(
    total_probability = sum(probability_weight),
    .groups = "drop"
  )

# HC50: mean of 100 subsamplings
hc50_from_subsampling <- sheet3 %>%
  dplyr::group_by(mutation_point) %>%
  dplyr::summarise(
    mean_probability = mean(probability_weight, na.rm = TRUE),
    .groups = "drop"
  )


# HC50: final donor-level distribution used in the analysis
hc50_final <- mutation_density_complete %>%
  dplyr::filter(sample_id == "HC50") %>%
  dplyr::select(
    mutation_point,
    HC50_probability = weight
  )


# Compare
hc50_check <- hc50_from_subsampling %>%
  dplyr::left_join(
    hc50_final,
    by = "mutation_point"
  ) %>%
  dplyr::mutate(
    difference = mean_probability - HC50_probability
  )

hc50_check

# 3.4 sheet4
sheet4 <- mutation_density_complete %>%
  dplyr::select(
    sample_id,
    mutation_point,
    probability_weight = weight
  ) %>%
  dplyr::arrange(sample_id, mutation_point)

sheet4 <- dplyr::bind_rows(
  
  mutation_density_complete %>%
    dplyr::transmute(
      sample_id,
      mutation_point,
      probability_weight = weight
    ),
  
  mutation_density_total %>%
    dplyr::transmute(
      sample_id = "Overall",
      mutation_point,
      probability_weight = weight
    )
  
) %>%
  dplyr::arrange(sample_id, mutation_point)

# save results
#sheet1 %>% write_tsv("~team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/res_final/XCI_timing_summary.tsv")
#sheet2 %>% write_tsv("~team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/res_final/XCI_timing_pairwise_XCI.tsv")
#sheet3 %>% write_tsv("~team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/res_final/XCI_timing_HC50_subsampling.tsv")
#sheet4 %>% write_tsv("~team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/res_final/XCI_timing_timing_probability.tsv")

