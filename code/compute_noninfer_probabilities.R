compute_noninfer_probabilities <- function(suscept_df, antibiotics, reference = "MEM", delta = 0.20) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(binom)
  
  # Add draw ID
  suscept_df <- suscept_df %>%
    group_by(Antibiotic) %>%
    mutate(draw = row_number()) %>%
    ungroup()
  
  # Pivot to wide format
  suscept_mat <- suscept_df %>%
    pivot_wider(names_from = Antibiotic, values_from = Suscept) %>%
    arrange(draw) %>%
    dplyr::select(-draw) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.data.frame()
  
  # Reference draws
  ref_draws <- suscept_mat[[reference]]
  
  # Calculate non-inferiority probabilities
  results <- antibiotics[antibiotics != reference] %>%
    map_df(function(abx) {
      diff <- suscept_mat[[abx]] - ref_draws
      in_margin <- as.integer(diff >= -delta)
      
      p_hat <- mean(in_margin)
      bin_ci <- binom.confint(sum(in_margin), length(in_margin), method = "wilson")
      
      tibble(
        Antibiotic = abx,
        `P(Non-inferior to MEM)` = p_hat,
        `95% CI Lower` = bin_ci$lower,
        `95% CI Upper` = bin_ci$upper
      )
    }) %>%
    arrange(desc(`P(Non-inferior to MEM)`))
  
  return(results)
}


