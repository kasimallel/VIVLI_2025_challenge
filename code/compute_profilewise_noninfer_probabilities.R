compute_profilewise_noninfer_probabilities <- function(final_output, reference = "MEM", delta = 0.20) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(binom)
  
  profiles <- unique(final_output$ProfileID)
  antibiotics <- unique(final_output$Antibiotic)
  antibiotics <- antibiotics[antibiotics != reference]
  
  results <- map_dfr(profiles, function(pid) {
    df_profile <- final_output %>%
      filter(ProfileID == pid) %>%
      dplyr::select(Antibiotic, ProfileID, ProfileLabel, Draw, Suscept)
    
    if (!all(c(reference, antibiotics) %in% df_profile$Antibiotic)) return(NULL)
    
    # Reshape to wide format
    wide <- df_profile %>%
      pivot_wider(names_from = Antibiotic, values_from = Suscept) %>%
      arrange(Draw)
    
    if (!reference %in% colnames(wide)) return(NULL)
    
    ref_draws <- wide[[reference]]
    
    map_dfr(antibiotics, function(abx) {
      if (!abx %in% colnames(wide)) return(NULL)
      
      diff <- wide[[abx]] - ref_draws
      in_margin <- as.integer(diff >= -delta)
      
      p_hat <- mean(in_margin, na.rm = TRUE)
      bin_ci <- binom.confint(sum(in_margin, na.rm = TRUE), length(in_margin), method = "wilson")
      
      tibble(
        ProfileID = pid,
        ProfileLabel = df_profile$ProfileLabel[1],
        Antibiotic = abx,
        `P(Non-inferior)` = p_hat,
        `95% CI Lower` = bin_ci$lower,
        `95% CI Upper` = bin_ci$upper
      )
    })
  })
  
  return(results)
}

