run_antibiotic_model <- function(data_input, antibiotics,
                                 model_path = "multivariate_probit_model.stan",
                                 ref_abx = "MEM") {
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(abind)
  library(MetBrewer)
  library(ggplot2)
  
  # === Detect dataset name as string ===
  dataset_name <- deparse(substitute(data_input))
  
  # === Error check ===
  if (!(ref_abx %in% antibiotics)) {
    stop(paste0("Reference antibiotic '", ref_abx, "' is not in the antibiotic list."))
  }
  
  set.seed(123)
  rand_tag <- paste0("_", sample(1000:9999, 1))
  
  # === Step 1: Prepare Data ===
  data <- data_input %>%
    dplyr::select(all_of(antibiotics), Gender, Agegroup, Department_group) %>%
    na.omit() %>%
    mutate(across(c(Gender, Agegroup, Department_group), as.factor))
  
  X <- model.matrix(~ Gender + Agegroup + Department_group, data = data)
  Y <- as.matrix(data[, antibiotics])
  N <- nrow(data); D <- ncol(Y); P <- ncol(X)
  
  # === Step 2: Fit Model ===
  mod <- cmdstan_model(model_path)
  fit <- mod$sample(
    data = list(N = N, D = D, P = P, X = X, Y = Y),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000,
    adapt_delta = 0.95, max_treedepth = 15, seed = 123
  )
  
  draws_df <- as_draws_df(fit$draws())
  
  # === Step 3: Posterior vs Observed ===
  observed_df <- data_input %>%
    dplyr::select(all_of(antibiotics)) %>%
    summarise(across(everything(), ~ mean(. == 0, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Antibiotic", values_to = "Observed")
  
  suscept_df <- purrr::map_dfr(1:length(antibiotics), function(i) {
    tibble(
      Antibiotic = antibiotics[i],
      Suscept = 1 - pnorm(draws_df[[paste0("beta[1,", i, "]")]])
    )
  })
  
  estimated_df <- suscept_df %>%
    group_by(Antibiotic) %>%
    summarise(
      Median = mean(Suscept),
      Lower = quantile(Suscept, 0.025),
      Upper = quantile(Suscept, 0.975)
    )
  
  comparison_df <- left_join(estimated_df, observed_df, by = "Antibiotic") %>%
    arrange(Median) %>%
    mutate(Antibiotic = factor(Antibiotic, levels = Antibiotic))
  
  p1 <- ggplot(comparison_df) +
    geom_errorbarh(
      aes(xmin = Lower * 100, xmax = Upper * 100, y = Antibiotic),
      height = 0.25,
      color = "black"
    ) +
    geom_point(
      aes(x = Median * 100, y = Antibiotic, shape = "Posterior", fill = "Posterior"),
      color = "black",
      size = 3
    ) +
    geom_point(
      aes(x = Observed * 100, y = Antibiotic, shape = "Observed", fill = "Observed"),
      color = "black",
      size = 2.5
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("Posterior" = 21, "Observed" = 8)
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("Posterior" = "gold", "Observed" = "black")
    ) +
    labs(
      x = "Posterior probability susceptibility (%)",
      y = NULL,
      title = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 12)
    )
  
  
  ggsave(paste0(dataset_name, "_comparison", rand_tag, ".png"), p1, width = 7, height = 5)
  
  # === Step 4: Violin Plot ===
  intercepts <- draws_df %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    filter(str_detect(param, "^beta\\[1,")) %>%
    mutate(index = as.integer(str_extract(param, "(?<=,)[0-9]+")),
           Antibiotic = antibiotics[index],
           prob_suscept = 1 - pnorm(value))

  
  mem_median <- intercepts %>%
    filter(Antibiotic == ref_abx) %>%
    summarise(median_prob = median(prob_suscept)) %>%
    pull()
  
  ni_threshold <- mem_median - 0.20
  ordered_abx <- intercepts %>%
    group_by(Antibiotic) %>%
    summarise(med = median(prob_suscept)) %>%
    arrange(desc(med)) %>%
    pull(Antibiotic)
  
  intercepts$Antibiotic <- factor(intercepts$Antibiotic, levels = ordered_abx)
  palette <- MetBrewer::met.brewer("Hiroshige", length(ordered_abx))
  
  p2 <- ggplot(intercepts, aes(x = prob_suscept * 100, y = Antibiotic, fill = Antibiotic)) +
    geom_violin(alpha = 0.6, trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
    geom_vline(xintercept = ni_threshold * 100, linetype = "dashed", color = "red", linewidth = 0.6) +
    annotate("text",
             x = ni_threshold * 100 + 1.5,
             y = length(ordered_abx) + 0.45,
             label = "Non-inferiority area",
             hjust = 0, color = "red", size = 3, fontface = "italic") +
    scale_fill_manual(values = palette) +
    scale_x_continuous(limits = c(0, 100)) +  # <-- THIS ENFORCES HARD LIMIT
    labs(x = "Posterior probability susceptibility (%)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold")
    )
  
  ggsave(paste0(dataset_name, "_violin", rand_tag, ".png"), p2, width = 7, height = 5)
  
  # === Step 5: Individual Risk Profiles ===
  covariate_grid <- expand.grid(
    Gender = levels(data$Gender),
    Agegroup = levels(data$Agegroup),
    Department_group = levels(data$Department_group)
  )
  
  X_new <- model.matrix(~ Gender + Agegroup + Department_group, data = covariate_grid)[, -1]
  draws <- fit$draws(variables = c("alpha", "beta", "L"), format = "draws_matrix")
  
  n_draws <- nrow(draws)
  alpha_draws <- draws[, paste0("alpha[", 1:D, "]")]
  beta_draws <- array(NA, dim = c(n_draws, ncol(X_new), D))
  for (d in 1:D) {
    for (p in 1:ncol(X_new)) {
      beta_draws[, p, d] <- draws[, paste0("beta[", p, ",", d, "]")]
    }
  }
  
  get_L_matrix <- function(row) {
    idxs <- grep("^L\\[", colnames(draws))
    dim_L <- sqrt(length(idxs))
    matrix(draws[row, idxs], ncol = dim_L, byrow = TRUE)
  }
  
  n_profiles <- nrow(X_new)
  pred_probs_array <- array(NA, dim = c(n_draws, n_profiles, D))
  for (s in 1:n_draws) {
    alpha_s <- alpha_draws[s, ]
    beta_s <- beta_draws[s, , ]
    L_s <- get_L_matrix(s)
    mu <- sweep(X_new %*% beta_s, 2, alpha_s, "+")
    Z_noise <- matrix(rnorm(n_profiles * D), n_profiles, D) %*% t(L_s)
    Z_sim <- mu + Z_noise
    pred_probs_array[s, , ] <- 1 - pnorm(Z_sim)
  }
  
  pred_summary <- apply(pred_probs_array, c(2, 3), function(x) {
    c(mean = mean(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  })
  
  pred_summary_long <- do.call(rbind, lapply(1:n_profiles, function(i) {
    tibble(
      ProfileID = i,
      Antibiotic = antibiotics,
      Mean = pred_summary[1, i, ],
      Lwr = pred_summary[2, i, ],
      Upr = pred_summary[3, i, ]
    )
  }))
  
  final_output <- pred_summary_long %>%
    left_join(covariate_grid %>% mutate(ProfileID = row_number()), by = "ProfileID") %>%
    mutate(ProfileLabel = paste(Gender, Agegroup, Department_group, sep = " | "))
  
  # === Step 6: Faceted Profile Plot ===
  ref_thresholds <- final_output %>%
    filter(Antibiotic == ref_abx) %>%
    dplyr::select(ProfileID, RefMean = Mean) %>%
    mutate(NI_threshold = pmax(0, RefMean - 0.20))
  
  final_output_plot <- final_output %>%
    left_join(ref_thresholds, by = "ProfileID")
  
  ordered_abx <- final_output_plot %>%
    group_by(Antibiotic) %>%
    summarise(avg = mean(Mean, na.rm = TRUE)) %>%
    arrange(desc(avg)) %>%
    pull(Antibiotic)
  
  final_output_plot$Antibiotic <- factor(final_output_plot$Antibiotic, levels = ordered_abx)
  
  hline_data <- final_output_plot %>%
    distinct(ProfileID, ProfileLabel, NI_threshold) %>%
    mutate(x_start = 0.5, x_end = length(ordered_abx) + 0.5)
  
  p_profiles <- ggplot(final_output_plot, aes(x = Antibiotic, y = Mean, ymin = Lwr, ymax = Upr)) +
    geom_pointrange(size = 0.3, fatten = 1.5) +
    geom_segment(data = hline_data,
                 aes(x = x_start, xend = x_end,
                     y = NI_threshold, yend = NI_threshold),
                 linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~ ProfileLabel, ncol = 2) +
    labs(y = "Predicted susceptibility probability", x = NULL) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(size = 10),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(dataset_name, "_profileplot", rand_tag, ".png"), p_profiles, width = 10, height = 12)
  
  # === Step 7: Save data tables ===
  write_csv(final_output, paste0(dataset_name, "_profiles", rand_tag, ".csv"))
  write_csv(estimated_df, paste0(dataset_name, "_summary", rand_tag, ".csv"))
  
  return(list(
    fit = fit,
    posterior_summary = estimated_df,
    profile_predictions = final_output
  ))
}

