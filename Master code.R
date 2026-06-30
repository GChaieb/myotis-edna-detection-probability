

                                        
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(patchwork)

bunker_qpcr <- read.csv(
  "bunker_qpcr.csv",
  header = TRUE
)

hibernacula_qpcr <- read.csv(
  "hibernacula_qpcr.csv",
  header = TRUE
)

maternityroost_qpcr <- read.csv(
  "maternityroost_qpcr.csv",
  header = TRUE
)


                               ### 1- CONDITIONAL HIERARCHICAL DETECTION MODEL ###


# PREPARING qPCR COUNTS


prepare_counts <- function(dat) {
  
  pcr_cols <- grep("^pcr", names(dat), value = TRUE)
  
  dat_clean <- dat %>%
    select(all_of(pcr_cols))
  
  dat_clean[pcr_cols] <- lapply(dat_clean[pcr_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  y <- rowSums(dat_clean, na.rm = TRUE)
  K <- rowSums(!is.na(dat_clean))
  
  out <- data.frame(
    sample_id = seq_along(y),
    y = y,
    K = K
  )
  
  return(out)
}


# CONDITIONAL MODEL BY GRID POSTERIOR


fit_conditional_detection <- function(count_data,
                                      system_name,
                                      grid_size = 1000) {
  
  theta_grid <- seq(0.001, 0.999, length.out = grid_size)
  p_grid     <- seq(0.001, 0.999, length.out = grid_size)
  
  log_post <- matrix(NA, nrow = grid_size, ncol = grid_size)
  
  y <- count_data$y
  K <- count_data$K
  
  for (i in seq_along(theta_grid)) {
    
    theta <- theta_grid[i]
    
    for (j in seq_along(p_grid)) {
      
      p <- p_grid[j]
      
      log_lik <- 0
      
      for (s in seq_along(y)) {
        
        if (y[s] == 0) {
          
          prob_y <- (1 - theta) + theta * dbinom(0, size = K[s], prob = p)
          
        } else {
          
          prob_y <- theta * dbinom(y[s], size = K[s], prob = p)
        }
        
        log_lik <- log_lik + log(prob_y)
      }
      
      log_post[i, j] <- log_lik
    }
  }
  
  log_post <- log_post - max(log_post)
  post <- exp(log_post)
  post <- post / sum(post)
  
  theta_marginal <- rowSums(post)
  p_marginal     <- colSums(post)
  
  weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    approx(cw, x, xout = probs, rule = 2)$y
  }
  
  theta_q <- weighted_quantile(theta_grid, theta_marginal)
  p_q     <- weighted_quantile(p_grid, p_marginal)
  
  result <- data.frame(
    system = system_name,
    theta_median = theta_q[2],
    theta_lower = theta_q[1],
    theta_upper = theta_q[3],
    p_median = p_q[2],
    p_lower = p_q[1],
    p_upper = p_q[3]
  )
  
  posterior_grid <- expand.grid(
    theta = theta_grid,
    p = p_grid
  )
  
  posterior_grid$posterior <- as.vector(post)
  posterior_grid$system <- system_name
  
  return(list(
    result = result,
    posterior_grid = posterior_grid,
    count_data = count_data
  ))
}


#  BUNKER DATA


bunker_counts <- prepare_counts(bunker_qpcr)


#  HIBERNACULA / WINTER DATA


winter_counts <- prepare_counts(hibernacula_qpcr)


#   BUNKER MODEL


bunker_fit <- fit_conditional_detection(
  count_data = bunker_counts,
  system_name = "Bunker",
  grid_size = 1000
)


#  HIBERNACULA / WINTER MODEL


winter_fit <- fit_conditional_detection(
  count_data = winter_counts,
  system_name = "Hibernacula",
  grid_size = 1000
)


# FINAL MODELLED THETA AND p TABLE


theta_p_modelled <- bind_rows(
  bunker_fit$result,
  winter_fit$result
)

theta_p_modelled_round <- theta_p_modelled %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(theta_p_modelled_round)


# MATERNITY ROOST DIRECT ESTIMATES


maternity_pcr_cols <- grep("^pcr", names(maternityroost_qpcr), value = TRUE)

maternity_pcr <- maternityroost_qpcr %>%
  select(all_of(maternity_pcr_cols))

maternity_pcr[maternity_pcr_cols] <- lapply(maternity_pcr[maternity_pcr_cols], function(x) {
  as.numeric(as.character(x))
})

maternity_sample_positive <- rowSums(maternity_pcr, na.rm = TRUE) >= 2

maternity_theta <- mean(maternity_sample_positive)
maternity_p <- sum(as.matrix(maternity_pcr), na.rm = TRUE) / sum(!is.na(as.matrix(maternity_pcr)))

maternity_results <- data.frame(
  system = "Maternity roost",
  theta_median = maternity_theta,
  theta_lower = NA,
  theta_upper = NA,
  p_median = maternity_p,
  p_lower = NA,
  p_upper = NA
)


# FINAL TABLE 


theta_p_final <- bind_rows(
  maternity_results,
  theta_p_modelled
)

theta_p_final_round <- theta_p_final %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(theta_p_final_round)




                                              ### 2- SENSITIVITY ANALYSIS  ###


# PANEL A: theta and p comparison

theta_p_final <- data.frame(
  system = c("Maternity roost", "Bunker", "Hibernacula"),
  
  theta = c(1.000, 0.226, 0.184),
  theta_lower = c(0.664, 0.113, 0.045),
  theta_upper = c(1.000, 0.383, 0.441),
  
  p = c(0.972, 0.293, 0.365),
  p_lower = c(0.903, 0.179, 0.147),
  p_upper = c(0.997, 0.419, 0.609)
)

theta_p_final <- theta_p_final %>%
  mutate(
    system = factor(
      system,
      levels = c("Maternity roost", "Bunker", "Hibernacula")
    )
  )

panel_a_data <- theta_p_final %>%
  pivot_longer(
    cols = c(theta, p),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    lower = ifelse(parameter == "theta", theta_lower, p_lower),
    upper = ifelse(parameter == "theta", theta_upper, p_upper),
    
    parameter = factor(parameter, levels = c("p", "theta")),
    
    x_base = as.numeric(system),
    x_pos = ifelse(parameter == "p", x_base - 0.06, x_base + 0.06),
    
    trial_label = recode(
      as.character(system),
      "Maternity roost" = "High density trial",
      "Bunker" = "Low density trial",
      "Hibernacula" = "Winter trial"
    ),
    
    trial_label = factor(
      trial_label,
      levels = c("High density trial", "Low density trial", "Winter trial")
    )
  )

system_cols <- c(
  "High density trial" = "#0072F0",
  "Low density trial" = "#00C853",
  "Winter trial" = "#E83ED6"
)

# PLOT PANEL A
panel_a <- ggplot(panel_a_data) +
  
  geom_errorbar(
    aes(
      x = x_pos,
      ymin = lower,
      ymax = upper,
      color = trial_label
    ),
    width = 0.06,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      x = x_pos,
      y = estimate,
      shape = parameter,
      color = trial_label
    ),
    size = 4.2
  ) +
  
  scale_x_continuous(
    breaks = 1:3,
    labels = c("High density trial", "Low density trial", "Winter trial")
  ) +
  
  scale_shape_manual(
    values = c(
      "p" = 17,
      "theta" = 16
    ),
    labels = c(
      expression(italic(p)),
      expression(theta)
    )
  ) +
  
  scale_color_manual(
    values = system_cols,
    labels = c(
      "High density trial",
      "Low density trial",
      "Winter trial"
    )
  ) +
  
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        color = "black",
        size = 4
      )
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 16,
        size = 4
      )
    )
  ) +
  
  coord_cartesian(ylim = c(0, 1.05)) +
  
  labs(
    x = NULL,
    y = NULL,
    shape = NULL,
    color = NULL
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    
    legend.position = c(0.17, 0.39),
    legend.background = element_blank(),
    legend.box = "vertical",
    legend.text = element_text(size = 15),
    legend.spacing.y = unit(0.10, "cm")
  )


print(panel_a)

# PANEL B: qPCR replicate-level sensitivity


#  THETA AND p VALUES


theta_p_final <- data.frame(
  system = c("Maternity roost", "Bunker", "Hibernacula"),
  theta  = c(1.000, 0.226, 0.184),
  p      = c(0.972, 0.293, 0.365)
)

theta_p_final <- theta_p_final %>%
  mutate(
    system = factor(
      system,
      levels = c("Maternity roost", "Bunker", "Hibernacula")
    )
  )


qpcr_sensitivity <- expand.grid(
  system = levels(theta_p_final$system),
  k = 2:8
) %>%
  left_join(theta_p_final, by = "system") %>%
  mutate(
    sensitivity = 1 - (1 - p)^k - k * p * (1 - p)^(k - 1),
    system = factor(
      system,
      levels = c("Maternity roost", "Bunker", "Hibernacula")
    )
  )


# PLOT PANEL B

system_cols <- c(
  "Maternity roost" = "#0072F0",
  "Bunker" = "#00C853",
  "Hibernacula" = "#E83ED6"
)

panel_b <- ggplot(
  qpcr_sensitivity,
  aes(
    x = k,
    y = sensitivity,
    color = system,
    group = system
  )
) +
  geom_line(linewidth = 1.4) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = system_cols,
    breaks = c("Maternity roost", "Bunker", "Hibernacula"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = 2:8,
    limits = c(2, 8)
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    x = "Number of qPCR replicates",
    y = "Detection probability"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )



print(panel_b)

# PANEL C: air-sample / survey-level sensitivity



# THETA AND p VALUES


theta_p_final <- data.frame(
  system = c("Maternity roost", "Bunker", "Hibernacula"),
  theta  = c(1.000, 0.226, 0.184),
  p      = c(0.972, 0.293, 0.365)
)

theta_p_final <- theta_p_final %>%
  mutate(
    system = factor(
      system,
      levels = c("Maternity roost", "Bunker", "Hibernacula")
    )
  )


na_fixed <- 8

survey_sensitivity <- expand.grid(
  system = levels(theta_p_final$system),
  ns = 1:40
) %>%
  left_join(theta_p_final, by = "system") %>%
  mutate(
    qpcr_ge2 = 1 - (1 - p)^na_fixed -
      na_fixed * p * (1 - p)^(na_fixed - 1),
    
    sensitivity = 1 - (1 - theta * qpcr_ge2)^ns,
    
    system = factor(
      system,
      levels = c("Maternity roost", "Bunker", "Hibernacula")
    )
  )

# PLOT PANEL C

panel_c <- ggplot(
  survey_sensitivity,
  aes(
    x = ns,
    y = sensitivity,
    color = system,
    group = system
  )
) +
  geom_line(linewidth = 1.4) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = system_cols,
    breaks = c("Maternity roost", "Bunker", "Hibernacula"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq(0, 40, 5),
    limits = c(1, 40)
  ) +
  coord_cartesian(ylim = c(0, 1.05))  +
  
  labs(
    x = "Number of air samples",
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

panel_c 
###



add_panel_label <- function(plot, label) {
  plot +
    annotation_custom(
      grob = grobTree(
        rectGrob(
          x = unit(0.035, "npc"),
          y = unit(0.975, "npc"),
          width = unit(0.10, "npc"),
          height = unit(0.07, "npc"),
          gp = gpar(fill = "white", col = NA)
        ),
        textGrob(
          label,
          x = unit(0.055, "npc"),
          y = unit(0.965, "npc"),
          gp = gpar(fontsize = 18)
        )
      ),
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    )
}

panel_a <- add_panel_label(panel_a, "(A)")
panel_b <- add_panel_label(panel_b, "(B)")
panel_c <- add_panel_label(panel_c, "(C)")
###patchwork ##################################



final_plot <- panel_a / panel_b / panel_c

final_plot

ggsave(
  plot = final_plot,
  width = 8,
  height = 15,
  dpi = 800
)

#############################################





# 95% DETECTION THRESHOLDS


qpcr_95 <- qpcr_sensitivity %>%
  group_by(system) %>%
  summarise(
    qPCR_replicates_for_95 = ifelse(
      any(sensitivity >= 0.95),
      min(k[sensitivity >= 0.95]),
      NA
    ),
    .groups = "drop"
  )

survey_95 <- survey_sensitivity %>%
  group_by(system) %>%
  summarise(
    air_samples_for_95 = ifelse(
      any(sensitivity >= 0.95),
      min(ns[sensitivity >= 0.95]),
      NA
    ),
    .groups = "drop"
  )

sensitivity_95_summary <- qpcr_95 %>%
  left_join(survey_95, by = "system")

print(sensitivity_95_summary)
