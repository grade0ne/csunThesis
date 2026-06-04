library(tidyverse)

predict_logistic_alpha <- function(t, y0, r, alpha) {
  (r * y0) / (alpha * y0 + (r - alpha * y0) * exp(-r * t))
}

fit_theme <- theme_bw(9) + theme(
  panel.grid = element_blank(), strip.text = element_blank(),
  strip.background = element_blank(), legend.position = "right",
  legend.text = element_text(size = 7), legend.title = element_text(size = 7),
  axis.text = element_text(size = 7), axis.title = element_text(size = 10),
  panel.spacing = unit(0.3, "lines"))

t_grid <- seq(0, 510, by = 2)

### ROTIFER ###
rot_sum <- exp2_growth_summary_rotifer %>%
  mutate(treatment = paste(evolvedTemp, currentTemp, competition, sep = " | "),
         panel = paste(treatment, replicate, sep = " | Rep "),
         r2_text = sprintf("R² = %.2f", r2))

rot_fits <- rot_sum %>% rowwise() %>%
  mutate(d = list(tibble(day = t_grid/24,
                         count_pred = predict_logistic_alpha(t_grid, y0, r, alpha)))) %>%
  unnest(d) %>% ungroup()

rot_obs <- exp2_rotifer_growth %>%
  mutate(competition = if_else(competition == TRUE, "Present", "Absent"),
         treatment = paste(evolvedTemp, currentTemp, competition, sep = " | "),
         panel = paste(treatment, repNum, sep = " | Rep "),
         day = hour/24)

ggplot() +
  geom_point(data = rot_obs, aes(day, count, color = treatment), alpha = 0.7, size = 1.2) +
  geom_line(data = rot_fits, aes(day, count_pred, color = treatment), linewidth = 0.6) +
  geom_text(data = rot_sum, aes(max(rot_obs$day), max(rot_obs$count), label = r2_text),
            hjust = 1, vjust = 1.3, size = 2.5) +
  facet_wrap(~panel, ncol = 10) +
  scale_color_brewer(palette = "Dark2", name = "Treatment") +
  coord_cartesian(xlim = c(0, max(rot_obs$day)), ylim = c(0, max(rot_obs$count))) +
  labs(x = "Time (days)", y = expression("Rotifer count ml"^-1)) +
  fit_theme + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("figures/rotifer_model_fits_individual.pdf",
       width = 12, height = 6.5, units = "in", device = cairo_pdf)

### PROTIST ###
pro_sum <- exp2_growth_summary_protist %>%
  mutate(treatment = paste(compFactor, currentTemp, sep = " | "),
         panel = paste(treatment, replicate, sep = " | Rep "),
         r2_text = sprintf("R² = %.2f", r2))

pro_fits <- pro_sum %>% rowwise() %>%
  mutate(d = list(tibble(day = t_grid/24,
                         count_pred = predict_logistic_alpha(t_grid, y0, r, alpha)))) %>%
  unnest(d) %>% ungroup()

pro_obs <- exp2_protist_growth %>%
  mutate(compFactor = case_when(
    evolvedTemp == "25" & competition == TRUE ~ "25-evol-rotif",
    evolvedTemp == "30" & competition == TRUE ~ "30-evol-rotif",
    competition == FALSE ~ "no-comp"),
    treatment = paste(compFactor, currentTemp, sep = " | "),
    panel = paste(treatment, repNum, sep = " | Rep "),
    day = hour/24)

ggplot() +
  geom_point(data = pro_obs, aes(day, count, color = treatment), alpha = 0.7, size = 1.2) +
  geom_line(data = pro_fits, aes(day, count_pred, color = treatment), linewidth = 0.6) +
  geom_text(data = pro_sum, aes(max(pro_obs$day), max(pro_obs$count), label = r2_text),
            hjust = 1, vjust = 1.3, size = 2.5) +
  facet_wrap(~ panel, ncol = 10) +
  scale_color_brewer(palette = "Set2", name = "Treatment") +
  coord_cartesian(xlim = c(0, max(pro_obs$day)), ylim = c(0, max(pro_obs$count))) +
  labs(x = "Time (days)", y = expression("Protist count ml"^-1)) +
  fit_theme + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("figures/protist_model_fits_individual.pdf",
       width = 12, height = 6.5, units = "in", device = cairo_pdf)


################################################################################################
library(tidyverse)
library(growthrates)

# Resolve conflicts with MASS/stats that get loaded via growthrates
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag, .quiet = TRUE)

fit_theme <- theme_bw(9) + theme(
  panel.grid = element_blank(), strip.text = element_blank(),
  strip.background = element_blank(), legend.position = "right",
  legend.text = element_text(size = 7), legend.title = element_text(size = 7),
  axis.text = element_text(size = 7), axis.title = element_text(size = 10),
  panel.spacing = unit(0.3, "lines"))

### EXP1 — spline fits ###

group_vars <- c("Treatment", "Site", "Leaf", "Clone", "ID")

day_grid <- seq(min(exp1_growthcurves$Day, na.rm = TRUE),
                max(exp1_growthcurves$Day, na.rm = TRUE),
                length.out = 200)

split_name <- function(nm) {
  parts <- str_split_fixed(nm, ":", length(group_vars))
  colnames(parts) <- group_vars
  as_tibble(parts)
}

# Per-fit pieces: spline, tangent, R²
exp1_pieces <- imap(exp1_growth_models@fits, function(fit, nm) {
  
  meta <- split_name(nm)
  ss   <- fit@fit
  
  spline_df <- bind_cols(meta, tibble(
    Day = day_grid,
    Count1_pred = exp(predict(ss, x = day_grid)$y)
  ))
  
  d1        <- predict(ss, x = day_grid, deriv = 1)$y
  i_max     <- which.max(d1)
  t_tan     <- day_grid[i_max]
  log_y_tan <- predict(ss, x = t_tan)$y
  slope     <- d1[i_max]
  
  tan_df <- bind_cols(meta, tibble(
    Day = day_grid,
    Count1_pred = exp(log_y_tan + slope * (day_grid - t_tan))
  ))
  
  r2 <- as.numeric(fit@rsquared)
  r2_df <- bind_cols(meta, tibble(
    r2_text = if (is.finite(r2)) sprintf("R² = %.2f", r2) else "R² = NA"
  ))
  
  list(spline = spline_df, tangent = tan_df, r2 = r2_df)
})

exp1_spline  <- map_dfr(exp1_pieces, "spline")
exp1_tangent <- map_dfr(exp1_pieces, "tangent")
exp1_r2lab   <- map_dfr(exp1_pieces, "r2")

# Observed data (keep Clone numeric so it sorts correctly)
exp1_obs <- exp1_growthcurves %>%
  mutate(across(all_of(c("Treatment", "Site", "Leaf", "ID")), as.character),
         Clone_num = suppressWarnings(as.numeric(as.character(Clone))))

# Build panel ordering: sort by Clone numerically, then Treatment, then ID.
# Within each Clone, assign per-replicate suffixes a, b, c, ... in that order
# to produce labels like "1a", "1b", ..., "2a", ...
panel_order <- exp1_obs %>%
  distinct(Clone, Clone_num, Treatment, ID) %>%
  arrange(Clone_num, Treatment, ID) %>%
  group_by(Clone) %>%
  mutate(rep_suffix = letters[row_number()],
         clone_label = paste0(Clone, rep_suffix),
         panel = paste(Clone, ID, sep = " | ")) %>%
  ungroup() %>%
  arrange(Clone_num, Treatment, ID) %>%
  mutate(panel = factor(panel, levels = panel))

# Attach panel factor + clone_label to every layer
attach_panel <- function(df) {
  df %>%
    mutate(panel_key = paste(Clone, ID, sep = " | ")) %>%
    left_join(panel_order %>% select(panel_key = panel, clone_label) %>%
                mutate(panel_key = as.character(panel_key)),
              by = "panel_key") %>%
    mutate(panel = factor(panel_key, levels = levels(panel_order$panel))) %>%
    select(-panel_key)
}

exp1_obs     <- attach_panel(exp1_obs)
exp1_spline  <- attach_panel(exp1_spline)
exp1_tangent <- attach_panel(exp1_tangent)
exp1_r2lab   <- attach_panel(exp1_r2lab)

# Shared axes
ymax_global <- max(exp1_obs$Count1, na.rm = TRUE)
xmax_global <- max(exp1_obs$Day,    na.rm = TRUE)

# Clip curves so a runaway one can't push the axis up
exp1_spline  <- exp1_spline  %>% filter(Count1_pred <= ymax_global, Count1_pred >= 0)
exp1_tangent <- exp1_tangent %>% filter(Count1_pred <= ymax_global)

ggplot() +
  geom_point(data = exp1_obs,
             aes(Day, Count1, color = Treatment),
             alpha = 0.7, size = 1.2) +
  geom_line(data = exp1_spline,
            aes(Day, Count1_pred, color = Treatment),
            linewidth = 0.6) +
  geom_line(data = exp1_tangent,
            aes(Day, Count1_pred, color = Treatment),
            linewidth = 0.5, linetype = "dashed") +
  # Clone label (e.g. "1a") in top-left
  geom_text(data = exp1_r2lab,
            aes(x = 0, y = ymax_global, label = clone_label),
            hjust = 0, vjust = 1.3, size = 2.3, fontface = "bold") +
  # R² in top-right
  geom_text(data = exp1_r2lab,
            aes(x = xmax_global, y = ymax_global, label = r2_text),
            hjust = 1, vjust = 1.3, size = 2.3) +
  facet_wrap(~ panel, ncol = 12) +
  scale_color_brewer(palette = "Dark2", name = "Treatment") +
  coord_cartesian(xlim = c(0, xmax_global), ylim = c(0, ymax_global)) +
  labs(x = "Time (days)", y = "Count") +
  fit_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("figures/exp1_spline_fits_individual.pdf",
       width = 12, height = 6.5, units = "in", device = cairo_pdf)
