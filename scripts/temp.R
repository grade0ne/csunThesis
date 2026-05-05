library(tidyverse)

get_deriv_df <- function(fit, id, n = 200) {
  x <- seq(min(fit@obs$time), max(fit@obs$time), length.out = n)
  d1 <- predict(fit@fit, x = x, deriv = 1)
  
  tibble(
    id = id,
    time = x,
    deriv = d1$y,
    mumax_pkg = unname(fit@par["mumax"]),
    t_mumax_pkg = fit@xy[1]
  )
}

# assuming you already added t_mumax_pkg
near_zero <- exp1_growth_summary %>%
  filter(t_mumax < 0.01) %>%
  slice(1:3)

interior <- exp1_growth_summary %>%
  filter(t_mumax > 2) %>%
  slice(1:3)

subset_ids <- c(near_zero$ID, interior$ID)
subset_idx <- match(subset_ids, exp1_growth_summary$ID)

deriv_df <- map2_dfr(
  subset_idx,
  subset_ids,
  ~ get_deriv_df(exp1_growth_models[[.x]], .y)
)

ggplot(deriv_df, aes(time, deriv)) +
  geom_line() +
  geom_point(aes(x = t_mumax_pkg, y = mumax_pkg), color = "red", size = 2) +
  facet_wrap(~ id, scales = "free_y") +
  labs(
    title = "Derivative curves (dy/dt)",
    subtitle = "Red point = package μmax",
    x = "Time",
    y = "Growth rate (dy/dt)"
  ) +
  theme_bw()
