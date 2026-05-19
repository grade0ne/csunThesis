# =============================================================================
# SEM: Indirect effects of rotifer evolutionary history on protist growth
# =============================================================================
#
# Model structure:
#   Exogenous predictors : evolvedTemp, currentTemp (both effect-coded)
#   Mediators            : rotifer r and rotifer alpha (competition-present,
#                          replicate means), modeled as parallel mediators
#   Outcome              : delta_r_protist = protist r (competition) minus
#                          protist r (no competition), within currentTemp x
#                          replicate, averaged across replicates
#
# Effect coding rationale:
#   evolvedTemp and currentTemp are symmetric two-level treatments with no
#   natural reference level. Effect coding (-0.5 / +0.5) centers the
#   predictors at the grand mean so that (a) main-effect coefficients
#   represent contrasts relative to the overall mean rather than one
#   arbitrary baseline, and (b) any interaction term is interpretable as
#   the difference in the evolvedTemp slope between the two currentTemps.
#   This is preferred over dummy coding when neither level is a true control.
#
# =============================================================================

library(tidyverse)
library(lavaan)

# -----------------------------------------------------------------------------
# 1.  PREPARE ROTIFER MEDIATOR DATA
#     Unit: evolvedTemp x currentTemp x replicate (competition-present only)
# -----------------------------------------------------------------------------

rotifers_comp <- exp2_growth_summary_rotifer |>
  filter(competition == "Present") |>
  # Summarise to replicate means within each treatment cell
  group_by(evolvedTemp, currentTemp, repNum) |>
  summarise(
    rot_r     = mean(r_day,   na.rm = TRUE),
    rot_alpha = mean(alpha_day, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 2.  PREPARE PROTIST OUTCOME DATA
#     delta_r = r (competition) - r (no competition)
#     Averaged within currentTemp x replicate, separately for each evolvedTemp
#     context (protists are not evolved, but they are paired with a specific
#     rotifer evolutionary history replicate).
# -----------------------------------------------------------------------------

# Protist r under competition — keep evolvedTemp to match with rotifer partner
protist_comp <- exp2_growth_summary_protist |>
  filter(competition == TRUE) |>
  group_by(evolvedTemp, currentTemp, repNum) |>
  summarise(r_comp = mean(r_day, na.rm = TRUE), .groups = "drop")

# Protist r without competition — no evolvedTemp context; average over currentTemp x rep
protist_nocomp <- exp2_growth_summary_protist |>
  filter(competition == FALSE) |>
  group_by(currentTemp, repNum) |>
  summarise(r_nocomp = mean(r_day, na.rm = TRUE), .groups = "drop")

# Delta r: join on currentTemp x replicate, then subtract
protist_delta <- protist_comp |>
  left_join(protist_nocomp, by = c("currentTemp", "repNum")) |>
  mutate(delta_r_protist = r_comp - r_nocomp)

# -----------------------------------------------------------------------------
# 3.  JOIN INTO ONE ANALYSIS DATA FRAME
# -----------------------------------------------------------------------------

sem_data <- rotifers_comp |>
  left_join(
    protist_delta |> select(evolvedTemp, currentTemp, repNum, delta_r_protist),
    by = c("evolvedTemp", "currentTemp", "repNum")
  ) |>
  # Effect-code both temperature predictors: 25°C = -0.5, 30°C = +0.5
  mutate(
    evolvedTemp_e  = if_else(as.character(evolvedTemp)  == "30",  0.5, -0.5),
    currentTemp_e  = if_else(as.character(currentTemp)  == "30",  0.5, -0.5)
  )

# Quick check — should have one row per evolvedTemp x currentTemp x replicate
glimpse(sem_data)
cat("\nMissing values per column:\n")
colSums(is.na(sem_data))

# -----------------------------------------------------------------------------
# 4.  SPECIFY THE SEM
#
#   Measurement model : none (all observed variables — path model)
#
#   Structural paths  :
#     evolvedTemp_e  → rot_r        (a1)
#     currentTemp_e  → rot_r        (c1)
#     evolvedTemp_e  → rot_alpha    (a2)
#     currentTemp_e  → rot_alpha    (c2)
#     rot_r          → delta_r_protist  (b1)
#     rot_alpha      → delta_r_protist  (b2)
#     evolvedTemp_e  → delta_r_protist  (direct effect, c')
#     currentTemp_e  → delta_r_protist  (direct effect of currentTemp)
#
#   Indirect effects (defined via := operator for bootstrapping):
#     ind_r     = a1 * b1   (evolvedTemp → rot_r → delta_r_protist)
#     ind_alpha = a2 * b2   (evolvedTemp → rot_alpha → delta_r_protist)
#     ind_total = ind_r + ind_alpha  (total indirect)
#     total     = ind_total + c'    (total effect of evolvedTemp)
#
# -----------------------------------------------------------------------------

sem_model <- "
  # --- Mediator equations ---
  rot_r     ~ a1 * evolvedTemp_e + c1 * currentTemp_e
  rot_alpha ~ a2 * evolvedTemp_e + c2 * currentTemp_e

  # --- Outcome equation ---
  delta_r_protist ~ b1 * rot_r +
                    b2 * rot_alpha +
                    cp * evolvedTemp_e +
                    d  * currentTemp_e

  # --- Indirect and total effects ---
  ind_via_r     := a1 * b1
  ind_via_alpha := a2 * b2
  ind_total     := (a1 * b1) + (a2 * b2)
  total_evolv   := (a1 * b1) + (a2 * b2) + cp
"


sem_model_r_only <- "
  rot_r           ~ a1 * evolvedTemp_e + c1 * currentTemp_e
  delta_r_protist ~ b1 * rot_r + cp * evolvedTemp_e + d * currentTemp_e
  ind_via_r := a1 * b1
  total_evolv := a1 * b1 + cp
"

# -----------------------------------------------------------------------------
# 5.  FIT THE SEM WITH BOOTSTRAPPED CIs FOR INDIRECT EFFECTS
# -----------------------------------------------------------------------------

set.seed(42)

sem_fit <- sem(
  model   = sem_model_r_only,
  data    = sem_data,
  se      = "bootstrap",
  bootstrap = 2000,
  estimator = "ML"
)

# -----------------------------------------------------------------------------
# 6.  RESULTS
# -----------------------------------------------------------------------------

cat("\n===== MODEL SUMMARY =====\n")
summary(sem_fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

cat("\n===== FIT INDICES =====\n")
fitMeasures(sem_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea",
                       "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))

cat("\n===== ALL PARAMETER ESTIMATES WITH BOOTSTRAP CIs =====\n")
param_est <- parameterEstimates(
  sem_fit,
  boot.ci.type = "bca.simple",   # bias-corrected and accelerated
  level        = 0.95,
  standardized = TRUE
)
print(param_est, digits = 4)

# --- Focused table: indirect and total effects only ---
cat("\n===== INDIRECT & TOTAL EFFECTS (key paths) =====\n")
indirect_est <- param_est |>
  filter(op == ":=") |>
  select(label, est, se, z, pvalue, ci.lower, ci.upper, std.all) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))
print(indirect_est)

# -----------------------------------------------------------------------------
# 7.  OPTIONAL: STANDARDISED SOLUTION FOR REPORTING
# -----------------------------------------------------------------------------

cat("\n===== STANDARDISED SOLUTION =====\n")
standardizedSolution(sem_fit, type = "std.all") |>
  filter(op %in% c("~", ":=")) |>
  select(lhs, op, rhs, est.std, se, z, pvalue, ci.lower, ci.upper) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print()

# -----------------------------------------------------------------------------
# 8.  PLOT: Path diagram data (for use with semPlot or tidySEM)
# -----------------------------------------------------------------------------

# Install semPlot if needed: install.packages("semPlot")
library(semPlot)
semPaths(sem_fit,
          what       = "std",        # standardised coefficients
          whatLabels = "std",
          layout     = "tree2",
          residuals  = FALSE,
          intercepts = FALSE,
          edge.label.cex = 0.85,
          nCharNodes = 0,
          title      = FALSE)

# Alternatively with tidySEM:
# install.packages("tidySEM")
# library(tidySEM)
# graph_sem(sem_fit)
# 
# 
# 





library(lme4)

# Build matched dataset (you mostly have this already from the SEM prep)
matched_data <- sem_data  # evolvedTemp, currentTemp, repNum,
# rot_r, rot_alpha, delta_r_protist

# Primary model: does rot_r predict protist competitive suppression?
fit <- lmer(delta_r_protist ~ rot_r +
              evolvedTemp_e +
              currentTemp_e +
              (1 | repNum),   # replicate as random intercept
            data = matched_data,
            REML = FALSE)

# Compare to model without rot_r to get likelihood ratio test
fit_no_rotr <- lmer(delta_r_protist ~ evolvedTemp_e +
                      currentTemp_e +
                      (1 | repNum),
                    data = matched_data,
                    REML = FALSE)

anova(fit_no_rotr, fit)  # LRT: does rot_r explain protist outcomes?

library(ggplot2)

ggplot(matched_data,
       aes(x = rot_r, y = delta_r_protist,
           color = evolvedTemp, shape = currentTemp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE,
              aes(group = currentTemp),   # or group = 1 for pooled line
              color = "grey40") +
  labs(x = "Rotifer intrinsic growth rate (competition)",
       y = expression(Delta * "r  protist (competition \u2212 alone)"),
       color = "Evolved temp",
       shape = "Current temp") +
  theme_classic() +
  facet_wrap(~currentTemp)


fit_interaction <- lmer(delta_r_protist ~ rot_r * currentTemp_e +
                          evolvedTemp_e +
                          (1 | repNum),
                        data = matched_data,
                        REML = FALSE)

anova(fit, fit_interaction)  # is the interaction supported?
summary(fit_interaction)


# Fit a model with the three-way structure you care about
fit_check <- lmer(delta_r_protist ~ evolvedTemp_e * currentTemp_e + (1 | repNum),
                  data = matched_data, REML = FALSE)

# Get means by evolvedTemp within each currentTemp
emm <- emmeans(fit_check, ~ evolvedTemp_e | currentTemp_e)
contrast(emm, "pairwise")
