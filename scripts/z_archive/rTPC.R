library(rTPC)
library(tidyverse)
library(car)

data("bacteria_tpc")
library(ggplot2)
ggplot(bacteria_tpc) +
  geom_point(aes(temp, rate, col = phage))

# load in ggplot
library(ggplot2)
# subset for the first TPC curve
data('chlorella_tpc')
d <- subset(chlorella_tpc, curve_id == 1)
# get start values and fit model
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'gaussian_1987')
# fit model
mod <- nls.multstart::nls_multstart(rate~gaussian_1987(temp = temp,rmax, topt,a),
                                    data = d,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 10,
                                    start_upper = start_vals + 10,
                                    lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                                    upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)
# look at model fit
summary(mod)
# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)
# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw()