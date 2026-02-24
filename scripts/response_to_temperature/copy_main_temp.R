library(tidyverse)

#-------------- Data Handling --------------

popData <- read.csv("Data/tempres/tempRes2-30d.csv") %>%
  mutate(Day = Day - 3,
         across(c(Site, Leaf, Clone, Treatment, ID), as.factor),
         across(c(Day, Count), as.integer))

clone_key <- c(
  "1" = 1,  "2" = 2,  "17" = 3,  "18" = 4,  "19" = 5,
  "72" = 6, "73" = 7, "74" = 8, "39" = 9,  "40" = 10,
  "42" = 11,"75" = 12,"26" = 13,"25" = 14,"35" = 15,
  "36" = 16,"37" = 17
)

sizeData <- read.csv("Data/sizeData/diversitySizeData.csv") %>%
  separate(base, into = c("unk", "culture", "image", "object"), sep = "_", convert = TRUE) %>%
  separate(culture, into = c("siteLeaf", "clone", "replicate"), sep = "-", convert = TRUE) %>%
  mutate(site = as.integer(case_when(str_sub(siteLeaf, 1, 1) == "c" ~ 1,
                                     str_sub(siteLeaf, 1, 1) == "p" ~ 2)),
         leaf = str_sub(siteLeaf, start = 2),
         leaf = as.integer(case_when(leaf == "23" ~ 1,
                                     leaf == "7"  ~ 2,
                                     leaf == "36" & site == "1" ~ 3,
                                     leaf == "36" & site == "2" ~ 4,
                                     leaf == "24" ~ 5,
                                     leaf == "3"  ~ 6)),
         clone = as.integer(unname(clone_key[as.character(clone)])),
         temperature = as.factor(case_when(replicate %in% c(1, 2, 3) ~ "25C",
                                           replicate %in% c(4, 5, 6) ~ "30C"))
  ) %>%
  #select(
  #  temperature, site, leaf, clone, area_px) %>%
  mutate(across(c(temperature, site, leaf, clone), factor)) %>%
  filter_out(clone=="6")


#-------------- Growth Models --------------

library(growthrates)

popData$Count1 <- popData$Count + 1

models <- all_splines(Count1 ~ Day | Treatment + Site + Leaf + Clone + ID,
                      data = popData,
                      spar = 0.5)

# Visually inspect fit of models
plot_models <- function(models) {
  coefs <- as.data.frame(coef(models))
  model_ids <- c(1:length(coefs$y0))
  
  columns <- floor(sqrt(length(model_ids)))
  rows <- columns + 1
  par(mfrow = c(columns, rows), mar = c(1, 1, 1, 1))
  
  for (i in seq_along(model_ids)) {
    plot(models[[i]])
    
    this_r <- round(coefs$mumax[i], 6)
    this_model <- rownames(coefs)[i]
    
    mtext(this_model, side = 3, adj = 0.1, line = -1, cex = 0.5)
    mtext(paste("r:", this_r), side = 3, adj = 0.1, line = -2, cex = 0.5)
  }
}



#-------------- Linear Mixed Models -------------- 
library(car)
library(moments)
library(lme4)
library(lmerTest)

summaryData <- results(models) %>%
  mutate(across(c(Treatment, Site, Leaf, Clone), factor))

summaryData$Tr <- ifelse(summaryData$Treatment == "30C", 0.5, -0.5)


# with leaf
model_1 <- lmer(mumax ~ Tr + (1|Site:Leaf) + (0 + Tr|Site:Leaf) +
                  (1|Site:Leaf:Clone) + (0 + Tr|Site:Leaf:Clone), summaryData)
# clone only
model_2 <- lmer(mumax ~ Tr + (1|Site:Leaf:Clone) + (0 + Tr|Site:Leaf:Clone), summaryData)

vars <- c(
  clone25 = 3.028e-04,
  clone30 = 1.166e-03,
  leaf30 = 7.851e-08,
  leaf25 = 2.805e-05,
  resid = 1.947e-03)
total_var <- sum(vars)

var_df <- data.frame(
  source = c(names(vars), "total"),
  variance = c(vars, total_var),
  percent = c(vars / total_var * 100, 100),
  row.names = NULL)

library(emmeans)

fig1data_emm <- summary(emmeans(model_full, ~ Tr))



sizeData$Tr <- ifelse(sizeData$temperature == "30C", 0.5, -0.5)

lmer(mumax ~ Tr + (1|Site:Leaf) + (0 + Tr|Site:Leaf) +
       (1|Site:Leaf:Clone) + (0 + Tr|Site:Leaf:Clone), summaryData)



model_full_size <- lmer(area_px ~ Tr + (1 | site) + (1 | site:leaf) + (1 | site:leaf:clone), sizeData)

vars_size <- c(
  clone = 8345,
  leaf = 0,
  site = 5225,
  resid = 1208593)
total_var_size <- sum(vars_size)

var_df_size <- data.frame(
  source = c(names(vars_size), "total"),
  variance = c(vars_size, total_var_size),
  percent = c(vars_size / total_var_size * 100, 100),
  row.names = NULL)



# ------------ Data Vis ----------------

library(paletteer)


fig4data <- summaryData %>%
  group_by(Treatment, Leaf, Clone) %>%
  summarize(mean = mean(mumax),
            ci = 1.96 * (sd(mumax)/sqrt(length(mumax))))

ggplot(fig3data, aes(x = Leaf, y = mean, color = Treatment, group = Leaf, shape = Treatment)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.1), size = 3) +
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.1),
                aes(ymin = mean - ci, ymax = mean + ci), width = 0.3) +
  #geom_line(stat = "identity", position = position_dodge(width = 0.1)) +
  geom_point(data = fig4data, aes(x = Leaf, y = mean, color = Treatment), alpha = 0.5,
             position = position_jitter(width = 0.2)) +
  scale_color_manual(values = c("#1BB6AF", "#E9A17C")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  geom_text(data = data.frame(Leaf = c(0.9, 4), y = 0.39, label = c("Site 1", "Site 2")),
            aes(x = Leaf, y = y, label = label), inherit.aes = FALSE) +
  labs(x = "Leaf", y = "Maximum growth rate (μmax)") +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )


ggplot(fig4data, aes(x = Clone, y = mean, color = Treatment, shape = Treatment)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.1), size = 3) +
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.1),
                aes(ymin = mean - ci, ymax = mean + ci), width = 0.3) +
  #geom_point(data = summaryData, aes(x = Clone, y = mumax, color = Treatment), alpha = 0.5,
  #           position = position_jitter(width = 0.2)) +
  scale_color_manual(values = c("#1BB6AF", "#E9A17C")) +
  labs(x = "Clone", y = "Maximum growth rate (μmax)") +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )
                
# ------------ Extra Figs ----------------

fig1data <- summaryData %>%
  group_by(Treatment) %>%
  summarize(mean = mean(mumax),
            ci = 1.96 * (sd(mumax)/sqrt(length(mumax))))

ggplot(fig1data, aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_errorbar(stat = "identity", width = 0.8,
                aes(ymin = mean - ci, ymax = mean + ci)) +
  labs(x = "Temperature", y = "Maximum growth rate (μmax)") +
  scale_fill_manual(values = c("seagreen2", "slateblue2")) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )

fig2data <- summaryData %>%
  group_by(Treatment, Site) %>%
  summarize(mean = mean(mumax),
            ci = 1.96 * (sd(mumax)/sqrt(length(mumax))))

ggplot(fig2data, aes(x = Treatment, y = mean, color = Site, group = Site)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.1), size = 2.5) +
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.1),
                aes(ymin = mean - ci, ymax = mean + ci), width = 0.3) +
  geom_line(stat = "identity", position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("seagreen2", "slateblue2")) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )

fig3data <- summaryData %>%
  group_by(Treatment, Leaf) %>%
  summarize(mean = mean(mumax),
            ci = 1.96 * (sd(mumax)/sqrt(length(mumax))))

ggplot(fig3data, aes(x = Treatment, y = mean, color = Leaf, group = Leaf)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.1), size = 2.5) +
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.1),
                aes(ymin = mean - ci, ymax = mean + ci), width = 0.3) +
  geom_line(stat = "identity", position = position_dodge(width = 0.1)) +
  scale_color_paletteer_d("LaCroixColoR::PeachPear") +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )

                
                



