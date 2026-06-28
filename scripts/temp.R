library(lme4)
library(lmerTest)

# Full clone-level model (what the data can actually support)
m_full <- lmer(mumax ~ Tr + 
                 (1 | Site/Leaf/Clone) + 
                 (0 + Tr | Site:Leaf:Clone),
               data = exp1_growth_summary, REML = FALSE)

# Test: does leaf-level intercept variance contribute?
m_no_leaf <- lmer(mumax ~ Tr + 
                    (1 | Site) + 
                    (1 | Site:Leaf:Clone) + 
                    (0 + Tr | Site:Leaf:Clone),
                  data = exp1_growth_summary, REML = FALSE)

# Test: does site-level intercept variance contribute?
m_no_site <- lmer(mumax ~ Tr + 
                    (1 | Site:Leaf) + 
                    (1 | Site:Leaf:Clone) + 
                    (0 + Tr | Site:Leaf:Clone),
                  data = exp1_growth_summary, REML = FALSE)

# LRTs
anova(m_full, m_no_leaf)   # is leaf variance supported?
anova(m_full, m_no_site)   # is site variance supported?

m_full <- lmer(mumax ~ Tr + 
                 (1 | Site:Leaf:Clone) + 
                 (0 + Tr | Site:Leaf:Clone),
               data = exp1_growth_summary, REML = FALSE)
