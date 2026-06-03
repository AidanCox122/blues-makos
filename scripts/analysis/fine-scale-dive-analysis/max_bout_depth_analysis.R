# setup -------------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(DHARMa)
devtools::load_all('../../Documents/Ecology/Marine/tags2etuff')

source('scripts/wrangling/identify_mesoBouts.R')


# exploration -------------------------------------------------------------

# visualize max depth distribution
bout_duration %>% 
  ggplot() +
  geom_histogram(aes(x = (max.z)), color = 'black', fill = 'grey', binwidth = 25) +
  coord_cartesian() +
  theme_classic()

# max.z distribution by cluster
z.freqpoly.cluster <-
  bout_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(max.z), color = cluster, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(cluster, species) %>% summarize(median = median(max.z), .groups = 'drop'), aes(xintercept = log10(median), linetype = species), alpha = 0.5) +
  facet_grid(cluster~.) +
  coord_cartesian() +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  theme_classic()

# ggsave(
#   z.freqpoly.cluster,
#   file = 'products/figures/revisions/exploration/max_z_freqpoly_cluster.png',
#   height = 6,
#   width = 8,
#   dpi = 300)

# max.z distribution by species
z.freqpoly.species <- 
  bout_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(max.z), color = species, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(species, ptt) %>% summarize(median = median(max.z), .groups = 'drop'), aes(xintercept = log10(median), color = species), alpha = 0.5) +
  facet_grid(species~.) +
  coord_cartesian() +
  guides(color = 'none') +
  theme_classic()

# ggsave(
#   z.freqpoly.species,
#   file = 'products/figures/revisions/exploration/max_z_freqpoly_species.png',
#   height = 6,
#   width = 8,
#   dpi = 300)

# model max depth ---------------------------------------------------------

# full model
z.mod1 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

z.mod2 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

z.mod3 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species + (1|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

z.mod4 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

z.mod5 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species + (1|species:ptt) + (1|ptt:Date),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

z.mod6 <- 
  bout_duration %>% 
  glmmTMB(
    max.z ~ cluster + species + cluster:species + (species|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

# get AIC weights
vAIC <- AIC(z.mod1, z.mod2, z.mod3, z.mod4, z.mod5, z.mod6)[,2]
dAIC <- vAIC - min(vAIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw #d.mod4 is hands down best

drop1(d.mod4)

summary(d.mod4)

car::Anova(d.mod4, type = 'III')

sim1 <- DHARMa::simulateResiduals(d.mod4)
plot(sim1)

testDispersion(sim1) # no overdispersion

testCategorical(sim1, catPred = bout_duration$cluster) # heteroscedastic
testCategorical(sim1, catPred = bout_duration$species) # 
testCategorical(sim1, catPred = bout_duration$month) # heteroscedastic



