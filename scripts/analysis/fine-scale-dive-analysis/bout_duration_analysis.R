# setup -------------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(DHARMa)
devtools::load_all('../../Documents/Ecology/Marine/tags2etuff')

source('scripts/wrangling/identify_mesoBouts.R')


# exploration -------------------------------------------------------------

# visualize distribution
bout_duration %>% 
  ggplot() +
  geom_histogram(aes(x = (duration.min)), color = 'black', fill = 'grey', binwidth = 5) +
  annotate(geom = 'text', label = 'binwidth = 5 mins.', fontface = 'italic', color = 'grey', x = 750, y = 300, size = 5) +
  coord_cartesian() +
  theme_classic()

# distribution by cluster
duration.freqpoly.cluster <-
  bout_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(duration.min), color = cluster, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(cluster, species) %>% summarize(median = median(duration.min), .groups = 'drop'), aes(xintercept = log10(median), linetype = species), alpha = 0.5) +
  facet_grid(cluster~.) +
  coord_cartesian() +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  theme_classic()

# ggsave(
#   duration.freqpoly.cluster,
#   file = 'products/figures/revisions/exploration/duration_freqpoly_cluster.png',
#   height = 6,
#   width = 8,
#   dpi = 300)

# distribution by species
duration.freqpoly.species <- 
  bout_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(duration.min), color = species, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(species, ptt) %>% summarize(median = median(duration.min), .groups = 'drop'), aes(xintercept = log10(median), color = species), alpha = 0.5) +
  facet_grid(species~.) +
  coord_cartesian() +
  guides(color = 'none') +
  theme_classic()

# ggsave(
#   duration.freqpoly.species,
#   file = 'products/figures/revisions/exploration/duration_freqpoly_species.png',
#   height = 6,
#   width = 8,
#   dpi = 300)

# boxplots to visualize variation
bout_duration %>% 
  ggplot() + 
  geom_boxplot(aes(x = cluster, y = log10(duration.min), color = cluster)) +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  coord_cartesian() +
  guides(color = 'none') +
  theme_classic()

# modelling ---------------------------------------------------------------

# visualize correlations
bout_duration %>% 
  transmute(duration.min, cluster = as.numeric(cluster), species = as.numeric(factor(species)), month, lat) %>% 
  PerformanceAnalytics::chart.Correlation()

# remove outliers for test modelling
bout_duration <- 
  bout_duration %>% 
  filter(duration.min < 250) # this removes 36 records (DVM1-4; DVM2-14; DVM3-18)

# full model
d.mod1 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod2 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod3 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod4 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (species|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

# get AIC weights
vAIC <- AIC(d.mod1, d.mod2, d.mod3, d.mod4)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw #d.mod4 is hands down best

drop1(d.mod2)

summary(d.mod4)

car::Anova(d.mod2, type = 'III')

sim1 <- DHARMa::simulateResiduals(d.mod2)
plot(sim1)

plotQQunif(sim1)
testDispersion(sim1) # no overdispersion

testCategorical(sim1, catPred = bout_duration$cluster) # heteroscedastic
testCategorical(sim1, catPred = bout_duration$species) # 
testCategorical(sim1, catPred = bout_duration$month) # heteroscedastic

# add a dispersion parameter
d.mod7 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date),
    dispformula = ~cluster,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod8 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date),
    dispformula = ~species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod9 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date),
    dispformula = ~cluster+species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

d.mod10 <- 
  bout_duration %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date),
    dispformula = ~cluster*species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

vAIC <- AIC(d.mod7, d.mod8, d.mod9, d.mod10)[,2]
dAIC <- vAIC - min(vAIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw #d.mod9

sim2 <- simulateResiduals(d.mod7)
plot(sim2)  
plotQQunif(sim2)

testDispersion(sim2) # no overdispersion

testCategorical(sim2, catPred = sim2$fittedPredictedResponse)
testCategorical(sim2, catPred = bout_duration$cluster) # heteroscedastic
testCategorical(sim2, catPred = bout_duration$species) # heteroscedastic
testCategorical(sim2, catPred = bout_duration$month) # heteroscedastic
testQuantiles(sim2, predictor = bout_duration$lat) # heteroscedastic

sim2 %>% 
  recalculateResiduals(group = bout_duration$start) %>% 
  testTemporalAutocorrelation(time = unique(bout_duration$start)) # strong temporal autocorrelation

# test spatial autocorrelation

# create a location ID
test <- bout_duration %>% 
  mutate(latlon = paste(lat, lon)) 
# find unique locations
coords <- unique(
  test[, c("latlon", "lon", "lat")]
)
# recalculate residuals at each location ID
resids <- recalculateResiduals(sim2, group = test$latlon)

# remove any NA positions
keep <- complete.cases(coords$lon, coords$lat)

# get complete coords and corresponding resids
coords <- coords[keep, ]
resids$scaledResiduals <- resids$scaledResiduals[keep]

# test spatial autocorrelation
testSpatialAutocorrelation(resids, x = coords$lon, y = coords$lat) # no spatial autocorrelation

# cleanup
rm(test, coords, resids, keep)

# add autoregressive correlation structure
d.mod11 <- 
  bout_duration %>% 
  mutate(time = factor(as.numeric(start), levels = as.numeric(bout_duration$start) %>% unique() %>% sort())) %>% 
  glmmTMB(
    duration.min ~ cluster + species + cluster:species + (1|ptt) + (1|ptt:Date) + ar1(time + 0|ptt),
    dispformula = ~cluster+species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = Gamma(link = 'log'),
    REML = F,
    data = .)

AIC(d.mod9, d.mod11)

sim3 <- simulateResiduals(d.mod11)
plot(sim3)

testCategorical(sim3, catPred = sim3$fittedPredictedResponse)
testCategorical(sim3, catPred = bout_duration$cluster) # heteroscedastic
testCategorical(sim3, catPred = bout_duration$species) # heteroscedastic
testCategorical(sim3, catPred = bout_duration$month) # heteroscedastic
testQuantiles(sim3, predictor = bout_duration$lat) # heteroscedastic

sim3 %>% 
  recalculateResiduals(group = bout_duration$start) %>% 
  testTemporalAutocorrelation(time = unique(bout_duration$start)) # strong temporal autocorrelation


