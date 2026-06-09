
# setup -------------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(DHARMa)

source('scripts/wrangling/identify_mesoBouts.R')


# exploration -------------------------------------------------------------

# visualize distribution
bout_frequency %>% 
  ggplot() +
  geom_histogram(aes(x = (bout_freq)), color = 'black', fill = 'grey', bins = 100) +
  labs(x = 'Bout Frequency (bouts / hr.)') +
  coord_cartesian() +
  theme_classic()

# distribution by cluster
bout.freqpoly.cluster <- 
  bout_frequency %>% 
  ggplot() + 
  geom_freqpoly(aes(x = (bout_freq), color = cluster, group = ptt)) +
  geom_vline(data = bout_frequency %>% group_by(cluster, species) %>% summarize(median = median(bout_freq, na.rm =T), .groups = 'drop'), aes(xintercept = (median), linetype = species), alpha = 0.5) +
  facet_grid(cluster~.) +
  labs(x = 'Bout Frequency (bouts / hr.)') +
  coord_cartesian() +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  theme_classic()

ggsave(
  bout.freqpoly.cluster,
  file = 'products/figures/revisions/exploration/frequency_freqpoly_cluster.png',
  height = 6,
  width = 8,
  dpi = 300)

# distribution by species
bout.freqpoly.species <- 
  bout_frequency %>% 
  ggplot() + 
  geom_freqpoly(aes(x = (bout_freq), color = species, group = ptt)) +
  geom_vline(data = bout_frequency %>% group_by(species, ptt) %>% summarize(median = median(bout_freq, na.rm =T), .groups = 'drop'), aes(xintercept = (median), color = species), alpha = 0.5) +
  facet_grid(species~.) +
  labs(x = 'Bout Frequency (bouts / hr.)') +
  coord_cartesian() +
  theme_classic()

ggsave(
  bout.freqpoly.species,
  file = 'products/figures/revisions/exploration/frequency_freqpoly_species.png',
  height = 6,
  width = 8,
  dpi = 300)


# modelling ---------------------------------------------------------------

# full model
f.mod1 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + cluster:species,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

f.mod2 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + cluster:species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

f.mod3 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + cluster:species + (1|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)


f.mod4 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + cluster:species + (species|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

# get AIC weights
vAIC <- AIC(f.mod1, f.mod2, f.mod3, f.mod4)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw #f.mod4 is hands down best

drop1(f.mod2)

summary(f.mod2)

car::Anova(f.mod2, type = 'III')

sim1 <- DHARMa::simulateResiduals(f.mod3)
plot(sim1)

plotQQunif(sim1)
testDispersion(sim1) # no overdispersion

testCategorical(sim1, catPred = bout_frequency$cluster) # heteroscedastic
testCategorical(sim1, catPred = bout_frequency$species) # 
testCategorical(sim1, catPred = bout_frequency %>% mutate(month = lubridate::month(date)) %>% pull(month)) # heteroscedastic

## step 2
# remove non-significant predictors
f.mod7 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

# test AIC weights
vAIC <- AIC(f.mod2, f.mod7)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw


sim2 <- DHARMa::simulateResiduals(f.mod7)
plot(sim2)
plotQQunif(sim2)
testQuantiles(sim2)

testCategorical(sim2, catPred = bout_frequency$cluster) # heteroscedastic
testCategorical(sim2, catPred = bout_frequency$species) # 

## step 3
# add dispersion formulae
f.mod8 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    dispformula = ~cluster,
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

f.mod9 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    dispformula = ~species,
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

f.mod10 <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    dispformula = ~cluster + species,
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

# test AIC weights
vAIC <- AIC(f.mod7, f.mod8, f.mod9, f.mod10)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw

sim3 <- DHARMa::simulateResiduals(f.mod10)
plot(sim3)
plotQQunif(sim3)
testQuantiles(sim3)

testCategorical(sim3, catPred = bout_frequency$cluster) # heteroscedastic
testCategorical(sim3, catPred = bout_frequency$species) # 


# evaluation --------------------------------------------------------------

freq.mod <- 
  bout_frequency %>% 
  glmmTMB(
    bout_freq ~ cluster + species + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum'),
    dispformula = ~cluster + species,
    family = beta_family(link = 'logit'),
    REML = T,
    data = .)

# DHARMa evaluation
sim <- DHARMa::simulateResiduals(freq.mod)
plot(sim)

testQuantiles(sim)
testCategorical(sim, catPred = bout_frequency$cluster) 
testCategorical(sim, catPred = bout_frequency$species) # 
# testCategorical(sim, catPred = bout_frequency$month) # 
# testQuantiles(sim, predictor = bout_frequency$lat) # 

# check VIF
performance::check_collinearity(freq.mod, component = 'all')

# omnibus anova
car::Anova(freq.mod, type = 'II')

# emmeans
emmeans::emmeans(freq.mod, ~cluster, type = 'response')
emmeans::emmeans(freq.mod, ~species, type = 'response')


# re-do the selection with month / lat ------------------------------------

# full model
f.mod1 <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:species + cluster:month + cluster:lat + species:month + species:lat + month:lat,
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

f.mod2 <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:species + cluster:month + cluster:lat + species:month + species:lat + month:lat + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    family = beta_family(link = 'logit'),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
    REML = F,
    data = .)

f.mod3 <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:species + cluster:month + cluster:lat + species:month + species:lat + month:lat + (1|ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    family = beta_family(link = 'logit'),
    REML = F,
    data = .)

# get AIC weights
vAIC <- AIC(f.mod1, f.mod2, f.mod3)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw #f.mod4 is hands down best

drop1(f.mod2)

summary(f.mod2)

sim1 <- DHARMa::simulateResiduals(f.mod2)
plot(sim1)

plotQQunif(sim1)
testDispersion(sim1) # no overdispersion

## step 2
# remove non-signif. interactions
f.mod7 <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:lat + species:lat + month:lat + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    family = beta_family(link = 'logit'),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
    REML = F,
    data = .)

# test AIC weights
vAIC <- AIC(f.mod2, f.mod7)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw


sim2 <- DHARMa::simulateResiduals(f.mod7)
plot(sim2)
plotQQunif(sim2)
testQuantiles(sim2)

testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(cluster)) # heteroscedastic
testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(species)) # 
testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(month)) # heteroscedastic
testQuantiles(sim2, predictor = bout_frequency %>% filter(!is.na(lat)) %>% pull(lat))

## step 3
# assess a dispersion parameter
f.mod8 <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:lat + species:lat + month:lat + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    dispformula = ~ species,
    family = beta_family(link = 'logit'),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
    REML = F,
    data = .)

# test AIC weights
vAIC <- AIC(f.mod7, f.mod8)[,2]
dAIC <- vAIC - min(vAIC, na.rm = T)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2), na.rm = T)
AICw

sim3 <- DHARMa::simulateResiduals(f.mod8)
plot(sim3)
plotQQunif(sim3)
testQuantiles(sim3)

testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(cluster)) # heteroscedastic
testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(species)) # 
testCategorical(sim2, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(month)) # heteroscedastic
testQuantiles(sim2, predictor = bout_frequency %>% filter(!is.na(lat)) %>% pull(lat))


# evaluation --------------------------------------------------------------
freq.mod <- 
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  glmmTMB(
    bout_freq ~ cluster + species + month + lat + cluster:lat + species:lat + month:lat + (1|species:ptt),
    contrasts = list(cluster = 'contr.sum', species = 'contr.sum', month = 'contr.sum'),
    dispformula = ~ species,
    family = beta_family(link = 'logit'),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
    REML = T,
    data = .)

# DHARMa diagnostics
sim <- DHARMa::simulateResiduals(freq.mod)
plot(sim)
plotQQunif(sim)
testQuantiles(sim)

testCategorical(sim, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(cluster)) # heteroscedastic
testCategorical(sim, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(species)) # 
testCategorical(sim, catPred = bout_frequency %>% filter(!is.na(lat)) %>% pull(month)) # heteroscedastic
testQuantiles(sim, predictor = bout_frequency %>% filter(!is.na(lat)) %>% pull(lat))


# check VIF
performance::check_collinearity(freq.mod, component = 'all')

# omnibus anova results
car::Anova(freq.mod, type = 'III')

# emmeans
emmeans::emmeans(freq.mod, ~ cluster, type = 'response')
emmeans::emmeans(freq.mod, ~ species, type = 'response')
emmeans::emmeans(freq.mod, ~ month, type = 'response')

# visualize interactions
predict.data <-
  bout_frequency %>% 
  mutate(month = factor(month)) %>% 
  filter(!is.na(lat)) %>% 
  cbind(
    predict(freq.mod, newdata = ., type = 'response', se = T))

# cluster interaction
predict.data %>% 
  ggplot() +
  geom_point(aes(x = lat, y = bout_freq, color = cluster), alpha = 0.15) +
  geom_smooth(aes(x = lat, y = fit, color = cluster), se = F) +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  theme_classic()

# species interaction
predict.data %>% 
  ggplot() +
  geom_point(aes(x = lat, y = bout_freq, color = species), alpha = 0.15) +
  geom_smooth(aes(x = lat, y = fit, color = species), se = F) +
  theme_classic()

# month interaction
tibble(
  cluster = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3') %>% rep(each = 8) %>% rep(times = 2) %>% rep(times = 4),
  species = c('I.oxyrinchus', 'P.glauca') %>% rep(each = 40) %>% rep(times = 4),
  month = c(9,10,11,12,1,2,3,4) %>% rep(times = 10) %>% rep(times = 4) %>% factor(),
  lat = c(10.61, 31.60, 37.46, 43.53) %>% rep(each = 80)) %>% 
  cbind(
    predict(freq.mod, newdata = ., type = 'response', se = T, re.form = NA)) %>% 
  mutate(lat = factor(lat)) %>% 
  ggplot() +
  geom_point(aes(x = month, y = fit, color = lat), alpha = 0.15) +
  geom_smooth(aes(x = month, y = fit, color = lat, group = lat), se = F) +
  theme_classic()




