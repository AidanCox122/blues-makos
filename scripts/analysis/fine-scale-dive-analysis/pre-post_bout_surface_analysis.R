# setup -------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)

source('scripts/wrangling/identify_mesoBouts.R')


# exploration -------------------------------------------------------------

# visualize distribution
pre.dist <- 
  surface_use %>% 
  ggplot() +
  geom_histogram(aes(x = pre.prop.50), color = 'black', fill = 'grey') +
  labs(x = 'Pre-Bout Proportion Above 50m') +
  theme_classic()

post.dist <- 
  surface_use %>% 
  ggplot() +
  geom_histogram(aes(x = post.prop.50), color = 'black', fill = 'grey') +
  labs(x = 'Post-Bout Proportion Above 50m') +
  theme_classic()

base.dist <- 
  baseline_surface_use %>% 
  filter(!is.na(cluster)) %>% 
  ggplot() +
  geom_histogram(aes(x = base.prop.50), color = 'black', fill = 'grey') +
  labs(x = 'Baseline Proportion Above 50m') +
  theme_classic()

(pre.dist / post.dist)

base.dist

# distribution by cluster
## pre-stats 
pre.dist.cluster <-
  surface_use %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(pre.prop.50+0.001), color = cluster, group = ptt)) +
  geom_vline(data = surface_use %>% group_by(cluster, species) %>% summarize(median = median(pre.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), linetype = species)) +
  scale_color_manual(
    values = c("yellow",
               "yellow3",
               "#488E9EFF",
               "#404C8BFF",
               "#281A2CFF")) +
  facet_grid(cluster~.) +
  guides(color = 'none', linetype = 'none') +
  labs(x = 'log10(Pre-Bout Prop. Above 50m)') +
  theme_classic()

## post-stats 
post.dist.cluster <- 
  surface_use %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(post.prop.50+0.001), color = cluster, group = ptt)) +
  geom_vline(data = surface_use %>% group_by(cluster, species) %>% summarize(median = median(post.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), linetype = species)) +
  scale_color_manual(
    values = c("yellow",
               "yellow3",
               "#488E9EFF",
               "#404C8BFF",
               "#281A2CFF")) +
  facet_grid(cluster~.) +
  labs(x = 'log10(Post-Bout Prop. Above 50m)') +
  theme_classic()

base.dist.cluster <-
  baseline_surface_use %>%
  filter(!is.na(cluster)) %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(base.prop.50+0.001), color = cluster, group = ptt)) +
  geom_vline(data = baseline_surface_use %>% filter(!is.na(cluster)) %>% group_by(cluster, species) %>% summarize(median = median(base.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), linetype = species)) +
  scale_color_manual(
    values = c("yellow",
               "yellow3",
               "#488E9EFF",
               "#404C8BFF",
               "#281A2CFF")) +
  facet_grid(cluster~.) +
  labs(x = 'log10(Baseline Prop. Above 50m)') +
  theme_classic()

pre.dist.cluster + post.dist.cluster

base.dist.cluster

# distribution by ptts
## pre-stats
pre.dist.ptt <- 
  surface_use %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(pre.prop.50+0.001), color = species, group = ptt)) +
  geom_vline(data = surface_use %>% group_by(species, ptt) %>% summarize(median = median(pre.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), color = species)) +
  facet_grid(species~.) +
  guides(color = 'none') +
  labs(x = 'log10(Pre-Bout Prop. Above 50m)') +
  theme_classic()

## post-stats
post.dist.ptt <-
  surface_use %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(post.prop.50+0.001), color = species, group = ptt)) +
  geom_vline(data = surface_use %>% group_by(species, ptt) %>% summarize(median = median(post.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), color = species)) +
  facet_grid(species~.) +
  labs(x = 'log10(Post-Bout Prop. Above 50m)') +
  theme_classic()

base.dist.ptt <-
  baseline_surface_use %>% 
  filter(!is.na(cluster)) %>% 
  # mutate(cluster = factor(cluster, ordered = F)) %>% 
  ggplot() +
  geom_freqpoly(aes(x = log10(base.prop.50+0.001), color = species, group = ptt)) +
  geom_vline(data = baseline_surface_use %>% filter(!is.na(cluster)) %>% group_by(species, ptt) %>% summarize(median = median(base.prop.50, na.rm = T)), aes(xintercept = log10(median+0.001), color = species)) +
  facet_grid(species~.) +
  labs(x = 'log10(Baseline Prop. Above 50m)') +
  theme_classic()

pre.dist.ptt + post.dist.ptt

base.dist.ptt


## remake fig 6 (thorrold et al., 2014) ------------------------------------

# calculate weighted quantiles for boxplots
## simple way
pre_weighted_boxplots <-
  surface_use %>% 
  mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>%
  group_by(ptt, species, dn) %>% 
  summarize(
    ymin.ptt = min(pre.prop.50),
    lower.ptt = quantile(pre.prop.50, 0.25),
    median.ptt = median(pre.prop.50),
    upper.ptt = quantile(pre.prop.50, 0.75),
    ymax.ptt = max(pre.prop.50),
    n.ptt = n()) %>% 
  group_by(species, dn) %>% 
  summarize(ymin = min(ymin.ptt),
            lower = quantile(lower.ptt, 0.25),
            middle = median(median.ptt),
            upper = quantile(upper.ptt, 0.75),
            ymax = max(ymax.ptt),
            n = sum(n.ptt), .groups = 'drop')
## fancy way --
# pre_weighted_boxplots <-
#   surface_use %>% 
#   mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>% 
#   group_by(ptt) %>% 
#   mutate(weight = 1/n()) %>% 
#   ungroup() %>% 
#   group_by(species, dn) %>% 
#   summarize(
#     ymin   = Hmisc::wtd.quantile(pre.prop.50, weights = weight, probs = 0.00),
#     lower  = Hmisc::wtd.quantile(pre.prop.50, weights = weight, probs = 0.25),
#     middle = Hmisc::wtd.quantile(pre.prop.50, weights = weight, probs = 0.50),
#     upper  = Hmisc::wtd.quantile(pre.prop.50, weights = weight, probs = 0.75),
#     ymax   = Hmisc::wtd.quantile(pre.prop.50, weights = weight, probs = 1.00), .groups = 'drop')

post_weighted_boxplots <-
  surface_use %>% 
  mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>%
  group_by(ptt, species, dn) %>% 
  summarize(
    ymin.ptt = min(post.prop.50),
    lower.ptt = quantile(post.prop.50, 0.25),
    median.ptt = median(post.prop.50),
    upper.ptt = quantile(post.prop.50, 0.75),
    ymax.ptt = max(post.prop.50),
    n = n()) %>% 
  group_by(species, dn) %>% 
  summarize(ymin = min(ymin.ptt),
            lower = quantile(lower.ptt, 0.25),
            middle = median(median.ptt),
            upper = quantile(upper.ptt, 0.75),
            ymax = max(ymax.ptt), .groups = 'drop',
            n = sum(n)) 


base_weighted_boxplots <-
  baseline_surface_use %>%
  mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>% 
  filter(!is.na(dn)) %>% 
  group_by(ptt, species, dn) %>% 
  summarize(
    ymin.ptt = min(base.prop.50, na.rm = T),
    lower.ptt = quantile(base.prop.50, 0.25, na.rm = T),
    median.ptt = median(base.prop.50, na.rm = T),
    upper.ptt = quantile(base.prop.50, 0.75, na.rm = T),
    ymax.ptt = max(base.prop.50, na.rm = T),
    n = n()) %>% 
  group_by(species, dn) %>% 
  summarize(ymin = min(ymin.ptt),
            lower = quantile(lower.ptt, 0.25),
            middle = median(median.ptt),
            upper = quantile(upper.ptt, 0.75),
            ymax = max(ymax.ptt),
            n = sum(n), .groups = 'drop') 

surface_use %>% 
  mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>% 
  # calculate weights
  group_by(ptt) %>% 
  mutate(weight = 1/n()) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_violin(aes(x = -1, y = pre.prop.50), fill = 'coral4', alpha = 0.5, width = 0.75) +
  geom_violin(aes(x = 0, y = post.prop.50), fill = 'steelblue4', alpha = 0.5, width = 0.75) +
  geom_violin(data = baseline_surface_use %>% mutate(dn = if_else(dn == 'd', 'Daylight', 'Nighttime')) %>% filter(!is.na(dn)), aes(x = 1, y = base.prop.50), fill = 'grey25', alpha = 0.5, width = 0.75) +
  geom_boxplot(data = pre_weighted_boxplots, aes(x = -1, ymin = ymin, lower = lower,
                   middle = middle, upper = upper, ymax = ymax), stat = 'identity', color = 'black', fill = NA, width = 0.15) + 
  geom_boxplot(data = post_weighted_boxplots, aes(x = 0, ymin = ymin, lower = lower,
                                                 middle = middle, upper = upper, ymax = ymax), stat = 'identity', color = 'black', fill = NA, width = 0.15) +   
  geom_boxplot(data = base_weighted_boxplots, aes(x = 1, ymin = ymin, lower = lower,
                                                 middle = middle, upper = upper, ymax = ymax), stat = 'identity', color = 'black', fill = NA, width = 0.15) + 
  facet_grid(dn~species) + 
  labs(x = '', y = 'Proportion of Time') +
  theme_classic() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank())


