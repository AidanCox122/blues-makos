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


