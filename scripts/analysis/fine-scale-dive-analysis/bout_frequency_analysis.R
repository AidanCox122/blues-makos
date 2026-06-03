
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
