
# setup -------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)

source('scripts/wrangling/identify_mesoBouts.R')


# exploration -------------------------------------------------------------

# visualize distribution
temp_duration %>% 
  pivot_longer(cols = deg.14.min:deg.18.min, names_to = 'temp', values_to = 'duration') %>% 
  ggplot() +
  geom_histogram(aes(x = log10(duration)), color = 'black', fill = 'grey', bins = 100) +
  coord_cartesian() +
  facet_grid(temp ~.) +
  theme_classic()

# distribution by cluster
temp.freqpoly.cluster <- 
  temp_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(deg.18.min), color = cluster, group = ptt)) +
  geom_vline(data = temp_duration %>% group_by(cluster, species) %>% summarize(median = median(deg.18.min+0.001, na.rm =T), .groups = 'drop'), aes(xintercept = log10(median), linetype = species), alpha = 0.5) +
  facet_grid(cluster~.) +
  coord_cartesian() +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  theme_classic()

ggsave(
  temp.freqpoly.cluster,
  file = 'products/figures/revisions/exploration/18C_freqpoly_cluster.png',
  height = 6,
  width = 8,
  dpi = 300)

# distribution by species
temp.freqpoly.species <- 
  temp_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(deg.18.min), color = species, group = ptt)) +
  geom_vline(data = temp_duration %>% group_by(species, ptt) %>% summarize(median = median(deg.18.min+0.001, na.rm =T), .groups = 'drop'), aes(xintercept = log10(median), color = species), alpha = 0.5) +
  facet_grid(species~.) +
  coord_cartesian() +
  theme_classic()

ggsave(
  temp.freqpoly.species,
  file = 'products/figures/revisions/exploration/18C_freqpoly_species.png',
  height = 6,
  width = 8,
  dpi = 300)


