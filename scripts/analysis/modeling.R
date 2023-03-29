
# setup -------------------------------------------------------------------

library(foreign)
library(nnet)
library(reshape2)
library(mclogit)

# get cluster information
source('scripts/analysis/heirarchical_clustering.R')

# formatting ----------------------------------------------------------------

combo_mod <- 
  high_res %>% 
  # remove rare clusters
  filter(cluster <= 5 ) %>% 
  # set the reference level to the most common cluster
  mutate(
    clus2 = relevel(
      as.factor(cluster),
      ref = 1),
    # reference as blue sharks because we have more tags from them
    species = relevel(
      as.factor(species),
      ref = 'P.glauca'),
    ptt = factor(ptt),
    depth = bathy,
    l2 = lunar) %>% 
  # select response and predictors
  dplyr::select(c(Date,
                  kode,
                  species,
                  ptt,
                  cluster,
                  clus2,
                  latitude,
                  longitude,
                  depth,
                  sst,
                  sst_sd,
                  ssh,
                  ssh_sd,
                  eke,
                  ild.5,
                  n2,
                  chloro,
                  l2,
                  lunar)) %>% 
  # convert lunar to factor variable based on phase
  mutate(
    lunar = case_when(
      lunar >= 0.75 ~ 2,
      lunar < 0.75 ~ 1),
    lunar = relevel(
      # setting lunar to a factor with reference 1
      as.factor(lunar), ref = 1))


# cross correlation analysis ----------------------------------------------

c.cor <- cor(combo_mod[,c(9:18)], use = "complete.obs")
library(PerformanceAnalytics)
chart.Correlation(combo_mod[,c(9:18)])


# model construction ------------------------------------------------------

# load a function for model construction and AIC evaluation
source('scripts/analysis/get_weights.R')

get_weight(base = NULL,
           test = c('sst',
                    'sst_sd',
                    'ssh',
                    'ssh_sd',
                    'eke',
                    'n2',
                    'depth',
                    'ild.5',
                    'chloro',
                    'lunar',
                    'species'))
# ssh is first

get_weight(base = c('ssh'),
           test = c('sst_sd',
                    'ssh_sd',
                    'eke',
                    'n2',
                    'ild.5',
                    'lunar',
                    'species'))

# ssh_sd is second

get_weight(base = c('ssh', 'ssh_sd'),
           test = c('sst_sd',
                    'n2',
                    'ild.5',
                    'lunar',
                    'species'))
# lunar is third

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar'),
           test = c('sst_sd',
                    'n2',
                    'ild.5',
                    'species'))
# species is fourth

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species'),
           test = c('sst_sd',
                    'n2',
                    'ild.5'))


## interactions ------------------------------------------------------------

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2'),
           test = c('ssh:ssh_sd', 
                    'ssh:lunar',
                    'ssh:n2',
                    'ssh:species',
                    'ssh_sd:lunar',
                    'ssh_sd:n2',
                    'ssh_sd:species',
                    'lunar:n2',
                    'lunar:species',
                    'n2:species'))

# none of the interactions terms improve the base model

# predictor performance ---------------------------------------------------

# what is the significance level of each predictor?
m.mod <- 
  mblogit(
    formula = clus2 ~ ssh + ssh_sd + lunar + species,
    random = ~1|ptt,
    data = combo_mod)

# How much deviance is explained by each predictor
null <- 
  mblogit(
    formula = clus2 ~ 1,
    data = combo_mod)

m_random <-
  mblogit(
    formula = clus2 ~ 1,
    random = ~1|ptt,
    data = combo_mod)

m_ssh <-
  mblogit(
    formula = clus2 ~ ssh,
    data = combo_mod)

m_ssh_sd <-
  mblogit(
    formula = clus2 ~ ssh_sd,
    data = combo_mod)

m_lunar <-
  mblogit(
    formula = clus2 ~ lunar,
    data = combo_mod)

m_species <-
  mblogit(
    formula = clus2 ~ species,
    data = combo_mod)

1 - (m_ssh_sd$deviance/null$deviance)




get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species'),
           test = c('n2 + n2:ssh'))

