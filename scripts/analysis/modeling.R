
# setup -------------------------------------------------------------------

library(foreign)
library(nnet)
library(reshape2)
library(mclogit)

# get cluster information
source('scripts/analysis/heirarchical_clustering.R')

s <- summary(high_res)

# formatting ----------------------------------------------------------------

combo_mod <- 
  high_res %>% 
  # remove rare clusters
  filter(cluster <= 5 ) %>% 
  # reset cluster names
  mutate(
    depth = bathy,
    ptt = factor(ptt),
    species = factor(species,
                     levels = c('P.glauca', 'I.oxyrinchus')),
    clus2 = case_when(
      cluster == 1 ~ 'DVM 1',
      cluster == 2 ~ 'Epipelagic',
      cluster == 3 ~ 'DVM 2',
      cluster == 4 ~ 'DVM 3',
      cluster == 5 ~ 'DVM 4'),
    clus2 = factor(clus2, levels = c(
      'DVM 1',
      'Epipelagic',
      'DVM 2', 
      'DVM 3',
      'DVM 4')
    )) %>%
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
                  lunar)) %>% 
  # convert lunar to factor variable based on phase
  mutate(
    lunar = case_when(
      lunar >= 0.75 ~ 1,
      lunar < 0.75 ~ 0))


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

# n2 is fifth

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2'),
           test = c('sst_sd'))

# null model is best, moving on to interactions

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

# ssh:n2 interaction improves model

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2',
                    'ssh:n2'),
           test = c('ssh:ssh_sd', 
                    'ssh:lunar',
                    'ssh:species',
                    'ssh_sd:lunar',
                    'ssh_sd:n2',
                    'ssh_sd:species',
                    'lunar:n2',
                    'lunar:species',
                    'n2:species'))

# ssh:lunar is next improvement, but removes effect of ssh on epipelagic
# some argument for stopping here

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2',
                    'ssh:n2',
                    'ssh:lunar'),
           test = c('ssh:ssh_sd', 
                    'ssh:species',
                    'ssh_sd:lunar',
                    'ssh_sd:n2',
                    'ssh_sd:species',
                    'lunar:n2',
                    'lunar:species',
                    'n2:species'))

# ssh_sd:lunar is next improvement

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2',
                    'ssh:n2',
                    'ssh:lunar',
                    'ssh_sd:lunar'),
           test = c('ssh:ssh_sd', 
                    'ssh:species',
                    'ssh_sd:n2',
                    'ssh_sd:species',
                    'lunar:n2',
                    'lunar:species',
                    'n2:species'))

# ssh_sd:species is next

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2',
                    'ssh:n2',
                    'ssh:lunar',
                    'ssh_sd:lunar',
                    'ssh_sd:species'),
           test = c('ssh:ssh_sd', 
                    'ssh:species',
                    'ssh_sd:n2',
                    'lunar:n2',
                    'lunar:species',
                    'n2:species'))

# next is effect of n2 on species

get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'species',
                    'n2',
                    'ssh:n2',
                    'ssh:lunar',
                    'ssh_sd:lunar',
                    'ssh_sd:species',
                    'n2:species'),
           test = c('ssh:ssh_sd', 
                    'ssh:species',
                    'ssh_sd:n2',
                    'lunar:n2',
                    'lunar:species'))

# what happens if we remove species as a fixed effect
get_weight(base = c('ssh',
                    'ssh_sd',
                    'lunar',
                    'n2',
                    'ssh:n2',
                    'ssh:lunar',
                    'ssh_sd:lunar',
                    'ssh_sd:species',
                    'n2:species'),
           test = c(
             'species'
           ))

# species is best included as an interaction term

# predictor performance ---------------------------------------------------

# what is the significance level of each predictor?
m.mod <- 
  mblogit(
    formula = clus2 ~ ssh +
      ssh_sd +
      lunar +
      n2 +
      # species +
      ssh:n2 +
      ssh:lunar +
      ssh_sd:lunar +
      ssh_sd:species + 
      species:n2,
    random = ~1|ptt,
    data = combo_mod)

summary(m.mod)

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


1 - (m_ssh$deviance/null$deviance) # should create a loop to do this automatically someday


# vizualizing model results -----------------------------------------------

## Effect of SSH -----------------------------------------------------------

# vary ssh while holding other predictors constant
p_ssh <- 
  tibble(
    ssh = rep(seq(-1,1, 0.025),times = 34), 
    lunar = as.factor(rep(rep(c(0,1),each = 81), times = 17)),
    # hold ssh_sd at the mean value
    ssh_sd = rep(c(0.02605838), times = 2754),
    n2 = rep(c(1.472e-05), times = 2754),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 162))) %>% 
  # bind to another tibble with median values of n2
  rbind(
    tibble(
      ssh = rep(seq(-1,1, 0.025),times = 34), 
      lunar = as.factor(rep(rep(c(0,1),each = 81), times = 17)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(c(0.02605838), times = 2754),
      n2 = rep(c(3.096e-05), times = 2754),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 162)))) %>% 
  # bind to another tibble with high values of n2
  rbind(
    tibble(
      ssh = rep(seq(-1,1, 0.025),times = 34), 
      lunar = as.factor(rep(rep(c(0,1),each = 81), times = 17)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(c(0.02605838), times = 2754),
      n2 = rep(c(7.457e-05), times = 2754),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 162)))) %>% 
  left_join(high_res %>% 
              dplyr::select(ptt, species) %>%
              distinct() %>% 
              mutate(ptt = factor(ptt)),
            by = 'ptt')

# add predictions based on the best performing model
p_ssh <- 
  cbind(
    p_ssh,
    predict(m.mod,
            newdata = p_ssh,
            type = "response",
            se.fit = FALSE))
# condense predictions to a single column (value) dependent on cluster (variable)
lssh <- 
  reshape2::melt(p_ssh,
                 id.vars = c("ssh",
                             "lunar",
                             "ssh_sd",
                             'n2',
                             "ptt",
                             "species"),
                 measure.vars = c("DVM 1",
                                  "Epipelagic",
                                  "DVM 2",
                                  "DVM 3",
                                  "DVM 4")) %>% 
  mutate(
    lunar = factor(lunar))

rm(p_ssh)

# create a mean line by averaging all individuals 
opall <- 
  lssh %>% 
  group_by(species, 
           variable,
           ssh,
           lunar,
           n2) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()

# plot relationships for blues
opall %>% 
  # remove values of n2 if needed
  filter(n2 == 3.096e-05) %>% 
  filter(species == 'I.oxyrinchus') %>%
  # filter(species == 'P.glauca') %>%
  filter(variable != 'Epipelagic') %>%
  ggplot(aes(x = ssh,
             y = value,
             color = variable)) + 
  geom_ribbon(aes(x = ssh,
                  ymin = LL,
                  ymax = UL,
                  fill = variable),
              color = 'white',
              linewidth = 0.25,
              alpha = 0.25) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                # "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                               # "#FFFF5CFF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) +
  facet_wrap(n2~lunar) + 
  labs(fill = "Cluster") +
  guides(color = 'none') +
  ylab("Probability") +
  theme_linedraw()

opall %>% 
  group_by(species,
           variable,
           # we have to include lunar as a grouping variable to keep consistency 
           n2,
           lunar) %>%
  # find the maximum predicted probability
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(lunar == 0 & n2 == 0.00003096)

## Effect of SSH SD --------------------------------------------------------

p_ssh_sd <- 
  tibble(
    ssh = rep(c(-0.0500), times = 1666), 
    lunar = as.factor(rep(c(0), times = 1666)),
    # hold ssh_sd at the mean value
    ssh_sd = rep(seq(0.010611,0.035090, by = 0.00025), times = 17),
    n2 = rep(c(3.096e-05),times = 1666),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 98))) %>% 
  # repeat the same dataset with high lunar illumination
  rbind(
    tibble(
      ssh = rep(c(-0.0500), times = 1666), 
      lunar = as.factor(rep(c(1), times = 1666)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(seq(0.010611,0.035090, by = 0.00025), times = 17),
      n2 = rep(c(3.096e-05),times = 1666),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 98)))) %>% 
  left_join(high_res %>% 
              dplyr::select(ptt, species) %>%
              distinct() %>% 
              mutate(ptt = factor(ptt),
                     species = factor(species)),
            by = 'ptt') %>% 
  # add predicted probabilities based on model results
  cbind(
    predict(m.mod,
            newdata = .,
            type = "response",
            se.fit = FALSE)) %>% 
  # pivot the table into a longer format
  reshape2::melt(id.vars = c("ssh",
                             "lunar",
                             "ssh_sd",
                             'n2',
                             "ptt",
                             "species"),
                 measure.vars = c("DVM 1",
                                  "Epipelagic",
                                  "DVM 2",
                                  "DVM 3",
                                  "DVM 4")) %>% 
  # distill data from all individuals into an average, minimum and maximum predicted probability for each unique value of n2 
  group_by(species, 
           variable,
           ssh_sd,
           lunar) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()

# visualize the relationship

p_ssh_sd %>% 
  filter(lunar == 1) %>% 
  filter(variable != 'Epipelagic') %>%
  ggplot(aes(x = ssh_sd,
             y = value,
             color = variable)) + 
  geom_ribbon(aes(x = ssh_sd,
                  ymin = LL,
                  ymax = UL,
                  fill = variable),
              color = 'white',
              linewidth = 0.25,
              alpha = 0.25) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                # "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                               # "#FFFF5CFF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) +
  facet_grid(lunar~species) + 
  labs(fill = "Cluster") +
  guides(color = 'none') +
  ylab("Probability") +
  theme_linedraw()

# find the conditions under which each cluster reaches maximum predicted probability

p_ssh_sd %>% 
  group_by(species,
           variable,
           # we have to include lunar as a grouping variable to keep consistency 
           lunar) %>%
  # find the maximum predicted probability
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(lunar == 0)


## Effect of Lunar Illumination --------------------------------------------
p_lunar <- 
  tibble(
    ssh = rep(c(-0.0500), times = 17), 
    lunar = as.factor(rep(c(0), times = 17)),
    # hold ssh_sd at the mean value
    ssh_sd = rep(c(0.02605838), times = 17),
    n2 = rep(c(3.096e-05),times = 17),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 1))) %>% 
  # repeat the same dataset with high lunar illumination
  rbind(
    tibble(
      ssh = rep(c(-0.0500), times = 17), 
      lunar = as.factor(rep(c(1), times = 17)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(c(0.02605838), times = 17),
      n2 = rep(c(3.096e-05),times = 17),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 1)))) %>% 
  left_join(high_res %>% 
              dplyr::select(ptt, species) %>%
              distinct() %>% 
              mutate(ptt = factor(ptt),
                     species = factor(species)),
            by = 'ptt') %>% 
  # add predicted probabilities based on model results
  cbind(
    predict(m.mod,
            newdata = .,
            type = "response",
            se.fit = FALSE)) %>% 
  # pivot the table into a longer format
  reshape2::melt(id.vars = c("ssh",
                             "lunar",
                             "ssh_sd",
                             'n2',
                             "ptt",
                             "species"),
                 measure.vars = c("DVM 1",
                                  "Epipelagic",
                                  "DVM 2",
                                  "DVM 3",
                                  "DVM 4")) %>% 
  # distill data from all individuals into an average, minimum and maximum predicted probability for each unique value of n2 
  group_by(species, 
           variable,
           lunar) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()

p_lunar %>% 
  filter(variable != 'Epipelagic') %>%
  ggplot(aes(x = lunar,
             y = value,
             color = variable)) + 
  geom_pointrange(aes(x = lunar,
                  ymin = LL,
                  ymax = UL,
                  color = variable),
                  position = position_dodge(width = 0.75)) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                # "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                               # "#FFFF5CFF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) +
  labs(x = 'lunar illumination') +
  facet_wrap(~species) + 
  labs(color = "Cluster") +
  ylab("Probability") +
  theme_linedraw()



## Effect of N2 ----
p_n2 <- 
  tibble(
    ssh = rep(c(-0.0500), times = 2040), 
    lunar = as.factor(rep(c(0), times = 2040)),
    # hold ssh_sd at the mean value
    ssh_sd = rep(c(0.02605838), times = 2040),
    n2 = rep(seq(1.472e-05,7.457e-05, by = 0.05e-05),times = 17),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 120))) %>% 
  # repeat the same dataset with high lunar illumination
  rbind(
    tibble(
      ssh = rep(c(-0.0500), times = 2040), 
      lunar = as.factor(rep(c(1), times = 2040)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(c(0.02605838), times = 2040),
      n2 = rep(seq(1.472e-05,7.457e-05, by = 0.05e-05),times = 17),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 120)))) %>% 
  left_join(high_res %>% 
              dplyr::select(ptt, species) %>%
              distinct() %>% 
              mutate(ptt = factor(ptt),
                     species = factor(species)),
            by = 'ptt') %>% 
  # add predicted probabilities based on model results
  cbind(
    predict(m.mod,
            newdata = .,
            type = "response",
            se.fit = FALSE)) %>% 
  # pivot the table into a longer format
  reshape2::melt(id.vars = c("ssh",
                             "lunar",
                             "ssh_sd",
                             'n2',
                             "ptt",
                             "species"),
                 measure.vars = c("DVM 1",
                                  "Epipelagic",
                                  "DVM 2",
                                  "DVM 3",
                                  "DVM 4")) %>% 
  # distill data from all individuals into an average, minimum and maximum predicted probability for each unique value of n2 
  group_by(species, 
           variable,
           n2,
           lunar) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()
  
# visualize the relationship

p_n2 %>% 
  filter(lunar == 0) %>%
  filter(variable != 'Epipelagic') %>%
  ggplot(aes(x = n2,
             y = value,
             color = variable)) + 
  geom_ribbon(aes(x = n2,
                  ymin = LL,
                  ymax = UL,
                  fill = variable),
              color = 'white',
              linewidth = 0.25,
              alpha = 0.25) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                # "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                               # "#FFFF5CFF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) +
  facet_grid(lunar~species) + 
  labs(x = 'Buoyancy Frequency (s-1)', fill = "Cluster") +
  guides(color = 'none') +
  ylab("Probability") +
  theme_linedraw()

# find the conditions under which each cluster reaches maximum predicted probability

p_n2 %>% 
  group_by(species,
           variable,
           # we have to include lunar as a grouping variable to keep consistency 
           lunar) %>%
  # find the maximum predicted probability
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(lunar == 0)


