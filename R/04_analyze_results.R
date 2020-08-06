## This script compiles the results from the local, regional and interaction scale analyses and presents
# the mean sums of squares for each term

library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)


# load the results from adonis 

local <- readRDS("Results/local.rds") %>% 
  purrr::map_df(~as.data.frame(.x)) %>%  
  mutate(run = rep(1:1000, each = 5))

regional <- readRDS("Results/regional.rds") %>% 
  purrr::map_df(~as.data.frame(.x)) %>%  
  mutate(run = rep(1:1000, each = 6))

interaction <- readRDS("Results/interaction.rds") %>% 
  purrr::map_df(~as.data.frame(.x)) %>%  
  mutate(run = rep(1:1000, each = 17))

#### local results ####

total_local <- local %>%
  filter(term == 'Total') %>% 
  rename(SS_total = SumOfSqs) %>% 
  select(SS_total, run)

# compute percentages for figure
local_percentage <- local %>% 
  filter(term != 'Total') %>% 
  left_join(total_local) %>% 
  mutate(per_ss = (SumOfSqs/SS_total)*100) %>% 
  mutate(per_ss = ifelse(term == 'Residual', 100-per_ss, per_ss)) %>% 
  mutate(term = ifelse(term == 'Residual', "Full", term))

# Figure A1#

nice_name_local <- c("Canopy Cover",  "Water Volume", "Detritus", "Full Model")

term_name_local <- data.frame(term = as.character(unique(local_percentage$term)), nice_name = nice_name_local)

A1_hist <- local_percentage %>% 
  left_join(term_name_local) %>% 
  ggplot(aes(x = per_ss)) + geom_histogram() + facet_wrap(~nice_name, scales = 'free') + 
  theme_cowplot() + ylab("Number of runs") + xlab("Percentage of total sums of squares") +
  theme(strip.background = element_blank())

ggsave(A1_hist, filename = "Figures/A1_hist.jpeg")

#Table 1 fine grained analysis

mean_ss_local <- local %>% 
  group_by(term) %>% 
  summarise(mean_ss = mean(SumOfSqs))

total_ss_local <- filter(mean_ss_local, term == 'Total')

# Percentage sums of squares

# We calcualte first the mean sums of squares and then calculate the proportion based on the mean total sums of squares
# we do this procedure in this order as some of the distributions were skewed and 
#allowed us to take advantage of the central limit theorem to ensure that the division occurs on normally distributed means

mean_ss_local %>% 
  filter(!term %in% c("Residual", "Total")) %>% 
  mutate(percentage_ss = 100*(mean_ss/total_ss_local$mean_ss))

full_local <- 
  sum(filter(mean_ss_local, !term %in% c("Residual", "Total"))$mean_ss)/total_ss_local$mean_ss

  
#### regional ####  
  
total_regional <- regional %>%
  filter(term == 'Total') %>% 
  rename(SS_total = SumOfSqs) %>% 
  select(SS_total, run)

# compute percentages for figure

regional_percentage <- regional %>% 
  filter(term != 'Total') %>% 
  left_join(total_regional) %>% 
  mutate(per_ss = (SumOfSqs/SS_total)*100) %>% 
  mutate(per_ss = ifelse(term == 'Residual', 100-per_ss, per_ss)) %>% 
  mutate(term = ifelse(term == 'Residual', "Full", term))

# Figure A2#

nice_name_regional <- c("Mean Diurnal\nRange", "Temperature\nSeasonality", "Precipitation of\nDriest Quarter", 
                        "Precipitation\nSeasonality", "Full Model")

term_name_regional <- data.frame(term = as.character(unique(regional_percentage$term)), nice_name = nice_name_regional)

A2_hist <- regional_percentage %>% 
  left_join(term_name_regional) %>% 
  ggplot(aes(x = per_ss)) + geom_histogram() + facet_wrap(~nice_name, scales = 'free') + 
  theme_cowplot() + ylab("Number of runs") + xlab("Percentage of total sums of squares") +
  theme(strip.background = element_blank())

ggsave(A2_hist, filename = "Figures/A2_hist.jpeg")

#Table 1 fine grained analysis

mean_ss_regional <- regional %>% 
  group_by(term) %>% 
  summarise(mean_ss = mean(SumOfSqs))

total_ss_regional <- filter(mean_ss_regional, term == 'Total')

# Percentage sums of squares

# We calcualte first the mean sums of squares and then calculate the proportion based on the mean total sums of squares
# we do this procedure in this order as some of the distributions were skewed and 
#allowed us to take advantage of the central limit theorem to ensure that the division occurs on normally distributed means


mean_ss_regional %>% 
  filter(!term %in% c("Residual", "Total")) %>% 
  mutate(percentage_ss = 100*(mean_ss/total_ss_regional$mean_ss))
  
full_regional <- 
  sum(filter(mean_ss_regional, !term %in% c("Residual", "Total"))$mean_ss)/total_ss_regional$mean_ss
  

#### interaction ####  

total_interaction <- interaction %>%
  filter(term == 'Total') %>% 
  rename(SS_total = SumOfSqs) %>% 
  select(SS_total, run)

# compute percentages for figure

interaction_percentage <- interaction %>% 
  filter(term != 'Total') %>% 
  left_join(total_interaction) %>% 
  mutate(per_ss = (SumOfSqs/SS_total)*100) %>% 
  mutate(per_ss = ifelse(term == 'Residual', 100-per_ss, per_ss)) %>% 
  mutate(term = ifelse(term == 'Residual', "Full", term))

nice_name_interaction <- c("Detritus", "Canopy Cover", "Water Volume", "Mean Diurnal\nRange", "Temperature\nSeasonality", "Precipitation of\nDriest Quarter", 
              "Precipitation\nSeasonality", "Detritus:\nMean Diurnal\nRange", "Detritus:\nTemperature\nSeasonality", "Detritus:\nPrecipitation of\nDriest Quarter", 
              "Detritus:\nPrecipitation\nSeasonality", "Water Volume:\nMean Diurnal\nRange", "Water Volume:\nTemperature\nSeasonality", "Water Volume:\nPrecipitation of\nDriest Quarter", 
              "Water Volume:\nPrecipitation\nSeasonality", "Full Model")

term_name_interaction <- data.frame(term = as.character(unique(interaction_percentage$term)), nice_name = nice_name_interaction, 
                        colour = c(rep("local", 3), rep('reg', 4), rep('int', 8),"full"))

A9_int_hist <- interaction_percentage %>% 
  left_join(term_name_interaction) %>% 
  mutate(nice_name = factor(nice_name, levels = nice_name_interaction)) %>% 
  ggplot(aes(per_ss, fill = colour)) + geom_histogram() + facet_wrap(~nice_name, scales = 'free') +
  theme_cowplot() + ylab("Number of runs") + xlab("Percentage of total sums of squares") +
  theme(strip.background = element_blank(), 
        legend.position = 'none') + scale_fill_manual(values = c("black","#445F94","#D0603D",  "#C0F56E"))

ggsave(A9_int_hist, filename = "Figures/A9_int_hist.jpeg", width = 10, height =10)


## Table 2##

## total ss ##

mean_ss_interaction <- interaction %>% 
  group_by(term) %>% 
  summarise(mean_ss = mean(SumOfSqs))

total_ss_interaction <- filter(mean_ss_interaction, term == 'Total')

# Percentage sums of squares

# We calcualte first the mean sums of squares and then calculate the proportion based on the mean total sums of squares
# we do this procedure in this order as some of the distributions were skewed and 
#allowed us to take advantage of the central limit theorem to ensure that the division occurs on normally distributed means


mean_ss_interaction %>% 
  filter(!term %in% c("Residual", "Total")) %>% 
  mutate(percentage_ss = 100*(mean_ss/total_ss_interaction$mean_ss))

full_interaction <- 
  sum(filter(mean_ss_interaction, !term %in% c("Residual", "Total"))$mean_ss)/total_ss_interaction$mean_ss

## partition ##

local_interaction <- 
  sum(filter(mean_ss_interaction, term %in% c("final_canopy", "log(total_detritus)", "log(actual_water)"))$mean_ss)/total_ss_interaction$mean_ss

regional_interaction <- 
  sum(filter(mean_ss_interaction, term %in% c("BC4", "BC2", "BC15", "BC17"))$mean_ss)/total_ss_interaction$mean_ss

interaction_interaction <- 
  sum(filter(mean_ss_interaction, term %in% c("log(actual_water):BC15", "log(actual_water):BC17", "log(actual_water):BC4", "log(actual_water):BC2",
                                              "log(total_detritus):BC15", "log(total_detritus):BC17", "log(total_detritus):BC4", "log(total_detritus):BC2"
                                              ))$mean_ss)/total_ss_interaction$mean_ss


