library(vegan)
library(SYNCSA)
library(broom)
library(tidyr)
source("R/functions.R")

###### interaction #######

## load in data
abundance_interaction <- readRDS("Data/abundance_interaction.rds")
region_organized_interaction <- readRDS("Data/region_organized_interaction.rds")
trait_axis_interaction <- readRDS("Data/trait_axis_interaction.rds")
all_environment_interaction <- readRDS("Data/all_environment_interaction.rds")


# run adonis with original data #

cwm <- Fuzzit4a(abundance_interaction, trait_axis_interaction) %>% 
  dplyr::select(-region.id) 

sampled_dist_matrix_pwt <- dist(cwm, method = "euclidean")

output_adonis_interaction <- adonis2(sampled_dist_matrix_pwt ~  log(total_detritus) + final_canopy + log(actual_water) + BC2 + BC4 + BC17 + BC15 +
                               log(total_detritus):BC2 + log(total_detritus):BC4 + log(total_detritus):BC17 + log(total_detritus):BC15 +
                               log(actual_water):BC2 + log(actual_water):BC4 + log(actual_water):BC17 + log(actual_water):BC15,
                             data = all_environment_interaction, by = 'term', permutations = 2)

interaction_adonis <- tidy(output_adonis_interaction)

## permute traits ##

f_table <- data.frame()

for(p in 1:999){
  
  p_abun_mat <- permut.row.matrix(data = abundance_interaction, strata = region_organized_interaction, seqpermutation = NULL)$permut.matrix

  p_cwm <- Fuzzit4a(p_abun_mat, trait_axis_interaction) %>% 
    dplyr::select(-region.id) 
  
  p_sampled_dist_matrix_pwt <- dist(p_cwm, method = "euclidean")
  
  p_output_adonis <- adonis2(p_sampled_dist_matrix_pwt ~  log(total_detritus) + final_canopy + log(actual_water) + BC2 + BC4 + BC17 + BC15 +
                               log(total_detritus):BC2 + log(total_detritus):BC4 + log(total_detritus):BC17 + log(total_detritus):BC15 +
                               log(actual_water):BC2 + log(actual_water):BC4 + log(actual_water):BC17 + log(actual_water):BC15,
                             data = all_environment_interaction, by = 'term', permutations = 2)
  Fs = p_output_adonis$F
  
  f_table <- rbind(f_table, data.frame(Fs, permutation = p))

}

original_interaction_adonis <- data.frame(term = interaction_adonis$term, F_orig = output_adonis_interaction$F)

spread_f_table <- f_table %>% 
  mutate(term = rep(interaction_adonis$term, 999)) %>% 
  spread(key = term, value = Fs) %>% 
  select(interaction_adonis$term)

## greater than

ps <- sapply(1:15, FUN = function(x) ((sum(output_adonis_interaction$F[x] > spread_f_table[,x]) +1)/1000))

all_output <- data.frame(term = interaction_adonis$term[1:length(ps)], 
                              SS = output_adonis_all$SumOfSqs[1:length(ps)], P = ps) %>% 
  mutate(P= ifelse(P> 0.5, 1-P, P)) %>% 
  mutate(P = 2*P)

all_output %>%View()

##### Local scale#######

cwm <- Fuzzit4a(abundance_interaction, trait_axis_interaction) %>% 
  dplyr::select(-region.id) 

sampled_dist_matrix_pwt <- dist(cwm, method = "euclidean")

output_adonis_local <- adonis2(sampled_dist_matrix_pwt ~  log(total_detritus) + final_canopy + log(actual_water),
                             data = all_environment_interaction, by = 'margin', permutations = 2)

local_adonis <- tidy(output_adonis_local)

f_table <- data.frame()

for(p in 1:999){
  
  p_abun_mat <- permut.row.matrix(data = abundance_interaction, strata = region_organized_interaction, seqpermutation = NULL)$permut.matrix

  p_cwm <- Fuzzit4a(p_abun_mat, trait_axis_interaction) %>% 
    dplyr::select(-region.id) 
  
  p_sampled_dist_matrix_pwt <- dist(p_cwm, method = "euclidean")
  
  p_output_adonis <- adonis2(p_sampled_dist_matrix_pwt ~  log(total_detritus) + final_canopy + log(actual_water),
                             data = all_environment_interaction, by = 'margin', permutations = 2)
  Fs = p_output_adonis$F
  
  f_table <- rbind(f_table, data.frame(Fs, permutation = p))
  
}

original_local_adonis <- data.frame(term = local_adonis$term, F_orig = output_adonis_local$F)

spread_f_table <- f_table %>% 
  mutate(term = rep(local_adonis$term, 999)) %>% 
  spread(key = term, value = Fs) %>% 
  select(local_adonis$term)

## greater than

ps <- sapply(1:3, FUN = function(x) ((sum(output_adonis_local$F[x] > spread_f_table[,x]) +1)/1000))

all_output_local <- data.frame(term = local_adonis$term[1:length(ps)], 
                         SS = output_adonis_local$SumOfSqs[1:length(ps)], P = ps) %>% 
  mutate(P= ifelse(P> 0.5, 1-P, P)) %>% 
  mutate(P = 2*P)

all_output_local %>%View()

###### region ##########

abundance_region <- readRDS("Data/abundance_region.rds") 
trait_axis_region <- readRDS("Data/trait_axis_region.rds")
all_environment_region <- readRDS("Data/all_environment_region.rds")


cwm <- Fuzzit4a(abundance_region, trait_axis_region)%>% 
  dplyr::select(-region.id) %>% 
  as.matrix()


sampled_dist_matrix_pwt <- dist(cwm, method = "euclidean")

output_adonis_region <- adonis2(sampled_dist_matrix_pwt ~ BC2 + BC4 + BC17 + BC15,
                             data = all_environment_region, by = 'margin', permutations = 2)

region_adonis <- tidy(output_adonis_region)

f_table <- data.frame()

for(p in 1:999){
  
  p_abun_mat <- permut.row.matrix(abundance_region)$permut.matrix
  
  p_cwm <- Fuzzit4a(p_abun_mat, trait_axis_region) %>% 
    dplyr::select(-region.id) 
  
  p_sampled_dist_matrix_pwt <- dist(p_cwm, method = "euclidean")
  
  p_output_adonis <- adonis2(p_sampled_dist_matrix_pwt ~ BC2 + BC4 + BC17 + BC15,
                             data = all_environment_region, by = 'margin', permutations = 2)
  
  Fs = p_output_adonis$F
  
  f_table <- rbind(f_table, data.frame(Fs, permutation = p))
  
}


original_regional_adonis <- data.frame(term = region_adonis$term, F_orig = output_adonis_region$F)

spread_f_table <- f_table %>% 
  mutate(term = rep(region_adonis$term, 999)) %>% 
  spread(key = term, value = Fs) %>% 
  select(region_adonis$term)

## greater than

ps <- sapply(1:4, FUN = function(x) ((sum(output_adonis_region$F[x] > spread_f_table[,x]) +1)/1000))

all_output_region <- data.frame(term = region_adonis$term[1:length(ps)], 
                               SS = output_adonis_region$SumOfSqs[1:length(ps)], P = ps) %>% 
  mutate(P= ifelse(P> 0.5, 1-P, P)) %>% 
  mutate(P = 2*P)

all_output_region %>%View()

