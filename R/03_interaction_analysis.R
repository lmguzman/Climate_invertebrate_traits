library(vegan)
library(broom)

## Load cwm and local variables
# the data has already been sub sampled

interaction_cwm_response <- readRDS("Data/interaction_cwm_response.rds")
interaction_var_explanatory <- readRDS("Data/interaction_var_explanatory.rds")

trait_permute_adonis <- function(x){
  
  #create distance matrix
  sampled_dist_matrix_pwt <- dist(interaction_cwm_response[[x]], method = "euclidean")
  
  #create a new variable for the regions to add caribean central america and south america 
  
  output_adonis <- adonis2(sampled_dist_matrix_pwt ~  log(total_detritus) + final_canopy + log(actual_water) + BC2 + BC4 + BC17 + BC15 +
                                 log(total_detritus):BC2 + log(total_detritus):BC4 + log(total_detritus):BC17 + log(total_detritus):BC15 +
                                 log(actual_water):BC2 + log(actual_water):BC4 + log(actual_water):BC17 + log(actual_water):BC15,
                               data = interaction_var_explanatory[[x]], by = 'term', permutations = 2)
  
  all_output <- tidy(output_adonis)[,c("term", "SumOfSqs")]
  
  return(list(all_output))
}

all_output_list <- lapply(1:1000, FUN = trait_permute_adonis)

saveRDS(all_output_list, file = "Results/interaction.rds")
