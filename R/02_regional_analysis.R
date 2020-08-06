library(vegan)
library(broom)

## Load cwm and local variables
# the data has already been sub sampled

regional_cwm_response <- readRDS("Data/regional_cwm_response.rds")
regional_var_explanatory <- readRDS("Data/regional_var_explanatory.rds")

trait_permute_adonis <- function(x){
  
  #create distance matrix
  sampled_dist_matrix_pwt <- dist(regional_cwm_response[[x]], method = "euclidean")
  
  #create a new variable for the regions to add caribean central america and south america 
  
  output_adonis <- adonis2(sampled_dist_matrix_pwt ~  
                             BC2 + BC4 + BC17+ BC15,  
                           data = regional_var_explanatory[[x]], by = 'margin', permutations = 2)
  
  all_output <- tidy(output_adonis)[,c("term", "SumOfSqs")]
  
  return(list(all_output))
}

all_output_list <- lapply(1:1000, FUN = trait_permute_adonis)

saveRDS(all_output_list, file = "Results/regional.rds")
