### This script runs the permanova, only collecting the sums of squares.
## the data has already been sub-sampled and it is loaded here as CWM and local variables
## the raw data is not presented as per agreement with the data contributors

#Local scale analysis 

library(vegan)
library(broom)

## Load cwm and local variables
# the data has already been sub sampled

local_cwm_response <- readRDS("Data/local_cwm_response.rds")
local_var_explanatory <- readRDS("Data/local_var_explanatory.rds")

trait_permute_adonis <- function(x){
  
  #create distance matrix
  sampled_dist_matrix_pwt <- dist(local_cwm_response[[x]], method = "euclidean")
  
  #run adonis
  output_adonis <- adonis2(sampled_dist_matrix_pwt ~  final_canopy + log(actual_water) + log(total_detritus),
                           data = local_var_explanatory[[x]], by = 'margin', strata = sampled_bromeliads_environment$visit_id, permutations = 2)
  
  #save only sums of squares
  all_output <- tidy(output_adonis)[,c("term", "SumOfSqs")] 
  
  return(list(all_output))
}

all_output_list <- lapply(1:1000, FUN = trait_permute_adonis)

# save results 
saveRDS(all_output_list, file = "Results/local.rds")
