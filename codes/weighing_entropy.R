#Entropy weighing

# Function to calculate entropy
calc_entropy <- function(x){

  m <- length(!is.na(x))
  probabilities = x/sum(x, na.rm = T)

  entropy = -1*(sum(probabilities*log(probabilities), na.rm = T)/log(m))
  entropy_weight = 1-entropy
  return(entropy_weight)
}


entropy_weighting <- function(data, target_variable){
  namelist <- colnames(testing)
  
  weights_values = c()
    for(k in 1:length(target_variable)){
      weights_values[[k]] = calc_entropy(
        testing[,namelist[grep(target_variable[[k]],namelist)]]
      )}
  
  weights_values <- as_tibble(unlist(weights_values))
  (weights_values <- weights_values/sum(weights_values))
  return(weights_values)
}