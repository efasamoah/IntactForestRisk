

topsis <- function(X, W, is_benefit, std = FALSE) {
  # Normalize the decision matrix
  if(std){
    X_normalized <- t(apply(X, 1, function(x) x / sqrt(sum(x^2))))
  } else { X_normalized <- X}
  
  # Weighted normalized decision matrix
  X_weighted <- X_normalized * W
  
  # Determine ideal and negative-ideal solutions
  ideal_best <- ifelse(is_benefit, apply(X_weighted, 2, max), apply(X_weighted, 2, min))
  ideal_worst <- ifelse(is_benefit, apply(X_weighted, 2, min), apply(X_weighted, 2, max))
  
  # Calculate distances from ideal and negative-ideal solutions
  d_best <- sqrt(rowSums((sweep(X_weighted, 2, ideal_best, `-`))^2))
  d_worst <- sqrt(rowSums((sweep(X_weighted, 2, ideal_worst, `-`))^2))
  
  # Calculate the relative closeness to the ideal solution
  closeness <- d_worst / (d_best + d_worst)
  
  # Rank the alternatives based on their closeness to the ideal solution
  # ranks <- order(closeness, decreasing = TRUE)/length(closeness)
  
  # Return the ranked alternatives
  return(closeness)
}



