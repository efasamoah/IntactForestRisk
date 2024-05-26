
weightedPCA <- function(x, focal, scale, center, npc) {
  x = na.omit(x[focal])
  pca.out = prcomp(x, scale = scale, center = center)
  biplot(pca.out)
  
  dat_ = as.data.frame(t(pca.out$rotation))
  dat_ = dat_[1:npc,] #subset to the required PCs
  eigen = (pca.out$sdev^2)[1:npc]
  
  n_var = length(focal)
  weights_values <- lapply(1:n_var, FUN = function(x) { colSums(abs(dat_[x])*eigen[x]) })
  std_weights <- as_tibble(unlist(weights_values)); std_weights <- std_weights/sum(std_weights, na.rm = T)
  return(std_weights)
}
