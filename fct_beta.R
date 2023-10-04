fct_beta <- function(phyloseq = ps, transformation = FALSE, distance = "bray"){
  
  #source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R") 
  
  if (transformation != FALSE) {
    phyloseq %>%
      microbiome::transform(transformation) -> ps_transformed
  } else {
    ps_transformed <- phyloseq
  }
  
  ps_transformed %>%
    phyloseq::distance(method = distance) -> dis 
  
  ps_transformed %>%
    ordinate(method = "PCoA",
             distance = dis) -> ord
  
  return(list("distance_matrix" = dis, "ordination" = ord))
  
}