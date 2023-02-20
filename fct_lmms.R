fct_lmms <- function(ps, time_col = "Day_numeric"){
  
  library(timeOmics)
  library(lmms)
  
  time <- as.data.frame(as.matrix(sample_data(ps)))[, time_col]
  
  otu <- as.data.frame(t(otu_table(ps)))
  meta <- as.data.frame(as.matrix(sample_data(ps)[, c("sample_name", time_col)]))
  
  
  lmms.output <- lmms::lmmSpline(data = otu, time = time,
                                 sampleID = rownames(otu), deri = FALSE,
                                 basis = "p-spline", numCores = 4, timePredict = 1:10,
                                 keepModels = TRUE)
  
  df_fitted_values <- as.data.frame(slot(lmms.output, 'predSpline'))
  
  #col_df <- colnames(summary(lmms.output@models[[2]])$coefficients)
  
  
  models_used <- lmms.output %>% slot("modelsUsed") %>% as.data.frame()
  rownames(models_used) <- rownames(df_fitted_values)
  
  list_models <- lmms.output %>% slot("models")
  names(list_models) <- rownames(df_fitted_values)
  
  index0_ASV <- models_used %>% filter(models_used == 0) %>% rownames()
  index0_df <- data.frame(matrix(nrow = length(index0_ASV), ncol = 4))
  rownames(index0_df) <- index0_ASV
  colnames(index0_df) <- c("Estimate", "Std.Error", "t.val", "p.val")
  
  for (i in index0_ASV){
    index0_df[i, ] <- summary(list_models[[i]])$coefficients[2,]
  }
  
  index1_ASV <- models_used %>% filter(models_used >=1) %>% rownames()
  index1_df <- data.frame(matrix(nrow = length(index1_ASV), ncol = 5))
  rownames(index1_df) <- index1_ASV
  colnames(index1_df) <- c("Estimate", "Std.Error", "DF", "t.val", "p.val")
  
  for (i in index1_ASV){
    index1_df[i, ] <- summary(list_models[[i]])$tTable[2,]
  }
  
  
  
  res_df <- bind_rows(index0_df, index1_df)
  res_df$ASV <- rownames(res_df)
  res_df$model_type <- models_used
  res_df <- res_df %>% arrange(ASV)
  
  results_lmms <- list(fit = df_fitted_values, models = res_df)
  
  return(results_lmms)
  
}