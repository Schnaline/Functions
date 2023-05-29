fct_lmms <- function(ps, time_col = "Day_numeric", group = "Sample_type_per_Donor", n_tax = 15) {
  
  library(timeOmics)
  library(lmms)
  library(dplyr)
  
  
  #if (group != "none"){
  names <- sample_data(ps)[, group] %>% as.matrix() %>% as.vector %>% unique()
  ps_list <- list()
  otu_list <- list()
  
  for (i in names){
    prune <- sample_data(ps)[, group] == i
    colnames(prune) <- "a"
    
    prune_samples(prune[, "a"], ps) -> ps_list[[i]]
    
    subset_taxa(ps_list[[i]], rownames(tax_table(ps_list[[i]])) %in% microbiome::top_taxa(ps_list[[i]], n = n_tax)) -> ps_list[[i]]
    
    ps_list[[i]] %>%
      otu_table() %>% as.data.frame() -> otu_list[[i]] }
  
  # Define a function to rename the row names of a data frame
  rename_row_names <- function(df, name) {
    rownames(df) <- paste(name, rownames(df), sep = "_")
    return(df)}
  
  # Use lapply and a for loop to rename the row names of each data frame in the list
  otu_list <- lapply(names(otu_list), function(name) {
    rename_row_names(otu_list[[name]], name)})
  
  names(otu_list) <- names
  
  list_lmms <- list()
  
  for (i in names){
    time <- as.numeric(as.data.frame(as.matrix(sample_data(ps_list[[i]])))[, time_col])
    otu <- as.data.frame(t(otu_list[[i]]))
    
    list_lmms[[i]] <- lmms::lmmSpline(data = otu, time = time,
                                      sampleID = rownames(otu), deri = FALSE,
                                      basis = "p-spline", numCores = 4, timePredict = 1:length(unique(time)),
                                      keepModels = TRUE)
  }
  
  df_fitted_values <- data.frame(matrix(ncol = max(sample_data(ps)[,time_col])))
  
  list_fitted_values <- lapply(list_lmms, function(x) x@predSpline)
  df_fitted_values <- bind_rows(list_fitted_values)
  df_fitted_values$ID <- rownames(df_fitted_values)
  
  list_models_used <- lapply(list_lmms, function(x) x@modelsUsed) 
  vec_models_used <- unlist(list_models_used)
  df_models_used <- data.frame(vec_models_used, row.names = rownames(df_fitted_values))
  
  
  list_models <- unlist(lapply(list_lmms, function(x) x@models), recursive=FALSE)
  names(list_models) <- rownames(df_fitted_values)
  
  index0_ASV <- df_models_used %>% filter(df_models_used == 0) %>% rownames()
  index0_df <- data.frame(matrix(nrow = length(index0_ASV), ncol = 4))
  rownames(index0_df) <- index0_ASV
  colnames(index0_df) <- c("Estimate", "Std.Error", "t.val", "p.val")
  
  for (i in index0_ASV){
    index0_df[i, ] <- summary(list_models[[i]])$coefficients[2,]
  }
  
  index1_ASV <- df_models_used %>% filter(df_models_used >=1) %>% rownames()
  index1_df <- data.frame(matrix(nrow = length(index1_ASV), ncol = 5))
  rownames(index1_df) <- index1_ASV
  colnames(index1_df) <- c("Estimate", "Std.Error", "DF", "t.val", "p.val")
  
  for (i in index1_ASV){
    index1_df[i, ] <- summary(list_models[[i]])$tTable[2,]
  }
  
  # } else {
  #   
  #   time <- as.data.frame(as.matrix(sample_data(ps)))[, time_col]
  #   meta <- as.data.frame(as.matrix(sample_data(ps)[, c("sample_name", group, time_col)]))
  #   
  #   otu_ex <- as.data.frame(t(otu_table(ps)))
  #   
  #   lmms.output <- lmms::lmmSpline(data = otu, time = time,
  #                                  sampleID = rownames(otu), deri = FALSE,
  #                                  basis = "p-spline", numCores = 4, timePredict = 1:length(unique(time)),
  #                                  keepModels = TRUE)
  #   
  #   df_fitted_values <- as.data.frame(slot(lmms.output, 'predSpline'))
  #   
  #   
  #   df_models_used <- lmms.output %>% slot("modelsUsed") %>% as.data.frame()
  #   rownames(df_models_used) <- rownames(df_fitted_values)
  #   
  #   list_models <- lmms.output %>% slot("models")
  #   names(list_models) <- rownames(df_fitted_values)
  #   
  #   index0_ASV <- df_models_used %>% filter(df_models_used == 0) %>% rownames()
  #   index0_df <- data.frame(matrix(nrow = length(index0_ASV), ncol = 4))
  #   rownames(index0_df) <- index0_ASV
  #   colnames(index0_df) <- c("Estimate", "Std.Error", "t.val", "p.val")
  #   
  #   for (i in index0_ASV){
  #     index0_df[i, ] <- summary(list_models[[i]])$coefficients[2,]
  #   }
  #   
  #   index1_ASV <- df_models_used %>% filter(df_models_used >=1) %>% rownames()
  #   index1_df <- data.frame(matrix(nrow = length(index1_ASV), ncol = 5))
  #   rownames(index1_df) <- index1_ASV
  #   colnames(index1_df) <- c("Estimate", "Std.Error", "DF", "t.val", "p.val")
  #   
  #   for (i in index1_ASV){
  #     index1_df[i, ] <- summary(list_models[[i]])$tTable[2,]
  #   }
  # }
  
  
  res_df <- bind_rows(index0_df, index1_df)
  res_df$ASV <- rownames(res_df)
  res_df$model_type <- df_models_used
  res_df <- res_df %>% arrange(ASV)
  
  results_lmms <- list(fit = df_fitted_values, models = res_df)
  
}