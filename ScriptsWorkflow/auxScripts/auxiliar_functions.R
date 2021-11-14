
pause = function()
{
  if (interactive())
  {
    invisible(readline(prompt = "Press <Enter> to continue..."))
  }
  else
  {
    cat("Press <Enter> to continue...")
    invisible(readLines(file("stdin"), 1))
  }
}


imputation_method = function(data_wo_im,parameters){
  
  if (parameters[1] == "knn") {
    
    parameters_A = c(10, 0.5, 0.8, 1500, 362436069)
    
    for (i in 1:length(parameters[2:length(parameters)])){
      if ((!is.na(parameters[i+1])) & (parameters[i+1] != "")) {
        parameters_A[i] = as.numeric(parameters[i+1])
      }
    }
    print(parameters_A)
    data_imp = impute(data_wo_im,fun = parameters[1],k = parameters_A[1],
                      rowmax = parameters_A[2], colmax = parameters_A[3],
                      maxp = parameters_A[4], rng.seed = parameters_A[5])
  }
  
  if (parameters[1] == "MinDet") {
    q = 0.01
    if ((!is.na(parameters[2])) & (parameters[2] != "")){
      q = as.numeric(parameters[2])
    }
    print(q)
    data_imp = impute(data_wo_im, fun = parameters[1], q = q)
  }
  
  if (parameters[1] == "MinProb"){
    parameters_A = c(0.01, 1)
    for (i in 1:length(parameters[2:length(parameters)])){
      if ((!is.na(parameters[i+1])) & (parameters[i+1] != "")){
        parameters_A[i] = as.numeric(parameters[i+1])
      }
    }
    print(parameters_A)
    data_imp = impute(data_wo_im, fun = parameters[1], q = parameters_A[1],
                      tune.sigma = parameters_A[2])
  }
  
  if (parameters[1] == "man"){
    parameters_A = c(1,0.3)
    for (i in 1:length(parameters[2:length(parameters)])){
      if ((!is.na(parameters[i+1])) & (parameters[i+1] != "")){
        parameters_A[i] = as.numeric(parameters[i+1])
      }
    }
    print(parameters_A)
    data_imp = impute(data_wo_im,fun = parameters[1],
                      shift = parameters_A[1], scale = parameters_A[2])
  }
  
  else{
    data_imp = impute(data_wo_im,fun = parameters[1])
  }
  return (data_imp)
  }
  

condition_var_matrix = function(data_imp,group,replicates){
  intensities_id = get_df_wide(data_imp)[,c(colnames(data_imp))]
  conditions = unique(group)
  var_data_frame = data.frame()
  var_data_frame[1,conditions] = 1:length(conditions)
  for (j in 1:nrow(intensities_id)){
    row_intensities = intensities_id[j,]
    start = 0
    var_row = c()
    for (i in replicates){
      variance = var(as.numeric(row_intensities[(start+1):(start+i)]))
      var_row = c(var_row,variance)
      start = start + i
    }
    var_data_frame[j,] = var_row
  }
  return(var_data_frame)
}



equal_variance_assumption = function(data_imp,group,replicates){
  if (any(replicates < 4)){
    cat("WARNING : low replicate number in at least one condition, equal variances assumption is FALSE.\n")
    #cat("Recomendation : Use Welch method for comparisons.\n")
  }
  variance_matrix = condition_var_matrix(data_imp,group,replicates)
  combinat = combinations(unique(group),2)
  data_ratios = data.frame()
  ratio_names = c()
  for (i in 1:nrow(combinat)){
    ratio_names = c(ratio_names,paste(combinat[i,1],combinat[i,2],sep = "/"))
  }
  data_ratios[1:nrow(variance_matrix),ratio_names] = 1:length(ratio_names)
  for (i in 1:nrow(combinat)){
    data_ratios[,i] = variance_matrix[,combinat[i,1]] / variance_matrix[,combinat[i,2]]
  }
  return(data_ratios)
  }

#contrast_generator

test_compare = function(user_input,data_imp,group, variance, alpha, gamma){
  data_diff = data_imp
  intensities_id = get_df_wide(data_diff)[,c(colnames(data_diff))]
  if(any(user_input == "ANOVA")){
    data_diff = anova_func(data_diff,intensities_id , group , variance, alpha, gamma)
    
  }
  else{
  data_diff = t_test(data_diff,intensities_id,group,variance,alpha,gamma)
  }
  return(data_diff)
  
}
  
  
# rowData(se)[,"nome"] <- vetor
