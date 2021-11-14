get_some_columns = function(data_frame,conditionS_2){
  number_contrasts = nrow(combinations(length(unique(conditionS_2)),2))
  stationary_df = data_frame[,1:(ncol(data_frame)-(48*number_contrasts + 6))]
  to_format_df = data_frame[,(ncol(data_frame)-(48*number_contrasts + 5)):(ncol(data_frame))]
  for (increment in 0:5){
    col = 1:(2*number_contrasts + 1) + ((2*number_contrasts + 1) + (6*number_contrasts))*increment #0-5
    adj = (2*number_contrasts + 2):(2*number_contrasts + 7) + 
      ((2*number_contrasts + 1) + (6*number_contrasts))*increment + (increment*6)
    stationary_df[,colnames(to_format_df[,col])] <- to_format_df[,col]
    stationary_df[,colnames(to_format_df[,adj])] <- to_format_df[,adj]
  }
  return(stationary_df)
}


output_tables = function(table,important_cols,names,conditionS_2){
  comb = combinations(unique(conditionS_2),2)
  contrast_vector = c()
  for (i in 1:nrow(comb)){
    contrast_vector = c(contrast_vector,paste(comb[i,1],comb[i,2],sep = "_vs_"))
  }
  l = list()
  mega_list_names = c("BS_BS","CS_CS","BH_BH","BY_BY","Bonf_Bonf","S_S")
  for (i in 0:5){
    slicE_1 = important_cols[,1][(1:nrow(comb)) + (i*nrow(comb))]
    slicE = important_cols[,2][(1:nrow(comb))+(i*nrow(comb))]
    protein_names = table[!is.na(table[,slicE[1]]),names]
    d_f = data.frame(matrix(nrow = length(contrast_vector)))
    rownames(d_f) = contrast_vector
    for(j in protein_names){
      protein_j_contrasts_diff = as.numeric(table[table[,names] == j,slicE_1])
      protein_j_contrasts_p_value = as.numeric(table[table[,names] == j,slicE])
      d_f[,paste(j,"diff",sep="_")] = protein_j_contrasts_diff
      d_f[,paste(j,"p_adj",sep="_")] = protein_j_contrasts_p_value
    }
    d_f = d_f[,-1]
    l[[mega_list_names[i+1]]] = d_f
  }
  return(l)
}


