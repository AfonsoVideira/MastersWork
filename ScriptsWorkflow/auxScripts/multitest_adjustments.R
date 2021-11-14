
anova_adjusment_function = function(data_diff,alpha,gamma,vector_new_names){
  interest_cols = grep("p.val",vector_new_names,value = TRUE)
  interest = rowData(data_diff)[,interest_cols]
  interest = as.data.frame(interest)
  full_list_p_value = c()
  for (i in colnames(interest)){
    p_value_vector = interest[,i]
    p_value_vector = p_value_vector[!is.na(p_value_vector)]
    full_list_p_value = c(full_list_p_value,p_value_vector)
  }
  names_df = c("Binom_SGoF","Conser_SGoF","BH","BY","Bonferroni","Sidak")
  df = anova_adjusment_function_aux(full_list_p_value,alpha,gamma)
  for (i in 1:ncol(df)){
    adj_p_all = df[,i]
    index = 1
    for (j in 1:ncol(interest)){
      new_name = paste(interest_cols[j],names_df[i],sep = ".adj_")
      info = interest[,j]
      for (k in 1:length(info)){
        if (!is.na(info[k])){
          info[k] = adj_p_all[index]
          index = index + 1
        }
      }
      rowData(data_diff)[,new_name] <- info
    }
  }
  return(data_diff)
}


anova_adjusment_function_aux = function(p_value_list,alpha,gamma){
  ord = order(order(p_value_list))
  Binom_SGoF = sgof::Binomial.SGoF(p_value_list,alpha,gamma)$Adjusted.pvalues[ord]
  Conser_SGoF = sgof::SGoF(p_value_list,alpha, gamma)$Adjusted.pvalues[ord]
  BH = sgof::BH(p_value_list,alpha)$Adjusted.pvalues[ord]
  BY = sgof::BY(p_value_list,alpha)$Adjusted.pvalues[ord]
  bonf = p.adjust(p_value_list,"bonferroni")
  sidak = 1 - (1 - p_value_list)**length(p_value_list)
  df = data.frame(Binom_SGoF,Conser_SGoF,BH,BY,bonf,sidak)
  return(df)
}

adjustment_function = function(vs,data_diff, p_value_list, alpha, gamma, df1){
  ord = order(order(p_value_list))
  Binom_SGoF = sgof::Binomial.SGoF(p_value_list,alpha,gamma)$Adjusted.pvalues[ord]
  Conser_SGoF = sgof::SGoF(p_value_list,alpha, gamma)$Adjusted.pvalues[ord]
  BH = sgof::BH(p_value_list,alpha)$Adjusted.pvalues[ord]
  BY = sgof::BY(p_value_list,alpha)$Adjusted.pvalues[ord]
  bonf = p.adjust(p_value_list,"bonferroni")
  sidak = 1 - (1 - p_value_list)**length(p_value_list)
  names = c("Binom_SGoF","Conser_SGoF","BH","BY","Bonferroni","Sidak")
  df = data.frame(Binom_SGoF,Conser_SGoF,BH,BY,bonf,sidak)
  if (vs == "ANOVA_p.adj"){
    for (i in 1:length(names)){
      col_name = paste(vs,names[i],sep = "_")
      rowData(data_diff)[,col_name] <- df[,i]
      vector_new_names = c()
      #format each contrast given adjusted p-values
      for (j in colnames(df1)){
        dup = df1[,j]
        dup[df[,i] > alpha] <- NA
        new_name = paste(col_name,j,sep="__")
        vector_new_names = c(vector_new_names,new_name)
        rowData(data_diff)[,new_name] <- dup
      }
      data_diff = anova_adjusment_function(data_diff,alpha,gamma,vector_new_names)
    }
  }
  
  else{
    for (i in 1:length(names)){
      col_name = paste(vs,names[i],sep = "_")
      rowData(data_diff)[,col_name] <- df[,i]
      
    }
  }
  return(data_diff)
}

