

vulcano_plots = function(Diferential_test,data_diff_wide,nameS,
                         conditionS_2,alpha,fold_change_cutof,simplified)
  {
  
if (Diferential_test != "ANOVA"){
  diff_columns <- grep("diff",colnames(data_diff_wide), value = FALSE)[1]
  vulcano_df = data_diff_wide[,diff_columns:length(colnames(data_diff_wide))]
  for (i in colnames(vulcano_df)){
    print(i)
    if (grepl(".diff",i))
    {
      contrast = i
    }
    else
    {
      d_f = data.frame(contrast = as.numeric(vulcano_df[,contrast]) , p_value = as.numeric(vulcano_df[,i]))
      rownames(d_f) = nameS
      cat(contrast,i,":\n")
      print(EnhancedVolcano(toptable = d_f, rownames(d_f),x = "contrast",y = "p_value",
                            title = paste(i),pCutoff = alpha , 
                           FCcutoff = fold_change_cutof))
      
      cat(rownames(d_f[d_f[,2]<= alpha,]))
      pause()
    }
  }
  return(NULL)
}
if (Diferential_test == "ANOVA"){
  n1 = nrow(combinations(unique(conditionS_2),2))
  Q=grep("_diff",colnames(data_diff_wide),value = TRUE)
  Q1 = c()
  index = 1:n1
  for (i in 1:6){
    Q1 = c(Q1,rep(Q[index],6))
    index = index + n1
  }
  if(simplified == TRUE){
    Q1 = Q
    }
  Q2 = grep("p.val.adj_",colnames(data_diff_wide),value = TRUE)
  for (j in 1:length(Q1)){
    d_f = data_diff_wide[,c(Q1[j],Q2[j])]
    rownames(d_f) = nameS
    d_f = d_f[!is.na(d_f[,1]),]
    cat(Q1[j],Q2[j],":\n")
    print(EnhancedVolcano(toptable = d_f, rownames(d_f), x = Q1[j], y = Q2[j],
                          title = paste(Q2[j]),pCutoff = alpha,
                          FCcutoff = fold_change_cutof))
    
    cat(rownames(d_f[d_f[,2]<= alpha,]))
    pause()
    enD = readline(prompt = "End loop ? [y/n] ")
    if (enD == "y"){
      break
    }
  }
  return(data.frame(Q1,Q2))
}
}
