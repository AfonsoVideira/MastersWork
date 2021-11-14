
############################## ANOVA ##########################################

anova_contrast = function(data_diff,intensities_id, group, variance){
  contrastS = combinations(unique(group),2)
  for (i in 1:nrow(contrastS)){
    c1 = contrastS[i,1]
    c2 = contrastS[i,2]
    cat("Contrasts:\n",c1,"_vs_",c2,"\n")
    diff_c1_c2 = c()
    p_value_list = c()
    for (j in 1:nrow(intensities_id)){
      intensities = as.numeric(intensities_id[j,])
      data_for_contrast = data.frame(intensities = intensities , group = group)
      t_test = t.test(x = data_for_contrast[data_for_contrast$group == c1,1],
                      y = data_for_contrast[data_for_contrast$group == c2,1],var.equal = variance)
      diff_c1_c2 = c(diff_c1_c2,as.numeric(t_test$estimate[1] - t_test$estimate[2]))
      p_value_list = c(p_value_list,as.numeric(t_test$p.value))
    }
    vs1 = paste(c1,"_vs_",c2,"_diff",sep = "")
    vs2 = paste(c1,"_vs_",c2,"_p.val",sep = "")
    rowData(data_diff)[,vs1] <- diff_c1_c2
    rowData(data_diff)[,vs2] <- p_value_list
  }
  return(as.data.frame(rowData(data_diff)[,(length(rowData(data_diff))-(nrow(contrastS)*2 - 1)):length(rowData(data_diff))]))
}

anova_func = function(data_diff,intensities_id, group, variance, alpha, gamma){
  anova_p_values = c()
  for (i in 1:nrow(intensities_id)){
    intensities = as.numeric(intensities_id[i,])
    data_for_diftest = data.frame(intensities = intensities, group = group)
    anovA = oneway.test(intensities ~ group,data_for_diftest,var.equal = variance)
    p_value = as.numeric(anovA$p.value)
    anova_p_values = c(anova_p_values,p_value)
  }
  rowData(data_diff)$ANOVA_p.val = anova_p_values
  #ajustar ANOVA
  vs = "ANOVA_p.adj"
  df1 = anova_contrast(data_diff,intensities_id,group,variance)
  data_diff = adjustment_function(vs,data_diff,anova_p_values,alpha,gamma,df1)
  return(data_diff)
}

####################### t-test ####################################################

t_test = function(data_diff,intensities_id, group, variance, alpha, gamma){
  contrastS = combinations(unique(group),2)
  for (i in 1:nrow(contrastS)){
    c1 = contrastS[i,1]
    c2 = contrastS[i,2]
    cat(c1,"_vs_",c2,"\n")
    diff_c1_c2 = c()
    p_value_list = c()
    for (j in 1:nrow(intensities_id)){
      intensities = as.numeric(intensities_id[j,])
      data_for_contrast = data.frame(intensities = intensities , group = group)
      t_test = t.test(x = data_for_contrast[data_for_contrast$group == c1,1],
                      y = data_for_contrast[data_for_contrast$group == c2,1],var.equal = variance)
      diff_c1_c2 = c(diff_c1_c2,as.numeric(t_test$estimate[1] - t_test$estimate[2]))
      p_value_list = c(p_value_list,as.numeric(t_test$p.value))
    }
    vs1 = paste(c1,"_vs_",c2,"_diff",sep = "")
    vs2 = paste(c1,"_vs_",c2,"_p.val",sep = "")
    rowData(data_diff)[,vs1] <- diff_c1_c2
    rowData(data_diff)[,vs2] <- p_value_list
    #Ajustar p_values
    vs = paste(c1,"_vs_",c2,"_p.adj",sep = "")
    data_diff = adjustment_function(vs,data_diff,p_value_list,alpha,gamma,NULL)
  }
  return(data_diff)
}






#an_t = pairwise.t.test(data_for_diftest$intensities,data_for_diftest$group,p.adjust.method = "none",paired = FALSE)

#oneway.test() 
#v1 E V2 , outra coluna quociente, V1/V2 <- boxplot
 
#?pairwise.t.test()

