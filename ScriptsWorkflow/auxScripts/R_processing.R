#library(proDA)
library(SummarizedExperiment)
library(DEP)
library(dplyr)
source("input.R")
source("tecnicas_ajustamento.R")



data_se <- make_se(data_unique, LFQ_columns, experimental_design)

plot_frequency(data_se)

data_filt <- filter_missval(data_se, thr = 0)
data_filt2 <- filter_missval(data_se, thr = 1)


plot_numbers(data_filt)
plot_coverage(data_filt)

data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)


plot_missval(data_filt)
plot_detect(data_filt)

############################### Missing Values imputation ##############################

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

data_imp_man <- impute(data_norm, fun = "man", shift = 1, scale = 0.3)

data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

plot_imputation(data_norm, data_imp_knn)



################### Differential enrichment analysis #################################

data_diff_all_contrasts <- test_diff(data_imp, type = "all")
p_value_list <- rowData(data_diff_all_contrasts)
p_value_list <- p_value_list$C_vs_F_p.val

n_reject_sgof <- Numero_testes_rejeitados_Sgof(p_value_list,0.05,0.05)
a <- sort(rowData(data_diff_all_contrasts)$C_vs_F_p.adj)

a <- (a[n_reject_sgof] + a[n_reject_sgof+1])/2

dep <- add_rejections(data_diff_all_contrasts, alpha = a, lfc = 0)
dep2 <- add_rejections(data_diff_all_contrasts, alpha = 0.03, lfc = 0)
################ Result Visualization #########################################

plot_pca(dep, x = 1, y = 2, n = 3, point_size = 4)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Blues")


plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

plot_volcano(dep, contrast = "C_vs_F", label_size = 3, add_names = TRUE)

#plot_single(dep, proteins = c("Ckb", "Tpi"))
#plot_cond(dep,plot = FALSE)

data_results <- get_results(dep)

data_results %>% filter(significant) %>% nrow()

colnames(data_results)

df_wide <- get_df_wide(dep)


####################################################################################################
#df_long <- get_df_long(dep)
#data_results_1 <- LFQ(data,experimental_design,fun ="MinProb",type = "all",alpha = 0.2,lfc = 0)
#data_results_1$results
#colnames(data)[14]
