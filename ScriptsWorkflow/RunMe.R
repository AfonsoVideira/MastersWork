######################## Directory Set up ##############################
#directory_name = readline(prompt = "Please enter directory name for file output: ")
#pause()
#dir.create(directory_name)

######################### EXPERIMENT SET UP #############################

file_direct <- file.path("proteinGroups.txt")

full_data <- read.delim(file_direct, stringsAsFactors = FALSE)
View(full_data)
pause()

conditionS <- readline(prompt = "Enter the names of the conditions (separated by , ) : ")
pause()
conditionS <- strsplit(conditionS,",")[[1]]

replicateS <- readline(prompt = "Enter the number of replicates per condition (separated by , ) : ")
pause()
replicateS <- lapply(strsplit(replicateS,","),as.numeric)[[1]]

conditionS_2 <- c()
replicateS_2 <- c()
for (i in 1:length(conditionS)){
  conditionS_2 <- c(conditionS_2, rep(conditionS[i],replicateS[i]))
  replicateS_2 <- c(replicateS_2,1:replicateS[i])
}

labeL <- readline(prompt = "Input the labels as they were set in MaxQuant: ")
pause()
labeL <- strsplit(labeL,",")[[1]]
#labeL <- paste0(conditionS_2,"",replicateS_2)
replicateS_2 <- as.character(replicateS_2)

experimental_design <- data.frame(label = labeL,
                                  condition = conditionS_2,
                                  replicate = replicateS_2,
                                  stringsAsFactors = FALSE)


Gene_names__Protein_ID <- c(grep("gene.name",colnames(full_data),ignore.case = TRUE,value = TRUE),
                            grep("fasta",colnames(full_data),ignore.case = TRUE,value = TRUE),
                            grep("protein.id",colnames(full_data),ignore.case = TRUE,value = TRUE))

idSite_Reverse_Contaminant <- c(grep("only.identified",colnames(full_data),ignore.case = TRUE,value = TRUE),
                                grep("reverse",colnames(full_data),ignore.case = TRUE,value = TRUE),
                                grep("contaminant",colnames(full_data),ignore.case = TRUE,value = TRUE))
q_score_PEP <- c(grep("Q.value",colnames(full_data),ignore.case = TRUE,value = TRUE),
                 grep("PEP",colnames(full_data),value = TRUE,fixed = TRUE),
                 grep("score",colnames(full_data),ignore.case = TRUE,value = TRUE))

LFQ_names <- grep("LFQ",colnames(full_data),ignore.case = TRUE,value = TRUE)

interesting_features <- c(Gene_names__Protein_ID,q_score_PEP,idSite_Reverse_Contaminant,LFQ_names)

data <- full_data[,interesting_features]

######################### CLEANING UNWANTED DATA #############################

cat("Potencialy unwanted proteins","\n")
for (i in idSite_Reverse_Contaminant){
  cat("Viewing table with proteins :",i,"\n")
  View(full_data[full_data[,i] =="+",])
  pause()
}

cat("Filter out [Only identified by site, Reverse, Contaminants] ?\nInput 1,2,3 respectivly (e.g Reverse and Contaminants should be 2,3)")
clean_data <- readline(prompt = "Input filters : ")
pause()
clean_data <- lapply(strsplit(clean_data,",")[[1]], as.numeric)
  
cat("Number of proteins before filtering: ",nrow(data),"\n")
for (i in clean_data){
  cat("Cleaning :",idSite_Reverse_Contaminant[i],"\n")
  data <- data[data[,idSite_Reverse_Contaminant[i]] !="+",]
}
cat("Number of proteins after filtering: ",nrow(data),"\n")


######################################################################

data_unique <- make_unique(data, Gene_names__Protein_ID[1], 
                           Gene_names__Protein_ID[2], delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

data_lfq <- data[,LFQ_columns]
data_lfq[data_lfq == 0] <- NA

######################### LOG2 TRANSFORM #############################
cat("Viewing histogram of each sample","\n")

display = ceiling(sqrt(sum(replicateS)))
par(mfrow = c(display,display))

for (i in (colnames(data_lfq))){
  hist(data_lfq[,i], xlab = i,main = paste("Histogram of ",i))
}
pause()

cat("Viewing histograms of the Log2 transformed data","\n")
par(mfrow = c(display,display))
for(i in (colnames(assay(data_se)))){
  hist(assay(data_se)[,i],xlab = i, main = paste("Histogram of log2",i))
}
pause()

log_fit <- readline(prompt = "Keep the Log2 transformed data ? [y/n] : ")
pause()

if (log_fit == "n") {
  assay(data_se) <- 2**(assay(data_se))
}

View(assay(data_se))
par(mfrow = c(1,1))

######################## MISSING VALUES ################################
cat("Showing the overlap protein identifications","\n")
print(plot_frequency(data_se))
pause()

filter_miss <- as.numeric(readline(prompt = "How much missing values are allowed in at least one condition ? "))
pause()
data_filt <- filter_missval(data_se, thr = filter_miss)
View(assay(data_filt))
cat("Viewing coverage of proteins per sample","\n")
print(plot_numbers(data_filt))
pause()

cat("Imputation of missing values","\n")
#cat("Viewing correlation of missing values with each sample","\n")
#plot_missval(data_filt)
#pause()
#cat("Viewing Density and Cumulative sum plots of intensities of proteins with and without missing values","\n")
#plot_detect(data_filt)
#pause()

cat("Check the imputation documentation (right panel) to select an impuation method","\n")
print(?MSnbase::impute)
pause()

imputation <- readline(prompt = "What imputation method to use? ")
pause()
imputation <- strsplit(imputation,",")[[1]]

data_imp = imputation_method(data_filt,imputation)
View(assay(data_imp))
#print(plot_imputation(data_filt, data_imp))
pause()

####################### EQUAL VARIANCE ? #############################

cat("Checking equal variance assumption.\n")
data_ratios <- equal_variance_assumption(data_imp,conditionS_2,replicateS)
rownames(data_ratios) <- rownames(data_imp)
View(data_ratios)
boxplot(data_ratios, main = "Variance ratios between conditions",
        xlab = "Condition(s)",ylab = "Ratio value")
pause()
boxplot(data_ratios, main = "Variance ratios between conditions",
        xlab = "Condition(s)",ylab = "Ratio value",ylim = c(0,10))
variance <- as.logical(readline(prompt = "Should the variances between groups be considered equal? [TRUE/FALSE] : "))
pause()

####################### Diferential test #############################
cat("### Diferencial expression analysis ###","\n")
Diferential_test <- "T-test or Welch"
if (length(unique(conditionS_2)) > 2){
  Diferential_test <- "ANOVA"
}

significance_level <- as.numeric(readline(prompt = "Input the significance level: "))
pause()

alpha = significance_level
gamma = significance_level

data_diff = test_compare(Diferential_test, data_imp,conditionS_2,variance,alpha,gamma)
data_diff_wide = get_df_wide(data_diff)

simplified = FALSE

if (Diferential_test == "ANOVA"){
  data_diff_wide_FULL_ANOVA = data_diff_wide
  data_diff_wide = get_some_columns(data_diff_wide,conditionS_2)
  simplified = TRUE
}
View(data_diff_wide)
pause()

############################ Vulcano plots ##########################################
cat("### Vulcano plots generation ###","\n")
View(data_diff_wide[,c("name",Gene_names__Protein_ID)])
label_names = readline(prompt = "Input the name of the column which contains the protein labels : ")
pause()

label_names_table = data_diff_wide[,label_names]
fold_change_cutof = as.numeric(readline(prompt = "Input fold change cut off : "))


colnames_inter = vulcano_plots(Diferential_test,data_diff_wide,label_names_table,conditionS_2,
              alpha,fold_change_cutof,simplified)


cat("### Visualizing tables of pairwise comparisons ###","\n")
if (is.null(colnames_inter)){
  cat("Showing data frame of tests","\n")
  View(data_diff_wide[,c(label_names,grep("_vs_",colnames(data_diff_wide),value = TRUE))])
  pause()
}

if(!is.null(colnames_inter)) {
  if (!dir.exists("outTables")){
  dir.create("outTables")
  }
  table_names = c("BS_BS","CS_CS","BH_BH","BY_BY","Bonf_Bonf",
                  "S_S")
  table_out = output_tables(data_diff_wide,colnames_inter,label_names,conditionS_2)
  for (i in table_names){
    cat("Viewing table of tests by : ",i)
    View(table_out[[i]])
    write.csv(table_out[[i]],file = paste("outTables/",i,".csv",sep = ""))
    pause()
  }
}



