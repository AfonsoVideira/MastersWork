list.of.packages <- c("SummarizedExperiment", "DEP","arrangements","sgof","EnhancedVolcano")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,dependencies = T)

library(SummarizedExperiment)
library(DEP)
library(arrangements)
library(sgof)
library(EnhancedVolcano)
source("auxScripts/auxiliar_functions.R")
source("auxScripts/anova_t_test.R")
source("auxScripts/multitest_adjustments.R")
source("auxScripts/get_anova.R")
source("auxScripts/Vulcanoplot.R")
