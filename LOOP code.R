#==================================Mendelian Rrandomization Analysis=========================
#Load required packages
library(TwoSampleMR)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(stringr)

#Set the file path
PATH<-'local path'
setwd(PATH)

#Extract instrumental variables for alcohol consumption
EXP<-extract_instruments(outcomes = "ieu-b-73", p1 = 5e-08,clump = T,r2=0.01,kb=10000, p2 = 5e-08)
EXP$id.exposure<-"Alcohol consumption"
EXP$exposure<-"Alcohol consumption"
head(EXP,6)

#Extract strong SNP with calculate the F value of SNP
EXP$R2<-EXP$beta.exposure*EXP$beta.exposure*2*(EXP$eaf.exposure)*(1-EXP$eaf.exposure)
EXP$Fvalue<-(EXP$samplesize.exposure-2)*EXP$R2/(1-EXP$R2)
EXP<- subset(EXP, Fvalue >= 10)
write.csv(EXP, "F value.csv")

# Read outcome GWAS IDs and names list
gwas_data <- read.csv("gwas_ids.csv", header = FALSE, stringsAsFactors = FALSE)
gwas_ids <- gwas_data$V1  # The first column contains GWAS IDs
gwas_names <- gwas_data$V2  # The second column contains GWAS names

#Check for confounding factors
#Extract the name of all SNP with confounding factors (BMI, smoking, liver disease and diabetes)
data_BMI <- fread("BMI.tsv") # Confounding factor SNP file path_BMI
data_BMI$rs_id <- sub("^(rs[0-9]+).*", "\\1", data_BMI$riskAllele)
data_Smoking <- fread("Smoking.tsv")# Confounding factor SNP file path_Smoking
data_Smoking$rs_id <- sub("^(rs[0-9]+).*", "\\1", data_Smoking$riskAllele)
data_liver <- fread("liver.tsv")# Confounding factor SNP file path_liver
data_liver$rs_id <- sub("^(rs[0-9]+).*", "\\1", data_liver$riskAllele)
data_diabetes <- fread("diabetes.tsv")# Confounding factor SNP file path_diabetes
data_diabetes$rs_id <- sub("^(rs[0-9]+).*", "\\1", data_diabetes$riskAllele)

#Extract all SNP numbers starting with rs from the rs_id column of the matched_rows_BMI data frame and store them in a new data frame 
snp_BMI <- data.frame(SNP = str_extract(data_BMI$rs_id, "rs\\d+"))
snp_Smoking <- data.frame(SNP = str_extract(data_Smoking$rs_id, "rs\\d+"))
snp_liver <- data.frame(SNP = str_extract(data_liver$rs_id, "rs\\d+"))
snp_diabetes <- data.frame(SNP = str_extract(data_diabetes$rs_id, "rs\\d+"))

# Merge four confounding SNP data frames
all_snp <- bind_rows(snp_BMI, snp_Smoking, snp_liver, snp_diabetes)

# Loop through each outcome (metabolites)
for (i in 1:length(gwas_ids)){
  # Loop through each outcome
  folder_name <- as.character(i)
  dir.create(folder_name)
  
  # Extract outcome data
  OUT <- extract_outcome_data(snps=EXP$SNP, outcomes=gwas_ids[i], proxies=T, maf_threshold = 0.01)
  OUT$id.outcome <- gwas_ids[i]
  OUT$outcome <- gwas_names[i]
  OUT <- OUT[!duplicated(OUT$SNP),]
  
  # Harmonize SNP alleles between exposure and outcome
  data_h <- harmonise_data(exposure_dat=EXP, outcome_dat=OUT, action=2)
  data_h <- data_h %>% subset(data_h$mr_keep == TRUE)
  write.csv(data_h, file=paste0(folder_name, "/SNP.csv"))
  
  # Remove confounding SNPs from data_h
  data_h_SNP <- anti_join(data_h, all_snp, by = "SNP")
  
  # Conduct Mendelian Randomization (MR) analysis
  data_h_SNP_mr <- data_h_SNP[!duplicated(data_h_SNP$SNP),]
  mr <- mr(data_h_SNP_mr)
  mr_OR <- generate_odds_ratios(mr)# Convert to odds ratios
  write.csv(mr_OR, file=paste0(folder_name, "/mr_OR.csv"))

  # Conduct heterogeneity test
  H <- mr_heterogeneity(data_h_SNP_mr)
  write.csv(H, file=paste0(folder_name, "/H.csv"))

  # Conduct pleiotropy test
  ple <- mr_pleiotropy_test(data_h_SNP_mr)
  write.csv(ple, file=paste0(folder_name, "/ple.csv"))
  
  #Use MR-PRESSO to detect outliers and horizontal pleiotropy
  library(MRPRESSO)
  set.seed(1234)#Set random seed
  results<-mr_presso(BetaOutcome = "beta.outcome",BetaExposure="beta.exposure", 
                     SdOutcome = "se.outcome", OUTLIERtest = T, DISTORTIONtest = T,SdExposure = "se.exposure",data=data_h_SNP_mr)
  #Extract main results of ME_PRESSO
  main_results <- results$`Main MR results`
  presso_results <- results$`MR-PRESSO results`
  pval_global <- presso_results$`Global Test`$Pvalue
  # Check if p-value is greater than 0.05
  if (pval_global > 0.05){
    # Extract main results of ME_PRESSO
    causal_estimate <- main_results$`Causal Estimate`[1]
    corrected_beta <- main_results$`Causal Estimate`[2]
    sd <- main_results$Sd[1]
    corrected_se <- main_results$Sd[2]
    t_stat <- main_results$`T-stat`[1]
    corrected_t_stat <- main_results$`T-stat`[2]
    pval_main <- main_results$`P-value`[1]
    corrected_pval_main <- main_results$`P-value`[2]
    rssobs <- presso_results$`Global Test`$RSSobs
    pval_global <- presso_results$`Global Test`$Pvalue
    # Combined into a matrix
    result_matrix <- matrix(c(causal_estimate, sd, t_stat, pval_main, 
                              corrected_beta, corrected_se, corrected_t_stat, corrected_pval_main, rssobs, pval_global), 
                            nrow = 1, byrow = TRUE)
    # Add column names
    colnames(result_matrix) <- c("Causal Estimate", "Sd", "T-stat", "P-value", 
                                 "Corrected beta", "Corrected Se", "Corrected T-stat", "Corrected P-value", 
                                 "RSSobs", "Global P-value")
    df <- as.data.frame(result_matrix) #Convert matrix to data frame
    df <- df %>%
      mutate(
        OR = exp(corrected_beta),  # calculate corrected_OR value
        low_CI = exp(corrected_beta - 1.96 * corrected_se),  # Calculate the lower 95% CI
        up_CI = exp(corrected_beta + 1.96 * corrected_se)    # Calculate the upper 95% CI
      )
    write.csv(df, file=paste0(folder_name, "/presso.csv"))
  } else {# Filter outliers if p-value ≤ 0.05
    outlier_test <- as.data.frame(results$`MR-PRESSO results`$`Outlier Test`)# Extract Outlier Test Results
    significant_rows <- which(outlier_test$Pvalue < 0.05)# Filter the row numbers with a P value less than 0.05
    # Check if there is a P value less than 0.05
    if (length(significant_rows) > 0){
      data_h_SNP1 <- data_h_SNP[-significant_rows, ]
    } else {
      data_h_SNP1 <- data_h_SNP
    } 
    # Re-run MR analysis
    data_h_SNP_mr <- data_h_SNP1
    mr <- mr(data_h_SNP_mr)
    mr_OR <- generate_odds_ratios(mr)
    write.csv(mr_OR, file=paste0(folder_name, "/mr_OR.csv"))
    H <- mr_heterogeneity(data_h_SNP_mr)#
    write.csv(H, file=paste0(folder_name, "/H.csv"))
    ple <- mr_pleiotropy_test(data_h_SNP_mr)
    write.csv(ple, file=paste0(folder_name, "/ple.csv"))
    
    library(MRPRESSO)
    set.seed(1234)#Set random seed
    results<-mr_presso(BetaOutcome = "beta.outcome",BetaExposure="beta.exposure", 
                       SdOutcome = "se.outcome", OUTLIERtest = T, DISTORTIONtest = T,SdExposure = "se.exposure",data=data_h_SNP_mr)
    #Extract main results of ME_PRESSO
    causal_estimate <- main_results$`Causal Estimate`[1]
    corrected_beta <- main_results$`Causal Estimate`[2]
    sd <- main_results$Sd[1]
    corrected_se <- main_results$Sd[2]
    t_stat <- main_results$`T-stat`[1]
    corrected_t_stat <- main_results$`T-stat`[2]
    pval_main <- main_results$`P-value`[1]
    corrected_pval_main <- main_results$`P-value`[2]
    rssobs <- presso_results$`Global Test`$RSSobs
    pval_global <- presso_results$`Global Test`$Pvalue
    
    result_matrix <- matrix(c(causal_estimate, sd, t_stat, pval_main, 
                              corrected_beta, corrected_se, corrected_t_stat, corrected_pval_main, rssobs, pval_global), 
                            nrow = 1, byrow = TRUE)
    
    colnames(result_matrix) <- c("Causal Estimate", "Sd", "T-stat", "P-value", 
                                 "Corrected beta", "Corrected Se", "Corrected T-stat", "Corrected P-value", 
                                 "RSSobs", "Global P-value")
    df <- as.data.frame(result_matrix)
    df <- df %>%
      mutate(
        OR = exp(corrected_beta),  # calculate corrected_OR value
        low_CI = exp(corrected_beta - 1.96 * corrected_se),  # Calculate the lower 95% CI
        up_CI = exp(corrected_beta + 1.96 * corrected_se)    # Calculate the upper 95% CI
      )
    write.csv(df, file=paste0(folder_name, "/presso.csv"))
  }
  #Visualize Leave-one-out analysis
  mr_outcome_loo <- mr_leaveoneout(data_h_SNP_mr)
  p3 <- mr_leaveoneout_plot(mr_outcome_loo) # Leave-one-out plot
  pdf(file = paste0(folder_name, "/leave.pdf"), width = 10, height = 10)
  print(p3[[1]])
  dev.off()
}

