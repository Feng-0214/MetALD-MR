Mendelian Randomization (MR) Analysis for Multiple Outcomes
===
To streamline the MR analysis process and avoid repetitive tasks, we have optimized the MR analysis code. This repository demonstrates a batch-processing approach to MR analysis, using the causal relationship between alcohol consumption and metabolites as an example.

Project Background
---
With the introduction of MASLD (Metabolic dysfunction-associated steatotic liver disease) and MetALD (Metabolic dysfunction-associated alcoholic liver disease) as distinct subtypes of fatty liver disease, accurately classifying and differentiating these conditions has become a significant challenge. Currently, clinical differentiation between MASLD and MetALD primarily relies on self-reported alcohol consumption. However, this approach is fraught with limitations due to recall bias, social desirability bias, and stigma associated with alcohol consumption, making it inherently unreliable. Therefore, there is an urgent need for objective and quantifiable biomarkers to enable more precise differentiation between these diseases.

The rapid advancement of lipidomics offers new possibilities for identifying such biomarkers. Using liver MRI data from the UK Biobank, this study analyzed the lipidomic profiles of MASLD and MetALD patients, uncovering significant differences in plasma lipidomics between the two groups. Notably, HDL-related markers were significantly elevated in MetALD and ALD patients, indicating distinct pathophysiological differences in lipoprotein metabolism between alcohol-related fatty liver disease and MASLD.

To further validate whether these metabolites can serve as potential biomarkers for distinguishing MASLD from MetALD, we conducted a Mendelian Randomization (MR) analysis at the genetic level.

Study Design
---
Data on alcohol consumption and metabolites were obtained from two large-scale genome-wide association studies (GWAS), focusing on individuals of European ancestry. Alcohol consumption was defined based on the average weekly intake of alcohol. All data are accessible through the IEU OpenGWAS summary database.

* The GWAS ID for alcohol consumption is `ieu-b-73`.<br>
* GWAS IDs for metabolites can be found in the `gwas_ids.csv` file within the "`Data`" folder.

Below is the flowchart illustrating the study design.

![Figure S1_00](https://github.com/user-attachments/assets/d294b1e9-e457-491a-86d3-be78ecaa64ce)

How to Reproduce This Analysis
---
The main MR analysis can be reproduced using the code in LOOP code.R by following these steps:

1. Ensure Required R Packages Are Installed and Loaded<br>
Install and load the necessary R packages, which are used for data processing, extracting instrumental variables, and conducting MR analyses.

2. Extract Instrumental Variables for Alcohol Consumption<br>
Extract SNPs associated with alcohol consumption from the GWAS data and calculate the F-statistic to evaluate the strength of the SNPs as instrumental variables.

3. Load Metabolite GWAS IDs<br>
Read the CSV file (`gwas_ids.csv`) containing the GWAS IDs and names of metabolites to retrieve all target GWAS IDs and corresponding names.

4. Account for Confounding Factors<br>
Download SNP information related to confounding factors (e.g., smoking, BMI, abnormal liver function, and diabetes) from the `GWAS Catalog`. As shown in the `Confounder` folder:
* Read the SNP data related to these confounders.<br>
* Combine them into a single dataframe to exclude SNPs associated with confounding factors.

5. Construct a Batch Loop for Metabolite Analysis
Create a batch loop to process data for each metabolite and perform MR analysis. Save the results for each metabolite in corresponding folders. The loop executes the following steps for each GWAS ID:
* Extract relevant outcome data.
* Harmonize the SNP alleles for consistency across datasets.
* Exclude SNPs associated with confounding factors.
* Conduct Mendelian randomization analysis, along with heterogeneity and pleiotropy tests.
* Perform MR-PRESSO analysis to detect outliers. If no outliers are identified, directly obtain the MR-PRESSO results. If outliers are detected (P < 0.05), remove them, rerun the MR analysis, and calculate the corrected MR results.
  
6. Perform `Leave-One-Out` Analysis and Visualization
Conduct "leave-one-out" analyses on the MR results and generate visualizations.

7. Combine and Organize Results
Use the code in Code-results to merge and organize results from all files, ensuring easy access to the final outputs.



