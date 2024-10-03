# Prolonged COVID Samples Analysis

## Description:
This script performs an analysis on prolonged SARS-CoV-2 infection samples. It processes raw mutation data, calculates statistical models, generates heatmaps, visualizes iSNV variability by subject, and produces daily pairwise genetic diversity metrics. The end goal is to determine how minor variants change over time and gain insights into prolonged COVID cases.

## Dependencies:
- R packages: `readxl`, `ggplot2`, `plotly`, `zoo`, `dplyr`, `ggrepel`, `tidyverse`, `ggridges`

## Usage:
1. Ensure all dependencies are installed.
2. Set the `prefix`, `suffix`, and other user-defined variables accordingly.
3. Make sure the provided datasets (`mutations_minor_variant` and `prolonged-covid-samples`) are available in the paths specified.
4. The external R script (`long-term_functions_v2.92.R`) should be available in the specified directory.

## Main Steps:

1. **Loading Libraries and Data:** Loads necessary libraries and raw data.
  
2. **Data Initialization:** Establishes directories and transforms data, including DNA to RNA conversions and data subsetting.
   
3. **Heatmap Generation:** Creates heatmaps for visualization of mutations.
    
4. **Statistical Analysis:** Determines pairwise genetic diversity, conducts regression analysis on mutations, and assesses statistics from the linear regression.
    
5. **Data Summarization:** Produces a general data summary table showcasing key metrics like the number of subjects, total samples, and more.
    
6. **Minor Variant Detection:** Analyzes the influence of different minor variant cutoffs on the accumulation of iSNVs over time.

## Outputs:

- `long-term_functions_v2.92.R`: This script contains functions used in the main code.
- A directory named based on the `dir_data_name`, `prefix`, and `suffix` variables. This directory contains several outputs like heatmaps, regression statistics, and sample summaries.
- `regression-stats.csv`: Contains statistics from the regression analysis.
- `sample-summary-table.csv`: A table summarizing various data points.

