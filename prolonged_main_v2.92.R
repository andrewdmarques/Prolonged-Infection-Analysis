################################################################################
# Load Libraries
################################################################################

library('readxl')
library('ggplot2')
library('plotly')
library('zoo')
library('dplyr')
library('ggrepel') # For timeline (new)
library('tidyverse')
library('ggridges')

################################################################################
# User Defined Variables
################################################################################

prefix <- '2023-10-30'
suffix <- 'v5'
file_met3 <- paste0(prefix,'_prolonged-covid-samples_',suffix,'.csv')
min_change <- 0.10
max_mut <- 15
dir_data_name <- 'Data'

## Read in the data 
mut1 <- read.csv('2023-06-12_mutations_minor_variant_v3.csv')
met1 <- readxl::read_xlsx('2023-09-05_prolonged-covid-samples_v80.xlsx')

# Specify where the functions are.
source('D:/Desktop/OneDrive/Bushman-Lab-Cloud/03_Coding-Scripts/30_Minor-Variants-Heatmap/2023-10-30_long-term_functions_v2.92.R')

################################################################################
# Initialize
################################################################################

# Make the directory for data
dir_data <- paste0(dir_data_name,'_',prefix,'_',suffix)
if (!file.exists(dir_data)) dir.create(dir_data)

## Clean the data.
# Convert the DNA sequenced T to the RNA genome U.
mut1$REF <- gsub('T','U',mut1$REF)
mut1$ALT <- gsub('T','U',mut1$ALT)

# Make a lookup table to assign correct subject information to for each VSP.
# Assign the data to be found in the table.
lookup_id <- met1$study_id
# Give the data names.
names(lookup_id) <- met1$VSP
# Search a vector of names and outputs the associated data.
mut1$subject <- lookup_id[mut1$VSP]

# met0 <- subset(met1, sample_included_in_study %in% c("yes")) # We will not use this because we are using the statistical outliers
met0 <- subset(met1, sample_included_in_study_statistically %in% c("yes"))

# Get main heatmap
subj_all <- unique(met0$study_id)[1:length(unique(met0$study_id))]
mut6_single <- get_heatmap(met0,mut1,prefix,subj_all,dir_data)

# Get the main mutations data frame mut6.
met0_no_reinfection <- subset(met0, met0$study_id != 663) # 663 is our likely reinfection
met1_no_reinfection <- subset(met1,met1$study_id != 663)  # 663 is our likely reinfection
mut6_with_reinfection <- get_mut6(met0,mut1,prefix,subj_all,dir_data)
mut6 <- subset(mut6_with_reinfection,!grepl('663',mut6_with_reinfection$group))

# Sort the data frame so that the plots have the correct order
mut6 <- mut6[order(mut6$group, decreasing = FALSE),]

# Generate the iSNV variability plots by subject.
get_variability_plots(mut6)

# Determine the pairwise genetic diversity
dis5 <- get_daily_pairwise_genetic_diversity(mut6_with_reinfection)

################################################################################
# Statistical Model
################################################################################

# Function to save the summary plots for the analysis.
results <- get_raw_time_plots(mut6,dir_data,min_change,max_mut,10)

# Function for logistic regression for mutations for a given subject.
results2 <- get_regression_plots(mut6,prefix,suffix,min_change,max_mut)

# Function for logistic regression for mutations for a all subjects together.
stat1 <- get_regression_plots_large_groups(mut6,prefix,suffix,min_change,max_mut)

write.csv(stat1,paste0(dir_data,'/',prefix,'_regression-stats_',suffix,'.csv'),row.names = F)

# Assess the statistics from the linear regression.
stat2 <- read.csv(paste0(dir_data,'/',prefix,'_regression-stats_',suffix,'.csv'))
# Get just the data for the slope.
stat3 <- subset(stat2,stat2$stat == 'slope')
# Get the gene and the position.
stat3$gene <- sub("_.*", "", stat3$mutation)
stat3$POS <- sub(".*_", "", stat3$mutation)
stat3$POS <- as.numeric(gsub("[^0-9]", "", stat3$POS))
stat3 <- stat3[order(stat3$POS, decreasing = FALSE),]
stat3$gene <- factor(stat3$gene,levels = unique(stat3$gene))

# Choose only those values which are nonzero
stat4 <- subset(stat3, stat3$pred != 0)

# Rate of change by genes.
stat4_gene_p_value <- get_stat_plots(prefix,dir_data,stat4)

################################################################################
# General Data Summary (sample summary table)
################################################################################

# Make a blank data frame
col <- c('subjects','subjects_with_data','total_samples','total_samples_located','total_samples_used',
         'samples_per_subject','timepoints_per_subject','subjects_with_more_than_3_timepoints',
         'number_days_covered')
row <- c('value','min','max')
sum1 <- data.frame(matrix(NA, nrow = length(row), ncol = length(col)))
colnames(sum1) <- col
rownames(sum1) <- row

# Get the total number of subjects.
sum1$subjects[1] <- length(unique(met1$study_id))
# Get the number of subjects with data.
sum1$subjects_with_data[1] <- length(unique(met0$study_id))
# Get the total number of samples biobanked.
sum1$total_samples[1] <- nrow(met1)
# Get the total number of samples located.
sum1$total_samples_located[1] <- sum(grepl("located", met1$located_for_prolonged_covid_project, ignore.case = TRUE))
# Get the total number of samples used.
sum1$total_samples_used[1] <- nrow(met0)
# Get the value for median, min, and max number of samples used.
summary_df1 <- aggregate(VSP ~ study_id, data = met0, FUN = function(x) length(unique(x)))
summary_stats1 <- summary(summary_df1$VSP)
sum1$samples_per_subject[1] <- median(summary_stats1)
sum1$samples_per_subject[2] <- min(summary_stats1)
sum1$samples_per_subject[3] <- max(summary_stats1)
# Get the value for median, min, and max number of samples used.
summary_df2 <- aggregate(days_since_symptom_onset ~ study_id, data = met0, FUN = function(x) length(unique(x)))
summary_stats2 <- summary(summary_df2$days_since_symptom_onset)
sum1$timepoints_per_subject[1] <- median(summary_stats2)
sum1$timepoints_per_subject[2] <- min(summary_stats2)
sum1$timepoints_per_subject[3] <- max(summary_stats2)
# Get the number of subjects with more than 3 timepoints.
sum1$subjects_with_more_than_3_timepoints[1] <- sum(summary_df2 > 2, na.rm = TRUE)
# Get the number of days covered for each of the subjects.
met2 <- subset(met0,met0$study_id %in% subj_all)
mut2 <- subset(mut1, mut1$VSP %in% met2$VSP)
summary_df3 <- met2 %>%
  group_by(study_id) %>%
  summarize(max_days = max(days_since_symptom_onset) - min(days_since_symptom_onset))

# Format the table for display
sum2 <- data.frame(t(sum1))
# Replace NA with ''
sum2[is.na(sum2)] <- ''
# Replace underscores in row names with spaces
rownames(sum2) <- gsub("_", " ", rownames(sum2))
write.csv(sum2,paste0(dir_data,'/',prefix,'_sample-summary-table_',suffix,'.csv'))

################################################################################
# Getting minor variant detection (isnv-by-day)
################################################################################

# Geenrate a blank list that can record the plots.
plot_list <- list()

# Specify cutoff values for iSNV detection.
cutoffs <- c(0.01,0.03,0.05,0.1,0.15,0.2,0.25)

# Iterate through the cutoff values to assess affects of iSNV accumulation.
for (kk in 1:length(cutoffs)) {
  minor_variant_cutoff <- cutoffs[kk]
  # Get only the minor variants:
  mut7 <- subset(mut6,mut6$value >= minor_variant_cutoff)
  mut7 <- subset(mut7,mut7$value <= (1-minor_variant_cutoff))
  mut7$dayss <- as.numeric(gsub("[^0-9]", "", mut7$dim2))
  mut7$group <- gsub(' ','',mut7$group) # Remove the space from the group column.
  
  subj_all <- unique(mut7$group)
  for(ii in 1:length(subj_all)){
    gg <- subj_all[ii]
    temp1 <- subset(mut7,mut7$group == gg)
    
    # Get all the mutation date combinations.
    mm <- unique(paste0(temp1$dim1,'___',as.character(temp1$dayss)))
    
    # Make a blank data frame
    col <- c('mutation','dayss', 'value')
    # temp2 <- data.frame(matrix(NA, nrow = length(mm), ncol = length(col)))
    temp2 <- data.frame(matrix(NA, nrow = nrow(temp1), ncol = length(col)))
    colnames(temp2) <- col
    
    temp2$mutation <- temp1$dim1
    temp2$dayss <- temp1$dayss
    temp2$value <- temp1$value
    temp3 <- temp2
    temp3$remove <- FALSE
    # Remove mutations that don't have multiple detections of a mutation for a given day.
    for(ll in 1:nrow(temp3)){
      tt <- as.character(temp3$mutation[ll])
      uu <- temp3$dayss[ll]
      temp3a <- subset(temp3,temp3$mutation == tt)
      temp3b <- subset(temp3a,temp3a$dayss == uu)
      if(nrow(temp3b) < 2){temp3$remove[ll] <- TRUE} # It must have 2 or more occurances of a minor variant for that day for it to be counted.
    }
    
    temp3c <- subset(temp3,temp3$remove == FALSE) # Remove the mutations that did not have multiple occurances for a given day.
    temp3d <- temp3c[ , which(names(temp3c) %in% c('mutation','dayss'))] # Keep only the necessary columns
    temp3e <- temp3d[!duplicated(temp3d),] # Remove the duplicates so we have a count of how many iSNV there are 
    
    # Report the number of minor variants by day.
    temp4 <- data.frame(table(temp3e$dayss))
    colnames(temp4) <- c('dayss','minor_variants')
    temp4$dayss <- as.numeric(as.character(temp4$dayss))
    temp4$group <- gg
    
    # If make sure that days that had 0 minor variants are also recorded.
    tt <- subset(mut6,grepl(gg,mut6$group))
    dayss_all <- unique(as.character(tt$dim2))
    dayss_all2 <- unique(as.numeric(sapply(dayss_all, function(x) gsub("[^0-9]", "", x))))

    missing_days <- setdiff(dayss_all2, temp4$dayss)
    
    if(length(missing_days) > 0){
      # Create a data frame for the missing days
      missing_df <- data.frame(dayss = missing_days,
                               minor_variants = 0,
                               group = gg)
      
      # Combine the original data with the missing days data
      temp4 <- rbind(temp4, missing_df)
      
      # Order the dataframe by dayss for a better visualization
      temp4 <- temp4[order(temp4$dayss), ]
    }
    
    # Record the information for the group data.
    if(ii == 1){
      mut8 = temp4
    } else {
      mut8 <- rbind(mut8,temp4)
    }
    
    # Record the rates for each of the subjects.
    rate_temp <- as.numeric(coef(lm(minor_variants ~ dayss, data = temp4))[2])
    # Record the rates
    col <- c('cutoff', 'subject','rate')
    rate1 <- data.frame(matrix(NA, nrow = 1, ncol = length(col)))
    colnames(rate1) <- col
    rate1$cutoff[1] <- cutoffs[kk]
    rate1$subject[1] <- gg
    rate1$rate[1] <- rate_temp
    if(kk == 1 & ii == 1){
      rate_all <- rate1
    }else{
      rate_all <- rbind(rate_all,rate1)
    }
  }  
  
  # Determine the y limit based on the lowest cutoff value.
  if(kk == 1){ylimit <- max(mut8$minor_variants) + 10}
  
  # Perform linear regression
  model <- lm(minor_variants ~ dayss, data = mut8)
  
  # Define color palette
  color_palette <- colorRampPalette(c("red", "#FFA630","#4DA1A9", "#2E5077"))(length(unique(mut8$group)))
  
  # Create the plot
  p <- ggplot(mut8, aes(x = dayss, y = minor_variants, color = group)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "solid") +
    labs(x = "Days Since Onset", y = "iSNV Counts") +
    scale_color_manual(values = color_palette) +
    theme_bw() +  # this sets up the theme with white background
    theme(
      panel.grid.major = element_blank(),  # this removes major gridlines
      panel.grid.minor = element_blank(),  # this removes minor gridlines
      panel.background = element_blank(),  # this removes background color/lines
    ) +
    ggtitle(paste0("Cutoff = ",as.character(minor_variant_cutoff))) +
    ylim(0, ylimit) + 
    guides(color = "none") +
    theme(text = element_text(size = 12),  # this changes global text size
          axis.title = element_text(size = 12),  # axis labels
          axis.text = element_text(size = 12),  # axis tick labels
          plot.title = element_text(size = 12),  # plot title
          legend.title = element_text(size = 12),  # legend title
          legend.text = element_text(size = 12))  # legend items
  
  # Display the plot
  print(p)
  
  # Add the plot to the list
  plot_list[[kk]] <- p
}

# What is the mean rate of emergence?
aa <- subset(rate_all,rate_all$cutoff == 0.01)
mean(aa$rate)
min(aa$rate)
max(aa$rate)

# Add a legend.
# Create the color palette
color_palette <- colorRampPalette(c("red", "#FFA630","#4DA1A9", "#2E5077"))(length(unique(mut8$group)))
# Create a data frame for the legend
legend_data <- data.frame(
  name = unique(mut8$group),
  color = color_palette
)
# Create the legend plot
p <- ggplot(legend_data) +
  geom_tile(aes(x = 1, y = name, fill = I(color)), width = 0.3) +  # Adjust width here
  geom_text(aes(x = 1, y = name, label = name), hjust = 0.5, color = "white") +
  xlim(0.5, 1.5) +  # Add space on the left and right
  theme_void() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  scale_fill_identity()

# Add the plot to the list
plot_list[[kk + 1]] <- p

# Convert dimensions from mm to inches
width_in_inches <- 170 / 25.4
height_in_inches <- 300 / 25.4

# Save all the plots in a single PDF with the specified dimensions
pdf(paste0(dir_data,'/',prefix,'_isnv-by-day_',suffix,'.pdf'), width = width_in_inches, height = height_in_inches, onefile = TRUE)
gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
dev.off()

################################################################################
# Plot the patient status timeline
################################################################################

# User defined variables 
day_earliest <- -30
max_cell <- 10

# Open the file containing the prolonged infection patient status.
file_pt <- '2023-10-11_prolonged-infection-patient-status.csv'
pt1 <- read.csv(file_pt)

# Clean the data.
pt1$date <- as.Date(pt1$date)
pt1$date_end <- as.Date(pt1$date_end) 

# Generate the immune timeline plot.
get_immune_plot(pt1,day_earliest,max_cell)
