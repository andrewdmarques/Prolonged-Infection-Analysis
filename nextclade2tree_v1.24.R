################################################################################
# Load libraries
################################################################################

library(ape)
library(ggplot2)
library(pdftools) # To combine pdf files together.

################################################################################
# User defined parameters
################################################################################

dir_working <- '/media/andrewdmarques/Data01/Bioinformatics/39_Nextstrain/RTT'
prefix <- '2023-10-24'
suffix <- 'v1'
subjs <- c('486','637', '640','641')
cols <- c("#FF0000", "#FFA630", "#4DA1A9", "#2E5077")
nextclade_dir <- '/home/andrewdmarques/miniconda3/bin/nextclade'
nextclade_dataset <- '/media/andrewdmarques/Data01/Bioinformatics/39_Nextstrain/nextclade-sars-cov-2-dataset/sars-cov-2'

################################################################################
# Functions
################################################################################
# Function to determine what the file names of the subjects are.
get_files <- function(subjs){
  file_df <- data.frame(subj = character(), file = character(), stringsAsFactors = FALSE)
  
  for (subj in subjs) {
    file_pattern <- subj
    matching_files <- list.files(pattern = file_pattern, ignore.case = TRUE)
    
    if (length(matching_files) > 0) {
      file_df <- rbind(file_df, data.frame(subj = subj, file = matching_files[1], stringsAsFactors = FALSE))
    }
  }
  
  return(file_df)
}

# Function to combine pdf files together.
stitch_subject_pdfs <- function(directory, output_file) {
  # Find PDF files containing "subject" in their name
  pdf_files <- list.files(directory, pattern = "subject", full.names = TRUE)
  if (length(pdf_files) == 0) {
    print("No PDF files containing 'subject' in their name found in the directory.")
    return(NULL)
  }
  # Combine them with the other one
  pdf_combine(pdf_files, output = output_file)
  # Combine the pdf files.
  pdf_length(output_file)
  
  print("PDF files stitched successfully!")
}


################################################################################
# Initialize
################################################################################

# Set the working directory.
setwd(dir_working)

# Create directory if it doesn't exist
data_dir <- paste0(prefix,'_Data_',suffix)
if (!file.exists(data_dir)) {dir.create(data_dir)}

################################################################################
# Run
################################################################################

# Iterate through the samples.
for(ii in 1:length(subjs)){
  subj <- subjs[ii]
  dir_out <- paste0(data_dir,'/',prefix,'_',subj,'_Data_',suffix)
  if(!file.exists(dir_out)){dir.create(dir_out)}
  dir_nextclade_out <- paste0(dir_out,'/Nextclade-Output')
  file_nextclade_aligned <- paste0(dir_nextclade_out,'/nextclade.aligned.fasta')
  file_nextclade_aligned_new <- paste0(dir_out,'/',prefix,'_',subj,'_nextclade-aligned_',suffix,'.fasta')
  files <- list.files(dir_working)
  file_met1 <- files[grepl(subj,files) & grepl('metadata_v',files)]
  
  # Determine the metadata and sequence files from nextclade
  sub_m1 <- files[grepl(subj,files) & grepl('Prolonged',files) & grepl('metadata',files)]
  sub_s1 <- files[grepl(subj,files) & grepl('Prolonged',files) & grepl('sequence',files)]
  
  # unzip the nextclade files.
  system(paste("unxz", sub_m1))
  system(paste("unxz", sub_s1))
  
  # Determine which sequence should be the root.
  met1 <- read.csv(file_met1,sep = '\t')
  root_acc <- met1$strain[1] # The root should be the first sequence in this metadata file
  vsps <- unlist(regmatches(met1$strain, gregexpr("VSP\\d+", met1$strain))) # The prolong infected samples should be the remaining sequences in this file.
  
  # Save the medata for all subjects
  if(ii == 1){
    met1_temp <- met1[2:nrow(met1),]
    met1_temp$subject <- subj
    met_sam_all <- met1_temp
  }else{
    met1_temp <- met1[2:nrow(met1),]
    met1_temp$subject <- subj
    met_sam_all <- rbind(met_sam_all,met1_temp)
  }
  
  # Add dates to the end of the sequences from their metadata sheet.
  files_update <- list.files(dir_working) # Pull the files that are present
  sub_m2 <- files_update[grepl(subj,files_update) & grepl('Prolonged',files_update) & grepl('metadata',files_update)] # Locate the subsamples metadata
  sub_s2 <- files_update[grepl(subj,files_update) & grepl('Prolonged',files_update) & grepl('sequence',files_update)] # Locate the subsampled sequence data
  sub_seq1 <- read.table(sub_s2,header = F) # Open the sequence file
  sub_met1 <- read.csv(sub_m2,sep = '\t') # Open the metadata file
  root_header <- 'none'
  for(jj in 1:nrow(sub_seq1)){
    header <- sub_seq1$V1[jj]
    if(grepl('>',header)){
      temp <- subset(sub_met1,sub_met1$strain == gsub('>','',header))
      if(nrow(temp) > 0){
        # Write the new header
        header_new <- paste0(header,'__',as.character(temp$date[1]))
        # Determine if it is the root
        if(grepl(root_acc,header)){
          print(paste0("Root located as: ",header))
          header_new <- paste0(header,'__root__',as.character(temp$date[1]))
          root_header <- gsub('/','_',gsub('>','',header_new))
        }
        
        # Determine if it is one of hte prolong infected samples.
        if (any(sapply(vsps, function(item) grepl(item, header)))) {
          print(paste0("Prolonged infection sample located as: ",header))
          header_new <- paste0(header,'_--_',subj,'_-_prolonged-infection__',as.character(temp$date[1]))
        }
        # Assign the new header.
        sub_seq1$V1[jj] <- header_new
      }else{
        print(paste0("No matching date for: ",header))
      }
    }
  }
  
  # Save the new file.
  file_seq_renamed <- paste0(dir_out,'/',prefix,'_',subj,'_sequences-renamed_',suffix,'.fasta')
  write.table(sub_seq1,file_seq_renamed,row.names = F,quote = F,col.names = F)

  # Align the sequences using nextclade.
  system(paste0(nextclade_dir,' run --input-dataset ',nextclade_dataset,' --output-all=',dir_nextclade_out,'/ ', file_seq_renamed))

  # Move the aligned file to a unique location and file name.
  file.copy(file_nextclade_aligned,file_nextclade_aligned_new)

  # Generate a tree using iqtree only if one is not yet made.
  dir_iqtree_out <- paste0(dir_out,'/IQTree-Output')
  if(!file.exists(dir_iqtree_out)){dir.create(dir_iqtree_out)}
  files_iqtree <- list.files(paste0(dir_iqtree_out,'/'))
  if(length(files_iqtree)< 1){
    if(root_header == 'none'){
      system(paste0('iqtree -s ',file_nextclade_aligned_new,' -m TEST -bb 1000 -alrt 1000 -nt 6'))
      print('root was not found')
    }else{
      system(paste0('iqtree -s ',file_nextclade_aligned_new,' -m TEST -bb 1000 -alrt 1000 -nt 6 -o ',root_header))
    }
    # Save the data in iqtree directory.
    system(paste0('mv ',dir_out,'/*.fasta.* ',dir_iqtree_out,'/'))
  }

  ##############################################################################
  
  # Generate root to tip plot.
  files_iqtree <- list.files(paste0(dir_iqtree_out,'/'))
  file_tree <- files_iqtree[grepl('.treefile',files_iqtree)]
  file_tree <- paste0(dir_iqtree_out,'/',file_tree[1])
  
  # Read the tree file
  tree <- read.tree(file_tree)
  
  # Pick an outgroup taxon (e.g. the taxon with the smallest tip label)
  if(ii == 1){
    outgroup <- min(tree$tip.label)
  }else{
    outgroup <- root_header
  }
  
  # Root the tree at the outgroup
  tree2 <- root(tree, outgroup)
  
  # Check the resulting tree
  plot(tree2, main="Outgroup Rooted Tree")
  
  # Get the distances from the root to each node
  dis1 <- cophenetic(tree2)
  dis2 <- data.frame(t(subset(dis1,colnames(dis1) == outgroup)))
  # Make a data frame with all the information.
  dis3 <- dis2
  colnames(dis3) <- c('distance')
  # Get the sample name.
  dis3$name <- rownames(dis3)
  # Remove the rows with unknown dates.
  dis3 <- subset(dis3,!grepl('_NA',dis3$name))
  # Get the date.
  dis3$name <- gsub('__root','_-_root',dis3$name)
  dis3$name <- gsub('__prolonged','_-_prolonged',dis3$name)
  dis3$date <- as.Date(sub("^.*?__", "", dis3$name))
  # dis3 <- subset(dis3,dis3$date > as.Date('2021-12-01'))
  # Convert the distance to number of mutations from the root.
  dis3$distance <- dis3$distance*29903 # 29903 is the genome length.
  #plot(dis3$date,dis3$distance)
  
  # Remove sequence that appeared before the root.
  if(ii != 1){
    rt <- subset(dis3,grepl('root',dis3$name))
    # If there is a root found then choose it as the date otherwise choose the earliest time as the root date.
    if(length(rt[,1]) > 0){
      rt_date <- rt$date[1]
    }else{
      rt_date <- min(dis3$date)
    }
    dis3 <- subset(dis3,dis3$date >= rt$date[1])
  }else{
    rt_date <- min(dis3$date)
  }

  # Modify the distances of prolonged infection samples to have the distance from the earliest time.
  dis_pro1 <- data.frame(dis1[grep("_-_prolonged", rownames(dis1)), ])
  # Determine what the earliest subject timepoint would be.
  subj_min <- subset(met1_temp,met1_temp$date == min(met1_temp$date)) 
  subj_min <- subj_min[1,]
  subj_vsp <- sub(".*-([A-Z]+[0-9]+)/.*", "\\1", subj_min$strain) # This recovers the VSP for the earliest timepoint. 
  subj_isolates <- rownames(dis_pro1)
  subj_iso <- subj_isolates[grep(subj_vsp, subj_isolates)] # This returns the isolate name for the earliest collected sample for the subject.
  dis_pro2 <- data.frame(dis_pro1[ , which(names(dis_pro1) %in% c(gsub('-','.',subj_iso)))])
  dis_pro2[,1] <- dis_pro2[,1]*29903
  rownames(dis_pro2) <- rownames(dis_pro1)
  
  # Combine the data ising the RTT for the background and the Tip-to-Tip for the prolonged infection.
  lookup_tab_dis <- dis_pro2[,1]
  names(lookup_tab_dis) <- rownames(dis_pro2)
  lookup_result <- lookup_tab_dis[dis3$name]
  # If lookup_result has NA (not found in the lookup table), then keep the existing dis3$distance value, 
  # otherwise use the value from lookup_result.
  dis3$distance <- ifelse(is.na(lookup_result), dis3$distance, lookup_result)
  
  # Convert the date into days since "root" where "root" is the earliest occurrence in the US for the background samples and the first timepoint for the prologned infected individuals.
  start_root <- as.Date(min(dis3$date))
  start_subj <- as.Date(subj_min$date)
  dis3$day <- 0
  dis3$day <- ifelse(grepl("_-_prolonged", dis3$name),
                     as.numeric(difftime(as.Date(dis3$date), start_subj, units="days")),
                     as.numeric(difftime(as.Date(dis3$date), start_root, units="days")))
  
  # Subset data for samples with "prolonged"
  prolonged <- subset(dis3,grepl("_-_prolonged",dis3$name))
  
  # Subset data for samples without "prolonged"
  non_prolonged <- subset(dis3,!grepl("_-_prolonged",dis3$name) & !grepl("_-_root",dis3$name) )
  
  # Get the linear regression model for non-prolonged infection.
  non_prolonged_model <- lm(distance ~ date, data = non_prolonged)
  conf_intervals <- confint(non_prolonged_model, level = 0.95)
  non_prolonged_lower <- round(30 * conf_intervals[2,1], 2)
  non_prolonged_upper <- round(30 * conf_intervals[2,2], 2)
  non_prolonged_slope <- round(30 * coef(lm(distance ~ date, data = non_prolonged))[2], 2)
  
  # Get the linear regression model for prologned infection.
  prolonged_model <- lm(distance ~ date, data = prolonged)
  conf_intervals <- confint(prolonged_model, level = 0.95)
  prolonged_lower <- round(30 * conf_intervals[2,1], 2)
  prolonged_upper <- round(30 * conf_intervals[2,2], 2)
  prolonged_slope <- round(30 * coef(lm(distance ~ date, data = prolonged))[2], 2)
  
  # Plot the data with separate regression lines for each subset
  ggplot() +
    geom_point(data = non_prolonged, aes(x = date, y = distance),
               color = "black", size = 2) +
    geom_point(data = prolonged, aes(x = date, y = distance),
               color = cols[ii], size = 2) +
    geom_smooth(data = non_prolonged, aes(x = date, y = distance),
                method = "lm", color = "black") +
    geom_smooth(data = prolonged, aes(x = date, y = distance),
                method = "lm", fill = cols[ii],color = cols[ii]) +
    scale_color_identity() +
    scale_x_date(date_labels = "%b %Y") +
    labs(title = paste0("Subject ",  subj,"\nBackground = ",
                        non_prolonged_slope," (",non_prolonged_lower,", ",non_prolonged_upper,")\nProlonged Infection = ",prolonged_slope," (",prolonged_lower,", ",prolonged_upper,")"),
         x = "Date", y = "Distance (Number of Mutations)") +
    theme_bw()+
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14))
  ggsave(paste0(data_dir,'/',prefix,'_',suffix,'_root-to-tip_subject_',subj,'.pdf'),width = 10, height = 6)
  
  # Save the data dataframe.
  write.csv(dis3,paste0(dir_out,'/',prefix,'_',subj,'_rtt-data_',suffix,'.csv'),row.names = F)
  
  # Save a dataframe of all the plots.
  dis4 <- dis3
  rownames(dis4) <- NULL
  dis4$subject <- subj
  if(ii == 1){
    rtt1 <- dis4
  }else{
    rtt1 <- rbind(rtt1,dis4)
  }
  write.csv(rtt1,paste0(data_dir,'/',prefix,'_rtt-data_',suffix,'.csv'),row.names = F,quote = F)
}

# Combine the root-to-tip pdf files.
stitch_subject_pdfs(data_dir,paste0(data_dir,'/',prefix,'_',suffix,'_root-to-tip_summary.pdf'))



