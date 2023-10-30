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
# Define Functions
################################################################################

# Function to get the mutations in a standard format
get_mutation <- function(mut_df){
  
  mut_df$protein_product <- NA
  # Assign the protein product
  pro1 <- seq(1:29903)
  col <- c('POS', 'protein_product')
  pro2 <- data.frame(matrix(NA, nrow = length(pro1), ncol = length(col)))
  colnames(pro2) <- col
  pro2$POS <- pro1
  pro2$protein_product <- 'intergenic'
  
  for(jj in 1:length(pro2$POS)){
    if(pro2$POS[jj] >= 266 && pro2$POS[jj] < 805){pro2$protein_product[jj] <- "nsp1"}
    if(pro2$POS[jj] >= 806 && pro2$POS[jj] < 2719){pro2$protein_product[jj] <- "nsp2"}
    if(pro2$POS[jj] >= 2720 && pro2$POS[jj] < 8554){pro2$protein_product[jj] <- "nsp3"}
    if(pro2$POS[jj] >= 8555 && pro2$POS[jj] < 10054){pro2$protein_product[jj] <- "nsp4"}
    if(pro2$POS[jj] >= 10055 && pro2$POS[jj] < 10972){pro2$protein_product[jj] <- "nsp5 (3CL-PRO)"}
    if(pro2$POS[jj] >= 10973 && pro2$POS[jj] < 11842){pro2$protein_product[jj] <- "nsp6"}
    if(pro2$POS[jj] >= 11843 && pro2$POS[jj] < 12091){pro2$protein_product[jj] <- "nsp7"}
    if(pro2$POS[jj] >= 12092 && pro2$POS[jj] < 12685){pro2$protein_product[jj] <- "nsp8"}
    if(pro2$POS[jj] >= 12686 && pro2$POS[jj] < 13024){pro2$protein_product[jj] <- "nsp9"}
    if(pro2$POS[jj] >= 13025 && pro2$POS[jj] < 13441){pro2$protein_product[jj] <- "nsp10"}
    if(pro2$POS[jj] >= 13442 && pro2$POS[jj] < 16236){pro2$protein_product[jj] <- "nsp12 (RdRp)"}
    if(pro2$POS[jj] >= 16237 && pro2$POS[jj] < 18039){pro2$protein_product[jj] <- "nsp13 (Hel)"}
    if(pro2$POS[jj] >= 18040 && pro2$POS[jj] < 19620){pro2$protein_product[jj] <- "nsp14 (ExoN)"}
    if(pro2$POS[jj] >= 19621 && pro2$POS[jj] < 20658){pro2$protein_product[jj] <- "nsp15 (EndoU)"}
    if(pro2$POS[jj] >= 20659 && pro2$POS[jj] < 21552){pro2$protein_product[jj] <- "nsp16 (2'-O-MT)"}
    if(pro2$POS[jj] >= 21563 && pro2$POS[jj] < 25381){pro2$protein_product[jj] <- "spike"}
    # if(pro2$POS[jj] >= 21599 && pro2$POS[jj] < 23617){pro2$protein_product[jj] <- "S1"}
    # if(pro2$POS[jj] >= 23618 && pro2$POS[jj] < 25381){pro2$protein_product[jj] <- "S2"}
    if(pro2$POS[jj] >= 25393 && pro2$POS[jj] < 26217){pro2$protein_product[jj] <- "ORF3a"}
    if(pro2$POS[jj] >= 26245 && pro2$POS[jj] < 26469){pro2$protein_product[jj] <- "envelope"}
    if(pro2$POS[jj] >= 26523 && pro2$POS[jj] < 27188){pro2$protein_product[jj] <- "membrane"}
    if(pro2$POS[jj] >= 27202 && pro2$POS[jj] < 27384){pro2$protein_product[jj] <- "ORF6"}
    if(pro2$POS[jj] >= 27439 && pro2$POS[jj] < 27756){pro2$protein_product[jj] <- "ORF7a"}
    if(pro2$POS[jj] >= 27756 && pro2$POS[jj] < 27884){pro2$protein_product[jj] <- "ORF7b"}
    if(pro2$POS[jj] >= 27939 && pro2$POS[jj] < 28256){pro2$protein_product[jj] <- "ORF8"}
    if(pro2$POS[jj] >= 28274 && pro2$POS[jj] < 29530){pro2$protein_product[jj] <- "nucleocapsid"}
    if(pro2$POS[jj] >= 29558 && pro2$POS[jj] < 29671){pro2$protein_product[jj] <- "ORF10"}
  }
  
  
  # Assign the data to be found in the table.
  pro_tab <- pro2$protein_product
  # Give the data names.
  names(pro_tab) <- pro2$POS
  # Search a vector of names and outputs the associated data.
  mut_df$protein_product <- pro_tab[mut_df$POS]
  
  # Assign the nucleotide mutation for all of the groups.
  mut_df$nt_mutation <- paste0(mut_df$REF,mut_df$POS,mut_df$ALT)
  
  
  # Determine the aa mutation
  mut_df$mutation <- paste0(mut_df$protein_product,'_',mut_df$type,'_',mut_df$REF, as.character(mut_df$POS),mut_df$ALT)
  # Clean up the mutations.
  for(ii in 1:length(mut_df$VSP)){
    if(grepl('intergenic',mut_df$mutation[ii])){mut_df$mutation[ii] <- paste0(mut_df$protein_product[ii],'_',mut_df$REF[ii], as.character(mut_df$POS[ii]),mut_df$ALT[ii])}
    if(grepl('del',mut_df$mutation[ii])){mut_df$mutation[ii] <- paste0(mut_df$protein_product[ii],'_del_',as.character(nchar(mut_df$ALT[ii])-3),'_',as.character(mut_df$POS[ii]))}
    if(grepl('ins',mut_df$mutation[ii])){mut_df$mutation[ii] <- paste0(mut_df$protein_product[ii],'_ins_',as.character(nchar(mut_df$ALT[ii])-3),'_',as.character(mut_df$POS[ii]))}
  }
  return(mut_df)
}

# Function to make a dataframe that is a x b into one that is a*b long (square data frame to list data frame).
linearize <- function(data_frame){
  r <- rownames(data_frame)
  c <- colnames(data_frame)
  
  # Make a blank data frame.
  col <- c('dim1','dim2', 'value')
  linear <- data.frame(matrix(NA, nrow = length(r)*length(c), ncol = length(col)))
  colnames(linear ) <- col
  
  x <- 1
  y <- 1
  for(i in 1:length(linear$dim1)){
    linear$dim1[i] <- r[x]
    linear$dim2[i] <- c[y]
    linear$value[i] <- data_frame[x,y]
    y <- y + 1
    if(y > length(c)){
      y <- 1
      x <- x + 1
    }
  }
  return(linear)
}

# Check if the position for a given sample had coverage.
check_coverage <- function(mut5,mut3_no_cov,sam1){
  mut_t <- mut5
  col_t <- colnames(mut5) 
  # mut3_no_cov <- mut3_no_cov[1:0, ] # This will remove the zero coverage for the plot.
  
  # Use lookup table to assign the name for the no coverage data frame.
  # Assign the data to be found in the table.
  lookup_tab <- sam1$name
  # Give the data names.
  names(lookup_tab) <- sam1$VSP
  # Search a vector of names and outputs the associated data.
  mut3_no_cov$name <- lookup_tab[mut3_no_cov$VSP]
  
  # Determine what positions are used in the mutation table
  # Split each value in "dim1" by the underscore character
  parts <- strsplit(mut_t$dim1, "_")
  # Extract the last element of each list and remove all non-numeric characters
  mut_t$alias <- paste0(mut_t$dim2,'_',as.character(sapply(parts, function(x) gsub("[^0-9]", "", x[[length(x)]]))))
  
  # Use a lookup table to determine if the position had no coverage.
  mut3_no_cov$alias <- paste0(mut3_no_cov$name,'_',as.character(mut3_no_cov$POS))
  # Assign the data to be found in the table.
  lookup_cov <- rep(-1,length(mut3_no_cov$alias))
  # Give the data names.
  names(lookup_cov) <- mut3_no_cov$alias
  # Search a vector of names and outputs the associated data.
  mut3_no_cov$name <- lookup_cov[mut3_no_cov$VSP]
  mut_t$temp <- lookup_cov[mut_t$alias]
  # Check which rows have NA in the "temp" column
  na_rows <- is.na(mut_t$temp)
  # Overwrite the values in "value" with the values in "temp" for rows where "temp" is not NA
  mut_t$value <- ifelse(!na_rows, mut_t$temp, mut_t$value)
  
  # Curate the data so that it is the same format that was put in.
  mut_t <- mut_t[ , which(names(mut_t) %in% col_t)]
}

# Function to get a dataframe that has the sample type and numbers for converting in other functions.
get_type_num <- function(met2){
  samp_type_unique <- unique(met2$sample_type_tube)
  # Make a blank data frame
  col <- c('sample_type', 'num')
  type_num <- data.frame(matrix(NA, nrow = length(samp_type_unique), ncol = length(col)))
  colnames(type_num) <- col
  type_num$sample_type <- samp_type_unique
  type_num$num <- seq(1,length(samp_type_unique))
  return(type_num)
}

# Function that retrieves the sample type and adds it to the data frame to be plotted.
get_sample_type <- function(met2,mut4,mut6,sam1,type_num,prop_samp){
  # # Go from the day data to the VSP data using a lookup table.
  mut6$name <- paste0(mut6$group,'_',mut6$dim2)
  # Make a lookup table to assign data to data frames more efficiently.
  # Assign the data to be found in the table.
  lookup_tab <- sam1$VSP
  # Give the data names.
  names(lookup_tab ) <- sam1$name
  # Search a vector of names and outputs the associated data.
  mut6$VSP <- lookup_tab[mut6$name]
  
  # Use a lookup table to get the sample type for each of the VSP
  samp <- unique(mut6$name)
  # Get xx number of rows to display the sampe data, for example, 2% would be .02
  xx <- max(1,round(length(mut4[,1])*prop_samp)) # Doing the max here makes sure that there is at least a value of 1
  samp_type_unique <- unique(met2$sample_type_tube)
  
  # Record the factors and temprarily make them characters.
  df_levels <- levels(mut6$dim1)
  mut6$dim1 <- as.character(mut6$dim1)
  samp_type_rep <- paste0('samp_type','-',seq(1,xx))
  blanks <- paste0('blank-',seq(1,xx))
  df_levels <- c(samp_type_rep,blanks,df_levels)
  
  # Get a conversion for the sample types.
  # Assign the data to be found in the table.
  lookup_num <- type_num$num
  # Give the data names.
  names(lookup_num) <- type_num$sample_type
  
  for(ii in 1:length(samp)){
    col <- colnames(mut6)
    temp <- data.frame(matrix(NA, nrow = xx, ncol = length(col)))
    colnames(temp) <- col
    temp2 <- subset(mut6, mut6$name == samp[ii])
    samp_vsp <- temp2$VSP[1]
    samp_type <- subset(met2,met2$VSP == samp_vsp)
    samp_type <- samp_type$sample_type_tube[1]
    
    # Populate the data frame with sample types.
    temp$dim1 <- paste0(rep('samp_type',xx),'-',seq(1,xx))
    temp$dim2 <- temp2$dim2[1]
    temp$group <- temp2$group[1]
    temp$value <- (lookup_num[samp_type] + 2)*-1
    
    # Add the blank lines.
    temp_b <- temp
    temp_b$dim1 <- paste0('blank-',seq(1,xx))
    temp_b$value <- -1
    
    # Add it to the existing data frame.
    mut6 <- rbind(mut6,temp,temp_b)
  }
  
  # Return the dim1 column back to a factor with the correct levels.
  mut6$dim1 <- factor(mut6$dim1,levels = df_levels)
  
  
  # Remove any columns that shouldn't be there.
  mut6 <- mut6[ , which(names(mut6) %in% c('dim1','dim2','group','value'))]
  return(mut6)
}

# Function to get the gradient legend for the mutation heatmap.
get_gradient_legend <- function() {
  # Create a dummy data frame
  df <- data.frame(x = c(0, 1), y = c(0, 1))
  
  # Plot
  plot_legend <- ggplot(df, aes(x = x, y = y)) + 
    geom_tile(aes(fill = x), height = 1, width = 1) + 
    scale_fill_gradientn(colors = c("lightblue", "darkblue"), 
                         name = "iSNV Proportion", 
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = c("0", "0.25", "0.5", "0.75", "1"),
                         guide = guide_colorbar(ticks = FALSE)) +  # Removing ticks
    theme_minimal() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.width = unit(1, "cm")) +
    labs(x = NULL, y = NULL)
  # Save the plot to a PDF
  ggsave(filename = paste0(dir_data,'/', prefix,'_heatmap-legend_',suffix,'.pdf'), plot = plot_legend, device = "pdf", width = 7, height = 5)
  
}

# Function to get the sample type legend for the mutation heatmap.
get_legend <- function(colors, labels){
  # Check that the inputs are of the correct length
  if (length(colors) != length(labels)) {
    stop("Error: the lists 'colors' and 'labels' must have the same length.")
  }
  
  # Create a blank plot
  plot <- ggplot()
  
  # Add a dummy layer to the plot
  plot <- plot + 
    geom_point(aes(x = 1, y = 1, color = labels, fill = labels), 
               size = 10, shape = 15, stroke = 1)
  
  # Add a legend to the plot
  plot <- plot + 
    scale_color_manual(values = colors, name = "Sample Type", labels = labels) +
    scale_fill_manual(values = colors, name = "Sample Type", labels = labels)
  
  # Remove the axis labels and tick marks
  plot <- plot + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  plot
  return(plot)
}


# Function to make the heatmap.
get_heatmap <- function(met0,mut1,prefix,subj_all,dir_data){
  # get_heatmap(met0,mut1,prefix3,subj3,dir_data)
  # Make a blank that all the plots are saved to.
  plot_list <- list()
  
  # Iterate through each of the subjects.
  for(ss in 1:length(subj_all)){
    # Record how long it takes for each sample. 
    start_time_seconds <- as.numeric(Sys.time())
    
    subj <- subj_all[ss]
    print(paste0('Beginning ', subj),quote = F)
    
    # Look at only the samples are included in the project.
    # met2 <- subset(met1, sample_included_in_study %in% c("yes", "poor coverage"))
    met2 <- subset(met0,met0$study_id %in% subj)
    mut2 <- subset(mut1, mut1$VSP %in% met2$VSP)
    
    # Extract the mutations that are present.
    mut3 <- get_mutation(mut2)
    dir.create(paste0(dir_data,'/'), showWarnings = FALSE)
    write.csv(mut3,paste0(dir_data,'/',prefix,'_minor-variant-examples_',suffix,'.csv'),row.names = F)
    
    # Separate the positions with no coverage from those with coverage.
    mut3_no_cov <- subset(mut3,mut3$ALT == 'N')
    mut3 <- subset(mut3,mut3$ALT != 'N')
    
    # Get the sample data frame set up.
    col <- c('VSP', 'passage','replicate','host','type')
    sam1 <- data.frame(matrix(NA, nrow = length(unique(met2$VSP)), ncol = length(col)))
    colnames(sam1) <- col
    sam1$VSP <- met2$VSP
    sam1$passage <- met2$days_since_symptom_onset
    sam1$host <- met2$study_id
    sam1 <- sam1[order(sam1$passage, decreasing = FALSE),]
    
    # Format the sample sheet.
    # sam1$name <- paste0(sam1$VSP,'_', sam1$host,'_Day ',as.character(sprintf("%02d",as.numeric(sam1$passage))))
    
    # Combine the values in the "host" and "passage" columns XXX
    sam1$name <- paste(sam1$host,'_',sam1$VSP,' Day ',as.character(formatC(as.numeric(sam1$passage), width = 3, format = "d", flag = "0")))
    sam1$name <- paste(sam1$host,'_','Day ',as.character(formatC(as.numeric(sam1$passage), width = 3, format = "d", flag = "0")))
    
    # Add suffixes to the "name" column for duplicate rows
    sam1$name <- ave(sam1$name, sam1$name, FUN = function(x) {
      if(length(x) > 1) {
        suffix <- paste(".", letters[1:length(x)], sep="")
        return(paste(x, suffix, sep=""))
      } else {
        return(x)
      }
    })
    
    # Populate the data frame to be plotted.
    col <- sam1$name
    row <- unique(mut3$mutation)
    mut4 <- data.frame(matrix(0, nrow = length(row), ncol = length(col)))
    colnames(mut4) <- col
    rownames(mut4) <- row
    for(ii in 1:length(col)){
      for(jj in 1:length(row)){
        t1 <- col[ii]
        hold <- subset(sam1,sam1$name == t1)
        t1 <- hold$VSP[1]
        t2 <- row[jj]
        temp <- subset(mut3,mut3$VSP == t1)
        temp <- subset(temp,temp$mutation == t2)
        if(length(temp$VSP > 0)){
          mut4[jj,ii] <- temp$percentAlt[1]
        }
      }
    }
    
    # Make columns where there is no mutation (there was no consensus sequence made) equal to -1 so it is clear that there is no coverage.
    mut4[, apply(mut4, 2, function(x) all(x == 0))] <- -1
    
    # Remove the rows where the mutation proption is 1 for all timepoints.
    mut4 <- mut4[rowSums(mut4 >= 0.90) < ncol(mut4), ]
    
    mut5 <- linearize(mut4)
    
    # Determine if the specified position for each of the samples had no coverage.
    mut5 <- check_coverage(mut5,mut3_no_cov,sam1)
    
    # Specify the order of the variables, displayed in the order of positions.
    index <- mut3[ , which(names(mut3) %in% c('POS','mutation'))]
    index <- index[order(index$POS, decreasing = FALSE),]
    index <- index[!duplicated(index),]
    # index$POS[length(index$POS)-1] <- index$POS[(length(index$POS)-2)] + 2
    # index$POS[length(index$POS)] <- index$POS[(length(index$POS)-2)] + 1
    index <- index[order(index$POS, decreasing = TRUE),]
    levels_factor <- index$mutation
    mut5$dim1 <- factor(mut5$dim1,levels = levels_factor)
    
    # Format the names correctly.
    mut6 <- mut5
    mut6$group <- gsub("_.*","",mut6$dim2)
    mut6$dim2 <- gsub(".*_","",mut6$dim2)
    
    # # Add the sample type to the plot.
    # type_num <- get_type_num(met2)
    # prop_samp <- 0.02
    # mut6 <- get_sample_type(met2,mut4,mut6,sam1,type_num,prop_samp)
    
    # Make the days in the mutation table a factor.
    factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", unique(mut6$dim2), perl = TRUE))
    factor_v <- sort(as.character(unique(mut6$dim2)))
    factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", factor_v, perl = TRUE))
    mut6$dim2 <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", mut6$dim2, perl = TRUE))
    mut6$dim2 <- factor(mut6$dim2,levels = factor_v)
    # mut6$dim2 <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", mut6$dim2, perl = TRUE))
    
    # Make the heatmap
    p_cols <- c("#693a01","#9e5600",   "#ffd6a6","#bfaf9b",'white','lightblue', 'darkblue')
    p_resc <- c(-5,-4,-3,-2,-1,0, 1)
    
    # For those plots without sample type.
    p_cols <- c('lightblue', 'darkblue')
    p_resc <- c(0, 1)
    
    # If there are no "no coverage" mutations, then just us light blue or dark blue.
    if(min(mut6$value) < 0){
      p_cols <- c('white','lightblue', 'darkblue')
      p_resc <- c(-1,0, 1)
    }
    
    heatmap_plot <- ggplot(data = mut6, aes(x=dim2, y=dim1)) +
      facet_grid(~group, scales = "free") +
      geom_tile(aes(fill = value)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
      xlab('') + ylab('') + labs(fill = "Proportion") + ggtitle("") +
      scale_fill_gradientn(
        # colours=c('grey','white', '#c2e7fc'),
        #colours=c("#693a01","#9e5600","#ffae4f","#ffd6a6","#bfaf9b",'white','lightblue', 'darkblue'),
        colours=p_cols,
        # colours=c('lightblue', 'darkblue'),
        values=scales::rescale(p_resc),
        guide="colorbar")+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank() 
      )
    
    # Add the plot to the list
    plot_list[[ss]] <- heatmap_plot
    
    # Report back the progress going through this loop.
    end_time_seconds <- as.numeric(Sys.time())
    print(paste0('Completed ', subj,'     Seconds: ',as.character(round(end_time_seconds - start_time_seconds),2)),quote = F)
    print('',quote = F)
  }
  
  # Save all the plots in a single PDF
  pdf(paste0(dir_data,'/', prefix,'_mutation-heatmap_',suffix,'.pdf'), width = ceiling(ss)*1.8, height = ceiling(ss)*2)
  gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(ss/2))
  dev.off()
  desired_height <- 100
  pdf(paste0(dir_data,'/', prefix,'_mutation-heatmap-paper_',suffix,'.pdf'), width = 375/25.4, height = desired_height/25.4)  # 25.4 mm per inch
  gridExtra::grid.arrange(grobs = plot_list, ncol = length(plot_list))
  dev.off()
  return(mut6)
}

# Function to get the mut6 data frame that has all the data sorted.
get_mut6 <- function(met0,mut1,prefix,subj,dir_data){
  
  # prefix <- prefix3
  # subj <- subj3
  
  # Look at only the samples are included in the project.
  # met2 <- subset(met1, sample_included_in_study %in% c("yes", "poor coverage"))
  met2 <- subset(met0,met0$study_id %in% subj)
  mut2 <- subset(mut1, mut1$VSP %in% met2$VSP)
  
  # Extract the mutations that are present.
  mut3 <- get_mutation(mut2)
  dir.create(paste0(dir_data,'/'), showWarnings = FALSE)
  write.csv(mut3,paste0(dir_data,'/',prefix,'_minor-variant-examples_',suffix,'.csv'),row.names = F)
  
  # Separate the positions with no coverage from those with coverage.
  mut3_no_cov <- subset(mut3,mut3$ALT == 'N')
  mut3 <- subset(mut3,mut3$ALT != 'N')
  
  # Get the sample data frame set up.
  col <- c('VSP', 'passage','replicate','host','type')
  sam1 <- data.frame(matrix(NA, nrow = length(unique(met2$VSP)), ncol = length(col)))
  colnames(sam1) <- col
  sam1$VSP <- met2$VSP
  sam1$passage <- met2$days_since_symptom_onset
  sam1$host <- met2$study_id
  sam1 <- sam1[order(sam1$passage, decreasing = FALSE),]
  
  # Format the sample sheet.
  # sam1$name <- paste0(sam1$VSP,'_', sam1$host,'_Day ',as.character(sprintf("%02d",as.numeric(sam1$passage))))
  
  # Combine the values in the "host" and "passage" columns XXX
  sam1$name <- paste(sam1$host,'_',sam1$VSP,' Day ',as.character(formatC(as.numeric(sam1$passage), width = 3, format = "d", flag = "0")))
  sam1$name <- paste(sam1$host,'_','Day ',as.character(formatC(as.numeric(sam1$passage), width = 3, format = "d", flag = "0")))
  
  # Add suffixes to the "name" column for duplicate rows
  sam1$name <- ave(sam1$name, sam1$name, FUN = function(x) {
    if(length(x) > 1) {
      suffix <- paste(".", letters[1:length(x)], sep="")
      return(paste(x, suffix, sep=""))
    } else {
      return(x)
    }
  })
  
  # Populate the data frame to be plotted.
  col <- sam1$name
  row <- unique(mut3$mutation)
  mut4 <- data.frame(matrix(0, nrow = length(row), ncol = length(col)))
  colnames(mut4) <- col
  rownames(mut4) <- row
  for(ii in 1:length(col)){
    for(jj in 1:length(row)){
      t1 <- col[ii]
      hold <- subset(sam1,sam1$name == t1)
      t1 <- hold$VSP[1]
      t2 <- row[jj]
      temp <- subset(mut3,mut3$VSP == t1)
      temp <- subset(temp,temp$mutation == t2)
      if(length(temp$VSP > 0)){
        mut4[jj,ii] <- temp$percentAlt[1]
      }
    }
  }
  
  # Make columns where there is no mutation (there was no consensus sequence made) equal to -1 so it is clear that there is no coverage.
  mut4[, apply(mut4, 2, function(x) all(x == 0))] <- -1
  mut5 <- linearize(mut4)
  
  # Determine if the specified position for each of the samples had no coverage.
  mut5 <- check_coverage(mut5,mut3_no_cov,sam1)
  
  # Specify the order of the variables, displayed in the order of positions.
  index <- mut3[ , which(names(mut3) %in% c('POS','mutation'))]
  index <- index[order(index$POS, decreasing = FALSE),]
  index <- index[!duplicated(index),]
  # index$POS[length(index$POS)-1] <- index$POS[(length(index$POS)-2)] + 2
  # index$POS[length(index$POS)] <- index$POS[(length(index$POS)-2)] + 1
  index <- index[order(index$POS, decreasing = TRUE),]
  levels_factor <- index$mutation
  mut5$dim1 <- factor(mut5$dim1,levels = levels_factor)
  
  # Format the names correctly.
  mut6 <- mut5
  mut6$group <- gsub("_.*","",mut6$dim2)
  mut6$dim2 <- gsub(".*_","",mut6$dim2)
  
  # Add the sample type to the plot.
  type_num <- get_type_num(met2)
  prop_samp <- 0.02
  # mut6 <- get_sample_type(met2,mut4,mut6,sam1,type_num,prop_samp)
  
  # Make the days in the mutation table a factor.
  factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", unique(mut6$dim2), perl = TRUE))
  factor_v <- sort(as.character(unique(mut6$dim2)))
  factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", factor_v, perl = TRUE))
  mut6$dim2 <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", mut6$dim2, perl = TRUE))
  mut6$dim2 <- factor(mut6$dim2,levels = factor_v)
  # mut6$dim2 <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", mut6$dim2, perl = TRUE))
  
  # Save the intermediate mutation file
  write.csv(mut4,paste0(dir_data,'/',prefix,'_mutation-table_',suffix,'.csv'))
  
  # Save the sample decoding file.
  sam2 <- sam1
  sam2$name <- gsub(".*_","",sam2$name)
  # Make the days in the mutation table a factor.
  factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", unique(sam2$name), perl = TRUE))
  factor_v <- sort(as.character(unique(sam2$name)))
  factor_v <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "", factor_v, perl = TRUE))
  sam2$name <- gsub('000','0',gsub("(?<=\\s)0+(?=[1-9])", "",sam2$name, perl = TRUE))
  sam2$name <- factor(sam2$name,levels = factor_v)
  write.csv(sam2,paste0(dir_data,'/',prefix,'_sample-vsp-replicates_',suffix,'.csv'),row.names = F)
  
  return(mut6)
}

# Function for logistic regression for mutations for a given subject.
# get_regression_plots(mut6,prefix,suffix)
get_regression_plots <- function(df,prefix,suffix,min_change,max_mut){
  # Add a new column called 'dayss' that contains only the numbers from the 'dim2' column
  df$dayss <- as.numeric(gsub("\\D", "", df$dim2))
  df$dayss <- df$dim2
  
  # Remove the sample type and the blank containing lines.
  df <- df[!grepl("samp_typ", df$dim1),]
  df <- df[!grepl("blank", df$dim1),]
  
  # Split the data frame by unique values in the 'group' column
  groups <- split(df, df$group)
  
  # Initialize an empty list to store the results
  results <- list()
  
  # Prepare to save all plots as a summary.
  plot_list <- list()
  
  # Loop through each group and perform the specified operations
  for (kk in 1:length(names(groups))) {
    name <- names(groups)[kk]
    # Extract the data for the current group
    group_df <- groups[[name]]
    
    # Create a temp data frame with the appropriate dimensions
    rr <- nrow(group_df)/length(unique(group_df$dayss))
    cc <- length(unique(group_df$dayss))
    temp_df <- data.frame(matrix(NA, rr, cc))
    colnames(temp_df) <- unique(group_df$dayss)
    
    # Populate the temp data frame with the appropriate values from the 'value' column
    for(j in 1:length(1:cc)){
      temp0 <- subset(group_df,group_df$dayss == unique(group_df$dayss)[j])
      for (i in 1:length(temp0$dim1)) {
        row <- temp0[i,]
        temp_df[i, j] <- row$value[1]
      }
    }
    
    # Define a custom function that determines the maximum range of values in a row, excluding values of -1
    max_range <- function(x) {
      x <- x[x != -1]  # Exclude values of -1
      return(max(x) - min(x))  # Calculate the maximum range of values
    }
    
    # Add a new column to the data frame called 'max_range' that applies the custom function to each row
    temp_df$max_range <- apply(temp_df, 1, max_range)
    
    # Replace -Inf with 0.
    temp_df$max_range <- ifelse(is.infinite(temp_df$max_range), 0, temp_df$max_range)
    
    # Add the mutations to the file.
    lev <- unique(as.character(df$dim1))
    rownames(temp_df) <- lev
    
    ## Determine the mutations that should be plotted.
    min_change <- min_change
    max_mut <- max_mut
    # Subset only those that show changes over time.
    temp_p <- subset(temp_df,temp_df$max_range>min_change)
    # Remove all but the top max_mut mutations
    # Sort a data frame by column content
    temp_p <- temp_p[order(temp_p$max_range, decreasing = TRUE),]
    temp_p <- head(temp_p,max_mut)
    
    # Make the mutations a factor in the order of the positions that they occur.
    temp_p$pos <- as.numeric(gsub("[^[:digit:]]", "", sub('^.*_','',rownames(temp_p))))
    temp_p <- temp_p[order(temp_p$pos, decreasing = FALSE),]
    mut_levels <- rownames(temp_p)
    
    # Remove the 'max_range' column from the data frame
    temp_p <- temp_p[, !names(temp_p) %in% c("max_range","pos")]
    # Modify the column names by removing all letters
    colnames(temp_p) <- gsub("[^[:digit:]]", "", colnames(temp_p))
    
    if(nrow(temp_p)>0){
      # Linearize the results
      temp_p2 <- linearize(temp_p)
      temp_p2$value[temp_p2$value == -1] <- NA
      temp_p2$dim1 <- factor(temp_p2$dim1,levels = mut_levels)
      temp_p2$dim2 <- as.numeric(temp_p2$dim2)
      
      subj <- gsub(' ','',name)
      muts <- as.character(unique(temp_p2$dim1))
      
      # Open the PDF device with 4 columns
      r_col <- 4
      r_row <- ceiling(length(muts)/r_col)
      pdf(paste0(dir_data,'/',prefix,'_regression_',subj,'_',suffix,'.pdf'), width=3*r_col, height=(3*r_row), onefile=TRUE)
      par(mfrow=c(ceiling(length(muts)/r_col), r_col))
      
      for(ii in 1:length(muts)){
        mutation <- muts[ii]
        
        temp_p3 <- temp_p2
        
        temp_p3 <- subset(temp_p3,temp_p3$dim1 == mutation)
        temp_p3 <- temp_p3[complete.cases(temp_p3$value),]
        temp_p3$dim2 <- as.integer(temp_p3$dim2)
        
        # Create example data
        set.seed(123)
        
        temp_p3 <- temp_p3[ , -which(names(temp_p3) %in% c('dim1'))]
        deconstructed_df <- capture.output(dput(temp_p3))
        paste(deconstructed_df, collapse = "")
        
        ##############################################################################
        # Add a polynomial plot.
        # Fit a logistic regression model
        data_df <- temp_p3
        colnames(data_df) <- c('x','y')
        # Fit Bayesian second order polynomial regression model
        # library(rstanarm)
        # library(brms)
        
        # Fit a linear regression model
        model <- lm(y ~ x, data = data_df)
        
        # Get a summary of the model
        summary(model)
        
        # Confidence intervals for the coefficients
        confint(model)
        
        # Generate predictions and confidence intervals for the predictions
        newdata <- data.frame(x = seq(min(data_df$x), max(data_df$x), length.out = 100))
        predictions <- predict(model, newdata = newdata, interval = "confidence")
        
        plot(data_df$x, data_df$y, pch = 1, xlab = "Days Since Onset", ylab = "iSNV Proportion",
             main = mutation, ylim = c(-0.005, 1.005), xlim = c(min(data_df$x), 1 + max(data_df$x)), xaxs = 'i', yaxs = 'i')
        lines(newdata$x, predictions[,"fit"], col = "red")
        polygon(c(newdata$x, rev(newdata$x)), 
                c(pmax(predictions[,"lwr"], 0), rev(pmin(predictions[,"upr"], 1))),
                col = rgb(0.5, 0.5, 0.5, alpha = 0.2), border = NA)
      }
      
      # Close the PDF device
      dev.off()
    }
  }
}

# Function for logistic regression for mutations for a all subjects together.
get_regression_plots_large_groups <- function(df,prefix,suffix,min_change,max_mut){
  # Add a new column called 'dayss' that contains only the numbers from the 'dim2' column
  df$dim1 <- as.character(df$dim1)
  df$dim2 <- as.character(df$dim2)
  df$dayss <- df$dim2
  # df$dayss <- gsub("\\D", "", df$dayss)
  # df$dayss <- as.numeric(df$dayss)
  
  # Subset to only have the subjects with 3 or more timepoints.
  df$dayss2 <- gsub("\\D", "", df$dayss) # Get the numerical values for the days
  tp1 <- data.frame(table(df$dayss2,df$group)) # Get the counts for how many day/sample combinations have sequence data
  tp2 <- subset(tp1,tp1$Freq > 0) # Only keep the sample/day that have sequence data
  tp3 <- data.frame(table(tp2$Var2)) # Get the number of days for each sample.
  tp4 <- subset(tp3,tp3$Freq >= 3) # Get the samples that have 3 or more days for them.
  
  # Subset the main data frame to have only those samples with 3 or more timepoints.
  df <- subset(df, df$group %in% tp4$Var1) # tp4$Var1 is a list of samples that have enough timepoints.
  
  # Remove the sample type and the blank containing lines. values < 0 | values > 1
  df <- df[df$value >= 0 & df$value <= 1, ] # Remove any values that are outside the range 0 to 1.
  df <- df[!grepl("samp_typ", df$dim1),]
  df <- df[!grepl("blank", df$dim1),]
  
  # Split the data frame by unique values in the 'group' column
  groups <- split(df, df$group)
  
  # Initialize an empty list to store the results
  results <- list()
  
  # Prepare to save all plots as a summary.
  plot_list <- list()
  
  # Initialize the list that will determine which mutations will be used.
  muts_interest1 <- NULL
  muts_change <- NULL
  
  # Loop through each group and determine the list of mutations that should be included.
  for (kk in 1:length(names(groups))) {
    name <- names(groups)[kk]
    # Extract the data for the current group
    group_df <- groups[[name]]
    
    # Create a temp data frame with the appropriate dimensions
    row_mut <- unique(group_df$dim1)
    rr <- length(unique(group_df$dim1))
    cc <- length(unique(group_df$dayss))
    temp_df <- data.frame(matrix(NA, rr, cc))
    colnames(temp_df) <- unique(group_df$dayss)
    rownames(temp_df) <- row_mut
    
    # Populate the temp data frame with the appropriate values from the 'value' column
    for(j in 1:length(1:cc)){
      temp0 <- subset(group_df,group_df$dayss == unique(group_df$dayss)[j])
      for(i in 1:rr){
        tt <- rownames(temp_df)[i]
        temp1 <- subset(temp0,temp0$dim1 == tt)
        if(nrow(temp1) > 0){
          temp_df[i, j] <- temp1$value[1]
        }
      }
    }
    
    # Define a custom function that determines the maximum range of values in a row, excluding values of -1
    max_range <- function(x) {
      x <- x[x != -1]  # Exclude values of -1
      return(max(x) - min(x))  # Calculate the maximum range of values
    }
    
    # Add a new column to the data frame called 'max_range' that applies the custom function to each row
    temp_df$max_range <- apply(temp_df, 1, max_range)
    
    # Replace -Inf with 0.
    temp_df$max_range <- ifelse(is.infinite(temp_df$max_range), 0, temp_df$max_range)
    
    ## Determine the mutations that should be plotted.
    min_change <- min_change
    max_mut <- max_mut
    # Subset only those that show changes over time.
    temp_p <- subset(temp_df,temp_df$max_range>min_change)
    
    # Add the mutations to the mutation data frame.
    muts_interest1 <- c(muts_interest1,rownames(temp_p))
    muts_change <- c(muts_change,temp_p$max_range)
  }
  
  # Select the mutations of interest.
  muts_interest2 <- data.frame(mutation = muts_interest1, max_change = muts_change)
  muts_interest3 <- unique(muts_interest2$mutation)
  
  # Set the mutations so that it is ordred by position number.
  extract_numbers <- function(input_string) {
    # Split the string by underscore
    split_string <- strsplit(input_string, "_")[[1]]
    # Extract thelast element (numbers)
    last_element <- split_string[length(split_string)]
    # Extract only the numbers
    numbers <- gsub("[^0-9]", "", last_element)
    return(numbers)
  }
  # Apply the function to each input string
  output_numbers <- as.numeric(lapply(muts_interest3, extract_numbers))
  # Create a dataframe
  muts_interest3b <- data.frame(muts = muts_interest3, position = output_numbers)
  # Order dataframe by position
  muts_interest3b <- muts_interest3b[order(muts_interest3b$position), ]
  # Extract ordered muts_interest3 list
  muts_interest3 <- muts_interest3b$muts
  
  timepoints <- as.numeric(unique(gsub("\\D", "", df$dayss)))
  
  # Inititalize the plots to be made.
  # Open the PDF device with 4 columns
  r_col <- 6
  r_row <- ceiling(length(muts_interest3)/r_col)
  subj <- gsub(' ','',name)
  subj_file <- gsub("\n", "", subj)
  # Convert width from mm to inches (1 inch = 25.4 mm) for the PDF
  pdf_width_in_inches <- 8.5
  r_col <- 3 # number of columns
  r_row <- ceiling(length(muts_interest3) / r_col) # number of rows, assuming muts_interest3 contains the mutations of interest
  
  # Create PDF
  pdf(file = paste0(dir_data, '/', prefix, '_regression-summary-compare-hosts_', suffix, '.pdf'), 
      width = pdf_width_in_inches, 
      height = 3 * r_row, 
      onefile = TRUE)
  
  # Set up the plotting area with 3 columns and a number of rows based on the mutations, and set the font size
  par(mfrow = c(r_row, r_col), cex = 0.8, cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8) # Set font size to 12, which is typically 1.2 * default size (10)
  
  
  # Make a pdf that contains a plot for each mutation, and within the plot predictions for all the samples.
  for(ll in 1:length(muts_interest3)){
    t0 <- muts_interest3[ll]
    print(paste0(as.character(ll),'/',as.character(length(muts_interest3)),' ', muts_interest3[ll]),quote = F)
    # Break this out to be just the two groups that could be there.
    group_large <- unique(df$group)
    for(kk in 1:length(group_large)) {
      
      tt <- group_large[kk]
      group_df <- subset(df,grepl(tt,df$group))
      group_df2 <- subset(group_df,group_df$dim1 == t0)
      mut_levels <- t0 # This is where the levels come from for factoring the mutations in the plot.
      
      temp_p2 <- group_df2
      temp_p2$dayss <- gsub("\\D", "", temp_p2$dayss)
      temp_p2 <- temp_p2[ , which(names(temp_p2) %in% c('dim1','value','dayss'))]
      temp_p2 <- temp_p2[, c('dim1','dayss','value')]
      colnames(temp_p2) <- c('dim1','dim2','value')
      
      
      temp_p2$value[temp_p2$value == -1] <- NA
      temp_p2$dim1 <- factor(temp_p2$dim1,levels = mut_levels)
      temp_p2$dim2 <- as.numeric(temp_p2$dim2)
      
      # Save all the data to a new data frame.
      subj <- tt
      temp_p4 <- temp_p2
      temp_p4$subj <- subj
      
      # Initialize for the stats.
      temp_p5 <- temp_p4
      
      # # Run the modelling. #####################################################
      # Format the data frame for the data.
      data_df <- temp_p4
      data_df <- data_df[ , which(names(data_df) %in% c('dim2','value'))]
      colnames(data_df) <- c('x','y')
      
      # Fit a linear regression model
      mm <- 'linear'
      if(mm == 'linear'){
        model <- lm(y ~ x, data = data_df)
        # Generate predictions and confidence intervals for the predictions
        newdata_range <- max(data_df$x) - min(data_df$x) + 1
        newdata <- data.frame(x = seq(min(data_df$x), max(data_df$x), length.out = newdata_range))
        predictions <- predict(model, newdata = newdata, interval = "confidence")
        
        # plot(data_df$x, data_df$y, pch = 1, xlab = "Days Since Onset", ylab = "Proportion of mutation",
        #      main = t0, ylim = c(-0.005, 1.005), xlim = c(min(data_df$x), 1 + max(data_df$x)), xaxs = 'i', yaxs = 'i')
        # lines(newdata$x, predictions[,"fit"], col = "red")
        # polygon(c(newdata$x, rev(newdata$x)), 
        #         c(pmax(predictions[,"lwr"], 0), rev(pmin(predictions[,"upr"], 1))),
        #         col = rgb(0.5, 0.5, 0.5, alpha = 0.2), border = NA)
        
        # Record the stats results (ONLY FOR LINEAR REGRESSION)
        # Obtain the 95% confidence interval for the slope
        slope_ci <- data.frame(confint(model, level = 0.95))
        colnames(slope_ci) <- c('lower','upper')
        rownames(slope_ci) <- c('intercept','slope')
        # Obtain the slope.
        slope_pred <- predict(model, interval = "confidence")
        slope <- slope_pred[, "fit"]
        slope_ci$pred <- c(coef(model)[1],coef(model)[2])
        # Add the slope and the p values
        p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
        slope_ci$p_val <- p_values
        # Format the data frame to be added to the other data frames.
        slope_ci$group <- tt
        slope_ci$mutation <- t0
        slope_ci$stat <- c('intercept','slope')
        rownames(slope_ci) <- NULL
        
        # Record the results.
        temp_p5_pred <- data.frame(pred = predictions[,"fit"], upper = predictions[,"upr"], lower = predictions[,"lwr"])
        temp_p5_pred$dim1 <- t0
        temp_p5_pred$dim2 <- seq(from = min(data_df$x),to = max(data_df$x), by = 1)
        temp_p5_pred$subj <- tt
      }else{
        data <- data_df
        newdata1 <- data.frame(x = seq(min(data$x[is.finite(data$x)]), max(data$x[is.finite(data$x)]))) # Extra lines are added to make sure that if there are NA then it is not used for this newdata calculation.
        # Define a function to calculate the predicted values for a GLM model
        glm_predict <- function(model, newdata) {
          return(predict(model, newdata=newdata, type="response"))
        }
        
        # Set the number of bootstrap samples
        B <- 100
        
        # Create a function to generate a bootstrap sample and fit a GLM model
        glm_bootstrap <- function(data) {
          index <- sample(nrow(data), replace=TRUE)
          # model <- glm(y ~ x, data=data[index,], family=binomial)
          model <- glm(y ~ poly(x, 1), data=data, family=binomial)
          return(model)
        }
        
        # Generate the bootstrap samples and fit the GLM models
        models <- lapply(1:B, function(i) glm_bootstrap(data))
        
        # Calculate the predicted values for each bootstrap sample
        pred_matrix <- sapply(models, function(model) glm_predict(model, newdata1))
        
        # Calculate the mean and confidence interval for each value of x
        pred_mean <- apply(pred_matrix, 1, mean)
        pred_lower <- apply(pred_matrix, 1, function(x) quantile(x, 0.025))
        pred_upper <- apply(pred_matrix, 1, function(x) quantile(x, 0.975))
        
        # Clip the predicted values and confidence interval between 0 and 1
        pred_mean <- pmax(pmin(pred_mean, 1), 0)
        pred_lower <- pmax(pmin(pred_lower, 1), 0)
        pred_upper <- pmax(pmin(pred_upper, 1), 0)
        
        # model <- glm(y ~ poly(x, 1), data=data, family=binomial)
        # model <- gam(y ~ s(x, k = 2), data = data, family = betar)
        model <- glm(y ~ x, data=data, family=binomial)
        
        newdata1 <- data.frame(x = seq(min(data$x),max(data$x))) # For higher resolution, might put, by = 0.2 here, but there are other complications that would need to be worked out 
        # Make predictions for the sequence of x values
        pred <- predict(model, newdata=newdata1, type="response")
        
        # Calculate the confidence interval
        se <- sqrt(pred * (1 - pred) / nrow(newdata1)) # Standard error
        z <- qnorm(0.975) # Z-value for 95% confidence interval
        lower <- pred - z * se # Lower bound
        upper <- pred + z * se # Upper bound
        
        # Clip the predicted values and confidence interval between 0 and 1
        pred <- pmax(pmin(pred, 1), 0)
        lower <- pmax(pmin(lower, 1), 0)
        upper <- pmax(pmin(upper, 1), 0)
        
        # Record the results.
        temp_p5_pred <- data.frame(pred = pred, upper = upper, lower = lower)
        temp_p5_pred$dim1 <- t0
        temp_p5_pred$dim2 <- seq(min(data$x),max(data$x))
        temp_p5_pred$subj <- tt
        
        # Make a blank data frame for now.
        # Make a blank data frame
        col <- c('item', 'value')
        slope_ci <- data.frame(matrix(NA, nrow = 1, ncol = length(col)))
        colnames(slope_ci) <- col
        
      }
      
      # Record the results
      if(kk == 1 && ll == 1){
        dat1 <- temp_p5
        dat2 <- temp_p5_pred
        stat1 <- slope_ci
      }else{
        dat1 <- rbind(dat1,temp_p5)
        dat2 <- rbind(dat2,temp_p5_pred)
        stat1 <- rbind(stat1,slope_ci)
      }
    }
    # Plot the data. ###########################################################
    temp_p6 <- subset(dat1,dat1$dim1 == muts_interest3[ll])
    temp_p6_pred <- subset(dat2,dat2$dim1 == muts_interest3[ll])
    # Define a color palette
    # color_palette <- rainbow(length(unique(temp_p6$subj)))
    # Create a color palette
    unique_subjects <- unique(temp_p6$subj)
    color_palette <- colorRampPalette(c("#FFA630","#D7E8BA","#4DA1A9", "#2E5077", "#611C35"))(length(unique(df$group)))
    color_palette <- colorRampPalette(c("red", "#FFA630","#4DA1A9", "#2E5077"))(length(unique(df$group)))
    
    
    # Create an empty plot
    plot(1, 1,xlim = c((min(temp_p6_pred$dim2)-1),(max(temp_p6_pred$dim2)+1)), ylim = c(-0.005,1.005), xlab="Days Since Onset", ylab="iSNV Proportion", 
         xaxs = 'i', yaxs = 'i', type = "n", main = paste0(gsub('_',' ',muts_interest3[ll])))
    # Loop over unique subjects
    for (i in unique(temp_p6$subj)) {
      
      # Subset the data for this subject
      data_subj <- subset(temp_p6, subj == i)
      data_subj_pred <- subset(temp_p6_pred, subj == i)
      
      # Determine the color for this subject
      color_subj <- color_palette[which(unique(temp_p6$subj) == i)]
      
      # Add points for this subject to the plot
      points(data_subj$dim2, data_subj$value, pch=1, col=color_subj)
      
      # Add prediction line for this subject
      lines(data_subj_pred$dim2, data_subj_pred$pred, col=color_subj)
      
      # Add shaded confidence interval
      polygon(c(data_subj_pred$dim2, rev(data_subj_pred$dim2)), c(data_subj_pred$lower, rev(data_subj_pred$upper)), 
              col=adjustcolor(color_subj, alpha=0.2), border=NA)
    }
    
    # Add a legend
    if(ll == 1){legend("topright", legend=unique(temp_p6$subj), fill=color_palette)}
    # legend("topright", legend=unique(temp_p6$subj), fill=color_palette)
  }
  dev.off()
  par(mfrow = c(1, 1), cex = 1, cex.main = 1, cex.lab = 1, cex.axis = 1)
  return(stat1)
}

# Function to save the summary plots for the analysis.
# results <- get_raw_time_plots(mut6,dir_data,min_change,10)
get_raw_time_plots <- function(df_summ, dir_data, min_change, max_mut, color_palette = c("#F2C57C","#240B36","#92140C","#748386","#EF6F6C")) {
  # Add a new column called 'dayss' that contains only the numbers from the 'dim2' column
  df_summ$dayss <- df_summ$dim2
  
  # Remove the sample type and the blank containing lines.
  df_summ <- df_summ[!grepl("samp_typ", df_summ$dim1),]
  df_summ <- df_summ[!grepl("blank", df_summ$dim1),]
  
  # Split the data frame by unique values in the 'group' column
  groups <- split(df_summ, df_summ$group)
  
  # Initialize an empty list to store the results
  results <- list()
  
  # Prepare to save all plots as a summary.
  plot_list <- list()
  
  # Generate a color palette with a gradient based on the number of unique groups
  color_func <- colorRampPalette(color_palette)
  group_colors <- color_func(length(unique(df_summ$group)))
  
  # Loop through each group and perform the specified operations
  for (kk in 1:length(names(groups))) {
    name <- names(groups)[kk]
    # Extract the data for the current group
    group_df <- groups[[name]]
    
    # Create a temp data frame with the appropriate dimensions
    rr <- nrow(group_df)/length(unique(group_df$dayss))
    cc <- length(unique(group_df$dayss))
    temp_df <- data.frame(matrix(NA, rr, cc))
    colnames(temp_df) <- unique(group_df$dayss)
    
    # Populate the temp data frame with the appropriate values from the 'value' column
    for(j in 1:length(1:cc)){
      temp0 <- subset(group_df, group_df$dayss == unique(group_df$dayss)[j])
      for (i in 1:length(temp0$dim1)) {
        row <- temp0[i,]
        temp_df[i, j] <- row$value[1]
      }
    }
    
    # Define a custom function that determines the maximum range of values in a row, excluding values of -1
    max_range <- function(x) {
      x <- x[x != -1]  # Exclude values of -1
      return(max(x) - min(x))  # Calculate the maximum range of values
    }
    
    # Add a new column to the data frame called 'max_range' that applies the custom function to each row
    temp_df$max_range <- apply(temp_df, 1, max_range)
    
    # Replace -Inf with 0.
    temp_df$max_range <- ifelse(is.infinite(temp_df$max_range), 0, temp_df$max_range)
    
    # Add the mutations to the file.
    lev <- unique(as.character(df_summ$dim1))
    rownames(temp_df) <- lev
    
    ## Determine the mutations that should be plotted.
    min_change <- min_change
    max_mut <- max_mut
    # Subset only those that show changes over time.
    temp_p <- subset(temp_df, temp_df$max_range > min_change)
    # Remove all but the top max_mut mutations
    # Sort a data frame by column content
    temp_p <- temp_p[order(temp_p$max_range, decreasing = TRUE),]
    temp_p <- head(temp_p, max_mut)
    
    # Make the mutations a factor in the order of the positions that they occur.
    temp_p$pos <- as.numeric(gsub("[^[:digit:]]", "", sub('^.*_','', rownames(temp_p))))
    temp_p <- temp_p[order(temp_p$pos, decreasing = FALSE),]
    mut_levels <- rownames(temp_p)
    
    # Remove the 'max_range' and 'pos' columns from the data frame
    temp_p <- temp_p[, !names(temp_p) %in% c("max_range", "pos")]
    # Modify the column names by removing all letters
    colnames(temp_p) <- gsub("[^[:digit:]]", "", colnames(temp_p))
    
    if(nrow(temp_p) > 0){
      # Linearize the results
      temp_p2 <- linearize(temp_p)
      temp_p2$value[temp_p2$value == -1] <- NA
      temp_p2$dim1 <- factor(temp_p2$dim1, levels = mut_levels)
      temp_p2$dim2 <- as.numeric(temp_p2$dim2)
      
      
      temp_p2_new <- temp_p2 %>%
        group_by(dim1, dim2) %>%
        summarise(mean_value = mean(value, na.rm = TRUE),
                  sd_value = sd(value, na.rm = TRUE),
                  n = n(), .groups = 'drop') %>% # use .groups = 'drop' to avoid the grouped output message
        filter(!is.na(mean_value))
      
      # Update the temp_p2_new data frame to no longer include the nt mutation unless it is a silent mutation
      # Convert factor to character for string manipulation
      temp_p2_new$dim1 <- as.character(temp_p2_new$dim1)
      # Remove everything after the last underscore if "silent" is not in the string, otherwise keep as is
      temp_p2_new$dim1 <- ifelse(grepl("silent|intergenic", temp_p2_new$dim1), temp_p2_new$dim1,
                                 sub("_[^_]+$", "", temp_p2_new$dim1))
      
      # Determine what the levels should be in their specified order
      temp_levels <- ifelse(grepl("silent|intergenic", levels(temp_p2$dim1)), levels(temp_p2$dim1),
                            sub("_[^_]+$", "", levels(temp_p2$dim1)))
      # Convert the character column back to a factor
      temp_p2_new$dim1 <- factor(temp_p2_new$dim1,levels = temp_levels)
      
      # Generate the color palette here, after we know how many unique levels there are in temp_p2_new$dim1
      color_func <- colorRampPalette(color_palette)
      group_colors <- color_func(length(unique(temp_p2_new$dim1)))
      
      p <- ggplot(temp_p2_new, aes(x = dim2, y = mean_value, group = dim1, color = dim1)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = mean_value - sd_value/sqrt(n), ymax = mean_value + sd_value/sqrt(n)), width = 0.2) +
        geom_line(size = 1) +
        scale_color_manual(values = group_colors) + # use the colors from group_colors
        labs(x = "Days Since Onset", y = "Proportion", title = paste0(name), color = "Mutation") +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + # set y limits from 0 to 1
        theme(text = element_text(size = 12), # this changes global text size
              axis.title = element_text(size = 12), # axis labels
              axis.text = element_text(size = 12), # axis tick labels
              plot.title = element_text(size = 12), # plot title
              legend.title = element_text(size = 12), # legend title
              legend.text = element_text(size = 12), # legend items
              legend.position = "right", # this moves the legend below the plot
              panel.grid.major = element_blank(), # this removes major grid lines
              panel.grid.minor = element_blank(), # this removes minor grid lines
              panel.background = element_blank()) # this removes the panel background
      
      
      print(p)
      # Add the plot to the list
      plot_list[[kk]] <- p
    }
  }
  
  # Convert dimensions from mm to inches
  width_in <- 170 / 25.4
  height_in <- 210 / 25.4
  
  # Save all the plots in a single PDF
  pdf(paste0(dir_data,'/',prefix,'_summary-plot_',suffix,'.pdf'), width = width_in, height = height_in)
  gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
  dev.off()
  
  return(results)
}


# Function to get plot that has mutation on the x axis, proportion on the y axis, grouped by mutation and colored by day.
get_variability_plots <- function(mut6){
  isnv_cutoff <- 0.03
  
  # Assuming mut6 is already defined
  mut9 <- mut6
  mut9$day <- as.character(sub("\\..*", "", mut9$dim2))
  
  v_group <- unique(mut9$group)
  
  pdf_filename <- paste0(dir_data,'/',prefix,'_isnv-variability-plot_',suffix,'.pdf')
  pdf(pdf_filename, width = 8.5, height = 11)
  
  for(ii in 1:length(v_group)){
    tt <- v_group[ii]
    mut9_temp <- subset(mut9, mut9$group == tt)
    
    mut9_temp <- mut9_temp %>%
      mutate(day = gsub('Day','', day) %>%
               gsub(' ','', .) %>%
               as.numeric() %>%
               factor(levels = unique(.))) %>%
      arrange(day) %>%
      group_by(dim1) %>%
      summarise(count_in_range = sum(value >= isnv_cutoff & value <= (1-isnv_cutoff))) %>%
      filter(count_in_range >= 2) %>%
      inner_join(mut9_temp, by = "dim1")
    mut9_temp <- mut9_temp[rev(order(mut9_temp$dim1)), ]
    
    
    max_plots_per_page <- 15
    unique_dims <- unique(mut9_temp$dim1)
    num_pages <- ceiling(length(unique_dims) / max_plots_per_page)
    
    plot_chunk <- function(data, group_name) {
      data$plot_title <- paste(group_name, gsub('_', ' ', data$dim1), sep = "- ")
      
      # Determine how many unique plots there are
      unique_plots <- length(unique(data$plot_title))
      
      # Ensure that the days are in the correct format and numerical.
      data$day <- factor(as.numeric(gsub('Day ','',as.character(data$day))),levels = rev(unique(as.numeric(gsub('Day ','',as.character(data$day))))))
      
      # get the plot title in position order.
      titles <- unique(data$plot_title)
      # Function to extract the substring after the last space
      extract_after_last_space <- function(string) {
        # Find the position of the last space
        last_space <- max(regexpr(" ", string))
        
        # If no space is found (i.e., regexpr returns -1), return NA
        if (last_space == -1) {
          return(NA)
        } else {
          # Extract everything after the last space
          substring <- substr(string, last_space + 1, nchar(string))
          return(substring)
        }
      }
      # Apply the function to the text column
      titles2 <- sub(".*\\s+", "", titles)
      titles3 <- as.numeric(gsub("\\D", "", titles2))
      titles_df <- data.frame(titles = unlist(titles), titles3 = unlist(titles3), stringsAsFactors = FALSE)
      
      # Sort a data frame by column content
      titles_df <- titles_df[order(titles_df$titles3, decreasing = FALSE),]
      data$plot_title <- factor(data$plot_title,levels = titles_df$titles)
      
      # Ensure that there are enough plots to fill the whole page.
      if (unique_plots < 15) {
        last_plot <- tail(data[data$plot_title == unique(data$plot_title)[unique_plots], ], 1)
        replicate_count <- 15 - unique_plots
        
        new_levels <- levels(data$plot_title) # Get current levels of the factor
        
        for (i in 1:replicate_count) {
          new_title <- paste(last_plot$plot_title, i + 1, sep = " ") # Change title for each duplication
          new_levels <- c(new_levels, new_title) # Add the new title to the level set
          last_plot$plot_title <- factor(new_title, levels = new_levels) # Update the factor with new level
          data <- rbind(data, last_plot) # Add duplicated plot to the original data
        }
        
        # After the loop, update the levels of the 'plot_title' column in the original 'data'
        data$plot_title <- factor(data$plot_title, levels = new_levels)
      }
      
      color_func <- colorRampPalette(c('orange', 'red', 'black'))
      color_vector <- color_func(length(unique(data$day)))
      colors <- setNames(color_vector, sort(unique(data$day)))
      
      p <- ggplot(data, aes(x = day, y = value, color = factor(day))) +
        geom_jitter(width = 0.4, height = 0) +
        facet_wrap(~ plot_title, ncol = 3, nrow = 5) + # Use the new plot_title column for facetting
        scale_y_continuous(limits = c(0, 1)) +
        scale_color_manual(values = colors, name = "Days") +
        labs(y = "Proportion") +
        theme_minimal() +
        theme(text = element_text(size = 10))
      
      return(p)
    }
    
    for(i in seq_len(num_pages)) {
      current_dims <- unique_dims[((i-1)*max_plots_per_page + 1):min(i*max_plots_per_page, length(unique_dims))]
      current_data <- mut9_temp %>% filter(dim1 %in% current_dims)
      
      print(plot_chunk(current_data, tt))
    }
    
    
  }
  dev.off()
}

# Function that determines the genetic diversity across days. xxx
get_daily_pairwise_genetic_diversity <- function(mut6){
  mut9 <- mut6
  mut9$day <- as.character(sub("\\..*", "", mut9$dim2))
  v_levels <- levels(mut9$dim1)
  v_group <- unique(mut9$group)
  
  sam2 <- read.csv(paste0(dir_data,'/',prefix,'_sample-vsp-replicates_',suffix,'.csv'))
  sam2$host <- as.character(sam2$host)
  
  # Determine which samples should be the ones to keep.
  # Assign the data to be found in the table.
  lookup_tab <- met0$mean_coverage
  # Give the data names.
  names(lookup_tab) <- met0$VSP
  # Search a vector of names and outputs the associated data.
  sam2$mean_coverage <- lookup_tab[sam2$VSP]
  sam2$include = FALSE
  sam2$day <- sub("\\..*", "", sam2$name)
  sam2$host_day <- paste0(as.character(sam2$host),sam2$day)
  sam2$host_day_rep <- paste0(as.character(sam2$host),sam2$name)
  sam2 <- sam2 %>%
    group_by(host_day) %>%
    mutate(include = mean_coverage %in% sort(mean_coverage, decreasing = TRUE)[1:1]) %>% #XXX this is the line that selects the n highest mean coverage.
    ungroup()
  sam3 <- subset(sam2,sam2$include == TRUE)
  sam3$host_day_rep <- gsub(' ','',sam3$host_day_rep)
  
  # Remove all of the samples that are not to be included in the pairwise analysis.
  mut9$host_day_rep <- gsub(' ','',paste0(mut9$group,mut9$dim2))
  mut9$host_day <- paste0(mut9$group,mut9$day)
  mut10 <- subset(mut9, mut9$host_day_rep %in% sam3$host_day_rep)
  
  for(ii in 1:length(v_group)){
    tt <- v_group[ii]
    mut10_temp <- subset(mut10,mut10$group == tt)
    
    # Remove the word "Day" from the day column. 
    mut10_temp$day <- gsub('Day','',mut10_temp$day)
    mut10_temp$day <- gsub(' ','',mut10_temp$day)
    mut10_temp <- mut10_temp %>%
      mutate(day = as.numeric(day)) %>%
      arrange(day) %>%
      mutate(day = factor(day, levels = unique(day)))
    
    # Filter to only show mutations that have 1 or more values that are above 0.
    filtered_data <- mut10_temp %>%
      group_by(dim1) %>%
      filter(!all(value <= 0) & !all(value == 1)) %>%
      ungroup()
    mut10_temp <- filtered_data
    
    # Change day to be a numeric factor.
    mut10_temp$day <- as.numeric(as.character(mut10_temp$day))
    
    # Ensure that the mutations remain in their correct order.
    mut10_temp$dim1 <- factor(mut10_temp$dim1,levels = v_levels)
    
    p_groups <- unique(as.character(mut10_temp$dim2))
    for(jj in 2:length(p_groups)){
      # Extract sequence 1 to compare to.
      p_1 <- subset(mut10_temp, mut10_temp$dim2 == p_groups[jj])
      p_curr_day <- p_1$day[1]
      
      # Determine what the next lowest number is
      next_lowest <- mut10_temp %>%
        filter(as.numeric(day) < p_curr_day) %>%
        summarise(next_lowest = max(as.numeric(day), na.rm = TRUE)) %>%
        pull(next_lowest)
      
      # Isolate the next lowest sample.
      p_2 <- subset(mut10_temp, mut10_temp$day == next_lowest)
      
      # Combine both samples into one data frame to compare the pairwise distance.
      p_both <- rbind(p_1,p_2)
      
      # Calculate the pairwise distance
      dis1 <- p_both %>%
        # Self join on dim1
        inner_join(p_both, by = "dim1", suffix = c("_x", "_y")) %>%
        # Filter out non-unique pairs
        filter(dim2_x != dim2_y) %>%
        # Compute the difference in value for each dim1 category
        mutate(value_difference = value_x - value_y) %>%
        # Select the relevant columns
        select(dim1, dim2_x, dim2_y, value_difference)
      
      # Summarize the pairwise distance
      dis2 <- dis1 %>%
        # Group by dim2_x and dim2_y
        group_by(dim2_x, dim2_y) %>%
        # Sum the absolute value of the differences
        summarize(total_abs_difference = sum(abs(value_difference))) %>%
        # Ensure resulting dataframe doesn't retain grouping
        ungroup()
      
      if(jj == 2){
        dis3 <- dis2[1,]
      }else{
        dis3 <- rbind(dis3,dis2[1,])
      }
    }
    colnames(dis3) <- c('previous','current','pairwise_distance') 
    # dis3 <- dis3[ , -which(names(dis3) %in% c('remove'))]
    
    dis3$group <- tt
    
    if(ii == 1){
      dis4 <- dis3
    }else{
      dis4 <- rbind(dis4,dis3)
    }
  }
  
  # Adjust to have the mutation distance per day.
  dis5 <- transform(dis4,
                    previous_day = as.numeric(gsub("[^0-9]", "", previous)),
                    current_day = as.numeric(gsub("[^0-9]", "", current)))
  
  
  dis5$difference_day <- dis5$current_day - dis5$previous_day
  dis5$pairwise_distance_adj <- dis5$pairwise_distance/dis5$difference_day 
  
  write.csv(dis5,paste0(dir_data,'/',prefix,'_pairwise-distance_',suffix,'.csv'),row.names = F,quote = F)
  
  return(dis5)
}

# Determine the number of SNP difference for each comparison.
distance_of_snp <- function(mut1, vsp1, vsp2) {
  # vsp1 <- "VSP0770"
  # vsp2 <- "VSP0772"
  vsp1 <- as.character(vsp1)
  vsp2 <- as.character(vsp2)
  a <- subset(mut1, VSP==vsp1)   # Look at only the vsp of interest
  a$mut <- paste(a$genes,a$type, sep = "_") # Combine the mutations
  a$diff <- TRUE # Make a column that will be determined if it is a unique difference
  b <- subset(mut1, VSP==vsp2)
  b$mut <- paste(b$genes,b$type, sep = "_")
  b$diff <- TRUE
  # Mark the number of mutations that are different.
  for(i in 1:length(a$diff)){
    for(j in 1:length(b$diff)){
      #print(paste(as.character(i), as.character(j),sep = "_"))
      if(a$mut[i] == b$mut[j]){
        a$diff[i] <- FALSE
        b$diff[j] <- FALSE
      }
    }
  }
  # Count the number of different mutations.
  a_diff <- subset(a, diff == TRUE)
  a_total <- length(a$diff)
  b_diff <- subset(b, diff == TRUE)
  b_total <- length(b$diff)
  dis <- length(a_diff$diff) + length(b_diff$diff)
  return(dis)
}

# Rate of change by genes.
get_stat_plots <- function(prefix,dir_data,stat4){
  # # Do non-synonymous mutations.
  # stat4 <- subset(stat4, !grepl('intergenic',stat4$mutation))
  # stat4 <- subset(stat4, !grepl('silent',stat4$mutation))
  # 
  # # Do synonymous only
  # stat4 <- subset(stat4, grepl('intergenic', stat4$mutation) | grepl('silent', stat4$mutation))
  
  #############################################
  # Plot the directionality of slope by gene.
  # stat4 <- subset(stat4,stat4$p_val < 0.05)
  # plot pred distributions for each gene
  stat4a <- stat4
  stat4a$pred <- stat4$pred * 30
  ggplot(stat4a, aes(x = pred, y = gene)) +
    geom_density_ridges(scale = 3) +
    labs(title = "Directionality of Slope by Gene",
         x = "Change in iSNV Proportion per Month",
         y = "Gene")
  # Save the ggplot to a PDF file
  ggsave(paste0(dir_data,'/',prefix,'_directionality-of-slope-by-gene_',suffix,'.pdf'), width = 11.6, height = 4.9, device = "pdf")
  
  #############################################
  # Plot the change in proportion and p values on a scatter plot
  # if you want to plot a scatter plot to investigate the relationship between pred and p_val for each gene
  # first, we may want to log transform the p-values for a more interpretable plot
  stat4 <- stat4 %>% mutate(log_p_val = -log10(p_val))
  
  # Make a data frame that makes it so that it condenses the nsp and ORF into "Other" category.
  stat5 <- stat4
  stat5$pred <- stat5$pred * 30
  levels_new <- c(levels(stat4$gene),'Other')
  stat5$gene <- factor(stat5$gene,levels = levels_new)
  # Identify rows to change based on conditions
  rows_to_change <- with(stat5,
                         grepl("nsp|ORF", gene) & 
                           !grepl("RdRp|nsp7|nsp8|3CL-Pro|Hel|ExoN", gene))
  
  # Change those rows' gene value to "Other"
  stat5$gene[rows_to_change] <- "Other"
  
  # Update the levels to only have the relevant ones.
  # Extract the original levels
  original_levels <- levels(stat5$gene)
  
  # Get unique gene values that are present in the data
  values_present <- as.character(unique(stat5$gene))
  
  # Get the levels that are present in the data in the order of original levels
  levels_present <- intersect(original_levels, values_present)
  stat5$gene <- factor(stat5$gene,levels = levels_present)
  
  # Get colors for genes.
  colfunc <- colorRampPalette(c('red','purple','green','orange','black'))
  p_col <- colfunc(length(unique(stat5$gene)))
  # Update the label_data dataframe to include the new mutation
  label_data <- stat5 %>% 
    filter(mutation %in% c('spike_K444N_G22894U', 'spike_E340D_A22582C', 'spike_del_9_22281', 'nsp13 (Hel)_L1220F_C17125U', 'nucleocapsid_P207L_C28893U', 'spike_V987F_G24521U')) %>%
    mutate(label_text = case_when(
      mutation == 'spike_K444N_G22894U' ~ "K444N",
      mutation == 'spike_E340D_A22582C' ~ "E340D",
      mutation == 'spike_del_9_22281' ~ "22281 Del 9",
      mutation == 'nsp13 (Hel)_L1220F_C17125U' ~ "L1220F",
      mutation == 'nucleocapsid_P207L_C28893U' ~ "",
      mutation == 'spike_V987F_G24521U' ~ "V987F"
    ))
  
  # Uncomment to include all labels.
  # label_data <- stat5 %>% 
  #   filter(mutation %in% mutation) 
  # label_data$label_text <- label_data$mutation
  
  # Make the rates positive only.
  stat5$pred <- abs(stat5$pred)
  label_data$pred <- abs(label_data$pred)
  
  # Create the plot
  plot <- ggplot(stat5, aes(x = pred, y = log_p_val, color = gene)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_text_repel(data = label_data, aes(label = label_text), size = 3, show.legend = FALSE) +  # use geom_text_repel() to avoid overlap
    theme_classic() +  # Black and White theme
    labs(title = "",
         x = "Change in Proportion per Month",
         y = "-log10(P-value)",
         color = "Coding Region") +  # Set legend title for color
    scale_color_manual(values = p_col) +  # Set color palette
    theme(legend.position = "bottom",
          text = element_text(size = 9),       # General text size
          axis.title = element_text(size = 9), # Axis title size
          axis.text = element_text(size = 9),  # Axis text size
          legend.text = element_text(size = 9),        # Legend text size
          legend.title = element_text(size = 9),       # Legend title size
          plot.title = element_text(size = 9),         # Plot title size
          legend.spacing.y = unit(-0.2, 'cm')) +       # Reduce spacing between legend items
    guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
    scale_x_continuous(limits = c(0, NA))
  
  
  # Save the plot
  ggsave(paste0(dir_data,'/',prefix,'_estimated-rate-of-change-by-gene_',suffix,'.pdf'), plot = plot, width = 140/25.4, height = 3, device = "pdf")
  
  #############################################
  # Plot the unique minor variants per gene per 1000 nucleotides.
  # Count the unique mutations per gene
  mutations_per_gene <- stat3 %>% 
    group_by(gene) %>% 
    summarise(n_unique_mutations = n_distinct(mutation)) 
  # Make a bar plot
  ggplot(mutations_per_gene, aes(x = gene, y = n_unique_mutations)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Number of unique mutations per gene (Not Normalized)",
         x = "Coding Region",
         y = "Number of iSNVs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(dir_data,'/',prefix,'_number-uique-variants-per-gene-NOT-NORMALIZED_',suffix,'.pdf'), device = "pdf", width = 11.6, height = 4.9)
  # Create vector with the provided information
  gene_info <- c(
    "266, 805, nsp1",
    "806, 2719, nsp2",
    "2720, 8554, nsp3",
    "8555, 10054, nsp4",
    "10055, 10972, nsp5 (3CL-PRO)",
    "10973, 11842, nsp6",
    "11843, 12091, nsp7",
    "12092, 12685, nsp8",
    "12686, 13024, nsp9",
    "13025, 13441, nsp10",
    "13442, 16236, nsp12 (RdRp)",
    "16237, 18039, nsp13 (Hel)",
    "18040, 19620, nsp14 (ExoN)",
    "19621, 20658, nsp15 (EndoU)",
    "20659, 21552, nsp16 (2'-O-MT)",
    "21563, 25381, spike",
    "25393, 26217, ORF3a",
    "26245, 26469, envelope",
    "26523, 27188, membrane",
    "27202, 27384, ORF6",
    "27439, 27756, ORF7a",
    "27756, 27884, ORF7b",
    "27939, 28256, ORF8",
    "28274, 29530, nucleocapsid",
    "29558, 29671, ORF10"
  )
  # Convert gene_info into a dataframe
  gene_len <- as.data.frame(do.call(rbind, strsplit(gene_info, ", ")), stringsAsFactors = FALSE)
  names(gene_len) <- c("start", "end", "gene")
  gene_len$start <- as.numeric(gene_len$start)
  gene_len$end <- as.numeric(gene_len$end)
  # Calculate gene length
  gene_len$length <- gene_len$end - gene_len$start + 1
  # Count unique mutations per gene in stat3
  mutations_per_gene <- stat3 %>%
    group_by(gene) %>%
    summarise(n_unique_mutations = n_distinct(mutation))
  # Merge the two dataframes
  gene_len <- left_join(gene_len, mutations_per_gene, by = "gene")
  # Calculate the number of unique mutations per 1000 positions
  gene_len$mutations_per_kbp <- (gene_len$n_unique_mutations / gene_len$length) * 1000
  # Make sure it keeps the same order.
  gene_len$gene <- factor(gene_len$gene,levels = unique(stat3$gene))
  gene_len <- gene_len[!is.na(gene_len$gene),]
  # Plot the number of unique mutations per gene per 1000 positions
  ggplot(gene_len, aes(x = gene, y = mutations_per_kbp)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Number of iSNVs per Coding Region per 1000 Nucleotides",
         x = "Coding Region",
         y = "Number of iSNVs per 1000 nucleotides") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(dir_data,'/',prefix,'_number-uique-variants-per-gene_',suffix,'.pdf'), device = "pdf", width = 11.6, height = 4.9)
  
  #############################################
  # Linear regression for number of unique mutations per gene.
  # Set up PDF file for output, convert width from mm to inches (1 inch = 25.4 mm)
  # Calculate the scaling factor for a 12pt font size based on the default size used by R in points
  # The default point size of text in R is generally 12 points, so if you're finding the text too large with cex = 1.5,
  # it suggests the default might be larger in your case. 
  # You may need to adjust the scaling factor according to your specific environment or the output on your device.
  
  scaling_factor <- 12 / 12 # desired font size divided by the default font size, adjust as needed
  
  # Variable to toggle gene labels on the plot
  label_genes <- FALSE  # Set this to FALSE if you don't want labels
  
  # Set up PDF file for output, convert width from mm to inches (1 inch = 25.4 mm)
  pdf(paste0(dir_data, '/', prefix, '_linear-regression-unique-variants-per-gene_', suffix, '.pdf'), 
      width = 170 / 25.4, 
      height = 95 / 25.4, 
      onefile = TRUE)
  
  # Scatter plot with adjusted font size
  plot(gene_len$length, gene_len$n_unique_mutations,
       ylim = c(0, (max(gene_len$n_unique_mutations) * 1.2)),
       xlab = 'Length (nt)',
       ylab = 'iSNV Counts',
       main = 'iSNV per Coding Region',
       xaxs = "i",
       yaxs = "i",
       col = "black", # Set point color to black
       cex.lab = scaling_factor, # Scale axis labels' text size 
       cex.axis = scaling_factor, # Scale axis notation's text size
       cex.main = scaling_factor) # Scale title text size
  
  if (label_genes) {
    # Specific genes to label
    specific_genes <- c("spike", "nsp3", "nsp12 (RdRp)", "nsp13 (Hel)", "nsp2", "nucleocapsid", 
                        "nsp15 (EndoU)", "envelope", "nsp14 (ExoN)")
    specific_genes_indices <- gene_len$gene %in% specific_genes
    
    # Add text labels with adjusted font size
    text(gene_len$length[specific_genes_indices], 
         gene_len$n_unique_mutations[specific_genes_indices],
         labels = gene_len$gene[specific_genes_indices],
         pos = 3,
         cex = scaling_factor) # Scale text size
  }
  
  # Fit linear regression model
  model <- lm(n_unique_mutations ~ length, data = gene_len)
  
  # Generate predictions and confidence intervals for the predictions
  newdata <- data.frame(length = seq(min(gene_len$length), max(gene_len$length), length.out = 100))
  predictions <- predict(model, newdata = newdata, interval = "confidence")
  
  # Add linear regression line
  lines(newdata$length, predictions[, "fit"], col = "black")
  
  # Add confidence interval shading (lighter and more translucent)
  polygon(c(newdata$length, rev(newdata$length)), 
          c(predictions[, "lwr"], rev(predictions[, "upr"])),
          col = rgb(1, 0, 0, alpha = 0.2), border = NA)
  
  # Close the PDF device
  dev.off()
  
  
  
  
  # Return data frame that has the significance of the genes.
  # Fit a linear regression model
  model <- lm(n_unique_mutations ~ length, data = gene_len)
  
  # Calculate standardized residuals
  standardized_resids <- rstandard(model)
  
  # Calculate p-values based on the t-distribution
  n <- nrow(gene_len)
  p <- length(coef(model))
  df <- n - p
  
  # Two-tailed p-values
  p_values <- 2 * (1 - pt(abs(standardized_resids), df))
  
  # Create a data frame with gene names and their corresponding p-values
  results_df <- data.frame(
    gene = gene_len$gene,
    p_value = p_values
  )
  
  # Get only significant genes
  significant_genes_df <- results_df[which(results_df$p_value < 0.05), ]
  
  
  # Adjust the p-values using BH procedure for FDR control
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # Add the adjusted p-values to the results data frame
  results_df$adjusted_p_value <- adjusted_p_values
  
  # If you want to get only significant genes based on adjusted p-values
  significant_genes_df <- results_df[which(results_df$adjusted_p_value < 0.05), ]
  
  results_df
  return(results_df)
}

get_immune_plot <- function(pt1,day_earliest,max_cell){
  library(dplyr)  # For data manipulation
  library(ggplot2)  # For creating plots
  library(RColorBrewer)  # For color palettes
  
  scale_alc <- 4  # Scale factor for ALC
  
  # Extract symptom onset dates
  symptom_onset_dates <- pt1 %>%
    filter(type %in% c('symptom onset', 'initially asymptomatic')) %>%
    select(patient, symptom_onset_date = date)
  
  # Merge symptom onset dates and calculate days since onset
  pt2 <- pt1 %>%
    left_join(symptom_onset_dates, by = "patient") %>%
    mutate(
      day = as.numeric(date - symptom_onset_date),
      day_end = as.numeric(date_end - symptom_onset_date)
    ) %>%
    filter(day >= day_earliest)
  
  # Prepare ALC data
  pt_alc <- pt2 %>%
    filter(type == 'alc_manual') %>%
    select(value, day, patient) %>%
    mutate(
      value = pmin(max_cell, value) * scale_alc,
      dataset = "pt_alc"
    )
  
  # Prepare ANC data
  pt_anc <- pt2 %>%
    filter(type == 'anc_manual') %>%
    group_by(patient, pt_day = paste0(patient, '_', day)) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = 'drop') %>%
    mutate(
      day = as.numeric(sub(".*_", "", pt_day)),
      value = pmin(max_cell, value),
      dataset = "pt_anc"
    )
  
  # Combine ANC and ALC data
  pt_anc <- pt_anc[ , -which(names(pt_anc) %in% c('pt_day'))]
  combined_data <- rbind(pt_anc, pt_alc)
  combined_data <- combined_data[order(combined_data$day),]
  combined_data <- combined_data[order(combined_data$dataset, decreasing = TRUE),]
  
  # Prepare treatment data
  pt_tx <- pt2 %>%
    filter(type %in% c('treatment', 'treatment covid')) %>%
    mutate(day_end = if_else(day_end - day < 3, day_end + 3, day_end))  # Ensure visibility for short treatments
  
  # Set y positions for treatments
  pt_subj <- unique(pt_tx$patient)
  pt_tx2 <- data.frame()  # Initialize empty data frame
  for (subj in pt_subj) {
    tx_data <- filter(pt_tx, patient == subj)
    num_tx <- unique(tx_data$note)
    y_positions <- data.frame(treatment = num_tx, ymin = seq(-length(num_tx) + 1, 0), ymax = seq(-length(num_tx), -1))
    tx_data <- left_join(tx_data, y_positions, by = c("note" = "treatment"))
    pt_tx2 <- rbind(pt_tx2, tx_data)  # Append data
  }
  pt_tx <- pt_tx2
  
  # Prepare colors
  num_notes <- n_distinct(pt_tx$note)
  note_colors <- brewer.pal(min(num_notes, brewer.pal.info["Dark2", "maxcolors"]), "Dark2")
  if (num_notes > brewer.pal.info["Dark2", "maxcolors"]) {
    note_colors <- colorRampPalette(brewer.pal(brewer.pal.info["Dark2", "maxcolors"], "Dark2"))(num_notes)
  }
  names(note_colors) <- unique(pt_tx$note)
  
  # Determine the replicates per timepoint
  rep1 <- subset(met1,met1$sample_included_in_study_statistically == 'yes')
  # Creating 'rep2' from 'rep1' where it summarizes the number of replicates for each timepoint of each sample.
  rep2 <- data.frame(table(rep1$study_id,rep1$days_since_symptom_onset))
  colnames(rep2) <- c('patient','day','rep')
  rep2 <- subset(rep2, rep2$rep != 0)
  rep2 <- rep2 %>%
    mutate(
      replicate = if_else(
        rep > 4, 
        ">4", 
        as.character(rep)  # Convert 'rep' to character for other cases
      )
    )
  
  # Defining shape based on 'replicate'
  rep2 <- rep2 %>%
    mutate(
      shape = case_when(
        replicate == "2" ~ 16, # Circle
        replicate == "3" ~ 17, # Triangle
        replicate == "4" ~ 15, # Square
        replicate == ">4" ~ 18, # Diamond
        TRUE ~ NA_real_ # NA for other cases or if 'replicate' is not a number or ">4"
      )
    )
  rep2$day <- as.numeric(as.character(rep2$day))
  rep2$patient <- as.character(rep2$patient)
  combined_data$patient <- as.character(combined_data$patient)
  pt_tx$patient <- as.character(pt_tx$patient)
  
  # Calculate min and max day for each patient
  line_data <- rep2 %>%
    group_by(patient) %>%
    summarise(
      start_day = min(day),
      end_day = max(day)
    )
  
  # Create plot
  plot_timeline <- ggplot() +
    geom_area(data = combined_data, aes(x = day, y = value, fill = dataset, group = interaction(dataset, patient)), alpha = 0.4, position = position_dodge(0.8), size = 0.25) +
    geom_line(data = combined_data, aes(x = day, y = value, color = dataset, group = interaction(dataset, patient)), size = 1) +
    geom_rect(data = pt_tx, aes(xmin = day, xmax = day_end, ymin = ymin, ymax = ymax, fill = note, group = patient), alpha = 1) +
    geom_vline(xintercept = 0, color = "black") +
    geom_point(data = rep2, aes(x = day, y = 7.5, shape = factor(shape)), size = 1.2, color = "black", show.legend = FALSE) +
    scale_shape_manual(values = c("16" = 16, "17" = 17, "15" = 15, "18" = 18)) +
    geom_segment(data = line_data, aes(x = start_day, xend = end_day, y = 7.5, yend = 7.5, group = patient), color = "black", size = 0.2) +  # black line for each patient
    scale_fill_manual(values = c("pt_anc" = "#fc7703", "pt_alc" = "#7b03fc", note_colors)) +
    scale_color_manual(values = c("pt_anc" = "#fc7703", "pt_alc" = "#7b03fc", note_colors)) +
    facet_wrap(~patient, scales = "free_y", ncol = 1) +
    coord_cartesian(ylim = c(-5, 10)) +
    scale_y_continuous(breaks = seq(0, 10, by = 2.5)) +
    labs(title = "Value Over Days by Patient", x = "Days Since Onset", y = "Thousands of Cells/ul", fill = "Dataset", color = "Dataset") +
    theme_minimal()
  
  plot_timeline
  
  ggsave(filename = paste0(dir_data,'/',prefix,'_immune-timeline_',suffix,'.pdf'), plot = plot_timeline, device = "pdf", width = 200/25.4, height = 170/25.4, units = "in")
  
}