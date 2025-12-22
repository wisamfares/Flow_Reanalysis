############################################################
# helpers.R
# A collection of utility functions for analysis and plotting
# ------------------------------------------------------------
# This file is meant to be sourced inside an R Markdown file:
#     source("flow_helpers.R")
############################################################


##################################
# 1. Analysis Pipeline Functions #
##################################

# Load two ligand experiments
loadTwoLigandData = function(data_dir){
  # Loads all samples from a given experiment into a single data frame
  # Cleans up data by selecting Norm_pSMAD2 and HA columns, and removing outliers
  # Outliers removed by filtering to 0.1 - 99.9 percentiles
  # Load Data
  csv_files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  mDox_files = csv_files[grepl("minusDox", csv_files)]
  pDox_files = csv_files[!grepl("minusDox", csv_files)]
  
  # Define old colnames
  old_HA_colname = "B1-A"
  old_pSMAD2_colname = "YG1-A"
  old_SMAD2_colname = "R2-A"
  old_DAPI_colname = "UV7-A"
  
  mD_data = mDox_files %>%
    set_names(basename(mDox_files)) %>%  # optional: name each element by file name
    map_dfr(~read_csv(.x, show_col_types = FALSE), .id = "Sample_Name") %>% # read each file and combine, keep filename
    mutate(Sample_Name = "Minus Dox Control",
           Cell_Line = "R3",
           TGFB1 = 0,
           GDF11 = 0) %>%
    mutate(Dox = "-dox")
  
  
  
  all_pD_data <- pDox_files %>%
    set_names(basename(pDox_files)) %>%  # optional: name each element by file name
    map_dfr(~read_csv(.x, show_col_types = FALSE), .id = "Sample_Name") %>% # read each file and combine, keep filename
    mutate(Sample_Name = str_replace(Sample_Name, ".csv", "")) %>% # Remove .csv from Sample_Name name
    mutate(Dox = "+dox") %>%
    
    # Extract Cell Name
    separate(Sample_Name, into = c("Cell_Line", "TGFB1", "GDF11"), sep = "_", fill = "right", remove = FALSE) %>%
    
    mutate(
      TGFB1 = as.numeric(str_remove(TGFB1, "T$")),
      GDF11 = as.numeric(str_remove(GDF11, "G$"))
    ) %>%
    
    # Clean up data frame
    relocate(Sample_Name, Cell_Line, TGFB1, GDF11)
  
  all_data = rbind(mD_data, all_pD_data) %>%
    # (2) Filter to only include columns of interest, B1.A and YG1.A for 5 laser, HA and pSMAD for 3 laser
    select("Sample_Name", "Cell_Line", "Dox", "TGFB1", "GDF11", contains("FSC-A"),
           starts_with(old_HA_colname), 
           starts_with(old_pSMAD2_colname),
           starts_with(old_SMAD2_colname),
           starts_with(old_DAPI_colname)) %>%
    rename_with(~ "HA", starts_with(old_HA_colname))%>%
    rename_with(~ "pSMAD2", starts_with(old_pSMAD2_colname)) %>%
    rename_with(~ "SMAD2", starts_with(old_SMAD2_colname)) %>%
    rename_with(~ "DAPI", starts_with(old_DAPI_colname)) %>%
    
    ## Log Transform HA and Norm_pSMAD2
    mutate(
      pSMAD2 = log10(pSMAD2),
      HA = log10(HA),
      SMAD2 = log10(SMAD2)
    ) %>%
    drop_na() %>% # Remove NAs
    
    ## Filter outliers 
    group_by(Cell_Line, TGFB1, GDF11) %>%
    
    # calculate the quantiles for Norm_pSMAD2 and HA
    mutate(
      p001_pSMAD2 = quantile(pSMAD2, 0.001),
      p999_pSMAD2 = quantile(pSMAD2, 0.999),
      p001_HA = quantile(HA, 0.001),
      p999_HA = quantile(HA, 0.999)
    ) %>%
    
    # REMOVE OUTLIERS
    # Keep only rows within the 0.1-99.9 percentile
    filter(pSMAD2 >= p001_pSMAD2,
           pSMAD2 <= p999_pSMAD2,
           HA >= p001_HA,
           HA <= p999_HA) %>%
    ungroup() %>%
    select(-p001_pSMAD2, -p999_pSMAD2, -p001_HA, -p999_HA) %>% # Remove percentile columns
    mutate(Norm_pSMAD2 = pSMAD2/SMAD2)
  
  return(all_data)
}

loadOldTwoLigandData <- function(data_dir){
  csv_files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  mDox_files = csv_files[grepl("minusDox", csv_files)]
  pDox_files = csv_files[!grepl("minusDox", csv_files)]
  
  # Define old colnames
  old_HA_colname = "HA"
  old_pSMAD2_colname = "pSMAD2"
  old_SMAD2_colname = "SMAD2"
  
  mD_data = mDox_files %>%
    set_names(basename(mDox_files)) %>%  # optional: name each element by file name
    map_dfr(~read_csv(.x, show_col_types = FALSE), .id = "Sample_Name") %>% # read each file and combine, keep filename
    mutate(Sample_Name = "Minus Dox Control",
           Cell_Line = "R3",
           TGFB1 = 0,
           GDF11 = 0) %>%
    mutate(Dox = "-dox")
  
  
  
  all_pD_data <- pDox_files %>%
    set_names(basename(pDox_files)) %>%  # optional: name each element by file name
    map_dfr(~read_csv(.x, show_col_types = FALSE), .id = "Sample_Name") %>% # read each file and combine, keep filename
    mutate(Sample_Name = str_replace(Sample_Name, ".csv", "")) %>% # Remove .csv from Sample_Name name
    mutate(Dox = "+dox") %>%
    
    # Extract Cell Name
    separate(Sample_Name, into = c("Cell_Line", "TGFB1", "GDF11"), sep = "_", fill = "right", remove = FALSE) %>%
    
    mutate(
      TGFB1 = as.numeric(str_remove(TGFB1, "T$")),
      GDF11 = as.numeric(str_remove(GDF11, "G$"))
    ) %>%
    
    # Clean up data frame
    relocate(Sample_Name, Cell_Line, TGFB1, GDF11)
  
  all_data = rbind(mD_data, all_pD_data) %>%
    # (2) Filter to only include columns of interest, B1.A and YG1.A for 5 laser, HA and pSMAD for 3 laser
    select("Sample_Name", "Cell_Line", "Dox", "TGFB1", "GDF11", contains("FSC-A"),
           starts_with(old_HA_colname), 
           starts_with(old_pSMAD2_colname),
           starts_with(old_SMAD2_colname)) %>%
    rename_with(~ "HA", starts_with(old_HA_colname))%>%
    rename_with(~ "pSMAD2", starts_with(old_pSMAD2_colname)) %>%
    rename_with(~ "SMAD2", starts_with(old_SMAD2_colname)) %>%
    
    ## Log Transform HA and Norm_pSMAD2
    mutate(
      pSMAD2 = log10(pSMAD2),
      HA = log10(HA),
      SMAD2 = log10(SMAD2)
    ) %>%
    drop_na() %>% # Remove NAs
    
    ## Filter outliers 
    group_by(Cell_Line, TGFB1, GDF11) %>%
    
    # calculate the quantiles for Norm_pSMAD2 and HA
    mutate(
      p001_pSMAD2 = quantile(pSMAD2, 0.001),
      p999_pSMAD2 = quantile(pSMAD2, 0.999),
      p001_HA = quantile(HA, 0.001),
      p999_HA = quantile(HA, 0.999)
    ) %>%
  
  # REMOVE OUTLIERS
  # Keep only rows within the 0.1-99.9 percentile
  filter(pSMAD2 >= p001_pSMAD2,
         pSMAD2 <= p999_pSMAD2,
         HA >= p001_HA,
         HA <= p999_HA) %>%
    ungroup() %>%
    select(-p001_pSMAD2, -p999_pSMAD2, -p001_HA, -p999_HA) %>% # Remove percentile columns
    mutate(Norm_pSMAD2 = pSMAD2/SMAD2)
  
  return(all_data)
  }










# Load all experiment data from a direcetory and its subdirectories
loadAllTwoLigandData <- function(data_dir){
  
  # Step 1: get all experiment folders (two levels deep)
  experiment_folders <- dir_ls(data_dir, recurse = 2, type = "directory")
  
  # Step 2: filter folders that contain CSVs
  experiment_folders <- experiment_folders[map_lgl(experiment_folders, ~ length(dir_ls(.x, glob = "*.csv")) > 0)]
  
  # Step 3: map loadData over all experiment folders
  all_data <- experiment_folders %>%
    set_names(experiment_folders) %>%
    map_dfr(~{
      path <- .x
      path_parts <- path_split(path)[[1]]  # split path into components
      
      experiment <- path_parts[length(path_parts)]
      
      
      loadTwoLigandData(path) %>%
        mutate(
          Experiment = experiment
        )
    })
  
  return(all_data)
}


loadAllOldTwoLigandData <- function(data_dir){
  
  # Step 1: get all experiment folders (two levels deep)
  experiment_folders <- dir_ls(data_dir, recurse = 2, type = "directory")
  
  # Step 2: filter folders that contain CSVs
  experiment_folders <- experiment_folders[map_lgl(experiment_folders, ~ length(dir_ls(.x, glob = "*.csv")) > 0)]
  
  # Step 3: map loadData over all experiment folders
  all_data <- experiment_folders %>%
    set_names(experiment_folders) %>%
    map_dfr(~{
      path <- .x
      path_parts <- path_split(path)[[1]]  # split path into components
      
      experiment <- path_parts[length(path_parts)]
      
      
      loadOldTwoLigandData(path) %>%
        mutate(
          Experiment = experiment
        )
    })
  
  return(all_data)
}



filterData = function(all_data, 
                      Norm_pSMAD2_perc_cutoff = 0,
                      HA_perc_cutoff = 0.95){
  # Uses perc_cutoffs to determine thresholds from Luc - No Cyto sample
  Norm_pSMAD2_thresh = all_data %>%
    filter(Dox == "-dox") %>%
    summarize(p_Norm_pSMAD2 = quantile(Norm_pSMAD2, Norm_pSMAD2_perc_cutoff)) %>%
    pull(p_Norm_pSMAD2)
  
  # Calculate threshold based on -dox, -cyto condition
  HA_thresh = all_data %>%
    filter(Dox == "-dox") %>%
    summarize(p_HA = quantile(HA, HA_perc_cutoff)) %>%
    pull(p_HA)
  
  # Use thresholds to filter data
  filtered_data = all_data %>%
    filter(Cell_Line == "R3") %>%
    filter(Norm_pSMAD2 > !!Norm_pSMAD2_thresh) %>%
    filter(HA > !!HA_thresh)
  
  return(filtered_data)
}

binHA = function(filtered_data, n_bins){
  # Bin HA data based off n_bins using previous HA_perc cutoff (new min) and max HA
  
  # Calculate min and max HA 
  filt_data_summ = filtered_data %>%
    summarize(min_HA = min(HA),
              max_HA = max(HA),
              # Adding 90p to prevent bins with low numbers of cells
              p_HA = quantile(HA, 0.95)) 
  
  # Pull min and max from summary
  min_HA = pull(filt_data_summ, min_HA)
  max_HA = pull(filt_data_summ, max_HA)
  p_HA = pull(filt_data_summ, p_HA)
  
  
  # Create bin edges 
  bins = seq(min_HA-1e-1, max_HA+1e-1, length.out = n_bins + 1) # Log spaced bin edges
  
  # Use bins to cut data 
  binned_data = filtered_data %>%
    mutate(HA_bin = cut(HA, breaks = bins, include.lowest = TRUE),
           binID = as.integer(HA_bin))
  
  return(binned_data)
}

binHA_deciles = function(filtered_data, n_bins = 10) {
  
  binned_data = filtered_data %>%
    mutate(
      binID = ntile(HA, n_bins),
      HA_bin = factor(binID)
    ) %>%
    ungroup()
  
  return(binned_data)
}


bootstrapBinnedData = function(binned_data, n_boot){
  boot_summary <- function(binned_data, n_boot) {
    # Bootstrap each bin
    boot_results <- replicate(n_boot, {
      sample_data <- binned_data[sample(nrow(binned_data), replace = TRUE), ]
      c(
        mean_Norm_pSMAD2 = mean(sample_data$Norm_pSMAD2),
        median_Norm_pSMAD2 = median(sample_data$Norm_pSMAD2)
      )
    }, simplify = FALSE)
    
    # Convert to data frame
    boot_df <- bind_rows(boot_results)
    
    # Compute mean and 95% CI for each stat
    tibble(
      mean_Norm_pSMAD2 = mean(boot_df$mean_Norm_pSMAD2, na.rm = TRUE),
      mean_Norm_pSMAD2_low = quantile(boot_df$mean_Norm_pSMAD2, 0.025, na.rm = TRUE),
      mean_Norm_pSMAD2_high = quantile(boot_df$mean_Norm_pSMAD2, 0.975, na.rm = TRUE),
      median_Norm_pSMAD2 = mean(boot_df$median_Norm_pSMAD2, na.rm = TRUE),
      median_Norm_pSMAD2_low = quantile(boot_df$median_Norm_pSMAD2, 0.025, na.rm = TRUE),
      median_Norm_pSMAD2_high = quantile(boot_df$median_Norm_pSMAD2, 0.975, na.rm = TRUE),
      n_cells = nrow(binned_data)
    )
  }
  
  binned_boot <- binned_data %>%
    filter(!is.na(HA_bin)) %>%
    group_by(Cell_Line, Dox, TGFB1, GDF11, HA_bin) %>%
    group_modify(~ boot_summary(.x, n_boot = 1000)) %>%
    ungroup()
  
  return(binned_boot)
}



#############################
# 2. Full Pipeline Analysis #
#############################

run_bootstrap_analysis_pipeline = function(all_exp_data, 
                                           pSMAD2_perc_cutoff = 0,
                                           HA_perc_cutoff = 0.9,
                                           n_bins = 10,
                                           n_boot = 1000){
  
  
  
  # (2) Filter Data based on pSMAD2 and HA threshold of -dox R3 
  filtered_data = all_exp_data %>%
    group_by(Experiment) %>%
    ## Apply filterData to each group
    group_modify(~ filterData(.x, # All samples for a given ligand and experiment
                              pSMAD2_perc_cutoff, 
                              HA_perc_cutoff)) %>%
    ungroup() %>%
    filter(Dox == "+dox")
  
  # (3) Bin Data
  binned_data = filtered_data %>%
    group_by(Experiment, TGFB1, GDF11) %>%
    group_modify(~ binHA(.x, n_bins)) %>%
    ungroup()
  
  
  # (4) Bootstrap Binned Data
  binned_boot = binned_data %>%
    group_by(Experiment) %>%
    group_modify(~bootstrapBinnedData(.x, n_boot)) %>%
    ungroup()%>% 
    group_by(Experiment, TGFB1, GDF11) %>%
    arrange(HA_bin) %>%
    mutate(binID = row_number()) %>%
    ungroup()
  
  results = list(all_exp_data = all_exp_data,
                 filtered_data = filtered_data,
                 binned_data = binned_data,
                 binned_boot = binned_boot)
}


###################
# TEST DECILES #
################

run_bootstrap_deciles_analysis_pipeline = function(all_exp_data, 
                                           pSMAD2_perc_cutoff = 0,
                                           HA_perc_cutoff = 0.9,
                                           n_bins = 10,
                                           n_boot = 1000){
  
  
  
  # (2) Filter Data based on pSMAD2 and HA threshold of Luc - No Cytokine 
  filtered_data = all_exp_data %>%
    group_by(Experiment) %>%
    ## Apply filterData to each group
    group_modify(~ filterData(.x, # All samples for a given ligand and experiment
                              pSMAD2_perc_cutoff, 
                              HA_perc_cutoff)) %>%
    ungroup()
  
  # (3) Bin Data
  binned_data = filtered_data %>%
    group_by(TGFB1, GDF11, Experiment, Replicate) %>%
    group_modify(~ binHA_deciles(.x, n_bins)) %>%
    ungroup()
  
  
  # (4) Bootstrap Binned Data
  binned_boot = binned_data %>%
    group_by(Experiment, Replicate) %>%
    group_modify(~bootstrapBinnedData(.x, n_boot)) %>%
    ungroup()%>% 
    group_by(Dox, TGFB1, GDF11, Experiment, Replicate) %>%
    arrange(HA_bin) %>%
    mutate(binID = row_number()) %>%
    ungroup()
  
  results = list(all_exp_data = all_exp_data,
                 filtered_data = filtered_data,
                 binned_data = binned_data,
                 binned_boot = binned_boot)
}



#################################
# 2. Summary and Utility Helpers #
#################################


###########################
# 3. ggplot Helper Themes #
###########################


############################################
# 4. Helper Functions for Group Operations #
############################################


#######################################
# 5. Miscellaneous Small Convenience  #
#######################################



############################################################
# End of helpers.R
############################################################
