# ------------------------------------------------------------------------------
# Script:        HelperFunctions.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Helper functions for the Interactive Analysis Suite 
# ------------------------------------------------------------------------------

################################################################################
# Helper function for reading numeric inputs
################################################################################
read_choice <- function(prompt_message) {
  choice <- suppressWarnings(as.numeric(trimws(readline(prompt = prompt_message))))
  return(choice)
}

################################################################################
# Source Selection and Variable Assignment
################################################################################

setup_variables <- function() {
  # Define allowed options
  regions <- c(unique(sinfo$region), "all")
  cancer_status <- c("cancer", "control", "all")
  panels <- c(unique(binfo$panel), "all")
  
  # Ask user for the data source selection
  cat("\nPlease select a source (1 = Region, 2 = Selected, 0 = Exit):\n")
  source_choice <- read_choice("")
  
  if(is.na(source_choice) || source_choice == 0) {
    cat("\nExiting script.\n")
    return(NULL)
  }
  
  if (source_choice == 1) {
    cat("\nPlease select a region (",
        paste0(seq_along(regions), " = ", regions, collapse = ", "),
        ", 0 = Exit)\n")
    region_choice <- read_choice("")
    if (is.na(region_choice) || region_choice == 0) {
      cat("No region selected. Exiting script.\n")
      return(NULL)
    }
    if (region_choice %in% seq_along(regions)) {
      selected_region <- if (regions[region_choice] != "all")
        regions[region_choice] else "all"
      cat("\nPlease select cancer status (",
          paste0(seq_along(cancer_status), " = ", 
                 cancer_status, collapse = ", "),
          ", 0 = Exit)\n")
      cancer_choice <- read_choice("")
      if (is.na(cancer_choice) || cancer_choice == 0) {
        cat("\nNo cancer status selected. Exiting script. \n")
        return(NULL)
      }
      if (cancer_choice %in% seq_along(cancer_status)) {
        selected_cancer <- if (cancer_status[cancer_choice] != "all")
          cancer_status[cancer_status] else "all"
        cat("\nPlease select a panel(",
        paste0(seq_along(panels), " = ", panels, collapse = ", "),
        ", 0 = Exit)\n")
        panel_choice <- read_choice("")
        if (is.na(panel_choice) || panel_choice == 0) {
          cat("\nNo panel selected. Exiting script.\n")
          return(NULL)
        }
        if (panel_choice %in% seq_along(panels)) {
          selected_panel <- if (panels[panel_choice] != "all")
            panels[panel_choice] else "all"
          
          # Construct variable name prefix
          variable_name_prefix <- paste0(
            selected_region,
            if (selected_cancer != "all") paste0(".", selected_cancer) else "",
            if (selected_panel != "all") paste0(".", selected_panel) else ""
          )
          ptx_name <- paste0(variable_name_prefix, ".ptx")
          sinfo_name <- paste0(variable_name_prefix, ".sinfo")
          binfo_name <- paste0(variable_name_prefix, ".binfo")
          
          # Retrieve objects from Global Environment
          if (!exists(ptx_name, envir = .GlobalEnv)) {
            cat("\nVariable not found for ptx:", ptx_name, "\n")
            return(NULL)
          }
          select.ptx_obj <- get(ptx_name, envir = .GlobalEnv)
          
          if (!exists(sinfo_name, envir = .GlobalEnv)) {
            cat("\nVariable not found for sinfo:", sinfo_name, "\n")
            return(NULL)
          }
          select.sinfo_obj <- get(sinfo_name, envir = .GlobalEnv)
          
          if (exists(binfo_name, envir = .GlobalEnv)) {
            select.binfo_obj <- get(binfo_name, envir = .GlobalEnv)
          } else { 
            cat("\nVariable not found for binfo. 
                Attempting to assign defaultbinfo...\n")
            if (exists("binfo", envir = .GlobalEnv)) {
              select.binfo_obj <- get("binfo", envir = .GlobalEnv)
          } else {
            cat("\nDefault binfo not found. Exiting script.\n")
            return(NULL)
            }
          }
          
          # Bind identified variables into new global variables
          assign("select.ptx", select.ptx_obj, envir = .GlobalEnv)
          assign("select.sinfo", select.sinfo_obj, envir = .GlobalEnv)
          assign("select.binfo", select.binfo_obj, envir = .GlobalEnv)
          cat("\nVariables assigned successfully:\n- select.ptx\n- select.sinfo\n- select.binfo\n")
          return(list(ptx = select.ptx_obj,
                      sinfo = select.sinfo_obj,
                      binfo = select.binfo_obj))
        }
      }
    }
  } else if (source_choice == 2) {
    cat("\nYou selected 'Selected'. Binding existing variables.\n")
    if (exists("select.ptx", envir = .GlobalEnv) &&
        exists("select.sinfo", envir = .GlobalEnv) &&
        exists("select.binfo", envir = .GlobalEnv)) {
      select.ptx_obj <- get("select.ptx", envir = .GlobalEnv)
      select.sinfo_obj <- get("select.sinfo", envir = .GlobalEnv)
      select.binfo_obj <- get("select.binfo", envir = .GlobalEnv)
      
      assign("select.ptx", select.ptx_obj, envir = .GlobalEnv)
      assign("select.sinfo", select.sinfo_obj, envir = .GlobalEnv)
      assign("select.binfo", select.binfo_obj, envir = .GlobalEnv)
      
      cat("\nVariables 'select.ptx', 'select.sinfo', and 'select.binfo' 
          assigned successfully.\n")
      return(list(ptx = select.ptx_obj,
                  sinfo = select.sinfo_obj,
                  binfo = select.binfo_obj))
      
    } else {
      cat("\nOne or more required variables not found. Exiting script.\n")
      return(NULL)
    }
  } else {
    cat("\nInvalid source choice. Exiting script.\n")
  }
}
################################################################################
# Trait selection function
################################################################################
choose_traits <- function(df) {
  # List columns with numbers
  cat("Available columns:\n")
  for (i in seq_along(names(df))) {
    cat(sprintf("[%d] %s\n", i, names(df)[i]))
  }
  
  # Ask user for input
  input <- readline(prompt = "Enter the numbers of the traits you want (comma or space separated): ")
  
  # Split input into numbers
  nums <- unlist(strsplit(input, "[, ]+"))
  nums <- as.integer(nums[nums != ""])
  
  # Validate
  if (any(is.na(nums)) || any(nums < 1 | nums > ncol(df))) {
    stop("Invalid selection. Please enter valid column numbers.")
  }
  
  # Return selected column names
  selected <- names(df)[nums]
  cat("You selected:\n")
  print(selected)
  
  return(selected)
}
################################################################################
# select non-continuous columns
################################################################################
list_noncontinous_columns <- function(df, max_unique = 10) {
  columns <- names(df)[sapply(df, function(x) {
    is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= max_unique)
  })]
  # List non-continous columns with numners
  paste(" ",seq_along(columns), "-", columns)
}

list_continous_columns <- function(df, min_unique = 11) {
  continuous <- names(df)[sapply(df, function(x) {
    is.numeric(x) && length(unique(x) >= min_unique)
  })]
  if (length(continuous) == 0) {
    cat("No comlums containing continuous values were found")
  } else {
    cat("Continuous columns:\n")
    for (i in seq_along(continuous)) {
      cat(i, "-", continuous[i], "\n")
    }
  }
  invisible(continuous)
}

################################################################################
# Function for selecting non-continuous columns
################################################################################
choose_noncontinuous_columns <- function(df, max_unique = 10) {
  nc_cols <- names(df)[sapply(df, function(x) {
    is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= max_unique)
  })]
  if (length(nc_cols) == 0){
    cat("No non-continous columns found. \n")
    return(character(0))
  }
  # List non-continous columns with numners
  cat("Available non-continuous columns:\n")
  for (i in seq_along(nc_cols)) {
    cat(sprintf("[%d] %s\n", i, nc_cols[i]))
  }
  # Ask for selection
  input <- readline(
    prompt = "Enter the numbers corresponding to the columns of choice:\n"
  )
  nums <- unlist(strsplit(input, "[, ]+"))
  nums <- as.integer(nums[nums != ""])
  
  # Validate choices
  if (any(is.na(nums)) || any(nums < 1 | nums > ncol(df))) {
    stop("Invalid selection. Please enter valid column numbers.")
  }
  selected <- nc_cols[nums]
  cat("You have selected:\n")
  return(selected)
  
}

################################################################################
# Function for selection of Trait and Level:
################################################################################
choose_filter <- function(df) {
  # List columns
  cat("Available columns:\n")
  for (i in seq_along(names(df))) {
    cat(sprintf("[%d] %s\n", i, names(df)[i]))
  }
  
  col_num <- as.integer(readline(prompt = "Enter the number of the trait of choice: "))
  if (is.na(col_num) || col_num < 1 || col_num > ncol(df)) {
    stop("Invalid selection for trait.")
  }
  filter_trait <- names(df)[col_num]
  
  # Coerce to character before listing
  unique_vals <- unique(as.character(df[[filter_trait]]))
  # Order display
  unique_vals <- unique_vals[order(tolower(unique_vals), na.last = TRUE)]
  
  cat(sprintf("\nUnique values in '%s':\n", filter_trait))
  for (i in seq_along(unique_vals)) {
    cat(sprintf("[%d] %s\n", i, unique_vals[i]))
  }
  # Select level of interest and use to filter
  val_num <- as.integer(readline(prompt = "Enter the number of desired level: "))
  if (is.na(val_num) || val_num < 1 || val_num > length(unique_vals)) {
    stop("Invalid selection for level.")
  }
  filter_level <- unique_vals[val_num]
  
  cat("\nYou selected:\n")
  cat("Trait:", filter_trait, "\n")
  cat("Level:", filter_level, "\n")
  
  return(list(filter_trait = filter_trait, filter_level = filter_level))
}

################################################################################
# Helper function for selection of specific cluster:
################################################################################
choose_cluster_filter <- function(df) {
  # Ask if to filter by cluster
  choice <- as.integer(readline(prompt = "Filter by cluster? (1 = yes, 2 = no): "))
  
  if (is.na(choice) || !(choice %in% c(1, 2))) {
    stop("Invalid choice. Please enter 1 or 2.")
  }
  
  if (choice == 2) {
    # No filtering
    return(list(filter_cluster = FALSE,
                cluster_col = NULL,
                cluster_value = NULL))
  }
  
  # Yes: filter by cluster
  filter_cluster <- TRUE
  
  # Find columns containing "cluster" (case-insensitive)
  cluster_cols <- grep("cluster", names(df), ignore.case = TRUE, value = TRUE)
  
  if (length(cluster_cols) == 0) {
    stop("No columns containing 'cluster' found in the data.")
  }
  
  cat("Cluster columns:\n")
  for (i in seq_along(cluster_cols)) {
    cat(sprintf("[%d] %s\n", i, cluster_cols[i]))
  }
  
  col_num <- as.integer(readline(prompt = "Enter the number of the cluster column: "))
  
  if (is.na(col_num) || col_num < 1 || col_num > length(cluster_cols)) {
    stop("Invalid selection for cluster column.")
  }
  
  cluster_col <- cluster_cols[col_num]
  
  # Show unique cluster values, sorted numerically by the number at the end
  unique_vals <- unique(df[[cluster_col]])
  
  # Extract numeric part from the end of each value
  num_part <- as.numeric(sub(".*_(\\d+)$", "\\1", unique_vals))
  
  # Order by numeric part
  ord <- order(num_part, na.last = TRUE)
  unique_vals <- unique_vals[ord]
  
  cat(sprintf("\nUnique values in '%s' (sorted numerically):\n", cluster_col))
  for (i in seq_along(unique_vals)) {
    cat(sprintf("[%d] %s\n", i, unique_vals[i]))
  }
  
  val_num <- as.integer(readline(prompt = "Enter the number of the cluster value: "))
  
  if (is.na(val_num) || val_num < 1 || val_num > length(unique_vals)) {
    stop("Invalid selection for cluster value.")
  }
  
  cluster_value <- unique_vals[val_num]
  
  cat("\nYou selected:\n")
  cat("Cluster column:", cluster_col, "\n")
  cat("Cluster value:", cluster_value, "\n")
  
  return(list(filter_cluster = filter_cluster,
              cluster_col = cluster_col,
              cluster_value = cluster_value))
}

################################################################################
# Global seed selector
################################################################################
set_global_seed <- function() {
  # Prompt for seed input
  seed_input <- readline(prompt = "\nPlease select a global seed for reproducibility: ")
  # Remove extra whitespace and convert to numeric
  seed_input <- trimws(seed_input)
  seed_numeric <- suppressWarnings(as.numeric(seed_input))
  
  # Check if the provided input is valid; if not, set default seed
  if (is.na(seed_numeric)) {
    cat("Invalid input. Using default global seed: 123\n")
    seed_numeric <- 123
  }
  
  # Store the seed in a global variable and set the seed for R's random number generator
  global_seed <<- seed_numeric
  set.seed(global_seed)
  
  cat("Global seed is set to:", global_seed, "\n")
}

################################################################################
# Clustering method selection function
###############################################################################
choose_clustering_method <- function() {
  cat("Choose a clustering method:\n")
  cat("1 = kmeans\n")
  cat("2 = HT-Kmeans\n")
  cat("3 = Fuzzy kmeans\n")
  cat("4 = DBSCAN\n")
  cat("5 = Hierarchical DBSCAN\n")
  cat("6 = GMM\n")
  cat("7 = Mclust (GMM)\n")
  cat("8 = kmeanspp\n")
  cat("9 = Spectral\n")
  cat("10 = Agglomerative Hierarchical\n") 
  
  choice <- as.integer(readline(prompt = "Enter your choice [1-10]: "))
  if (is.na(choice) || choice < 1 || choice > 9) {
    cat("Invalid response. Please enter a number between 1 and 10.\n")
    return(NULL)
  }
  methods <- c("kmeans", "HT-Kmeans", "Fuzzy kmeans",
               "DBSCAN", "Hierarchical DBSCAN", "GMM", "Mclust", "kmeanspp", 
               "Spectral", "Agglomerative Hierarchical" )
  chosen <- methods[choice]
  cat("You selected:", chosen, "\n\n")
  return(chosen)
}

################################################################################
# Dimensionality Reduction selector
################################################################################
choose_dr_method <- function() {
  cat("Choose a dimensionality reduction method:\n")
  cat("1. PCA\n")
  cat("2. Sparse PCA\n")
  cat("3: PLS-DA\n")
  cat("4: Sparse PLS-DA\n")
  cat("5: Kernel PCA\n")
  cat("6: UMAP\n")
  cat("7: t-SNE\n")
  
  
  choice <- as.integer(readline(prompt = "Enter your choice (1-7): "))
  if (is.na(choice) || choice < 1 || choice > 9) {
    cat("Invalid response. Please enter a number between 1 and 7.\n")
    return(NULL)
  }
  dr_methods <- c("PCA", "Sparse PCA", "PLS-DA", "Sparse PLS-DA", 
                  "Kernel PCA", "UMAP", "t-SNE")
  chosen <- dr_methods[choice]
  cat("You selected:", chosen, "\n\n")
  return(chosen)
}


################################################################################
# Custom readline wrapper
################################################################################

# Create a global variable for user inputs
user_inputs <<- character(0)

# Custom readline function
my_readline <- function(prompt = ""){
  input <- base::readline(prompt)
  timestamp <- Sys.time()
  
  # Capture context of calling function
  caller_context <- tryCatch({
    deparse(sys.call(sys.parent()))
  }, error = function(e) {
    "Global" # If no parent is found, the label context will be Global
  })
  
  # Create log_entry including context
  log_entry <- paste("Time:", as.character(timestamp),
                     "| Context:", caller_context,
                     "| Prompt:", prompt,
                     "| Input:", input)
  # Append log_entry to the global 'user_inputs' variable.
  user_inputs <<- c(user_inputs, log_entry)
  
  return(input)
}
# Override readline with custom readline function:
readline <- my_readline


################################################################################
# Function for user input extraction
################################################################################

extract_user_inputs <- function() {
  if (!exists("user_inputs")){
    stop("Global variable 'user_inputs' does not exist.")
  }
  n <- length(user_inputs)
  timestamps <- character(n)
  contexts <- character(n)
  prompts <- character(n)
  inputs <- character(n)
  
  # Loop over log entries and parse the parts
  for (i in seq_along(user_inputs)) {
    parts <- strsplit(user_inputs[i], " \\| ")[[1]]
    if (length(parts) >= 4) {
      timestamps[i] <- sub("^Time: ", "", parts[1])
      contexts[i]   <- sub("^Context: ", "", parts[2])
      prompts[i]    <- sub("^Prompt: ", "", parts[3])
      inputs[i]     <- sub("^Input: ", "", parts[4])
    } else{
      timestamps[i] <- NA
      contexts[i]   <- NA
      prompts[i]    <- NA
      inputs[i]     <- NA
    }
  }
  # Combine parsed results into a dataframe
  df_inputs <- data.frame(
    Timestamp = timestamps,
    Context   = contexts,
    Prompt    = prompts,
    Input     = inputs,
    stringAsFactors = FALSE
  )
  
  return(df_inputs)
}


