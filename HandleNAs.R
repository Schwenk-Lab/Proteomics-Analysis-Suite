# ------------------------------------------------------------------------------
# Script:        HandleNAs.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for handling of missing values (NAs)
# ------------------------------------------------------------------------------

################################################################################
# NAs / Imputation
################################################################################

################################################################################ 
# Impute NAs with median
################################################################################
impute_na_median <- function(x) {
  if (!is.numeric(x)) {
    warning("Median imputation is only applicable to numeric data. 
            Returning original vector.")
    return(x)
  }
  med_val <- median(x, na.rm = TRUE)
  x[is.na(x)] <- med_val
  return(x)
}

################################################################################ 
# Handle NAs
################################################################################
handle_NAs <- function(npx, sinfo = NULL, choices = NULL, seed = global_seed,
                       NA_handling_log <- NULL) {
  interactive <- is.null(choices)
  if (is.null("NA_handling_log")) {
    NA_handling_log <- data.frame(
      data_choice = integer(),
      option      = integer(),
      columns     = I(list()),
      stringsAsFactors = FALSE
    )
  }
  # Prompt for which data to process:
  if (interactive) {
    cat("\nWhich data would you like to process?\n")
    cat(" 1 = Abundance data (npx)\n")
    cat(" 2 = Clinical data (sinfo)\n")
    data_choice <- as.numeric(readline(prompt = "Select an option: "))
  } else {
    data_choice <- choices$data_choice
  }
  
  # For abundance Data
  if (!is.na(data_choice) && data_choice == 1) {
    if (interactive) {
      cat("\nCurrent number of NAs in Abundance data:", sum(is.na(npx)), "\n")
      cat("\nNA Handling Options:\n")
      cat(" 1 = minProb\n 2 = RandomForest\n 3 = kNN\n 4 = SVD\n 5 = minDet\n
          6 = complete.cases\n 7 = median\n 0 = Cancel\n")
      option <- as.numeric(readline(prompt = "Select an option: "))
    } else {
      option <- choices$option
    }
    if (is.na(option) || option == 0) 
      return(list(npx = npx, sinfo = sinfo))
    if (option == 1) { 
      cat("Imputing missing values in Abundance data using minProb method...")
      result.npx <- t(impute_na(t(npx), seed = seed))
      } else if (option == 2) {
        if (interactive) {
          maxiter <- as.numeric(readline("Enter maxiter (default = 5): "))
          ntree <- as.numeric(readline("Enter ntree (default = 100): "))
          if (is.na(maxiter)) maxiter <- 5
          if (is.na(ntree)) ntree <- 100
        } else{
          maxiter <- if (!is.null(choicees$maxiter)) choices$maxiter else 5
          ntree <- if (!is.null(choices$ntree)) choices$ntree else 100
        }
        cat("Imputing missing values in Abundance data using RF method...")
        result.npx <- t(impute_na(t(npx), method = "RF", maxiter = maxiter,
                                  ntree = ntree, seed = seed))
      } else if (option == 3) { 
        cat("Imputing missing values in Abundance data using kNN method...")
        result.npx <- t(impute_na(t(npx), method = "kNN", seed = seed))
        
      } else if (option == 4) {
        if (interactive) {
          n_pcs <- as.numeric(readline("Enter number of PCs (default = 3): "))
          if (is.na(n_pcs)) n_pcs <- 3
        } else {
          n_pcs <- if (!is.null(choices$n_pcs)) choices$n_pcs else 3
        }
        cat("Imputing missing values in Abundance data using SVD method...")
        result.npx <- t(impute_na(t(npx), method = "SVD",
                                  n_pcs = n_pcs, seed = seed))
        
      } else if (option == 5) {
        if (interactive) {
          q <- as.numeric(readline("Enter q (default = 0.001): "))
          if (is.na(q)) q <- 0.001
        } else {
          q <- if (!is.null(choices$q)) choices$q else 0.001
        }
        cat("Imputing missing values in Abundance data using minDet method...")
        result.npx <- t(impute_na(t(npx), method = "minDet", 
                                  q = q, seed = seed))
        
      } else if (option == 6) {
        cat("Filtering out rows with NA values using complete.cases...")
        result.npx <- npx[complete.cases(npx), ]
        cat("After filtering, number of samples remainig:", 
            nrow(result.npx),"\n")
        if (!is.null(sinfo)) {
          sinfo <- sinfo[rownames(sinfo) %in% rownames(result.npx), ,
                         drop = FALSE]
        }
        
      } else if (option == 7) {
        cat("Imputing missing values using median imputation...")
        imputed_matrix <- apply(as.matrix(npx), 2, impute_na_median)
        result.npx <- as.data.frame(imputed_matrix, stringsAsFactors = FALSE)
        rownames(result.npx) <- rownames(npx)
        
      } else {
        return(list(npx =npx, sinfo = sinfo))
      }
    # Log choices from interactive
    if(interactive) {
      new_entry <- data.frame(
        data_choice = data_choice,
        option = option,
        columns = I(list(NA)),
        maxiter = if (exists("maxiter")) maxiter else NA,
        ntree = if (exists("ntree")) ntree else NA,
        n_pcs = if (exists("n_pcs")) n_pcs else NA,
        q = if (exists("q")) q else NA,
        stringsAsFactors = FALSE
      )
      NA_handling_log <<- rbind(NA_handling_log, new_entry)
      
      # Prompt if to save used parameters as a named object
      logName <- readline("Enter variable name for saved parameters
                          (or press Enter to skip): ")
      if (nzchar(logName)) {
        assign(logName, new_entry, envir = .GlobalEnv)
      }
    }
    return(list(npx = result.npx, sinfo = sinfo))
  }
  
  # For Clinical Data
  else if (!is.na(data_choice) && data_choice == 2) {
    if (is.null(sinfo)) return(list(npx = npx, sinfo = sinfo))
  
    if (interactive) {
      cat("\nAvailable columns:\n")
      for (i in seq_along(names(sinfo))) cat(i, "-", names(sinfo)[i], "\n")
      columns_input <- readline(prompt = "Enter column numbers (comma-separated):
                                ")
      selected_indices <- as.integer(unlist(strsplit(columns_input, ",")))
      } else {
        selected_indices <- choices$columns
        }
    if  (interactive) {
      cat("\nNA Handling Options:\n")
      cat(" 1 - minProb\n 2 - RandomForest\n 3 - kNN\n 4 - SVD\n 5 - minDet\n
            6 - complete.cases\n 7 - median\n 0 - Cancel\n")
      option <- as.numeric(readline(prompt = "Select an option: "))
    } else {
      option <- choices$option
    }
  
    if (is.na(option) || option == 0) return(list(npx = npx, sinfo = sinfo))
  
    # Apply chosen method
    if (option == 1) {
      cat("Imputing missing values in Clinical data using minProb method...")
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, 
                                                      drop = FALSE])), 
                                    seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 2) {
      if (interactive) {
        maxiter <- as.numeric(readline("Enter maxiter (default = 5): "))
        ntree <- as.numeric(readline("Enter ntree (default = 100): "))
        if (is.na(maxiter)) maxiter <- 5
        if (is.na(ntree)) ntree <- 100
      } else{
        maxiter <- if (!is.null(choicees$maxiter)) choices$maxiter else 5
        ntree <- if (!is.null(choices$ntree)) choices$ntree else 100
      }
      cat("Imputing missing values in Clinical data using RF method...")
      sinfo[, selected_indices] <- impute_na(sinfo[, selected_indices,
                                                   drop = FALSE],
                                             method = "RF", maxiter = maxiter,
                                             ntree = ntree, seed = seed)
    } else if (option == 3) {
      cat("Imputing missing values in Clinical data using kNN method...")
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, 
                                                      drop = FALSE])), 
                                    method = "kNN", seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    
    } else if (option == 4) {
      if (interactive) {
        n_pcs <- as.numeric(readline("Enter number of PCs (default = 3): "))
        if (is.na(n_pcs)) n_pcs <- 3
      } else {
        n_pcs <- if (!is.null(choices$n_pcs)) choices$n_pcs else 3
      }
      cat("Imputing missing values in Clinical data using SVD method...")
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, 
                                                      drop = FALSE])), 
                                    method = "SVD", n_pcs = n_pcs, 
                                    seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 5) {
        if (interactive) {
          q <- as.numeric(readline("Enter q (default = 0.001): "))
          if (is.na(q)) q <- 0.001
        } else {
          q <- if (!is.null(choices$q)) choices$q else 0.001
        }
      cat("Imputing missing values in Clinical data using minDet method...")
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, 
                                                      drop = FALSE])), 
                                    method = "minDet", q = q, 
                                    seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    
    } else if (option == 6) {
      cat("Filtering out rows with NA values using complete.cases...")
      sinfo <- sinfo[complete.cases(sinfo[, selected_indices]), , drop = FALSE]
      common_samples <- intersect(rownames(sinfo), rownames(npx))
      npx <- npx[common_samples, , drop = FALSE]
      sinfo <- sinfo[common_samples, , drop = FALSE]
      cat("After filtering, number of samples remainig:", 
          nrow(npx),"\n")
  
    } else if (option == 7) {
      cat("Imputing missing values using median imputation...")
      imputed_matrix <- apply(as.matrix(sinfo[, selected_indices, 
                                              drop = FALSE])
                              , 2, impute_na_median)
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
   }
  # Log choices from the interactive mode
    if(interactive) {
      new_entry <- data.frame(
        data_choice = data_choice,
        option = option,
        columns = I(list(NA)),
        maxiter = if (exists("maxiter")) maxiter else NA,
        ntree = if (exists("ntree")) ntree else NA,
        n_pcs = if (exists("n_pcs")) n_pcs else NA,
        q = if (exists("q")) q else NA,
        stringsAsFactors = FALSE
    )
      NA_handling_log <<- rbind(NA_handling_log, new_entry)
  
      # Prompt if to save used parameters as a named object
      logName <- readline("Enter variable name for saved parameters
                          (or press Enter to skip): ")
      if (nzchar(logName)) {
        assign(logName, new_entry, envir = .GlobalEnv)
      }
    }
    return(list(npx = npx, sinfo = sinfo))
    }
    return(list(npx = npx, sinfo = sinfo))
  }
