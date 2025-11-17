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
handle_NAs <- function(ptx, sinfo = NULL, choices = NULL, seed = global_seed,
                       NA_handling_log = NULL) {
  interactive <- is.null(choices)
  
  # Initialize optional variables
  maxiter <- NA_integer_
  ntree   <- NA_integer_
  n_pcs   <- NA_integer_
  q       <- NA_real_
  selected_indices <- NA
  
  # Initialize log if not provided
  if (is.null(NA_handling_log)) {
    NA_handling_log <- data.frame(
      data_choice = integer(),
      option      = integer(),
      columns     = I(list()),
      maxiter     = integer(),
      ntree       = integer(),
      n_pcs       = integer(),
      q           = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Prompt for which data to process
  if (interactive) {
    cat("\nWhich data would you like to process?\n")
    cat(" 1 = Abundance data (ptx)\n")
    cat(" 2 = Clinical data (sinfo)\n")
    data_choice <- as.numeric(readline(prompt = "Select an option: "))
  } else {
    data_choice <- choices$data_choice
  }
  
  # Branch for Abundance data
  if (!is.na(data_choice) && data_choice == 1) {
    if (interactive) {
      cat("\nCurrent number of NAs in Abundance data:", sum(is.na(ptx)), "\n")
      cat("\nNA Handling Options:\n")
      cat(" 1 = minProb\n 2 = RandomForest\n 3 = kNN\n 4 = SVD\n 5 = minDet\n 6 = complete.cases\n 7 = median\n 0 = Cancel\n")
      option <- as.numeric(readline(prompt = "Select an option: "))
    } else {
      option <- choices$option
    }
    
    if (is.na(option) || option == 0) return(list(ptx = ptx, sinfo = sinfo))
    
    if (option == 1) {
      result.ptx <- t(impute_na(t(as.matrix(ptx)), seed = seed))
    } else if (option == 2) {
      maxiter <- if (!is.null(choices$maxiter)) choices$maxiter else 5
      ntree   <- if (!is.null(choices$ntree)) choices$ntree else 100
      result.ptx <- t(impute_na(t(as.matrix(ptx)), method = "RF",
                                maxiter = maxiter, ntree = ntree, seed = seed))
    } else if (option == 3) {
      result.ptx <- t(impute_na(t(as.matrix(ptx)), method = "kNN", seed = seed))
    } else if (option == 4) {
      n_pcs <- if (!is.null(choices$n_pcs)) choices$n_pcs else 3
      result.ptx <- t(impute_na(t(as.matrix(ptx)), method = "SVD",
                                n_pcs = n_pcs, seed = seed))
    } else if (option == 5) {
      q <- if (!is.null(choices$q)) choices$q else 0.001
      result.ptx <- t(impute_na(t(as.matrix(ptx)), method = "minDet",
                                q = q, seed = seed))
    } else if (option == 6) {
      result.ptx <- ptx[complete.cases(ptx), ]
      if (!is.null(sinfo)) {
        sinfo <- sinfo[rownames(sinfo) %in% rownames(result.ptx), , drop = FALSE]
      }
    } else if (option == 7) {
      imputed_matrix <- apply(as.matrix(ptx), 2, impute_na_median)
      result.ptx <- as.data.frame(imputed_matrix, stringsAsFactors = FALSE)
      rownames(result.ptx) <- rownames(ptx)
    } else {
      return(list(ptx = ptx, sinfo = sinfo))
    }
    
    # Logging for abundance
    if (interactive) {
      new_entry <- data.frame(
        data_choice = data_choice,
        option      = option,
        columns     = I(list(NA)),
        maxiter     = maxiter,
        ntree       = ntree,
        n_pcs       = n_pcs,
        q           = q,
        stringsAsFactors = FALSE
      )
      NA_handling_log <<- rbind(NA_handling_log, new_entry)
    }
    
    return(list(ptx = result.ptx, sinfo = sinfo))
  }
  
  # Branch for Clinical data
  else if (!is.na(data_choice) && data_choice == 2) {
    if (is.null(sinfo)) return(list(ptx = ptx, sinfo = sinfo))
    
    if (interactive) {
      cat("\nAvailable columns:\n")
      for (i in seq_along(names(sinfo))) cat(i, "-", names(sinfo)[i], "\n")
      columns_input <- readline(prompt = "Enter column numbers (comma-separated): ")
      selected_indices <- as.integer(unlist(strsplit(columns_input, ",")))
    } else {
      selected_indices <- choices$columns
    }
    
    if (interactive) {
      cat("\nNA Handling Options:\n")
      cat(" 1 - minProb\n 2 - RandomForest\n 3 - kNN\n 4 - SVD\n 5 - minDet\n 6 - complete.cases\n 7 - median\n 0 - Cancel\n")
      option <- as.numeric(readline(prompt = "Select an option: "))
    } else {
      option <- choices$option
    }
    
    if (is.na(option) || option == 0) return(list(ptx = ptx, sinfo = sinfo))
    
    if (option == 1) {
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, drop = FALSE])), seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 2) {
      maxiter <- if (!is.null(choices$maxiter)) choices$maxiter else 5
      ntree   <- if (!is.null(choices$ntree)) choices$ntree else 100
      sinfo[, selected_indices] <- impute_na(sinfo[, selected_indices, drop = FALSE],
                                             method = "RF", maxiter = maxiter,
                                             ntree = ntree, seed = seed)
    } else if (option == 3) {
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, drop = FALSE])),
                                    method = "kNN", seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 4) {
      n_pcs <- if (!is.null(choices$n_pcs)) choices$n_pcs else 3
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, drop = FALSE])),
                                    method = "SVD", n_pcs = n_pcs, seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 5) {
      q <- if (!is.null(choices$q)) choices$q else 0.001
      imputed_matrix <- t(impute_na(t(as.matrix(sinfo[, selected_indices, drop = FALSE])),
                                    method = "minDet", q = q, seed = seed))
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    } else if (option == 6) {
      sinfo <- sinfo[complete.cases(sinfo[, selected_indices]), , drop = FALSE]
      common_samples <- intersect(rownames(sinfo), rownames(ptx))
      ptx <- ptx[common_samples, , drop = FALSE]
      sinfo <- sinfo[common_samples, , drop = FALSE]
    } else if (option == 7) {
      imputed_matrix <- apply(as.matrix(sinfo[, selected_indices, drop = FALSE]),
                              2, impute_na_median)
      sinfo[, selected_indices] <- as.data.frame(imputed_matrix)
    }
    
    # Logging for clinical
    if (interactive) {
      new_entry <- data.frame(
        data_choice = data_choice,
        option      = option,
        columns     = I(list(selected_indices)),
        maxiter     = maxiter,
        ntree       = ntree,
        n_pcs       = n_pcs,
        q           = q,
        stringsAsFactors = FALSE
      )
      NA_handling_log <<- rbind(NA_handling_log, new_entry)
    }
    
    return(list(ptx = ptx, sinfo = sinfo))
  }
  
  # Default return if no branch taken
  return(list(ptx = ptx, sinfo = sinfo))
}
