# ------------------------------------------------------------------------------
# Script:        DimensionReduction.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for Dimensionality Reduction Workflow
# ------------------------------------------------------------------------------

################################################################################
# Menu:
################################################################################
dimReductionMenu <- function() {
  cat("Select a dimensionality reduction method:\n")
  cat("  1 = PCA\n")
  cat("  2 = Sparse PCA\n")
  cat("  3 = PLS-DA\n")
  cat("  4 = Sparse PLS-DA\n")
  cat("  5 = Kernel PCA\n")
  cat("  6 = UMAP\n")
  cat("  7 = t-SNE\n")
  cat("  0 = Exit Dimensionality Reduction\n")
  
  choice <- as.numeric(readline(prompt = "Enter choice number (0-7): "))
  unreduced.ptx <<- select.ptx
  if (is.na(choice)) {
    cat("Invalid input. Exiting.\n")
    return(invisible(NULL))
  } 
  if (choice == 1) {
    cat("\nYou chose PCA...\n")
    reducePCA(reproduce = NULL)
  } else if (choice == 2) {
    cat("\nYou chose Sparse PCA...\n")
    reduceSparsePCA(reproduce = NULL)
  } else if (choice == 3) {
    cat("\nYou chose PLS-DA...\n")
    reducePLSDA(reproduce = NULL)
  } else if (choice == 4) {
    cat("\nYou chose Sparse PLS-DA...\n")
    reduceSPLSDA(reproduce = NULL)
  } else if (choice == 5) {
    cat("\nYou chose Kernel PCA...\n")
    reduceKernelPCA(reproduce = NULL)
  } else if (choice == 6) {
    cat("\nYou chose UMAP...\n")
    reduceUMAP(reproduce = NULL)
  } else if (choice == 7) {
    cat("\nYou chose t-SNE...\n")
    reduceTSNE(reproduce = NULL)
  } else if (choice == 0) {
    break
  }
}

################################################################################
# PCA:
################################################################################
reducePCA <- function(reproduce = NULL) {
  cat("\nPerforming PCA...\n")
  
  # Run PCA with centering and scaling
  unreduced.ptx <<- select.ptx
  pca_result <- prcomp(select.ptx, center = FALSE, scale. = FALSE)
  
  # Calculate variance explained by each principal component:
  variances     <- pca_result$sdev^2
  prop_variance <- variances / sum(variances)
  cum_variance  <- cumsum(prop_variance)
  
  # Create dataframe to plot variance explanation
  variance_df <- data.frame(
    PC = 1:length(prop_variance),
    Variance = prop_variance * 100,       # percentage of variance for each PC
    CumVariance = cum_variance * 100      # cumulative percentage of variance
  )
  
  # Plot the variance explained by each PC and cumulative variance.
  p <- ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_point(aes(y = CumVariance), color = "red", size = 2) +
    geom_line(aes(y = CumVariance, group = 1), color = "red", linetype = "dashed") +
    labs(title = "Variance Explained by Each Principal Component",
         x = "Principal Component",
         y = "Percentage Variance Explained") +
    theme_minimal()
  
  print(p)
  
  if (is.null(reproduce)) { 
    reproduce <- data.frame(
      numPCs = NA
    )
    # Prompt for number of principal components to keep.
    num_components <- as.integer(readline(
      prompt = "\nEnter the number of principal components to keep (integer value): "
    ))
    reproduce$numPCs <- num_components
  } else {
      num_components <- reproduce$numPCs
    }
  
  
  # If the input is invalid, provide a default to 2 principal components
  if (is.na(num_components) || num_components < 1 || num_components > ncol(pca_result$x)) {
    cat("Invalid input. Using default of 2 components.\n")
    num_components <- 2
  }
  
  cat("\nNumber of components selected:", num_components, "\n")
  cat("Full cumulative variance vector:\n")
  print(cum_variance)
  
  # Extract the PCA scores for the selected components.
  pc_scores <- pca_result$x[, 1:num_components, drop = FALSE]
  pc_df     <- data.frame(pc_scores)
  
  # Create a PCA scores plot using ggplot2.
  # Uses the 'x_BC' column from select.sinfo (if available) for coloring. (study specific)
  if (exists("select.sinfo") && "x_BC" %in% colnames(select.sinfo)) {
    p2 <- ggplot(data.frame(pca_result$x), aes(x = PC1, y = PC2,
                                               color = as.factor(select.sinfo$x_BC))) +
      geom_point(size = 3) +
      stat_ellipse() +
      labs(title = "PCA Colored by x_BC",
           x = paste0("PC1: ", round((variances[1] / sum(variances)) * 100, 2), " %"),
           y = paste0("PC2: ", round((variances[2] / sum(variances)) * 100, 2), " %")) +
      theme_minimal()
    print(p2)
  } else {
    plot(pca_result$x[, 1:2], main = "PCA (PC1 vs PC2)")
  }
  
  # Save the PCA loadings into a new global variable "pca.loadings"
  pca_loadings <- pca_result$rotation[, 1:num_components, drop = FALSE]
  assign("pca.loadings", pca_loadings, envir = .GlobalEnv)
  
  # Update the global variable 'select.ptx' with the PCA scores.
  assign("select.ptx", pc_df, envir = .GlobalEnv)
  assign("reproduceDimRedPCA", reproduce, envir = .GlobalEnv) 
  cat("\nDimensionality reduction complete. 'select.ptx' now contains the PCA scores, ",
      "and loadings are saved in 'pca.loadings'.\n")
}

################################################################################
# Sparse PCA:
################################################################################

reduceSparsePCA <- function(reproduce = NULL) {
  if(is.null(reproduce)) {
    reproduce <- data.frame(
      ncompStart = NA,
      keepX = NA,
      ncompFinal = NA,
      stringsAsFactors = FALSE
    )
  }
  
  # Use grouping info (for coloring the scatter plot) if available. 
  # x_BC is a study specific variable and may not apply to all datasets
  # This can be replaced by a variable of interest using replace all
  grouping <- NULL
  if (exists("select.sinfo") && !is.null(select.sinfo)) {
    if ("x_BC" %in% colnames(select.sinfo)) {
      grouping <- as.factor(select.sinfo$x_BC)
    }
  }
  
  cat("\nPerforming Sparse PCA using mixOmics::spca...\n")
  
  # Prompt for the number of components.
  if (is.na(reproduce$ncompStart)) {
    ncomp_input <<- as.numeric(readline(prompt = "\nEnter the number of components for Sparse PCA: "))
    reproduce$ncompStart <- ncomp_input
    }
  else {
    ncomp_input <- reproduce$ncompStart
  }
  
  if (is.na(ncomp_input) || ncomp_input <= 0) {
    cat("Invalid input. Using default of 2 components.\n")
    ncomp_input <<- 2
  }
  
  # Compute maximum feasible components: min(nrow - 1, ncol)
  max_possible_comp <- min(nrow(select.ptx) - 1, ncol(select.ptx))
  if (ncomp_input > max_possible_comp) {
    cat("Requested number of components (", ncomp_input, 
        ") exceeds the maximum possible (", max_possible_comp, ").\n",
        "Setting number of components to ", max_possible_comp, ".\n", sep = "")
    ncomp_input <<- max_possible_comp
  }
  
  # Prompt for keepX values
  if (is.na(reproduce$keepX)) {
    keepX_input <- readline(prompt = "Enter the keepX value(s) (comma separated) for Sparse PCA [default 5]: ")
    reproduce$keepX <- keepX_input
  } else {
    keepX_input <- reproduce$keepX
  }
  
  if (nchar(trimws(keepX_input)) == 0) {
    keepX_vals <<- rep(5, ncomp_input)
    cat("Using default keepX = 5 for each component.\n")
  } else {
    keepX_split <- unlist(strsplit(keepX_input, split = ","))
    keepX_split <- as.numeric(trimws(keepX_split))
    if (length(keepX_split) == 1) {
      keepX_vals <<- rep(keepX_split, ncomp_input)
    } else if (length(keepX_split) == ncomp_input) {
      keepX_vals <<- keepX_split
    } else {
      cat("The number of keepX values provided does not match the number of components. Using the first value for all components.\n")
      keepX_vals <<- rep(keepX_split[1], ncomp_input)
    }
    cat("Using keepX =", paste(keepX_vals, collapse = ", "), "\n")
  }
  
  # Run Sparse PCA.
  sparsePCA_model <- mixOmics::spca(select.ptx, ncomp = ncomp_input, 
                                    keepX = keepX_vals, center = F, scale = F)
  
  # Check if the model returned scores.
  if (is.null(sparsePCA_model$x) || ncol(sparsePCA_model$x) < 1) {
    cat("Error: The Sparse PCA model did not return any components.\n")
    cat("This may be due to the chosen keepX values or data limitations.\n")
    cat("Try lowering the keepX values or reducing the number of components requested.\n")
    return(invisible(NULL))
  }
  
  avail_comps <- ncol(sparsePCA_model$x)
  if (ncomp_input > avail_comps) {
    cat("Requested number of components (", ncomp_input, 
        ") exceeds available components (", avail_comps, "). Using ", avail_comps, " components instead.\n", sep = "")
    ncomp_input <<- avail_comps
  }
  
  # Order components by explained variance 
  ########################################
  # Compute eigenvalues from the scores (x)
  eigenvals <- apply(sparsePCA_model$x, 2, function(x) sum(x^2))
  # Obtain the desired order (largest first)
  order_idx <- order(eigenvals, decreasing = TRUE)
  
  # Reorder scores and loadings accordingly.
  sparsePCA_model$x <- sparsePCA_model$x[, order_idx, drop = FALSE]
  sparsePCA_model$rotation <- sparsePCA_model$rotation[, order_idx, drop = FALSE]
  
  # Recalculate eigenvalues and proportion of variance explained.
  eigenvals <- eigenvals[order_idx]
  total_var <- sum(eigenvals)
  prop_expl <- eigenvals / total_var
  cum_prop <- cumsum(prop_expl)
  
  # Plot the results
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  
  # Barplot of percentage variance explained per component.
  barplot(prop_expl * 100,
          names.arg = paste("Comp", 1:ncomp_input),
          main = "Variance Explained per Component (%)",
          ylab = "Percent", ylim = c(0, max(prop_expl * 100) * 1.2))
  
  # Scatter plot for the first two components (if available).
  if (ncol(sparsePCA_model$x) >= 2) {
    if (!is.null(grouping)) {
      col_vector <- rainbow(length(levels(grouping)))
      point_colors <- col_vector[as.numeric(grouping)]
      plot(sparsePCA_model$x[, 1], sparsePCA_model$x[, 2],
           xlab = "Component 1", ylab = "Component 2",
           main = "Sparse PCA Scores (Comp1 vs Comp2)",
           pch = 16, col = point_colors)
      legend("topright", legend = levels(grouping), col = col_vector, pch = 16, cex = 0.8)
    } else {
      plot(sparsePCA_model$x[, 1], sparsePCA_model$x[, 2],
           xlab = "Component 1", ylab = "Component 2",
           main = "Sparse PCA Scores (Comp1 vs Comp2)",
           pch = 16)
      text(sparsePCA_model$x[, 1], sparsePCA_model$x[, 2],
           labels = rownames(sparsePCA_model$x),
           cex = 0.7, pos = 3)
    }
  }
  par(old_par)
  
  # List cumulative variance explained for each component.
  cat("\nCumulative variance explained by components:\n")
  for (i in 1:ncomp_input) {
    cat(i, ": ", round(cum_prop[i] * 100, 2), "%\n")
  }
  
  if (is.na(reproduce$ncompFinal)) {
    final_ncomp <<- as.numeric(readline(prompt = paste0("Enter number of components to include (1 to ", ncomp_input, "): ")))
    reproduce$ncompFinal <- final_ncomp
  } else {
    final_ncomp <- reproduce$ncompFinal
  }
  
  if (is.na(final_ncomp) || final_ncomp < 1 || final_ncomp > ncomp_input) {
    cat("Invalid input. Using all ", ncomp_input, " components.\n")
    final_ncomp <<- ncomp_input
  }
  
  # Subset the scores and loadings for the final specified number of components.
  sparsePCA_scores_final <- sparsePCA_model$x[, 1:final_ncomp, drop = FALSE]
  sparsePCA_df <- as.data.frame(sparsePCA_scores_final)
  sparsePCA_loadings <- sparsePCA_model$rotation[, 1:final_ncomp, drop = FALSE]
  
  # Update global variables.
  assign("sparsePCA.loadings", sparsePCA_loadings, envir = .GlobalEnv)
  assign("select.ptx", sparsePCA_df, envir = .GlobalEnv)
  assign("reproduceSPCA", reproduce, envir = .GlobalEnv)
  
  cat("\nDimensionality reduction complete. 'select.ptx' now contains the Sparse PCA scores for ",
      final_ncomp, " components, and loadings are saved in 'sparsePCA.loadings'.\n")
}

################################################################################
# PLS-DA:
################################################################################

reducePLSDA <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      group = NA,
      maxComp = NA,
      dist = NA,
      ncompFinal = NA,
      sortVar = NA,
      varMeasure = NA,
      stringsAsFactors = FALSE
    )
  }
  
  cat("\nPerforming PLS-DA tuning...\n")
  unreduced.ptx <<- select.ptx
  if (is.na(reproduce$group)) {
    # Interactive selection of Grouping Variable
    if (is.data.frame(select.sinfo)) {
      cat("Available columns in select.sinfo:\n")
      print(colnames(select.sinfo))
      
      group_input <- readline(prompt = "Enter the column(s) to use for grouping (comma separated, by names or indices): ")
      group_input <- unlist(strsplit(group_input, split = ","))
      group_input <<- trimws(group_input)
      
      # If the input can be converted to numbers, assume they are indices.
      if (all(!is.na(as.numeric(group_input)))) {
        sinfo_cols <<- as.numeric(group_input)
      } else {
        sinfo_cols <<- group_input
      }
      reproduce$group <- sinfo_cols
      # Create a temporary data frame using the selected columns.
      temp.plsda.sinfo <- select.sinfo[, sinfo_cols, drop = FALSE]
      # If more than one column is selected, combine them via interaction.
      if (ncol(temp.plsda.sinfo) > 1) {
        message("Multiple columns selected. Combining them into a single factor using interaction().")
        temp.plsda.sinfo <- interaction(temp.plsda.sinfo, drop = TRUE)
      } else {
        temp.plsda.sinfo <- as.factor(temp.plsda.sinfo[[1]])
      }
    } else {
      temp.plsda.sinfo <- select.sinfo
      warning("select.sinfo is not a data.frame. Using the provided vector directly.")
    }
  } else {
    sinfo_cols <- reproduce$group
    # Create a temporary data frame using the selected columns.
    temp.plsda.sinfo <- select.sinfo[, sinfo_cols, drop = FALSE]
    # If more than one column is selected, combine them via interaction.
    if (ncol(temp.plsda.sinfo) > 1) {
      message("Multiple columns selected. Combining them into a single factor using interaction().")
      temp.plsda.sinfo <- interaction(temp.plsda.sinfo, drop = TRUE)
    } else {
      temp.plsda.sinfo <- as.factor(temp.plsda.sinfo[[1]])
    }
  }
  
  cat("Using grouping variable with levels:\n")
  print(levels(temp.plsda.sinfo))
  
  if (is.na(reproduce$maxComp)) {
    # Maximum Component Selection for Tuning
    ncomp_input <- as.numeric(readline(prompt = "\nEnter the maximum number of components to consider for tuning: "))
    reproduce$maxComp <- ncomp_input
    if (is.na(ncomp_input) || ncomp_input <= 0) {
      cat("Invalid input. Using default of 2 components.\n")
      ncomp_input <- 2
    }
  } else {
      ncomp_input <- reproduce$maxComp
    }
  if (is.na(reproduce$dist)) {
    # Distance Measure Selection for tune.plsda
    cat("\nSelect the distance measure for tune.plsda:\n")
    cat("1: centroids.dist\n")
    cat("2: mahalanobis.dist\n")
    cat("3: max.dist\n")
    dist_choice <- as.numeric(readline(prompt = "Enter the distance choice (1-3): "))
    reproduce$dist <- dist_choice
  } else {
    dist_choice <- reproduce$dist
  }
  
  
  if (!is.na(dist_choice)) {
    if (dist_choice == 1) {
      selected_dist <- "centroids.dist"
    } else if (dist_choice == 2) {
      selected_dist <- "mahalanobis.dist"
    } else if (dist_choice == 3) {
      selected_dist <- "max.dist"
    } else {
      cat("Invalid selection, using default 'max.dist'\n")
      selected_dist <- "max.dist"
    }
  } else {
    cat("Invalid input, using default 'max.dist'\n")
    selected_dist <- "max.dist"
  }
  
  # Run tune.plsda
  cat("\nTuning PLS-DA model with tune.plsda using", selected_dist, "...\n")
  tune_res <- tune.plsda(select.ptx, temp.plsda.sinfo, ncomp = ncomp_input,
                         validation = "Mfold", folds = 10, dist = selected_dist,
                         progressBar = TRUE, seed = global_seed, nrepeat = 5, scale = FALSE)
  
  # Plot tuning results
  plot(tune_res)
  
  if (is.na(reproduce$ncompFinal)) {
    # Final Components Selection
    final_ncomp <- as.numeric(readline(prompt = paste0("Enter number of components to include (1 to ", ncomp_input, "): ")))
    reproduce$ncompFinal <- final_ncomp
    if (is.na(final_ncomp) || final_ncomp < 1 || final_ncomp > ncomp_input) {
      cat("Invalid input. Using all ", ncomp_input, " components.\n")
      final_ncomp <- ncomp_input
    }
  } else {
    final_ncomp <- reproduce$ncompFinal
  }
  
  # Run Final PLS-DA Model
  cat("\nRunning final PLS-DA model with", final_ncomp, "component(s)...\n")
  plsda_model <- plsda(select.ptx, temp.plsda.sinfo, ncomp = final_ncomp)
  
  if (is.na(reproduce$sortVar)) {
    # Prompt if to order components by variance explained
    order_choice <- toupper(readline(prompt = "\nDo you want to order components by variance explained? [Y/N]: "))
    reproduce$sortVar <- order_choice
  } else {
    order_choice <- reproduce$sortVar
  }
  
  if(order_choice == "Y"){
    if (is.na(varMeasure)) {
      # Prompt for choice of Variance Explained Measure 
      variance_choice <- toupper(readline(prompt = "\nChoose which variance explained measure to use [X/Y]:\nEnter 'X' for variance explained in the predictors (X), or 'Y' for variance explained in the response (Y): "))
      reproduce$varMeasure <- variance_choice
    } else {
      variance_choice <- reproduce$varMeasure
    }
     
    
    if (variance_choice == "Y") {
      if (!is.null(plsda_model$prop_expl_var$Y)) {
        var_expl <- plsda_model$prop_expl_var$Y
        cat("Using variance explained in Y.\n")
      } else {
        cat("Variance explained in Y is not available. Defaulting to variance explained in X.\n")
        var_expl <- plsda_model$prop_expl_var$X
      }
    } else {
      var_expl <- plsda_model$prop_expl_var$X
      cat("Using variance explained in X.\n")
    }
    
    # Convert to numeric.
    if (is.matrix(var_expl)) {
      var_expl <- as.numeric(var_expl[1, ])
    } else {
      var_expl <- as.numeric(var_expl)
    }
    
    # Sort components by chosen variance explained measure.
    order_idx <- order(var_expl, decreasing = TRUE)
    var_expl <- var_expl[order_idx]
    cum_var <- cumsum(var_expl)
    
    # Reorder the PLS-DA scores and loadings matrices.
    pls_scores_full <- plsda_model$variates$X[, order_idx, drop = FALSE]
    plsda_loadings <- plsda_model$loadings$X[, order_idx, drop = FALSE]
    
  } else {
    # Do not reorder; use default extraction order (X variance)
    var_expl <- plsda_model$prop_expl_var$X
    if (is.matrix(var_expl)) {
      var_expl <- as.numeric(var_expl[1, ])
    } else {
      var_expl <- as.numeric(var_expl)
    }
    cum_var <- cumsum(var_expl)
    pls_scores_full <- plsda_model$variates$X
    plsda_loadings <- plsda_model$loadings$X
  }
  
  # Plot Results
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  
  # Barplot for variance explained per component (in percentage)
  barplot(var_expl * 100,
          names.arg = paste("Comp", 1:final_ncomp),
          main = "Variance Explained per Component (%)",
          ylab = "Percent")
  
  # Scatter plot of the first two components.
  if (ncol(pls_scores_full) >= 2) {
    plot(pls_scores_full[, 1], pls_scores_full[, 2],
         xlab = "Component 1", ylab = "Component 2",
         main = "PLS-DA Scores (Comp1 vs Comp2)",
         pch = 16)
    text(pls_scores_full[, 1], pls_scores_full[, 2],
         labels = rownames(pls_scores_full),
         cex = 0.7, pos = 3)
  }
  par(old_par)
  
  # Show Cumulative Variance
  cat("\nCumulative variance explained by components (sorted by variance if ordered):\n")
  for (i in seq_along(var_expl)) {
    cat(i, ": ", round(cum_var[i] * 100, 2), "%\n")
  }
  
  # Update and assign to global environment
  pls_scores_final <- pls_scores_full[, 1:final_ncomp, drop = FALSE]
  pls_df <- as.data.frame(pls_scores_final)
  plsda_loadings <- plsda_loadings[, 1:final_ncomp, drop = FALSE]
  
  assign("plsda.loadings", plsda_loadings, envir = .GlobalEnv)
  assign("select.ptx", pls_df, envir = .GlobalEnv)
  assign("reproducePLSDA", reproduce, envir = .GlobalEnv)
  assign("temp.plsda.sinfo", temp.plsda.sinfo, envir =.GlobalEnv)
  assign("plsda.ncomp", final_ncomp, envir = .GlobalEnv)
  
  cat("\nDimensionality reduction complete. 'select.ptx' now contains the PLS-DA scores for", 
      final_ncomp, "component(s) (", if(order_choice == "Y") "sorted by variance" else "in extraction order", 
      "), and loadings are saved in 'plsda.loadings'.\n")
}


################################################################################
# Sparse PLS-DA:
################################################################################

reduceSPLSDA <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      group = NA,
      maxComp = NA,
      dist = NA,
      ncompFinal = NA,
      keepX = I(list(NULL)),   
      stringsAsFactors = FALSE
    )
  }
  
  cat("\nPerforming sPLS-DA tuning...\n")
  
  # Grouping variable selection
  if (is.na(reproduce$group)) {
    if (is.data.frame(select.sinfo)) {
      cat("Available columns in select.sinfo:\n")
      print(colnames(select.sinfo))
      
      group_input <- readline(prompt = "Enter the column(s) to use for grouping (comma separated, by names or indices): ")
      group_input <- trimws(unlist(strsplit(group_input, split = ",")))
      
      if (all(!is.na(as.numeric(group_input)))) {
        sinfo_cols <- as.numeric(group_input)
      } else {
        sinfo_cols <- group_input
      }
      reproduce$group <- sinfo_cols
      
      temp.select.sinfo <- select.sinfo[, sinfo_cols, drop = FALSE]
      if (ncol(temp.select.sinfo) > 1) {
        message("Multiple columns selected. Combining them into a single factor using interaction().")
        temp.select.sinfo <- interaction(temp.select.sinfo, drop = TRUE)
      } else {
        temp.select.sinfo <- as.factor(temp.select.sinfo[[1]])
      }
    } else {
      temp.select.sinfo <- select.sinfo
      warning("select.sinfo is not a data.frame. Using the provided vector directly.")
    }
  } else {
    sinfo_cols <- reproduce$group
    temp.select.sinfo <- select.sinfo[, sinfo_cols, drop = FALSE]
    if (ncol(temp.select.sinfo) > 1) {
      message("Multiple columns selected. Combining them into a single factor using interaction().")
      temp.select.sinfo <- interaction(temp.select.sinfo, drop = TRUE)
    } else {
      temp.select.sinfo <- as.factor(temp.select.sinfo[[1]])
    }
  }
  
  cat("Using grouping variable with levels:\n")
  print(levels(temp.select.sinfo))
  
  # Prompt for Max components ---
  if (is.na(reproduce$maxComp)) {
    ncomp_tuning <- as.numeric(readline(prompt = "\nEnter the maximum number of components to consider for tuning (default 5): "))
    if (is.na(ncomp_tuning) || ncomp_tuning <= 0) {
      cat("Invalid input. Using default of 5 components for tuning.\n")
      ncomp_tuning <- 5
    }
    reproduce$maxComp <- ncomp_tuning
  } else {
    ncomp_tuning <- reproduce$maxComp
  }
  
  # Prompt for Candidate keepX values
  if (is.null(reproduce$keepX[[1]])) {
    keepX_input <- readline(prompt = "Enter candidate keepX value(s) to test (comma separated, default: 5,10,15): ")
    if (nchar(trimws(keepX_input)) == 0) {
      candidate_keepX <- c(5, 10, 15)
      cat("Using default candidate keepX values: ", paste(candidate_keepX, collapse = ", "), "\n")
    } else {
      candidate_keepX <- as.numeric(unlist(strsplit(keepX_input, split = ",")))
      if (any(is.na(candidate_keepX))) {
        cat("Invalid input detected. Using default candidate keepX values: 5,10,15\n")
        candidate_keepX <- c(5, 10, 15)
      } else {
        cat("Using candidate keepX values: ", paste(candidate_keepX, collapse = ", "), "\n")
      }
    }
  } else {
    candidate_keepX <- reproduce$keepX[[1]]
  }
  
  # Prompt for Distance measure
  if (is.na(reproduce$dist)) {
    cat("\nSelect the distance measure for tuning sPLS-DA:\n")
    cat("1: centroids.dist\n2: mahalanobis.dist\n3: max.dist\n")
    dist_choice <- as.numeric(readline(prompt = "Enter the distance choice (1-3): "))
    reproduce$dist <- dist_choice
  } else {
    dist_choice <- reproduce$dist
  }
  
  if (!is.na(dist_choice)) {
    selected_dist <- switch(as.character(dist_choice),
                            "1" = "centroids.dist",
                            "2" = "mahalanobis.dist",
                            "3" = "max.dist",
                            "max.dist")
  } else {
    selected_dist <- "max.dist"
  }
  
  # Run tuning
  cat("\nTuning sPLS-DA model using tune.splsda with", ncomp_tuning, "components and distance =", selected_dist, "...\n")
  tune_res <- tune.splsda(select.ptx, temp.select.sinfo,
                          ncomp = ncomp_tuning,
                          test.keepX = candidate_keepX,
                          validation = "Mfold", folds = 10,
                          scale = FALSE, progressBar = TRUE,
                          dist = selected_dist, nrepeat = 10,
                          light.output = FALSE, seed = global_seed)
  
  plot(tune_res, sd = FALSE)
  
  cat("\nTuning results summary:\n")
  print(tune_res$choice.keepX)
  
  cat("\nAverage error rates (per component) across tuning repetitions:\n")
  avg_errors <- apply(tune_res$error.rate, 2, mean, na.rm = TRUE)
  print(avg_errors)
  
  chosen_ncomp <- which.min(avg_errors)
  cat("Tuned optimal number of components (minimizing average error): ", chosen_ncomp, "\n")
  
  # Prompt for Final number of components
  if (is.na(reproduce$ncompFinal)) {
    final_ncomp <- as.numeric(readline(prompt = paste0("Enter desired number of components (1 to ", ncomp_tuning, 
                                                       ") [default: ", chosen_ncomp, "]: ")))
    if (is.na(final_ncomp) || final_ncomp < 1 || final_ncomp > ncomp_tuning) {
      cat("Invalid input. Using tuned optimal number of components: ", chosen_ncomp, "\n")
      final_ncomp <- chosen_ncomp
    }
    reproduce$ncompFinal <- final_ncomp
  } else {
    final_ncomp <- reproduce$ncompFinal
  }
  
  # Final keepX vector
  ckx <- tune_res$choice.keepX
  if (is.list(ckx)) {
    final_keepX <- unlist(ckx, use.names = FALSE)[1:final_ncomp]
  } else if (is.matrix(ckx)) {
    final_keepX <- as.numeric(ckx)[1:final_ncomp]
  } else {
    final_keepX <- as.numeric(ckx)[1:final_ncomp]
  }
  cat("Using tuned keepX values for each component: ", paste(final_keepX, collapse = ", "), "\n")
  
  # Store final keepX vector for reproducibility
  reproduce$keepX[[1]] <- final_keepX
  
  # Run final sPLS-DA
  cat("\nRunning final sPLS-DA model with", final_ncomp, "component(s)...\n")
  splsda_model <- splsda(select.ptx, temp.select.sinfo,
                         ncomp = final_ncomp,
                         keepX = final_keepX)
  
  # Variance explained
  var_expl <- splsda_model$prop_expl_var$X 
  if (is.matrix(var_expl)) var_expl <- as.numeric(var_expl[1, ])
  cum_var <- cumsum(var_expl)
  
  # Plots
  spls_scores_full <- splsda_model$variates$X[, 1:final_ncomp, drop = FALSE]
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  barplot(var_expl * 100, names.arg = paste("Comp", 1:final_ncomp),
          main = "Variance Explained per Component (%)", ylab = "Percent")
  if (ncol(spls_scores_full) >= 2) {
    plot(spls_scores_full[, 1], spls_scores_full[, 2],
         xlab = "Component 1", ylab = "Component 2",
         main = "sPLS-DA Scores (Comp1 vs Comp2)", pch = 16)
    text(spls_scores_full[, 1], spls_scores_full[, 2],
         labels = rownames(spls_scores_full),
         cex = 0.7, pos = 3)
  }
  par(old_par)
  
  # Display cumulative variance explained per component
  cat("\nCumulative variance explained by components:\n")
  for (i in 1:final_ncomp) {
    cat(i, ": ", round(cum_var[i] * 100, 2), "%\n")
  }
  
  # Save results to global environment
  spls_scores_final <- splsda_model$variates$X[, 1:final_ncomp, drop = FALSE]
  spls_df <- as.data.frame(spls_scores_final)
  splsda_loadings <- splsda_model$loadings$X[, 1:final_ncomp, drop = FALSE]
  
  assign("splsda.loadings", splsda_loadings, envir = .GlobalEnv)
  assign("select.ptx", spls_df, envir = .GlobalEnv)
  assign("reproduceSPLSDA", reproduce, envir = .GlobalEnv)
  assign("temp.splsda.sinfo", temp.select.sinfo, envir = .GlobalEnv)
  assign("splsda.keepX", final_keepX, envir = .GlobalEnv)
  assign("splsda.ncomp", final_ncomp, envir = .GlobalEnv)
  
  cat("\nDimensionality reduction complete. 'select.ptx' now contains the sPLS-DA scores for", 
      final_ncomp, "component(s), and loadings are saved in 'splsda.loadings'.\n")
}


################################################################################
# Kernel PCA 
################################################################################
reduceKernelPCA <- function(reproduce = NULL) {
  # Check that select.ptx is available.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before kernel PCA.\n")
    return(invisible(NULL))
  }
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      nFeatures = NA,
      sigmaX = NA,
      keep = NA,
      stringsAsFactors = FALSE
    )
  
  repeat {
    # Prompt for the number of features (projected dimensions) for kernel PCA
    features <<- as.numeric(readline(prompt = "Enter the number of features for Kernel PCA [default: 20]: "))
    if (is.na(features) || features <= 0) {
      cat("Invalid input. Using default value of 20.\n")
      features <<- 20
    }
    reproduce$nFeatures <- features
    
    # Prompt for the sigma parameter for the RBF kernel.
    sigmaX <<- as.numeric(readline(prompt = "Enter the sigma parameter for the RBF kernel [default: 0.01]: "))
    if (is.na(sigmaX) || sigmaX <= 0) {
      cat("Invalid input. Using default sigma value of 0.01.\n")
      sigmaX <<- 0.01
    }
    reproduce$sigmaX <- sigmaX
    cat("\nComputing Kernel PCA with", features, "features and sigma =", sigmaX, "...\n")
    
    # Run kernel PCA using the rbfdot kernel.
    kpca_result <- kernlab::kpca(as.matrix(select.ptx),
                                 kernel = "rbfdot",
                                 kpar = list(sigma = sigmaX),
                                 features = features)
    
    # Extract the projected coordinates (Kernel PCA scores).
    pc_coords <- kernlab::rotated(kpca_result)
    
    # Extract eigenvalues from the KPCA result
    eig_vals <- kpca_result@eig
    prop_variance <- eig_vals / sum(eig_vals)
    cum_variance <- cumsum(prop_variance)
    
    # Create a data frame for plotting variance explained.
    variance_df <- data.frame(Component = factor(1:length(eig_vals)),
                              Variance = prop_variance * 100,       # percentage of variance for each component
                              CumVariance = cum_variance * 100)      # cumulative percentage
    
    # Plot the variance explained by each component with ggplot2.
    pvar <- ggplot2::ggplot(variance_df, ggplot2::aes(x = Component, y = Variance)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
      ggplot2::geom_point(ggplot2::aes(y = CumVariance), color = "red", size = 2) +
      ggplot2::geom_line(ggplot2::aes(y = CumVariance, group = 1), color = "red", linetype = "dashed") +
      ggplot2::labs(title = "Kernel PCA: Variance Explained by Each Component",
                    x = "Kernel PCA Component",
                    y = "Variance Explained (%)") +
      ggplot2::theme_minimal()
    
    print(pvar)
    
    # Prompt for how many components to keep.
    keep <- as.numeric(readline(prompt = "\nEnter the number of components to keep (integer value): "))
    if (is.na(keep) || keep < 1 || keep > ncol(pc_coords)) {
      cat("Invalid input. Using default value of 2 components.\n")
      keep <- 2
    }
    
    reproduce$keep <- keep
    # Subset the KPCA scores based on the number of selected components.
    pc_coords_subset <- pc_coords[, 1:keep, drop = FALSE]
    
    cat("\nUsing", keep, "components. The first few rows of these selected kernel PCA scores:\n")
    print(head(pc_coords_subset))
    
    # Build a data frame for plotting the first 2 components.
    # If only one component is selected, create a histogram instead.
    if (ncol(pc_coords_subset) >= 2) {
      df_dr <- data.frame(PC1 = pc_coords_subset[, 1], PC2 = pc_coords_subset[, 2])
    } else {
      df_dr <- data.frame(PC1 = pc_coords_subset[, 1])
    }
    
    # Basic plot:
    if (ncol(pc_coords_subset) >= 2) {
      plot(pc_coords_subset[, 1], pc_coords_subset[, 2],
           xlab = "Kernel PC1", ylab = "Kernel PC2",
           main = "Kernel PCA of Expression Matrix",
           pch = 19, col = "blue")
    } else {
      plot(pc_coords_subset[, 1], rep(1, nrow(pc_coords_subset)),
           xlab = "Kernel PC1", ylab = "Count",
           main = "Kernel PCA of Expression Matrix (1D)",
           pch = 19, col = "blue")
    }
    
    # Prepare a ggplot2 visualization.
    if (exists("select.sinfo") && !is.null(select.sinfo) &&
        "x_BC" %in% colnames(select.sinfo)) {
      if (ncol(pc_coords_subset) >= 2) {
        p <- ggplot(df_dr, aes(x = PC1, y = PC2,
                               color = as.factor(select.sinfo$x_BC))) +
          geom_point(size = 3) +
          labs(title = "Kernel PCA of Expression Matrix",
               x = "Kernel PC1", y = "Kernel PC2") +
          theme_minimal()
      } else {
        p <- ggplot(df_dr, aes(x = PC1)) +
          geom_histogram(fill = "blue", bins = 30) +
          labs(title = "Kernel PCA of Expression Matrix (1D)",
               x = "Kernel PC1", y = "Frequency") +
          theme_minimal()
      }
    } else {
      if (ncol(pc_coords_subset) >= 2) {
        p <- ggplot(df_dr, aes(x = PC1, y = PC2)) +
          geom_point(size = 3, color = "blue") +
          labs(title = "Kernel PCA of Expression Matrix",
               x = "Kernel PC1", y = "Kernel PC2") +
          theme_minimal()
      } else {
        p <- ggplot(df_dr, aes(x = PC1)) +
          geom_histogram(fill = "blue", bins = 30) +
          labs(title = "Kernel PCA of Expression Matrix (1D)",
               x = "Kernel PC1", y = "Frequency") +
          theme_minimal()
      }
    }
    print(p)
    
    # Prompt if to continue with chosen parameters.
    cont <- tolower(readline(prompt = "\nContinue with these parameters? (y/n): "))
    if (cont %in% c("y", "yes")) break
    cat("Let's try again with new parameters.\n\n")
  }
    assign("reproduceKernelPCA", reproduce, envir = .GlobalEnv)
  } else {
    sigmaX <- reproduce$sigmaX
    features <- reproduce$nFeatures
    keep <- reproduce$keep
    # Run kernel PCA using the rbfdot kernel.
    kpca_result <- kernlab::kpca(as.matrix(select.ptx),
                                 kernel = "rbfdot",
                                 kpar = list(sigma = sigmaX),
                                 features = features)
    
    # Extract the projected coordinates (Kernel PCA scores).
    pc_coords <- kernlab::rotated(kpca_result)
    pc_coords_subset <- pc_coords[, 1:keep, drop = FALSE]
  }
  
  # Bind selected scores to select.ptx and assign to global environment 
  assign("select.ptx", as.data.frame(pc_coords_subset), envir = .GlobalEnv)
  cat("\nKernel PCA complete. 'select.ptx' is now updated with the kernel PCA scores (", 
      keep, " components).\n")
}


################################################################################
# UMAP dimensionality reduction
################################################################################
reduceUMAP <- function(reproduce = NULL) {
  # Check that the select.ptx (abundance or expression matrix) is available.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before UMAP.\n")
    return(invisible(NULL))
  }
  
  if(is.null(reproduce)) {
    reproduce <- data.frame(
      ncomp = NA,
      minDist = NA,
      nNeighbors = NA,
      seed = NA,
      stringsAsFactors = FALSE
    )
    
  repeat {
    # Prompt for the number of UMAP components (output dimensions)
    ncomp <<- as.numeric(readline(prompt = "Enter the number of components for UMAP [default: 2]: "))
    if (is.na(ncomp) || ncomp <= 0) {
      cat("Invalid input. Using default value of 2.\n")
      ncomp <<- 2
    }
    reproduce$ncomp <- ncomp
    
    # Prompt for the min_dist parameter.
    min_dist <<- as.numeric(readline(prompt = "Enter the min_dist parameter for UMAP [default: 0.1]: "))
    if (is.na(min_dist) || min_dist <= 0) {
      cat("Invalid input. Using default min_dist value of 0.1.\n")
      min_dist <<- 0.1
    }
    reproduce$minDist <- min_dist
    
    # Prompt for the number of neighbors.
    n_neighbors <<- as.numeric(readline(prompt = "Enter the number of neighbors for UMAP [default: 15]: "))
    if (is.na(n_neighbors) || n_neighbors <= 0) {
      cat("Invalid input. Using default n_neighbors value of 15.\n")
      n_neighbors <<- 15
    }
    reproduce$nNeighbors <- n_neighbors
    
    # Prompt for a seed value to ensure reproducibility.
    seed_val <<- as.numeric(readline(prompt = "Enter seed value for reproducibility [default: 42]: "))
    if (is.na(seed_val)) {
      cat("Invalid input. Using default seed value of 42.\n")
      seed_val <<- 42
    }
    reproduce$seed <- seed_val
    
    cat("\nComputing UMAP with", ncomp, "components, min_dist =", min_dist,
        ", n_neighbors =", n_neighbors, "and seed =", seed_val, "...\n")
    
    # Create UMAP using uwot package
    umap_result <- uwot::umap(as.matrix(select.ptx),
                              n_components = ncomp,
                              n_neighbors = n_neighbors,
                              min_dist = min_dist,
                              metric = "euclidean",
                              seed = seed_val)
    
    cat("\nThe first few rows of the UMAP projection:\n")
    print(head(umap_result))
    
    # Build a data frame using the first two UMAP components.
    df_umap <- data.frame(UMAP1 = umap_result[, 1],
                          UMAP2 = umap_result[, 2])
    
    # Basic plot.
    plot(umap_result[, 1], umap_result[, 2],
         xlab = "UMAP1", ylab = "UMAP2",
         main = "UMAP Projection (Base Plot)",
         pch = 19, col = "blue")
    
    # Use ggplot2 for visualization. 
    # study specific variable CancerTime can be replaced with a variable of choice
    if (exists("select.sinfo") && !is.null(select.sinfo) &&
        "CancerTime" %in% colnames(select.sinfo)) {
      df_umap$Group <- as.factor(select.sinfo$CancerTime)
      p <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = Group)) +
        geom_point(size = 3) +
        labs(title = "UMAP Projection",
             x = "UMAP1", y = "UMAP2") +
        theme_minimal()
    } else {
      p <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
        geom_point(size = 3, color = "blue") +
        labs(title = "UMAP Projection",
             x = "UMAP1", y = "UMAP2") +
        theme_minimal()
    }
    print(p)
    
    # Prompt if to continue with the parameters or try new ones.
    cont <- tolower(readline(prompt = "\nContinue with these parameters? (y/n): "))
    if (cont %in% c("y", "yes")) break
    cat("Let's try again with new parameters.\n\n")
  }
    assign("reproduceUMAP", reproduce, envir = .GlobalEnv)
    
  } else{
    ncomp <- reproduce$ncomp
    min_dist <- reproduce$minDist
    n_neighbors <- reproduce$nNeighbors
    seed_val <- reproduce$seed
    
    umap_result <- uwot::umap(as.matrix(select.ptx),
                              n_components = ncomp,
                              n_neighbors = n_neighbors,
                              min_dist = min_dist,
                              metric = "euclidean",
                              seed = seed_val)
      
    }
  
  # Update select.ptx with the UMAP projection (assign to global environment)
  assign("select.ptx", as.data.frame(umap_result), envir = .GlobalEnv)
  cat("\nUMAP complete. 'select.ptx' is now updated with the UMAP coordinates.\n")
}

################################################################################ 
# t-SNE dimenionality reduction
################################################################################ 
reduceTSNE <- function(reproduce = NULL) {
  # Check that the select.ptx exists.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before t-SNE.\n")
    return(invisible(NULL))
  }
  
  if (is.null(reproduce)){
    reproduce <- data.frame(
      dims = NA,
      perplexity = NA,
      theta = NA,
      maxIter = NA,
      seed = NA,
      stringsAsFactors = FALSE
    )
  
  repeat {
    # Prompt for the output dimensions.
    dims2 <<- as.numeric(readline(prompt = "Enter the number of dimensions for t-SNE [default: 2]: "))
    if (is.na(dims2) || dims2 <= 0) {
      cat("Invalid input. Using default value of 2.\n")
      dims2 <<- 2
    }
    reproduce$dims <- dims2
    
    # Prompt for perplexity.
    perplexity <<- as.numeric(readline(prompt = "Enter perplexity for t-SNE [default: 30]: "))
    if (is.na(perplexity) || perplexity <= 0) {
      cat("Invalid input. Using default perplexity of 30.\n")
      perplexity <<- 30
    }
    reproduce$perplexity <- perplexity
    
    # Prompt for theta.
    theta <<- as.numeric(readline(prompt = "Enter theta for t-SNE (speed/accuracy trade-off) [default: 0.5]: "))
    if (is.na(theta) || theta < 0) {
      cat("Invalid input. Using default theta of 0.5.\n")
      theta <<- 0.5
    }
    reproduce$theta <- theta
    
    # Prompt for max iterations.
    max_iter <<- as.numeric(readline(prompt = "Enter max iterations for t-SNE [default: 1000]: "))
    if (is.na(max_iter) || max_iter <= 0) {
      cat("Invalid input. Using default max iterations of 1000.\n")
      max_iter <<- 1000
    }
    reproduce$maxIter <- max_iter
    
    # Prompt for a seed value.
    seed_val <<- as.numeric(readline(prompt = "Enter a seed value for reproducibility [default: 42]: "))
    if (is.na(seed_val)) {
      cat("Invalid input. Using default seed of 42.\n")
      seed_val <<- 42
    }
    reproduce$seed <- seed_val
    
    cat("\nComputing t-SNE with dims =", dims2,
        ", perplexity =", perplexity,
        ", theta =", theta,
        ", max_iter =", max_iter,
        "and seed =", seed_val, "...\n")
    
    set.seed(seed_val)
    # Run t-SNE.
    tsne_result <- Rtsne::Rtsne(as.matrix(select.ptx),
                                dims = dims2,
                                perplexity = perplexity,
                                theta = theta,
                                max_iter = max_iter,
                                verbose = TRUE,
                                check_duplicates = FALSE)
    
    # Extract the t-SNE coordinates.
    Y <- tsne_result$Y
    cat("\nThe first few rows of t-SNE coordinates:\n")
    print(head(Y))
    
    # Build a data frame for plotting the first two dimensions.
    df_tsne <- data.frame(TSNE1 = Y[, 1], TSNE2 = Y[, 2])
    
    # Basic plot.
    plot(Y[, 1], Y[, 2],
         xlab = "t-SNE 1", ylab = "t-SNE 2",
         main = "t-SNE Projection (Base Plot)",
         pch = 19, col = "blue")
    
    # Create a ggplot2 visualization.
    # Study specific variable x_BC can be replaced with variable of choice.
    if (exists("select.sinfo") && !is.null(select.sinfo) &&
        "x_BC" %in% colnames(select.sinfo)) {
      df_tsne$Group <- as.factor(select.sinfo$x_BC)
      p <- ggplot(df_tsne, aes(x = TSNE1, y = TSNE2, color = Group)) +
        geom_point(size = 3) +
        labs(title = "t-SNE Projection",
             x = "t-SNE 1", y = "t-SNE 2") +
        theme_minimal()
    } else {
      p <- ggplot(df_tsne, aes(x = TSNE1, y = TSNE2)) +
        geom_point(size = 3, color = "blue") +
        labs(title = "t-SNE Projection",
             x = "t-SNE 1", y = "t-SNE 2") +
        theme_minimal()
    }
    print(p)
    
    # Prompt if to continue with the chosen parameters.
    cont <- tolower(readline(prompt = "\nContinue with these parameters? (y/n): "))
    if (cont %in% c("y", "yes")) break
    cat("Let's try again with new parameters.\n\n")
  }
    assign("reproduceTSNE", reproduce, envir = .GlobalEnv)
  } else{
    dims2 <- reproduce$dims
    perplexity <- reproduce$perplexity
    theta <- reproduce$theta
    max_iter <- reproduce$maxIter
    set.seed(reproduce$seed)
    tsne_result <- Rtsne::Rtsne(as.matrix(select.ptx),
                                dims = dims2,
                                perplexity = perplexity,
                                theta = theta,
                                max_iter = max_iter,
                                verbose = TRUE,
                                check_duplicates = FALSE)
    
    # Extract the t-SNE coordinates.
    Y <- tsne_result$Y
  }
  
  # Update select.ptx with the t-SNE coordinates and assign to global environment
  assign("select.ptx", as.data.frame(Y), envir = .GlobalEnv)
  cat("\nt-SNE complete. 'select.ptx' is now updated with the t-SNE coordinates.\n")
}
