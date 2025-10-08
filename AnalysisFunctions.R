# ------------------------------------------------------------------------------
# Script:        AnalysisFunctions.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for Clustering and Analysis of obtained results
# ------------------------------------------------------------------------------

################################################################################
# Analysis Main Menu
################################################################################
analysis_main_menu <- function() {
  repeat {
    cat("\n=== Clustering and Analysis Menu: ===\n")
    cat("  1 = Clustering and Stability analysis \n")
    cat("  2 = Merge unstable clusters\n")
    cat("  3 = Validation of clustering results\n")
    cat("  4 = Visualize clustering results\n")
    cat("  5 = Cluster Similarity analysis\n")
    cat("  6 = Cluster-Trait correlation\n")
    cat("  7 = Mosaic plots\n")
    cat("  8 = Differential Abundance Analysis (categorical variables)\n")
    cat("  9 = Differential Abundance Analysis (continuous variables)\n")
    cat("  10 = Deviation Heatmap\n")
    cat("  0 = Exit\n")
    cat("-----------------------------------\n")
    
    choice <- as.integer(readline(prompt = "Enter your choice [0-10]: "))
    if (choice == 0) {
      cat("Exiting the program. Goodbye!\n")
      break
    } else if (choice == 1) {
      stability_analysis(reproduce = NULL)
      #stability_analysis()
    } else if (choice == 2) {
      merge_clusters(reproduce = NULL)  
    } else if (choice == 3) {
      validate_df <<- validateClusterIndices(expr = select.npx,
                                             clin = select.sinfo,
                                             parse_numbers = TRUE)
    } else if (choice == 4) {
      visualizationRunner()  
    } else if (choice == 5) {
      similarity_analysis()
    } else if (choice == 6) {
      trait_correlation()
    } else if (choice == 7) {
      interactive_mosaic_plot(select.sinfo)
    } else if (choice == 8) {
      res <- build_design_matrix_interactive(select.sinfo)
      run_differential_abundance(design_output = res,
                                 p_cutoff = 0.05, 
                                 logFC_cutoff = 0.3,
                                 adj_method = "bonferroni")
    } else if (choice == 9) {
      DEAcontinuous()
    } else if (choice == 10) {
      plotDeviationHeatmap()
    } else {
      break
    }
  }
}

################################################################################
# Clustering with stability analysis
################################################################################
stability_analysis <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      boot = NA,
      methods = I(list(NULL)),
      stringsAsFactors = FALSE)
    
    cat("### Stability Analysis ###\n")
    # Ask for the number of bootstrap replicates.
    B <- as.integer(readline(prompt = "Enter number of bootstrap replicates (B): "))
    if (is.na(B) || B <= 0) {
      cat("Invalid value for B. Exiting stability analysis.\n")
      return(invisible(NULL))
    }
    reproduce$boot <- B
    
    # Ensure that chosen_k exists.
    if (!exists("chosen_k")) {
      chosen_k <<- as.integer(readline(prompt = "Enter desired number of clusters (chosen_k): "))
    }
    
    # Make sure select.sinfo is defined as a global list.
    if (!exists("select.sinfo")) {
      select.sinfo <<- list()
    }
    
    # Define available methods.
    available_methods <- c("kmeans", "HT-Kmeans", "Fuzzy kmeans", "DBSCAN",
                           "Hierarchical DBSCAN", "GMM", "Mclust", "kmeanspp", "Spectral", "Agglomerative Hierarchical")
    
    # Display the available methods with a numbered list.
    cat("Select clustering method(s) by number:\n")
    for (i in seq_along(available_methods)) {
      cat(i, ". ", available_methods[i], "\n", sep = "")
    }
    cat("Tip: Enter a single number for one method OR comma-separated numbers for multiple methods.\n")
    
    method_input <- readline(prompt = "Enter your choice(s): ")
    chosen_nums <- as.numeric(trimws(unlist(strsplit(method_input, ","))))
    
    if (any(is.na(chosen_nums)) || length(chosen_nums) == 0) {
      cat("No valid method number(s) entered. Exiting.\n")
      return(invisible(NULL))
    }
    
    selected_methods <- available_methods[chosen_nums]
    reproduce$methods[[1]] <- selected_methods
  } else {
    B <- reproduce$boot
    selected_methods <- reproduce$methods[[1]]
  }
    
    # If only one method was selected:
    if (length(selected_methods) == 1) {
      chosen_method <- selected_methods[[1]]
      bootRes <- NULL
      
      if (chosen_method == "kmeans") {
        cat("Running stability analysis for kmeans...\n")
        bootRes <- clusterboot(select.npx, clustermethod = kmeansCBI,
                               k = chosen_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for kmeans:\n")
        print(bootRes$bootmean)
        select.sinfo$kmeans.clusters <<- paste0("kmeans_Cluster_", bootRes$result$result$cluster)
        
      } else if (chosen_method == "HT-Kmeans") {
        cat("Running stability analysis for HT-Kmeans...\n")
        bootRes <- clusterboot(select.npx, clustermethod = htkmeansWrapper,
                               k = HTKmeans_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for HT-Kmeans:\n")
        print(bootRes$bootmean)
        select.sinfo$HTKmeans.clusters <<- paste0("HTKmeans_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "Fuzzy kmeans") {
        cat("Running stability analysis for Fuzzy kmeans...\n")
        bootRes <- clusterboot(select.npx, clustermethod = fkmWrapper,
                               k = chosen_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for Fuzzy kmeans:\n")
        print(bootRes$bootmean)
        select.sinfo$FKM.clusters <<- paste0("FKM_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "DBSCAN") {
        cat("Running stability analysis for DBSCAN...\n")
        bootRes <- clusterboot(select.npx, clustermethod = dbscanWrapper,
                               k = chosen_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for DBSCAN:\n")
        print(bootRes$bootmean)
        select.sinfo$dbscan.clusters <<- paste0("dbscan_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "Hierarchical DBSCAN") {
        cat("Running stability analysis for Hierarchical DBSCAN...\n")
        bootRes <- clusterboot(select.npx, clustermethod = hdbscanWrapper,
                               k = chosen_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for Hierarchical DBSCAN:\n")
        print(bootRes$bootmean)
        select.sinfo$hdbscan.clusters <<- paste0("hdbscan_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "GMM") {
        cat("Running stability analysis for GMM...\n")
        bootRes <- clusterboot(select.npx, clustermethod = gmmWrapper_CR,
                               k = chosen_G, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for GMM (GMM):\n")
        print(bootRes$bootmean)
        select.sinfo$GMM.clusters <<- paste0("GMM_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "Mclust") {
        cat("Running stability analysis for Mclust...\n")
        bootRes <- clusterboot(select.npx, clustermethod = gmmWrapper,
                               k = chosen_G, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for Mclust (GMM):\n")
        print(bootRes$bootmean)
        select.sinfo$Mclust.clusters <<- paste0("Mclust_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "kmeanspp") {
        cat("Running stability analysis for kmeanspp...\n")
        bootRes <- clusterboot(select.npx, clustermethod = kmeansppWrapper,
                               k = chosen_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for kmeanspp (kmeanspp):\n")
        print(bootRes$bootmean)
        select.sinfo$kmeanspp.clusters <<- paste0("kmeanspp_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "Spectral") {
        cat("Running stability analysis for Spectral Clustering...\n")
        bootRes <- clusterboot(select.npx, clustermethod = speccWrapper,
                               k = spectral_k, B = B, bootmethod = "boot", seed = global_seed)
        cat("Bootstrap stability (bootmean) for Spectral Clustering:\n")
        print(bootRes$bootmean)
        select.sinfo$Spectral.clusters <<- paste0("Spectral_Cluster_", bootRes$result$partition)
        
      } else if (chosen_method == "Agglomerative Hierarchical") {
        cat("Running stability analysis for Agglomerative Hierarchical Clustering...\n")
        bootRes <- clusterboot(dist(select.npx, method = hclust_dist), clustermethod = hclustCBI,
                               k = hclust_k, B = B, bootmethod = "boot",
                               method = hclust_method, seed = global_seed)
        cat("Bootstrap stability (bootmean) for Agglomerative Hierarchical Clustering:\n")
        print(bootRes$bootmean)
        select.sinfo$hclust.clusters <<- paste0("hclust_Cluster_", bootRes$result$partition)
        
      } else {
        cat("Invalid clustering method selected.\n")
        return(invisible(NULL))
      }
      
    } else {
      # Process multiple selected methods.
      stabilities <- list()
      
      for (method in selected_methods) {
        cat("Running stability analysis for ", method, "...\n", sep = "")
        
        method_function <- switch(
          method,
          "kmeans"                   = kmeansCBI,
          "HT-Kmeans"                = htkmeansWrapper,
          "Fuzzy kmeans"             = fkmWrapper,
          "DBSCAN"                   = dbscanWrapper,
          "Hierarchical DBSCAN"      = hdbscanWrapper,
          "GMM"                      = gmmWrapper_CR,
          "Mclust"                   = gmmWrapper, 
          "kmeanspp"                 = kmeansppWrapper,
          "Spectral"                 = speccWrapper,
          "Agglomerative Hierarchical" = hclustCBI
        )
        
        # Decide k for the methods through saved variables
        current_k <- chosen_k
        if (method == "GMM") {
          current_k <- chosen_G
        } else if (method == "Mclust") { 
          current_k <- chosen_G_Mc
        } else if (method == "Spectral") {
          current_k <- spectral_k
        } else if (method == "HT-Kmeans") {
          current_k <- HTKmeans_K
        }
        
        if (method == "Agglomerative Hierarchical") {
          bootRes <- clusterboot(
            dist(select.npx, method = hclust_dist),
            clustermethod = hclustCBI,
            k            = hclust_k,
            method       = hclust_method,
            B            = B,
            bootmethod   = "boot",
            seed         = global_seed
          )
        } else {
          # Default call for all other methods
          bootRes <- clusterboot(
            select.npx,
            clustermethod = method_function,
            k            = current_k,
            B            = B,
            bootmethod   = "boot",
            seed         = global_seed
          )
        }
        
        stabilities[[method]] <- bootRes$bootmean
        
        
        # Save cluster assignments to select.sinfo based on the method
        if (method == "kmeans") {
          select.sinfo$kmeans.clusters <<- paste0("kmeans_Cluster_", bootRes$result$result$cluster)
        } else if (method == "HT-Kmeans") {
          select.sinfo$HTKmeans.clusters <<- paste0("HTKmeans_Cluster_", bootRes$result$partition)
        } else if (method == "Fuzzy kmeans") {
          select.sinfo$FKM.clusters <<- paste0("FKM_Cluster_", bootRes$result$partition)
        } else if (method == "DBSCAN") {
          select.sinfo$dbscan.clusters <<- paste0("dbscan_Cluster_", bootRes$result$partition)
        } else if (method == "Hierarchical DBSCAN") {
          select.sinfo$hdbscan.clusters <<- paste0("hdbscan_Cluster_", bootRes$result$partition)
        } else if (method == "GMM") {
          select.sinfo$GMM.clusters <<- paste0("GMM_Cluster_", bootRes$result$partition)
        } else if (method == "Mclust") {
          select.sinfo$Mclust.clusters <<- paste0("Mclust_Cluster_", bootRes$result$partition)
        } else if (method == "kmeanspp") {
          select.sinfo$kmeanspp.clusters <<- paste0("kmeanspp_Cluster_", bootRes$result$partition)
        } else if (method == "Spectral") {
          select.sinfo$Spectral.clusters <<- paste0("Spectral_Cluster_", bootRes$result$partition)
        } else if (method == "Agglomerative Hierarchical") {
          select.sinfo$hclust.clusters <<- paste0("hclust_Cluster_", bootRes$result$partition)
        }
      }
      
      # Construct the stability dataframe.
      max_clusters <- max(sapply(stabilities, length))
      stability.df <<- do.call(rbind, lapply(stabilities, function(v) {
        if (length(v) < max_clusters) {
          v <- c(v, rep(NA, max_clusters - length(v)))
        }
        c(v, Mean = mean(v, na.rm = TRUE))
      }))
      stability.df <<- as.data.frame(stability.df)
      stability.df$Method <<- rownames(stability.df)
      cluster_names <- paste0("Cluster", 1:max_clusters)
      colnames(stability.df)[1:max_clusters] <- cluster_names
      stability.df <<- stability.df[, c(ncol(stability.df), 1:(ncol(stability.df)-1))]
      assign("reproduceClustStability", reproduce, envir = .GlobalEnv)
      cat("Stability analysis complete. The stability dataframe (stability.df) is now available.\n")
      View(stability.df)
    }
} 


################################################################################
# Function to merge unstable clusters
################################################################################
merge_clusters <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      threshold = NA,
      method = NA,
      stringsAsFactors = FALSE
    )
    cat("### Merge Unstable Clusters ###\n")
    method_to_merge <- choose_clustering_method()
    reproduce$method <- method_to_merge
    
    threshold_input <- readline(prompt = "Enter stability threshold between 0 and 1 (clusters with bootmean < threshold will be renamed), or 'q' to go back: ")
    if (tolower(threshold_input) == "q") return(invisible(NULL))
    threshold <- as.numeric(threshold_input)
    if (is.na(threshold) || threshold < 0 || threshold > 1) {
      cat("Invalid threshold value. Aborting merge.\n")
      return(invisible(NULL))
    }
    reproduce$threshold <- threshold
  } else {
    method_to_merge <- reproduce$method
    threshold <- reproduce$threshold
  }
    if (!exists("stability.df")) {
      cat("Stability dataframe (stability.df) not found. Please run stability analysis (option 1 with multiple methods) first.\n")
      return(invisible(NULL))
    }
    
    method_key <- ifelse(method_to_merge == "GMM", "GMM", method_to_merge)
    sel_row <- stability.df[stability.df$Method == method_key, ]
    if (nrow(sel_row) == 0) {
      cat("Stability data for method", method_key, "not found.\n")
      return(invisible(NULL))
    }
    
    bootvalues <- as.numeric(sel_row[1, 2:(ncol(stability.df) - 1)])
    cat("Stability values for", method_to_merge, ":\n")
    print(bootvalues)
    
    unstable <- which(bootvalues < threshold)
    if (length(unstable) == 0) {
      cat("No clusters have stability below the threshold. No changes performed.\n")
      return(invisible(NULL))
    }
    cat("Clusters to be renamed as unstable:", paste(unstable, collapse = ", "), "\n")
    
    # Pull current assignments and coerce to character
    current_assignments_chr <- as.character(
      switch(method_to_merge,
             "kmeans" = select.sinfo$kmeans.clusters,
             "HT-Kmeans" = select.sinfo$HTKmeans.clusters,
             "Fuzzy kmeans" = select.sinfo$FKM.clusters,
             "DBSCAN" = select.sinfo$dbscan.clusters,
             "Hierarchical DBSCAN" = select.sinfo$hdbscan.clusters,
             "GMM" = select.sinfo$GMM.clusters,
             "Mclust" = select.sinfo$Mclust.clusters,
             "kmeanspp" = select.sinfo$kmeanspp.clusters,
             "Spectral" = select.sinfo$Spectral.clusters,
             "Agglomerative Hierarchical" = select.sinfo$hclust.clusters)
    )
    
    # Extract cluster number from label (e.g., "GMM_Cluster_4" -> 4)
    cluster_idx <- suppressWarnings(
      as.integer(sub(".*?(\\d+)$", "\\1", current_assignments_chr))
    )
    
    # Identify those with numeric suffix in unstable set
    to_mark <- !is.na(cluster_idx) & cluster_idx %in% unstable
    
    # Build the new uniform unstable label
    unstable_label <- paste0(method_key, "_Cluster_Unstable")
    
    # Overwrite unstable ones with the uniform label
    new_assignments <- ifelse(to_mark, unstable_label, current_assignments_chr)
    
    # Write back
    switch(method_to_merge,
           "kmeans" = { select.sinfo$kmeans.clusters <<- new_assignments },
           "HT-Kmeans" = { select.sinfo$HTKmeans.clusters <<- new_assignments },
           "Fuzzy kmeans" = { select.sinfo$FKM.clusters <<- new_assignments },
           "DBSCAN" = { select.sinfo$dbscan.clusters <<- new_assignments },
           "Hierarchical DBSCAN" = { select.sinfo$hdbscan.clusters <<- new_assignments },
           "GMM" = { select.sinfo$GMM.clusters <<- new_assignments },
           "Mclust" = { select.sinfo$Mclust.clusters <<- new_assignments},
           "kmeanspp" = { select.sinfo$kmeanspp.clusters <<- new_assignments },
           "Spectral" = { select.sinfo$Spectral.clusters <<- new_assignments },
           "Agglomerative Hierarchical" = { select.sinfo$hclust.clusters <<- new_assignments }
    )
    
    cat("Renaming completed. Updated cluster assignments for", method_to_merge, ":\n")
    print(new_assignments)
  }

################################################################################
# Calinski-Harabasz and Davies-Bouldin metrics for clustering results:
################################################################################

validateClusterIndices <- function(expr,
                                   clin,
                                   cluster_cols = grep("clusters$", names(clin), value = TRUE),
                                   parse_numbers = TRUE) {
  # Prepare expression matrix
  expr_mat <- if (is.data.frame(expr)) as.matrix(expr) else expr
  if (is.null(rownames(expr_mat)) || is.null(rownames(clin))) {
    warning("No rownames on expr or clin: assuming they are aligned by row order.")
  } else {
    # Align clinical to expression
    if (!all(rownames(expr_mat) %in% rownames(clin))) {
      stop("Not all samples in 'expr' are found in rownames of 'clin'.")
    }
    clin <- clin[rownames(expr_mat), , drop = FALSE]
  }
  
  # Iterate over each clustering column
  results <- lapply(cluster_cols, function(col) {
    # Raw labels
    raw <- as.character(clin[[col]])
    
    # Parse numeric suffix: e.g. "kmeans_Cluster_1" -> 1
    if (parse_numbers) {
      labels <- as.integer(sub(".*_(\\d+)$", "\\1", raw))
    } else {
      # Fallback to factor codes
      labels <- as.integer(factor(raw))
    }
    
    nclus <- length(unique(labels))
    if (nclus < 2) {
      return(data.frame(
        method             = col,
        clusters           = nclus,
        calinski_harabasz  = NA_real_,
        davies_bouldin     = NA_real_,
        stringsAsFactors   = FALSE
      ))
    }
    
    # Compute CH & DB indices
    idx <- clusterCrit::intCriteria(
      traj = expr_mat,
      part = labels,
      crit = c("Calinski_Harabasz", "Davies_Bouldin")
    )
    
    # Create a dataframe with the results
    data.frame(
      method             = col,
      clusters           = nclus,
      calinski_harabasz  = idx$calinski_harabasz,
      davies_bouldin     = idx$davies_bouldin,
      stringsAsFactors   = FALSE
    )
  })
  
  # Combine and return
  do.call(rbind, results)
}

################################################################################
# Function for cluster similarity analysis
################################################################################
similarity_analysis <- function() {
  cat("### Similarity Analysis ###\n\n")
  cat("Select an option:\n")
  cat("1. Overall (Adjusted Rand Similarity Index)\n")
  cat("2. Compare Methods (Jaccard Similarity %)\n")
  cat("q. Go back\n")
  
  # Prompt for analysis selection
  option <- readline(prompt = "Enter your choice: ")
  if(tolower(option) == "q") return(invisible(NULL))
  
  # Define available methods
  available_methods <- c("kmeans", "HT-Kmeans", "Fuzzy kmeans", 
                         "DBSCAN", "Hierarchical DBSCAN", "GMM", "Mclust",
                         "kmeanspp", "Spectral","Agglomerative Hierarchical")
  
  if(option == "1") {
    # OVERALL SIMILARITY 
    cat("\nSelect clustering method(s) to include for overall similarity analysis:\n")
    for(i in seq_along(available_methods)) {
      cat(i, ". ", available_methods[i], "\n", sep = "")
    }
    cat("Tip: Enter one or more numbers separated by commas (e.g., 1,3,5):\n")
    method_input <- readline(prompt = "Enter your choice(s): ")
    chosen_nums <- as.numeric(trimws(unlist(strsplit(method_input, ","))))
    
    if(any(is.na(chosen_nums)) || length(chosen_nums) == 0) {
      cat("No valid method number(s) entered. Exiting.\n")
      return(invisible(NULL))
    }
    
    selected_methods <- available_methods[chosen_nums]
    
    # Build a list of cluster assignments using the selected methods.
    clust_assignments <- list()
    for(method in selected_methods) {
      clust_assignments[[method]] <- switch(method,
                                            "kmeans" = select.sinfo$kmeans.clusters,
                                            "HT-Kmeans" = select.sinfo$HTKmeans.clusters,
                                            "Fuzzy kmeans" = select.sinfo$FKM.clusters,
                                            "DBSCAN" = select.sinfo$dbscan.clusters,
                                            "Hierarchical DBSCAN" = select.sinfo$hdbscan.clusters,
                                            "GMM" = select.sinfo$GMM.clusters,
                                            "Mclust" = select.sinfo$Mclust.clusters,
                                            "kmeanspp" = select.sinfo$kmeanspp.clusters,
                                            "Spectral" = select.sinfo$Spectral.clusters,
                                            "Agglomerative Hierarchical" = select.sinfo$hclust.clusters)
    }
    
    # Create similarity matrix based on the Adjusted Rand Index.
    nMethods <- length(selected_methods)
    sim_matrix <- matrix(NA, nrow = nMethods, ncol = nMethods,
                         dimnames = list(selected_methods, selected_methods))
    
    for(i in seq_along(selected_methods)) {
      for(j in seq_along(selected_methods)) {
        sim_matrix[i,j] <- adjustedRandIndex(clust_assignments[[i]], clust_assignments[[j]])
      }
    }
    
    cat("\nOverall Adjusted Rand Index Similarity Matrix:\n")
    print(sim_matrix)
    pheatmap(sim_matrix,
             main = "Overall Adjusted Rand Index Heatmap",
             cluster_rows = FALSE, cluster_cols = FALSE)
    
  } else if(option == "2") {
    # PAIRWISE COMPARISON 
    cat("\nSelect the first clustering method:\n")
    for(i in seq_along(available_methods)) {
      cat(i, ". ", available_methods[i], "\n", sep = "")
    }
    first_input <- readline(prompt = "Enter the number of the first method: ")
    first_num <- as.numeric(trimws(first_input))
    if(is.na(first_num) || first_num < 1 || first_num > length(available_methods)) {
      cat("Invalid choice for first method. Exiting.\n")
      return(invisible(NULL))
    }
    method1 <- available_methods[first_num]
    
    cat("\nSelect the second clustering method:\n")
    for(i in seq_along(available_methods)) {
      cat(i, ". ", available_methods[i], "\n", sep = "")
    }
    second_input <- readline(prompt = "Enter the number of the second method: ")
    second_num <- as.numeric(trimws(second_input))
    if(is.na(second_num) || second_num < 1 || second_num > length(available_methods)) {
      cat("Invalid choice for second method. Exiting.\n")
      return(invisible(NULL))
    }
    method2 <- available_methods[second_num]
    
    # Retrieve the cluster assignments for both methods.
    clusterMethod1 <- switch(method1,
                             "kmeans" = select.sinfo$kmeans.clusters,
                             "HT-Kmeans" = select.sinfo$HTKmeans.clusters,
                             "Fuzzy kmeans" = select.sinfo$FKM.clusters,
                             "DBSCAN" = select.sinfo$dbscan.clusters,
                             "Hierarchical DBSCAN" = select.sinfo$hdbscan.clusters,
                             "GMM" = select.sinfo$GMM.clusters,
                             "Mclust" = select.sinfo$Mclust.clusters,
                             "kmeanspp" = select.sinfo$kmeanspp.clusters,
                             "Spectral" = select.sinfo$Spectral.clusters,
                             "Agglomerative Hierarchical" = select.sinfo$hclust.clusters)
    
    clusterMethod2 <- switch(method2,
                             "kmeans" = select.sinfo$kmeans.clusters,
                             "HT-Kmeans" = select.sinfo$HTKmeans.clusters,
                             "Fuzzy kmeans" = select.sinfo$FKM.clusters,
                             "DBSCAN" = select.sinfo$dbscan.clusters,
                             "Hierarchical DBSCAN" = select.sinfo$hdbscan.clusters,
                             "GMM" = select.sinfo$GMM.clusters,
                             "Mclust" = select.sinfo$Mclust.clusters,
                             "kmeanspp" = select.sinfo$kmeanspp.clusters,
                             "Spectral" = select.sinfo$Spectral.clusters,
                             "Agglomerative Hierarchical" = select.sinfo$hclust.clusters)
    
    # Define a helper function to compute the Jaccard similarity matrix.
    compute_similarity_matrix <- function(clusterA, clusterB) {
      clustersA <- sort(unique(clusterA))
      clustersB <- sort(unique(clusterB))
      sim_matrix <- matrix(NA, nrow = length(clustersA), ncol = length(clustersB),
                           dimnames = list(as.character(clustersA), as.character(clustersB)))
      for(i in clustersA) {
        idxA <- which(clusterA == i)
        for(j in clustersB) {
          idxB <- which(clusterB == j)
          inter <- length(intersect(idxA, idxB))
          union_val <- length(union(idxA, idxB))
          if(union_val > 0) {
            sim_matrix[as.character(i), as.character(j)] <- (inter / union_val) * 100
          } else {
            sim_matrix[as.character(i), as.character(j)] <- NA
          }
        }
      }
      return(sim_matrix)
    }
    
    simMat <- compute_similarity_matrix(clusterMethod1, clusterMethod2)
    cat("\nJaccard Similarity (%) matrix between", method1, "and", method2, ":\n")
    print(simMat)
    pheatmap(simMat,
             main = paste("Jaccard Similarity (%) between", method1, "and", method2),
             cluster_rows = FALSE, cluster_cols = FALSE)
    
  } else {
    cat("Invalid option. Returning to main menu.\n")
    return(invisible(NULL))
  }
}

################################################################################
# Cluster-Trait correlation Function
################################################################################
trait_correlation <- function() {
  cat("### Trait Correlation Analysis ###\n\n")
  
  method_choice <- choose_clustering_method()
 
  # Map the chosen method to the proper cluster column name,
  # using only the valid names from select.sinfo.
  cluster_col <- switch(method_choice,
                        "kmeans" = "kmeans.clusters",
                        "HT-Kmeans" = "HTKmeans.clusters",
                        "Fuzzy kmeans" = "FKM.clusters",
                        "DBSCAN" = "dbscan.clusters",
                        "Hierarchical DBSCAN" = "hdbscan.clusters",
                        "GMM" = "GMM.clusters",
                        "Mclust" = "Mclust.clusters",
                        "kmeanspp" = "kmeanspp.clusters",
                        "Spectral" = "Spectral.clusters",
                        "Agglomerative Hierarchical" = "hclust.clusters")
  
  # Debug: Check if the cluster column exists in select.sinfo.
  if (!(cluster_col %in% colnames(select.sinfo))) {
    cat("The cluster column", cluster_col, "was not found in select.sinfo.\n")
    cat("Please ensure that the clustering results and trait data are present in select.sinfo.\n")
    return(invisible(NULL))
  }
  
  # Prompt for type of correlation (Within cluster Trait/Trait or Cluster/Trait)
  cat("Select Trait Correlation Option:\n")
  cat("1. Within Cluster Trait/Trait Correlation (using raw samples within one cluster)\n")
  cat("2. Cluster/Trait Correlation\n")
  cat("q. Go back\n")
  opt <- readline(prompt = "Enter your choice: ")
  if(tolower(opt) == "q") return(invisible(NULL))
  
  if(opt == "1") {
    # Within-cluster trait–trait correlation heatmap
    cat("Creating within-cluster trait–trait correlation heatmap...\n\n")
    
    # List and choose cluster
    clusters   <- select.sinfo[[cluster_col]]
    uniq_clust <- sort(unique(clusters))
    cat("Available clusters in '", cluster_col, "':\n", sep = "")
    for (i in seq_along(uniq_clust)) {
      cat(i, ": Cluster '", uniq_clust[i], "'\n", sep = "")
    }
    # Prompt for the cluster to focus on
    cl_choice <- as.integer(
      readline(prompt = "Enter the number of the cluster to focus on: ")
    )
    if (is.na(cl_choice) || !(cl_choice %in% seq_along(uniq_clust))) {
      stop("Invalid cluster selection.")
    }
    sel_cluster <- uniq_clust[cl_choice]
    cat("Selected cluster: '", sel_cluster, "'\n\n", sep = "")
    
    # Prompt for clinical-trait columns to include
    all_columns      <- colnames(select.sinfo)
    candidate_traits <- setdiff(all_columns, cluster_col)
    cat("Available clinical trait columns:\n")
    for (i in seq_along(candidate_traits)) {
      cat(i, ": ", candidate_traits[i], "\n", sep = "")
    }
    trait_selection_input <- trimws(
      readline(prompt = "Enter the numbers of the traits (comma-separated): ")
    )
    selected_indices <- as.numeric(unlist(strsplit(trait_selection_input, ",")))
    if (any(is.na(selected_indices)) ||
        any(selected_indices < 1) ||
        any(selected_indices > length(candidate_traits))) {
      stop("Invalid clinical trait column selection.")
    }
    selected_traits <- candidate_traits[selected_indices]
    cat("Selected clinical traits: ",
        paste(selected_traits, collapse = ", "), "\n\n", sep = "")
    
    # Subset to only chosen trait columns
    trait_data <- select.sinfo[, selected_traits, drop = FALSE]
    
    # One-hot encode characters/factors, pass numerics unchanged
    onehot_list <- lapply(names(trait_data), function(var) {
      vec <- trait_data[[var]]
      if (is.numeric(vec)) {
        df <- data.frame(vec)
        colnames(df) <- var
        return(df)
      }
      fac  <- as.factor(vec)
      lvls <- levels(fac)
      dummies <- sapply(lvls, function(lvl) as.numeric(fac == lvl))
      dummies <- as.data.frame(dummies)
      colnames(dummies) <- paste(var, lvls, sep = "_")
      return(dummies)
    })
    traits_onehot <- do.call(cbind, onehot_list)
    cat("One-hot encoded clinical trait data: ",
        nrow(traits_onehot), "×", ncol(traits_onehot), " features.\n")
    
    # Drop any dummy columns summing to zero
    nonEmptyCols  <- sapply(traits_onehot, function(x) sum(x, na.rm = TRUE) > 0)
    traits_onehot <- traits_onehot[, nonEmptyCols, drop = FALSE]
    cat("After dropping empty columns: ",
        nrow(traits_onehot), "×", ncol(traits_onehot), " features.\n\n")
    
    # Subset samples in selected cluster
    idx <- which(clusters == sel_cluster)
    if (length(idx) < 3) {
      stop("Need at least 3 samples in the cluster to compute correlations.")
    }
    traits_sub <- traits_onehot[idx, , drop = FALSE]
    
    # Compute Pearson correlation matrix
    corr_mat <- cor(traits_sub,
                    use    = "pairwise.complete.obs",
                    method = "pearson")
    
    # Prompt for variable name for saving plot
    var_name <- readline("Enter a name to save your plot object under: ")
    df_melt <- melt(corr_mat,
                    varnames   = c("Trait1","Trait2"),
                    value.name = "r")
    assign(var_name, df.melt, .GlobalEnv)
    
    p <- ggplot(df_melt, aes(x = Trait1, y = Trait2, fill = r)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
      scale_fill_gradient2(
        low      = "blue",
        mid      = "white",
        high     = "red",
        midpoint = 0,
        limits   = c(-1, 1),
        name     = paste0("Pearson r\n(cluster=", sel_cluster, ")")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank()
      ) +
      coord_fixed()
    
    print(p)
  }
  
  
  else if(opt == "2") {
    # Cluster vs. Clinical Trait Correlation Heatmap using One-Hot Encoding
    cat("Creating correlation heatmap between clusters and one-hot encoded clinical traits...\n")
    
    # Use the previously assigned cluster column
    cat("Using cluster column: ", cluster_col, "\n")
    
    # For clinical traits, remove the cluster column from the candidate list.
    all_columns <- colnames(select.sinfo)
    candidate_traits <- setdiff(all_columns, cluster_col)
    
    cat("Available clinical trait columns:\n")
    for (i in seq_along(candidate_traits)) {
      cat(i, ": ", candidate_traits[i], "\n", sep = "")
    }
    
    # Prompt for selection of traits (comma-separated numbers)
    trait_selection_input <- trimws(readline(prompt = "Enter the numbers corresponding to the clinical trait columns (separated by commas): "))
    selected_indices <- as.numeric(unlist(strsplit(trait_selection_input, ",")))
    if (any(is.na(selected_indices)) || any(selected_indices < 1) || any(selected_indices > length(candidate_traits))) {
      stop("Invalid clinical trait column selection.")
    }
    selected_traits <- candidate_traits[selected_indices]
    cat("Selected clinical traits: ", paste(selected_traits, collapse = ", "), "\n")
    
    # Subset the clinical trait data.
    trait_data <- select.sinfo[, selected_traits, drop = FALSE]
    
    # One-Hot Encode the clinical trait data
    onehot_list <- lapply(names(trait_data), function(var) {
      if (is.factor(trait_data[[var]])) {
        lvls <- levels(trait_data[[var]])
        # Create one column per level: indicator = 1 if the value equals the level, 0 otherwise.
        dummies <- sapply(lvls, function(lvl) as.numeric(trait_data[[var]] == lvl))
        dummies <- as.data.frame(dummies)
        colnames(dummies) <- paste(var, lvls, sep = "_")
        return(dummies)
      } else {
        # If already numeric, return the column.
        return(trait_data[var])
      }
    })
    traits_onehot <- do.call(cbind, onehot_list)
    cat("One-hot encoded clinical trait data. Original dimensions: ",
        nrow(traits_onehot), " samples x ", ncol(traits_onehot), " variables.\n")
    
    # Remove one-hot encoded columns with no group members 
    nonEmptyCols <- sapply(traits_onehot, function(col) sum(col, na.rm = TRUE)) > 0
    traits_onehot <- traits_onehot[, nonEmptyCols, drop = FALSE]
    cat("After removing empty columns, dimensions are: ",
        nrow(traits_onehot), " samples x ", ncol(traits_onehot), " variables.\n")
    
    # Compute Pearson Correlation Matrix
    # Retrieve the cluster assignments (using the already assigned cluster column)
    clusters <- select.sinfo[[cluster_col]]
    unique_clusters <- sort(unique(clusters))
    
    # Create an empty matrix for correlations.
    corr_matrix <- matrix(NA, nrow = length(unique_clusters), ncol = ncol(traits_onehot))
    rownames(corr_matrix) <- unique_clusters
    colnames(corr_matrix) <- colnames(traits_onehot)
    
    # For each unique cluster, create a binary indicator and compute Pearson correlation.
    for (i in seq_along(unique_clusters)) {
      clust <- unique_clusters[i]
      indicator <- as.numeric(clusters == clust)
      for (j in 1:ncol(traits_onehot)) {
        corr_matrix[i, j] <- cor(indicator, traits_onehot[[j]], use = "pairwise.complete.obs", method = "pearson")
      }
    }
    
    # Melt the correlation matrix into long format for plotting
    if(!requireNamespace("reshape2", quietly = TRUE))
      stop("Please install the reshape2 package.")
    melted_corr <- melt(corr_matrix, varnames = c("Cluster", "Trait"), value.name = "Correlation")
    
    # Plot the heatmap
    heatmap <- ggplot(melted_corr, aes(x = Trait, y = as.factor(Cluster), fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", Correlation)),
                color = "black", size = 3) +
      
      scale_fill_gradient2(
        low      = "blue",
        mid      = "white",
        high     = "red",
        midpoint = 0,
        limits   = c(-1, 1),
        space    = "Lab",
        name     = "Pearson\nCorrelation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10)
      ) +
      labs(
        title = "Heatmap: Clusters vs. One-Hot Encoded Clinical Traits (Pearson)",
        x     = "Clinical Trait (One-Hot Encoded)",
        y     = "Cluster"
      )
    var_name2 <- readline("Enter a name to save your plot object under: ")
    assign(var_name2, heatmap, .GlobalEnv )
    print(heatmap)
    
    
  }
  else {
    cat("Invalid option. Returning to main menu.\n")
    return(invisible(NULL))
  }
}

################################################################################
# Function for Deviation Heatmaps
################################################################################

plotDeviationHeatmap <- function() {
  # Check for required global objects (fit3 from differential abundance analysis)
  if (!exists("fit3", envir = .GlobalEnv))
    stop("Global object 'fit3' not found.")
  if (!exists("select.sinfo", envir = .GlobalEnv))
    stop("Global object 'select.sinfo' not found.")
  # At least one expression matrix should be present.
  if (!exists("unadjusted.npx", envir = .GlobalEnv) && !exists("unscaled.npx", envir = .GlobalEnv))
    stop("Neither 'unadjusted.npx' nor 'unscaled.npx' found in the global environment.")
  
  # Prompt for selection of contrast from fit3 
  contrast_options <- colnames(fit3$coefficients)
  if (length(contrast_options) == 0) stop("No contrasts found in fit3.")
  cat("Available contrasts in fit3:\n")
  for (i in seq_along(contrast_options)) {
    cat(i, ": ", contrast_options[i], "\n", sep = "")
  }
  contrast_choice <- readline(prompt = "Enter the number corresponding to the contrast to use: ")
  contrast_index <- as.numeric(contrast_choice)
  if (is.na(contrast_index) || contrast_index < 1 || contrast_index > length(contrast_options)) {
    stop("Invalid contrast selection.")
  }
  chosenContrast <- contrast_options[contrast_index]
  cat("Chosen contrast: ", chosenContrast, "\n")
  
  # Run topTable Using the Chosen Contrast
  top_df <- limma::topTable(fit3, coef = chosenContrast, number = Inf, 
                            adjust.method = "bonferroni", sort.by = "none")
  # Sort by adjusted p-value
  top_df <- top_df[order(top_df$adj.P.Val), ]
  cat("Displaying the first few rows of topTable:\n")
  print(head(top_df, 10))
  cat("Total proteins in topTable:", nrow(top_df), "\n")
  cat("Significant proteins in topTable:", length(rownames(top_df[top_df$adj.P.Val < 0.05,])), "\n")
  
  num_sig <- as.numeric(readline(prompt = "Enter how many significant proteins to include: "))
  if (is.na(num_sig) || num_sig < 1) stop("Invalid number entered.")
  sig_proteins <- rownames(top_df)[1:num_sig]
  cat("Top significant proteins:\n")
  print(sig_proteins)
  
  # One-Hot Encoding of Selected Variables in select.sinfo
  sinfo <- select.sinfo
  var_names <- colnames(sinfo)
  cat("Available columns in select.sinfo:\n")
  for (i in seq_along(var_names)) {
    cat(i, ": ", var_names[i], "\n", sep = "")
  }
  var_choice <- readline(prompt = "Enter the numbers (comma-separated) corresponding to the variables to include: ")
  var_indices <- as.numeric(unlist(strsplit(var_choice, ",")))
  if (any(is.na(var_indices)) || any(var_indices < 1) || any(var_indices > length(var_names))) {
    stop("Invalid selection of variables.")
  }
  chosen_vars <- var_names[var_indices]
  cat("Chosen variables:\n")
  print(chosen_vars)
  
  # Custom one-hot encoding that preserving original factor level labels.
  onehot_list <- lapply(chosen_vars, function(v) {
    # Convert to factor if not already.
    x <- sinfo[[v]]
    if (!is.factor(x)) x <- as.factor(x)
    levs <- levels(x)
    # Create a dummy column for each level.
    dummy_mat <- sapply(levs, function(l) as.integer(x == l))
    # Name the columns as "variable: level" to preserve informative labels.
    colnames(dummy_mat) <- paste(v, levs, sep = ": ")
    return(as.data.frame(dummy_mat))
  })
  # Combine dummy columns (side by side).
  onehot_df <- do.call(cbind, onehot_list)
  cat("One-hot encoded groups (column names):\n")
  print(colnames(onehot_df))
  
  # Choose Expression Matrix 
  cat("\nSelect the expression matrix to use for the heatmap:\n")
  cat(" 1 = unadjusted.npx\n")
  cat(" 2 = unscaled.npx\n")
  cat(" 3 = unreduced.npx\n")
  cat(" 4 = select.npx\n")
  matrix_choice <- as.integer(readline(prompt = "Enter your choice [1-4]: "))
  if (matrix_choice == 1) {
    if (!exists("unadjusted.npx", envir = .GlobalEnv))
      stop("Global object 'unadjusted.npx' not found.")
    expr_matrix <- unadjusted.npx
    cat("Using unadjusted.npx for the analysis.\n")
  } else if (matrix_choice == 2) {
    if (!exists("unscaled.npx", envir = .GlobalEnv))
      stop("Global object 'unscaled.npx' not found.")
    expr_matrix <- unscaled.npx
    cat("Using unscaled.npx for the analysis.\n")
  } else if (matrix_choice == 3) {
    if (!exists("unreduced.npx", envir = .GlobalEnv))
      stop("Global object 'unreduced.npx' not found.")
    expr_matrix <- unreduced.npx
    cat("Using unreduced.npx for the analysis.\n")
  } else if (matrix_choice == 4) {
    if (!exists("select.npx", envir = .GlobalEnv))
      stop("Global object 'select.npx' not found.")
    expr_matrix <- select.npx
    cat("Using select.npx for the analysis.\n")
  } else {
    stop("Invalid selection")
  }
  
  # Subset the Expression Matrix
  if (!all(sig_proteins %in% colnames(expr_matrix))) {
    warning("Some significant proteins not found in the selected expression matrix; using the intersection.")
    sig_proteins <- intersect(sig_proteins, colnames(expr_matrix))
  }
  # Subset the matrix to include only the significant proteins.
  prot_matrix <- expr_matrix[, sig_proteins, drop = FALSE]
  
  # Prompt for Transformation Method
  cat("\nSelect transformation method for the heatmap:\n")
  cat(" 1 - Raw Deviation (Group Mean - Overall Mean)\n")
  cat(" 2 - Z-score (Normalized Deviation: (Group Mean - Overall Mean) / Overall SD)\n")
  method_choice <- as.integer(readline(prompt = "Enter your choice (1 or 2): "))
  if (!(method_choice %in% c(1, 2))) {
    stop("Invalid method selection.")
  }
  
  # Compute the Transformation for Each One-Hot Group
  if (method_choice == 1) {
    # Compute the raw deviation (group mean - overall mean).
    deviation_matrix <- sapply(colnames(onehot_df), function(group) {
      group_idx <- which(onehot_df[[group]] == 1)
      apply(prot_matrix, 2, function(prot_expr) {
        overall_mean <- mean(prot_expr, na.rm = TRUE)
        group_mean <- mean(prot_expr[group_idx], na.rm = TRUE)
        group_mean - overall_mean
      })
    })
  } else if (method_choice == 2) {
    # Compute the z-score deviation (group mean - overall mean) / overall SD.
    deviation_matrix <- sapply(colnames(onehot_df), function(group) {
      group_idx <- which(onehot_df[[group]] == 1)
      apply(prot_matrix, 2, function(prot_expr) {
        overall_mean <- mean(prot_expr, na.rm = TRUE)
        overall_sd <- sd(prot_expr, na.rm = TRUE)
        group_mean <- mean(prot_expr[group_idx], na.rm = TRUE)
        if (overall_sd == 0) return(0)
        (group_mean - overall_mean) / overall_sd
      })
    })
  }
  deviation_matrix <- as.matrix(deviation_matrix)
  # Replace any non-finite values (NA, NaN, Inf) with 0.
  deviation_matrix[!is.finite(deviation_matrix)] <- 0
  
  # Prompt for if Columns Should be Clustered
  cat("Cluster columns in the heatmap?\n")
  cat(" 1 - Yes\n")
  cat(" 2 - No\n")
  cluster_choice <- as.integer(readline(prompt = "Enter your choice (1 for Yes, 2 for No): "))
  
  if (!is.na(cluster_choice) && cluster_choice == 1) {
    cluster_columns_flag <- TRUE
    cat("Clustering columns enabled.\n\n")
  } else {
    cluster_columns_flag <- FALSE
    cat("Clustering columns disabled.\n\n")
  }
  
  # Plot the Deviation Heatmap
  max_abs <- max(abs(deviation_matrix), na.rm = TRUE)
  custom_col <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))
  col_title <- readline("Enter plot title: ")
  rfntsiz <- as.integer(readline("Enter desired row font size: "))
  cfntsiz <- as.integer(readline("Enter desired column font size: "))
  heatmap_plot <- ComplexHeatmap::Heatmap(
    deviation_matrix,
    name = "Deviation",
    column_title = col_title,
    row_title = "Proteins",
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = rfntsiz),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = cfntsiz),
    col = custom_col,
    cluster_columns = cluster_columns_flag,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete"
  )
  save_choice <- readline(prompt ="Do you want to save the plot? y/n: \n")
  if (save_choice == "y") {
    var_name2 <- readline("Enter a name to save your plot object under: ")
    assign(var_name2, heatmap_plot, .GlobalEnv ) 
    ComplexHeatmap::draw(heatmap_plot, heatmap_legend_side = "right") 
  } else {
    cat("Plot not saved as a variable")
    ComplexHeatmap::draw(heatmap_plot, heatmap_legend_side = "right")
  }
  
  return(deviation_matrix)
}




