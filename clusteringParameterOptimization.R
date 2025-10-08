# ------------------------------------------------------------------------------
# Script:        clusteringParameterOptimization.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Parameter optimization functions for clustering analysis 
# ------------------------------------------------------------------------------

################################################################################
# Clustering Parameter Optimiztion
################################################################################

################################################################################
# Parameter Optimization Menu:
################################################################################
clustering_parameter_optimization_menu <- function() {
  cat("\n=== Clustering Parameter Optimization Menu: ===\n")
  cat("  1 = WSS\n")
  cat("  2 = Gap Statistic\n")
  cat("  3 = Silhouette\n")
  cat("  4 = NbClust Analysis\n")
  cat("  5 = Manual assignment of k\n")
  cat("  6 = Optimal Clusters GMM ClusterR (diagonal covariance)\n")
  cat("  7 = Optimal Clusters GMM Mclust\n")
  cat("  8 = Manual assignment of G\n")
  cat("  9 = EigenGap for SpectralClustering\n")
  cat("  10 = Optimal K for SpectralClustering\n")
  cat("  11 = Manual assignment of K for SpectralClustering\n")
  cat("  12 = Optimal K for k-means clustering\n")
  cat("  13 = Optimal K for Agglomerative Hierarchical Clustering\n")
  cat("  14 = Optimal lambda for HTKmeans\n")
  cat("  15 = Optimal eps for DBSCAN\n")
  cat("  16 = Optimal minPts Hierarchical DBSCAN\n")
  choice <- readline(prompt = "Enter your choice [1-16]: ")
  
  if (choice == "1") {
    run_wss_method()
  } else if (choice == "2") {
    run_gapStat_method()
  } else if (choice == "3") {
    run_silhouette_method()
  } else if (choice == "4") {
    run_nbclust_analysis()
  } else if (choice == "5") {
    run_choose_k_method()
  } else if (choice == "6") {
    optimalClustersGMM_Function(data = select.npx)
    chosen_G <<- as.integer(
      readline(prompt = "Enter choice of G for number of clusters: "))
  } else if (choice == "7") {
    optimal_G_Mclust(data = select.npx)
    chosen_G <<- as.integer(
      readline(prompt = "Enter choice of G for number of clusters: "))
    best_model <<- readline(prompt = "Enter name of best model: ")
  } else if (choice == "8") {
    run_choose_G_method()
  } else if (choice == "9") {
    kmax <- as.integer(readline(prompt = "Enter choice maximum number of clusters: "))
    spectral_eigengap(select.npx, n_eigs = kmax, return_all = TRUE)
  } else if (choice == "10") {
    optimizeSpectralK(data = select.npx)
    spectral_k <<- as.integer(
      readline(prompt = "Enter choice of k for number of clusters: "))
  } else if (choice == "11") {
    spectral_k <<- as.integer(
      readline(prompt = "Enter choice of k for number of clusters: "))
  } else if (choice == "12") {
    optimizeKmeansK(data = select.npx)
    chosen_k <<- as.integer(
      readline(prompt = "Enter choice of k for number of clusters: ")
    )
  } else if (choice == "13") {
    optimizeHclustK(data = select.npx)
    hclust_k <<- as.integer(
      readline(prompt = "Enter choice of k for number of clusters: ")
    )
  } else if (choice == "14") {
    optimal_lambda_HTKmeans(select.npx, k = chosen_k)
    
  } else if (choice == "15") {
    optimal_eps_DBSCAN(select.npx, plot_results = TRUE)
    run_dbscan_method()
  } else if (choice == "16") {
    optimal_minPts_hdbscan(data = select.npx, plot_results = TRUE)
    run_hdbscan_method()
  } else {
    cat("Invalid choice. Please try again.\n")
  }
}

################################################################################
# WSS elbow plot
################################################################################
run_wss_method <- function() {
  kmax <- as.integer(readline(prompt = "Enter kmax (maximum number of clusters): "))
  if (!is.na(kmax) && kmax > 0) {
    cat(sprintf("\nRunning WSS method with k.max = %d...\n", kmax))
    wss <- fviz_nbclust(select.npx, kmeans, k.max = kmax, method = "wss", nboot = 30)
    print(wss)
  } else {
    cat("Invalid kmax. Please enter a positive integer.\n")
  }
}

################################################################################
# Gap-statistic
################################################################################
run_gapStat_method <- function() {
  kmax <- as.integer(readline(prompt = "Enter kmax (maximum number of clusters): "))
  if (!is.na(kmax) && kmax > 0) {
    cat(sprintf("\nRunning gap_stat method with k.max = %d...\n", kmax))
    gapStat <-fviz_nbclust(select.npx, kmeans, k.max = kmax, 
                           nboot = 30,method = "gap_stat", maxSE = list(method = "Tibs2001SEmax", SE.factor = 1))
    print(gapStat)
  } else {
    cat("Invalid kmax. Please enter a positive integer.\n")
  }
}

################################################################################
# Silhouette score
################################################################################
run_silhouette_method <- function() {
  kmax <- as.integer(readline(prompt = "Enter kmax (maximum number of clusters): "))
  if (!is.na(kmax) && kmax > 0) {
    cat(sprintf("\nRunning silhouette method with k.max = %d...\n", kmax))
    silhouette <- fviz_nbclust(select.npx, kmeans, k.max = kmax, 
                               method = "silhouette", nboot = 30)
    print(silhouette)
  } else {
    cat("Invalid kmax. Please enter a positive integer.\n")
  }
}

################################################################################
# Interactive nbclust Analysis
################################################################################
run_nbclust_analysis <- function() {
  # List the options.
  cat("Select clustering method from the options below:\n")
  cat("  1. kmeans\n")
  cat("  2. ward.D\n")
  cat("  3. ward.D2\n")
  cat("  4. single\n")
  cat("  5. complete\n")
  cat("  6. average\n")
  cat("  7. mcquitty\n")
  cat("  8. median\n")
  cat("  9. centroid\n")
  
  # Prompt for choice of clusterig method
  method_choice <- readline(prompt = "Enter the number corresponding to the desired method: ")
  
  # Map the numeric choice to the corresponding method string.
  method_choices <- c("1" = "kmeans",
                      "2" = "ward.D",
                      "3" = "ward.D2",
                      "4" = "single",
                      "5" = "complete",
                      "6" = "average",
                      "7" = "mcquitty",
                      "8" = "median",
                      "9" = "centroid")
  
  if(!(method_choice %in% names(method_choices))) {
    cat("Invalid input. Defaulting to 'kmeans'.\n")
    chosen_method <- "kmeans"
  } else {
    chosen_method <- method_choices[[method_choice]]
  }
  
  # Prompt for minimum number of clusters
  min_nc_input <- readline(prompt = "Enter min.nc (minimum number of clusters; default = 2): ")
  if (min_nc_input == "") {
    min_nc <- 2
  } else {
    min_nc <- as.numeric(min_nc_input)
    if (is.na(min_nc) || min_nc < 1) {
      cat("Invalid input. Using default min.nc = 2.\n")
      min_nc <- 2
    }
  }
  
  # Prompt for maximum number of clusters
  max_nc_input <- readline(prompt = "Enter max.nc (maximum number of clusters; default = 20): ")
  if (max_nc_input == "") {
    max_nc <- 20
  } else {
    max_nc <- as.numeric(max_nc_input)
    if (is.na(max_nc) || max_nc < min_nc) {
      cat("Invalid input. Using default max.nc = 20.\n")
      max_nc <- 20
    }
  }
  
  # Check that the global variable select.npx exists.
  if (!exists("select.npx", envir = .GlobalEnv)) {
    cat("Error: 'select.npx' not found in the workspace.\n")
    return(invisible(NULL))
  }
  
  # Create temp variable for nbclust analysis
  data_for_nb <- select.npx
  
  cat("\nRunning NbClust Analysis with the following parameters:\n")
  cat("  Distance: 'euclidean'\n")
  cat("  Method: '", chosen_method, "'\n", sep = "")
  cat("  min.nc: ", min_nc, "\n", sep = "")
  cat("  max.nc: ", max_nc, "\n", sep = "")
  cat("  Index: 'all'\n\n")
  
  # Run NbClust
  res <- NbClust::NbClust(data = data_for_nb,
                          distance = "euclidean",
                          min.nc = min_nc,
                          max.nc = max_nc,
                          method = chosen_method,
                          index = "all")
  
  # Save the result to the global environment.
  assign("nbclust_result", res, envir = .GlobalEnv)
  
  cat("NbClust Analysis Completed.\n")
  cat("Best Number of Clusters Suggested by Each Index:\n")
  print(res$Best.nc)
  
  flush.console()  # Ensure that the console displays the output immediately.
  
  return(res)
}

################################################################################
# Manual selection of k
################################################################################
run_choose_k_method <- function() {
  # Prompt for selection of k (number of clusters)
  k <- as.integer(readline(prompt = "Enter your chosen k (number of clusters): "))
  if (!is.na(k) && k > 0) {
    chosen_k <<- k
    cat(sprintf("\nYou have chosen k = %d.\n", k))
    cat("You can now use this value of k in your clustering analysis.\n")
  } else {
    cat("Invalid input. Please enter a positive integer.\n")
  }
}

################################################################################
# Optimal G for GMM with ClusterR (diagonal covariance)
################################################################################
optimalClustersGMM_Function <- function(data) {
  # Convert data to a matrix if it is not already.
  data <- as.matrix(data)
  
  # Prompt for criterion for GMM clustering
  crit_input <- readline(prompt = "Enter criterion for GMM clustering (1 = BIC, 2 = AIC): ")
  crit_input <- as.numeric(crit_input)
  if (is.na(crit_input) || !(crit_input %in% c(1, 2))) {
    cat("Invalid input. Defaulting to BIC.\n")
    criterion <- "BIC"
  } else {
    criterion <- ifelse(crit_input == 1, "BIC", "AIC")
  }
  
  # Prompt for the Distance Metric
  dist_input <- readline(prompt = "Enter distance for GMM clustering (1 = eucl_dist, 2 = maha_dist): ")
  dist_input <- as.numeric(dist_input)
  if (is.na(dist_input) || !(dist_input %in% c(1, 2))) {
    cat("Invalid input. Defaulting to eucl_dist.\n")
    dist_mode <- "eucl_dist"
  } else {
    dist_mode <- ifelse(dist_input == 1, "eucl_dist", "maha_dist")
  }
  
  # Prompt for maximum number of clusters
  max_clusters_input <- readline(prompt = "Enter the maximum number of clusters to consider (e.g., 10): ")
  max_clusters <- as.numeric(max_clusters_input)
  if (is.na(max_clusters) || max_clusters < 1) {
    cat("Invalid input. Defaulting to 10 clusters.\n")
    max_clusters <- 10
  }
  
  # Set up a plotting layout (1x2)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  
  # Run Optimal_Clusters_GMM with chosen criterion
  cat("\nRunning Optimal_Clusters_GMM using", criterion, "and", dist_mode,
      "with max_clusters =", max_clusters, "...\n")
  criterion_values <- Optimal_Clusters_GMM(data, 
                                           max_clusters = max_clusters, 
                                           criterion = criterion, 
                                           dist_mode = dist_mode, 
                                           seed_mode = "random_subset",
                                           km_iter = 100, 
                                           em_iter = 100,
                                           var_floor = 1e-10,
                                           plot_data = TRUE, 
                                           seed = global_seed)
  
  optimal_clusters <- which.min(criterion_values)
  cat("Optimal number of clusters (based on", criterion, ") is:", optimal_clusters, "\n")
  
  chosen_G <<- optimal_clusters
  
  # Silhouette Analysis
  silhouette_values <- rep(NA, max_clusters)
  
  # Use Euclidean distance for silhouette analysis.
  dmat <- dist(data)
  
  # As silhouette analysis is undefined for a single cluster, start at k = 2.
  for (k in 2:max_clusters) {
    clusters <- tryCatch({
      tmp <- gmmWrapper_CR(x = data, k = k, full_covariance = FALSE)
      tmp$partition
    }, error = function(e) {
      cat("Error obtaining cluster assignments for", k, "clusters:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(clusters) && length(unique(clusters)) > 1) {
      sil <- silhouette(clusters, dmat)
      silhouette_values[k] <- mean(sil[, "sil_width"], na.rm = TRUE)
    } else {
      silhouette_values[k] <- NA
    }
  }
  
  if (all(is.na(silhouette_values[2:max_clusters]))) {
    cat("No valid silhouette values computed; skipping silhouette analysis.\n")
    silhouette_optimal <- NA
  } else {
    silhouette_optimal <- which.max(silhouette_values)
    cat("Optimal number of clusters (based on silhouette analysis) is:", silhouette_optimal, "\n")
    
    # Plot the silhouette values.
    plot(2:max_clusters, silhouette_values[2:max_clusters], type = "b",
         xlab = "Number of Clusters", ylab = "Average Silhouette Width",
         main = "Silhouette Analysis for GMM Clustering")
    abline(v = silhouette_optimal, col = "red", lty = 2)
  }
  
  # Reset graphics parameters to original.
  par(old_par)
  
  # Combine and Return the Results
  result_list <- list(
    criterion_optimal = optimal_clusters,
    silhouette_optimal = silhouette_optimal,
    criterion_values = criterion_values,
    silhouette_values = silhouette_values
  )
  
  return(result_list)
}

################################################################################
# Optimal G for GMM with Mclust
################################################################################
optimal_G_Mclust <- function(data) {
  # Prompt for Maximum Number of Clusters 
  max_clust_input <- readline(prompt = "Enter the maximum number of clusters to consider (e.g., 10): ")
  max_clust <- as.numeric(max_clust_input)
  if (is.na(max_clust) || max_clust < 1) {
    cat("Invalid input. Defaulting to 10 clusters.\n")
    max_clust <- 10
  }
  # Fit Mclust Model
  cat("Fitting Mclust model on data with G = 2:", max_clust, "...\n")
  
  model <<- mclust::Mclust(data, G = 2:max_clust,control = emControl(tol = 1e-4, itmax = 2000))
  
  # Retrieve the optimal number of clusters.
  optimal_G_Mc<<- model$G
  cat("Optimal number of clusters (G) from mclust is:", optimal_G, "\n")
  
  # Compute and Plot BIC Values
  cat("Computing BIC values using mclustBIC...\n")
  BIC_values <- mclust::mclustBIC(data, G = 2:max_clust)
  plot(BIC_values, main = "BIC Values for Different Numbers of Clusters")
  
  # Return a list with all results.
  return(list(model = model,
              optimal_G_Mc = optimal_G_Mc,
              BIC_values = BIC_values))
}

################################################################################
# Manual selection of G for GMM clustering
################################################################################
run_choose_G_method <- function() {
  # Prompt for G (number of clusters)
  G <- as.integer(readline(prompt = "Enter your chosen G (number of clusters): "))
  if (!is.na(G) && G > 0) {
    chosen_G <<- G  # Store the value globally
    cat(sprintf("\nYou have chosen G = %d.\n", G))
    cat("You can now use this value of G in your GMM clustering analysis.\n")
  } else {
    cat("Invalid input. Please enter a positive integer.\n")
  }
}


################################################################################
# Eigengap function for estimating k for spectral clustering
################################################################################
spectral_eigengap <- function(x,
                              n_eigs = k.max,
                              return_all = TRUE) {
  # Prepare data matrix
  data_mat <- if (is.data.frame(x)) as.matrix(x) else x
  set.seed(global_seed)
  sigma_val <- as.numeric(kernlab::sigest(data_mat, frac = 1)[2])
  # Build affinity
  data_mat <- as.matrix(x)       # samples in rows
  K        <- kernelMatrix(rbfdot(sigma = sigma_val), data_mat)
  
  # Create normalized Laplacian
  d       <- rowSums(K)
  Dinv12  <- diag(1 / sqrt(d))
  L       <- diag(nrow(K)) - Dinv12 %*% K %*% Dinv12
  
  # Get eigenvalues in increasing order
  ev      <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
  
  # Plot the eigenvalues
  idxs    <- seq_len(min(n_eigs, length(ev)))
  plot(idxs, ev[idxs], type = "b",
       xlab = "Index of eigenvalue (ascending)",
       ylab = "Eigenvalue",
       main = "Spectral Eigenvalues — look for the 'elbow'")
  
  # Compute gaps and find all maxima
  gaps    <- diff(ev[idxs])
  max_gap <- max(gaps)
  k_hats  <- which(gaps == max_gap)   # all tied indices
  
  # Print gaps to console
  cat("Gaps:\n")
  print(round(gaps, 4))
  
  if (length(k_hats) == 1) {
    cat("Largest gap between λ[", k_hats, "] and λ[", k_hats+1,
        "], suggesting k =", k_hats, "\n", sep = "")
  } else {
    cat("Multiple equal‐sized gaps found at positions:\n",
        paste0(k_hats, collapse = ", "), "\n")
    cat("Any of these indices could be a candidate k.\n")
  }
  
  if (return_all) {
    invisible(k_hats)
  } else {
    invisible(k_hats[1])
  }
}

################################################################################
# Estimate optimal K for spectral clustering (silhouette, gap-stat, Dunn)
################################################################################
optimizeSpectralK <- function(data) {
  
  # Prompt for maximum number of clusters
  k.max <- as.integer(
    readline(prompt = "Enter maximum number of clusters (k.max ≥ 2): ")
  )
  if (is.na(k.max) || k.max < 2) {
    stop("k.max must be an integer ≥ 2.")
  }
  
  # Prompt for bootstrap replicates
  nboot <- as.integer(
    readline(prompt = "Enter number of bootstrap replicates (nboot ≥ 1): ")
  )
  if (is.na(nboot) || nboot < 1) {
    stop("nboot must be an integer ≥ 1.")
  }
  
  # Prepare data matrix
  data_mat <- if (is.data.frame(data)) as.matrix(data) else data
  chosen_kernel <<- "rbfdot"
  # Estimate RBF sigma parameter and store it globally in kpar_list
  set.seed(global_seed)
  sigma_val <- as.numeric(kernlab::sigest(data_mat, frac = 1)[2])
  message("Using kernel parameter sigma = ", round(sigma_val, 4))
  kpar <- list(sigma = sigma_val)
  kpar_list <<- list(sigma = sigma_val)
  
  # Spectral clustering wrapper
  spectralFUN <- function(x, k) {
    x   <- as.matrix(x)
    res <- kernlab::specc(x, centers = k,
                          kernel = chosen_kernel,
                          kpar   = kpar)
    list(cluster = as.integer(res))
  }
  
  # Gap statistic (k ≥ 2)
  gap_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = spectralFUN,
    method     = "gap_stat",
    k.max      = k.max,
    nboot      = nboot
  )
  
  # Average silhouette width (k ≥ 2)
  sil_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = spectralFUN,
    method     = "silhouette",
    k.max      = k.max
  )
  
  # Dunn index (k ≥ 2)
  # Compute distance matrix 
  dist_mat <- dist(data_mat)
  
  # Loop over 2:k.max runing spectralFUN & extracting Dunn
  dunn_vals <- vapply(
    2:k.max,
    FUN   = function(k) {
      cls <- spectralFUN(data_mat, k)$cluster
      fpc::cluster.stats(dist_mat, cls)$dunn
    },
    numeric(1)
  )
  
  # Assemble data frame & find best k
  dunn_df   <- data.frame(k = 2:k.max, dunn = dunn_vals)
  k_dunn    <- dunn_df$k[which.max(dunn_df$dunn)]
  
  # Plot Dunn curve
  dunn_plot <- ggplot(dunn_df, aes(x = k, y = dunn)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue", size = 2) +
    geom_vline(xintercept = k_dunn, linetype = "dashed", color = "red") +
    labs(
      title = "Dunn Index vs. Number of Clusters",
      subtitle = paste0("optimal k = ", k_dunn),
      x     = "Number of clusters (k)",
      y     = "Dunn index"
    ) +
    theme_minimal()
  
  # Print diagnostics
  print(gap_plot)
  print(sil_plot)
  print(dunn_plot)
  
  # Return all three plots
  invisible(list(
    gap_statistic = gap_plot,
    silhouette    = sil_plot,
    dunn          = dunn_plot
  ))
}

################################################################################
# Estimate optimal k for k-means clustering:
################################################################################
optimizeKmeansK <- function(data) {
  
  # Prompt for maximum number of clusters
  k.max <- as.integer(
    readline(prompt = "Enter maximum number of clusters (k.max ≥ 2): ")
  )
  if (is.na(k.max) || k.max < 2) {
    stop("k.max must be an integer ≥ 2.")
  }
  
  # Prompt for number of bootstrap replicates (for Gap)
  nboot <- as.integer(
    readline(prompt = "Enter number of bootstrap replicates (nboot ≥ 1): ")
  )
  if (is.na(nboot) || nboot < 1) {
    stop("nboot must be an integer ≥ 1.")
  }
  
  # Prepare data matrix
  data_mat <- if (is.data.frame(data)) as.matrix(data) else data
  
  # Define a kmeans wrapper
  kmeansFUN <- function(x, k) {
    stats::kmeans(x, centers = k, nstart = 25)
  }
  
  # Gap Statistic
  message("▶ Computing Gap statistic (this may take some time)...")
  gap_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = kmeansFUN,
    method     = "gap_stat",
    k.max      = k.max,
    nboot      = nboot
  ) + ggtitle("Gap Statistic")
  
  # Average Silhouette Width
  sil_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = kmeansFUN,
    method     = "silhouette",
    k.max      = k.max
  ) + ggtitle("Average Silhouette")
  
  # Elbow Method (WSS)
  wss_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = kmeansFUN,
    method     = "wss",
    k.max      = k.max
  ) + ggtitle("Elbow Method (WSS)")
  
  # Print & return plots
  print(gap_plot)
  print(sil_plot)
  print(wss_plot)
  
  invisible(list(
    k.max        = k.max,
    nboot        = nboot,
    gap_plot     = gap_plot,
    silhouette   = sil_plot,
    wss_elbow    = wss_plot
  ))
}


################################################################################
# Estimate optimal k for hierarchical clustering (hclust):
################################################################################
optimizeHclustK <- function(data) {
  
  # Prompt for maximum number of clusters
  k.max <- as.integer(
    readline(prompt = "Enter maximum number of clusters (k.max ≥ 2): ")
  )
  if (is.na(k.max) || k.max < 2) {
    stop("k.max must be an integer ≥ 2.")
  }
  
  # Prompt for number of bootstrap replicates (for gap statistic)
  nboot <- as.integer(
    readline(prompt = "Enter number of bootstrap replicates (nboot ≥ 1): ")
  )
  if (is.na(nboot) || nboot < 1) {
    stop("nboot must be an integer ≥ 1.")
  }
  # Prompt for distance metric
  cat(" Distance methods for hclust:\n")
  cat(" 1. euclidian\n")
  cat(" 2. maximum\n")
  cat(" 3. manhattan\n")
  cat(" 4. canberra\n")
  cat(" 5. binary\n")
  cat(" 6. minkowski\n")
  choice <- readline(prompt = "Enter distance method of choice: ")
  
  if (choice == "1") {
    hclust_dist <<- "euclidian"
  } else if (choice == "2") {
    hclust_dist <<- "maximum"
  } else if (choice == "3") {
    hclust_dist <<- "manhattan"
  } else if (choice == "4") {
    hclust_dist <<- "canberra"
  } else if (choice == "5") {
    hclust_dist <<- "binary"
  } else if (choice == "6") {
    hclust_dist <<- "minkowski"
  }
  
  # Prompt for linkage method
  cat("Linkage criteria for hclust:\n")
  cat(" 1. ward.D\n")
  cat(" 2. ward.D2\n")
  cat(" 3. single\n")
  cat(" 4. complete\n")
  cat(" 5. average (UPGMA)\n")
  cat(" 6. mcquitty (WPGMA)\n")
  cat(" 7. median (WPGMC)\n")
  cat(" 8. centroid (UPGMC)\n")
  
  choice <- readline(prompt = "Enter linkage criteria of choice: ")
  
  if (choice == "1") {
    hclust_method <<- "ward.D"
  } else if (choice == "2") {
    hclust_method <<- "ward.D2"
  } else if (choice == "3") {
    hclust_method <<- "single"
  } else if (choice == "4") {
    hclust_method <<- "complete"
  } else if (choice == "5") {
    hclust_method <<- "average"
  } else if (choice == "6") {
    hclust_method <<- "mcquitty"
  } else if (choice == "7") {
    hclust_method <<- "median"
  } else if (choice == "8") {
    hclust_method <<- "centroid"
  }
  
  # Prepare data matrix
  data_mat <- if (is.data.frame(data)) as.matrix(data) else data
  
  # Hierarchical clustering wrapper for fviz_nbclust
  hclustFUN <- function(x, k) {
    x  <- as.matrix(x)
    hc <- hclust(dist(x, method = hclust_dist), method = hclust_method)
    list(cluster = as.integer(cutree(hc, k)))
  }
  
  # Gap Statistic (bootstrapped)
  message("▶ Computing Gap statistic (this may take some time)...")
  gap_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = hclustFUN,
    method     = "gap_stat",
    k.max      = k.max,
    nboot      = nboot
  ) + ggtitle("Gap Statistic (Hierarchical)")
  
  # Average Silhouette Width
  sil_plot <- fviz_nbclust(
    data_mat,
    FUNcluster = hclustFUN,
    method     = "silhouette",
    k.max      = k.max
  ) + ggtitle("Average Silhouette (Hierarchical)")
  
  # Elbow Method (Manual WSS for hierarchical cuts)
  n        <- nrow(data_mat)
  global_ctr <- colMeans(data_mat)
  wss_vals <- numeric(k.max)
  wss_vals[1] <- sum(rowSums((data_mat -
                                matrix(global_ctr, n, ncol(data_mat), byrow = TRUE))^2))
  
  # Precompute one hclust to avoid repeating the tree build
  hc_full <- hclust(dist(data_mat, method = hclust_dist), method = hclust_method)
  for (k in 2:k.max) {
    cls        <- cutree(hc_full, k)
    wss_vals[k] <- sum(sapply(unique(cls), function(clid) {
      pts <- data_mat[cls == clid, , drop = FALSE]
      ctr <- colMeans(pts)
      sum(rowSums((pts - 
                     matrix(ctr, nrow(pts), ncol(pts), byrow = TRUE))^2))
    }))
  }
  wss_df   <- data.frame(k = 1:k.max, wss = wss_vals)
  wss_plot <- ggplot(wss_df, aes(x = factor(k), y = wss, group = 1)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue", size = 2) +
    labs(
      title = "Elbow Method (WSS) — Hierarchical",
      x     = "Number of clusters k",
      y     = "Total within‐cluster sum of squares"
    ) +
    theme_minimal()
  
  # Print & return all three diagnostics
  print(gap_plot)
  print(sil_plot)
  print(wss_plot)
  
  invisible(list(
    k.max        = k.max,
    nboot        = nboot,
    gap_plot     = gap_plot,
    silhouette   = sil_plot,
    wss_elbow    = wss_plot
  ))
}

################################################################################
# Find optimal Lambda value (regularization parameter) for HTK-means:
################################################################################
optimal_lambda_HTKmeans <- function(data, k, plot_results = TRUE) {
  # Ensure the data is numeric.
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    data <- matrix(as.numeric(data), nrow = nrow(data))
    if (any(is.na(data))) {
      stop("Data conversion error: Some entries could not be converted to numeric values.")
    }
  }
  
  # Prompt for lambda range interactively.
  lambda_from <- as.numeric(readline(prompt = "Enter starting lambda value (From): "))
  lambda_to   <- as.numeric(readline(prompt = "Enter ending lambda value (To): "))
  lambda_by   <- as.numeric(readline(prompt = "Enter step size (By): "))
  
  if (any(is.na(c(lambda_from, lambda_to, lambda_by))) || lambda_by <= 0) {
    stop("Please provide valid numeric values for From, To, and By (with By > 0).")
  }
  lambda_seq <- seq(lambda_from, lambda_to, by = lambda_by)
  
  # Initialize storage for average silhouette and clustering results.
  silhouette_values <- numeric(length(lambda_seq))
  clustering_results <- vector("list", length(lambda_seq))
  
  # Pre-compute the distance matrix.
  d <- dist(data)
  
  for (i in seq_along(lambda_seq)) {
    this_lambda <- lambda_seq[i]
    
    # Run HTKmeans
    result <- HTKmeans(X = data, k = k, lambdas = this_lambda, standardize = FALSE)
    clustering_results[[i]] <- result
    
    # Extract cluster assignments
    cluster_assignments <- result$HTKmeans.out[[1]]$cluster
    
    # If fewer than 2 clusters, silhouette calculation is not defined.
    if (length(unique(cluster_assignments)) < 2) {
      silhouette_values[i] <- NA
    } else {
      sil_obj <- cluster::silhouette(cluster_assignments, d)
      if (is.matrix(sil_obj)) {
        silhouette_values[i] <- mean(sil_obj[, 3])
      } else {
        silhouette_values[i] <- NA
      }
    }
  }
  
  # Determine the lambda that gives the maximum average silhouette.
  best_idx <- which.max(silhouette_values)
  opt_lambda <- lambda_seq[best_idx]
  
  # Save the optimal lambda to the global environment
  opt_lambda_global <<- opt_lambda
  
  # Plotting of results
  if (plot_results) {
    # Create a data frame for plotting.
    df_plot <- data.frame(lambda = lambda_seq, silhouette = silhouette_values)
    
    # Build the ggplot.
    p <- ggplot(df_plot, aes(x = lambda, y = silhouette)) +
      geom_line() +                 
      geom_point() +                
      # Highlight the optimal lambda point.
      geom_point(data = data.frame(lambda = opt_lambda, 
                                   silhouette = silhouette_values[best_idx]),
                 aes(x = lambda, y = silhouette),
                 colour = "red", size = 3) +
      # Add a vertical dashed line at opt_lambda.
      geom_vline(xintercept = opt_lambda, colour = "red", linetype = "dashed") +
      labs(x = "Lambda", y = "Average Silhouette",
           title = "Optimal Lambda based on Average Silhouette") +
      scale_x_continuous(breaks = pretty(lambda_seq))
    
    # Print the ggplot.
    print(p)
  }
  
  
  # Return the computed values.
  return(list(optimal_lambda    = opt_lambda,
              silhouette_values = silhouette_values,
              lambda_seq        = lambda_seq,
              clustering_results = clustering_results))
}

################################################################################
# Optimal eps DBSCAN
################################################################################
optimal_eps_DBSCAN <- function(data, minPts = NULL, plot_results = TRUE) {
  # Ensure data variable (select.npx) is a numeric matrix.
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    data <- matrix(as.numeric(data), nrow = nrow(data))
    if (any(is.na(data))) {
      stop("Data conversion error: Some entries could not be converted to numeric values.")
    }
  }
  
  minPts <- readline(prompt = "Enter minPts value (common rule of thumb = features + 1)")
  # If minPts is not provided, use the common rule-of-thumb: 2 x number of features.
  if (is.null(minPts)) {
    minPts <- ncol(data) +1
    message("minPts not provided. Using default minPts = number of features + 1 = ", minPts)
  }
  
  # Compute the kNN distances (k = minPts) for all points.
  dists <- dbscan::kNNdist(data, k = minPts)
  
  # Sort the distances (ascending).
  sorted_dists <- sort(dists)
  n <- length(sorted_dists)
  
  # Automatic Elbow Detection using the geometric "maximal distance" method
  # Define the first and last points of the sorted distance curve.
  x1 <- 1
  y1 <- sorted_dists[1]
  x2 <- n
  y2 <- sorted_dists[n]
  x_ind <- seq_len(n)
  
  # Compute perpendicular distance from each point to the line
  numerator <- abs((x2 - x1) * (y1 - sorted_dists) - (x1 - x_ind) * (y2 - y1))
  denominator <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  distances_to_line <- numerator / denominator
  
  # The index with the maximum distance the elbow point.
  opt_index <- which.max(distances_to_line)
  opt_eps <- sorted_dists[opt_index]
  
  # Save opt_eps to global environment
  opt_eps_global <<- opt_eps
  cat("Optimal eps: ", opt_eps)
  
  # Plot the kNN distance plot and mark the optimal eps.
  if (plot_results) {
    
    
    # Plot the kNN distances using the built-in function.
    dbscan::kNNdistplot(data, k = minPts)
    title(main = paste("kNN Distance Plot (minPts =", minPts, ")"))
    
    # Add a horizontal red dashed line at the optimal eps.
    abline(h = opt_eps, col = "red", lty = 2)
  }
  
  # Return the optimal eps, chosen minPts, and the raw distances.
  return(list(opt_eps = opt_eps, minPts = minPts, distances = dists))
}

################################################################################
# Run DBSCAN to implement chosen parameters into pipeline:
################################################################################
run_dbscan_method <- function() {
  # Prompt for the DBSCAN eps parameter
  dbscan.eps <<- as.numeric(readline(prompt = "Enter eps value: "))
  assign("dbscan.eps", dbscan.eps, envir = .GlobalEnv)
  if (is.na(dbscan.eps) || dbscan.eps <= 0) {
    cat("Invalid eps value. Please enter a positive numeric value.\n")
    return(NULL)
  }
  
  # Prompt for the DBSCAN minPts parameter
  dbscan.minPts <<- as.integer(readline(prompt = "Enter minPts value: "))
  assign("dbscan.minPts", dbscan.minPts, envir = .GlobalEnv)
  if (is.na(dbscan.minPts) || dbscan.minPts <= 0) {
    cat("Invalid minPts value. Please enter a positive integer.\n")
    return(NULL)
  }
  
  cat(sprintf("Running DBSCAN with eps = %f and minPts = %d...\n", dbscan.eps, dbscan.minPts))
  
  if (!exists("select.npx")) {
    cat("Dataset 'select.npx' is not available in the environment.\n")
    return(NULL)
  }
  if (!exists("select.sinfo")) {
    cat("Dataset 'select.sinfo' is not available in the environment.\n")
    return(NULL)
  }
  
  # Run DBSCAN using specified eps and minPts values
  dbs <- dbscan::dbscan(select.npx, eps = dbscan.eps, minPts = dbscan.minPts)
  
  dbscan.clusters <- dbs$cluster
  names(dbscan.clusters) <- rownames(select.npx)
  dbscan.clusters <- as.factor(dbscan.clusters)
  
  # Adjust the factor levels by adding 1 so that cluster labels start at 1
  levels(dbscan.clusters) <- as.character(as.numeric(levels(dbscan.clusters)) + 1)
  
  # Calculate the highest (numeric) cluster value found
  numeric_clusters <- as.numeric(as.character(dbscan.clusters))
  max_cluster <- max(numeric_clusters, na.rm = TRUE)
  cat(sprintf("DBSCAN clustering completed. Highest cluster value found: %d\n", max_cluster))
}

################################################################################
# Optimal minPts for HDBSCAN (Hierarchical DBSCAN)
################################################################################
optimal_minPts_hdbscan <- function(data, minPts_seq = NULL, plot_results = TRUE) {
  # Ensure numeric matrix
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    data <- matrix(as.numeric(data), nrow = nrow(data))
    if (any(is.na(data))) {
      stop("Data conversion error: Some entries could not be converted to numeric values.")
    }
  }
  
  # Prompt for minPts sequence if minPts_seq not provided
  if (is.null(minPts_seq)) {
    cat("\nEnter a sequence of minPts values (e.g. 5:20 or c(5,10,15,20)):\n")
    input <- readline("minPts_seq = ")
    minPts_seq <- eval(parse(text = input))
    if (!is.numeric(minPts_seq)) stop("minPts_seq must evaluate to numeric values.")
  }
  
  silhouette_values <- numeric(length(minPts_seq))
  clustering_results <- vector("list", length(minPts_seq))
  
  for (i in seq_along(minPts_seq)) {
    cur_minPts <- minPts_seq[i]
    res <- dbscan::hdbscan(data, minPts = cur_minPts)
    clustering_results[[i]] <- res
    
    labels <- res$cluster
    valid_idx <- labels > 0   # exclude noise
    
    if (length(unique(labels[valid_idx])) < 2) {
      silhouette_values[i] <- NA
    } else {
      d_matrix <- dist(data[valid_idx, , drop = FALSE])
      sil <- cluster::silhouette(labels[valid_idx], d_matrix)
      silhouette_values[i] <- mean(sil[, 3])
    }
  }
  
  # Pick best
  if (all(is.na(silhouette_values))) {
    warning("No valid clustering (at least two clusters) was found for any candidate minPts value.")
    best_idx <- NA
    optimal_minPts <- NA
  } else {
    best_idx <- which.max(silhouette_values)
    optimal_minPts <- minPts_seq[best_idx]
  }
  
  opt_minPts_HDBSCAN <<- optimal_minPts
  
  if (plot_results) {
    df_plot <- data.frame(minPts = minPts_seq, silhouette = silhouette_values)
    
    p <- ggplot(df_plot, aes(x = minPts, y = silhouette)) +
      geom_line() +
      geom_point() +
      { if (!is.na(optimal_minPts)) 
        geom_point(data = data.frame(minPts = optimal_minPts, silhouette = silhouette_values[best_idx]),
                   aes(x = minPts, y = silhouette),
                   colour = "red", size = 3) } +
      labs(title = "Optimal minPts for HDBSCAN", x = "minPts", y = "Average Silhouette") +
      theme_minimal()
    
    print(p)
  }
  
  return(list(optimal_minPts = optimal_minPts,
              silhouette_values = silhouette_values,
              minPts_seq = minPts_seq,
              clustering_results = clustering_results))
}


################################################################################
# Run and test HDBSCAN to implement optimized parameters into pipeline
################################################################################
run_hdbscan_method <- function() {
  # Prompt for HDBSCAN parameters
  hdbscan.minPts <<- as.integer(readline(prompt = "Enter minPts value: "))
  assign("hdbscan.minPts", hdbscan.minPts, envir = .GlobalEnv)
  
  hdbscan.min_cluster_size <<- as.integer(readline(prompt = "Enter min_cluster_size value: "))
  assign("hdbscan.min_cluster_size", hdbscan.min_cluster_size, envir = .GlobalEnv)
  
  if (is.na(hdbscan.minPts) || hdbscan.minPts <= 0) {
    cat("Invalid minPts value. Please enter a positive integer.\n")
    return(NULL)
  }
  if (is.na(hdbscan.min_cluster_size) || hdbscan.min_cluster_size <= 0) {
    cat("Invalid min_cluster_size value. Please enter a positive integer.\n")
    return(NULL)
  }
  
  cat(sprintf("Running Hierarchical DBSCAN with minPts = %d and min_cluster_size = %d...\n", 
              hdbscan.minPts, hdbscan.min_cluster_size))
  
  # Check datasets
  if (!exists("select.npx")) {
    cat("Dataset 'select.npx' is not available in the environment.\n")
    return(NULL)
  }
  if (!exists("select.sinfo")) {
    cat("Dataset 'select.sinfo' is not available in the environment.\n")
    return(NULL)
  }
  
  # Run HDBSCAN
  hdbs <- dbscan::hdbscan(select.npx, 
                          minPts = hdbscan.minPts, 
                          min_cluster_size = hdbscan.min_cluster_size)
  
  hdbscan.clusters <- hdbs$cluster
  names(hdbscan.clusters) <- rownames(select.npx)
  hdbscan.clusters <- as.factor(hdbscan.clusters)
  
  # Adjust factor levels so clusters start at 1 (noise stays 0)
  levels(hdbscan.clusters) <- as.character(as.numeric(levels(hdbscan.clusters)))
  
  # Highest cluster label
  numeric_clusters <- as.numeric(as.character(hdbscan.clusters))
  max_cluster <- max(numeric_clusters, na.rm = TRUE)
  cat(sprintf("Hierarchical DBSCAN clustering completed. Highest cluster value found: %d\n", max_cluster))
  
  return(hdbscan.clusters)
}

