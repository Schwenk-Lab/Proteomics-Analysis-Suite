# ------------------------------------------------------------------------------
# Script:        VisualizationFunctions.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for visualization of clustering results
# ------------------------------------------------------------------------------

################################################################################
# Runner/ Visualization Menu
################################################################################

visualizationRunner <- function(reproduce = NULL) {
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      clusterMethod = NA,
      drMethod = NA,
      plotDim = NA,
      saveChoice = NA,
      plotVarName = NA,
      stringsAsFactors = FALSE
    )
  cat("### Visualization Runner ###\n\n")
  
  # Choose the clustering method
  clust_method <<- choose_clustering_method()
  if (is.null(clust_method)) return(invisible(NULL))
  assign("clust_method", clust_method, envir = .GlobalEnv)
  clusters <- NULL
  reproduce$clusterMethod <- clust_method
  
  
  # Choose the DR method.
  dr_method <- choose_dr_method()
  if (is.null(dr_method)) return(invisible(NULL))
  assign("dr_method", dr_method, envir = .GlobalEnv)
  cat("Proceeding with clustering method:", clust_method, 
      "and DR method:", dr_method, "\n\n")
  reproduce$drMethod <- dr_method
  } else {
    clust_method <- reproduce$clustMethod
    dr_method <- reproduce$drMethod
  }
  
  if (clust_method != "ALL") {
    clustersX <<- switch(clust_method,
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
    if (is.null(clustersX)) {
      cat("Error: Cluster assignments not found in select.sinfo.\n")
      return(invisible(NULL))
    }
  }
  # Visualize the data using the chosen DR method.
  if (dr_method == "PCA") {
    visualizePCA()
  } else if (dr_method == "Sparse PCA") {
    visualizeSparsePCA()
  } else if (dr_method == "PLS-DA") {
    visualizePLSDA(reproduce = reproduce)
  } else if (dr_method == "Sparse PLS-DA") {
    visualizeSPLSDA()
  } else if (dr_method == "Kernel PCA") {
    visualizeKernelPCA()
  } else if (dr_method == "UMAP") {
    visualizeUMAP(reproduce = reproduce)
  } else if (dr_method == "t-SNE") {
    visualizeTSNE()
  } else {
    cat("Selected DR method is not recognized.\n")
  } 
  if (is.na(reproduce$saveChoice)) {
    print(p)
    save_choice <- readline(prompt ="Do you want to save the plot? y/n: \n")
    reproduce$saveChoice <- save_choice
    if (save_choice == "y") {
      var_name2 <- readline("Enter a name to save your plot object under: ")
      reproduce$plotVarName <- var_name2
      assign(var_name2, p, .GlobalEnv ) 
      cat("\nVisualization complete. The plot has been saved globally as '", var_name2, "'.\n", sep = "")
    }
  } else {
    var_name2 <- reproduce$plotVarName
    assign(var_name2, p, .GlobalEnv ) 
    print(p)
  }
}

################################################################################
# PCA visualization function
################################################################################

visualizePCA <- function() {
  # Check that the DR result exists.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before PCA.\n")
    return(invisible(NULL))
  }
  
  cat("\nPerforming PCA visualization...\n")
  
  # If unreduced.ptx exists, but the select.ptx is same as unreduced.ptx, rerun PCA
  run_pca <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_pca <- TRUE
    }
  }
  
  if (run_pca) {
    # Re-run PCA without centering or scaling
    pca_result <- prcomp(select.ptx, center = FALSE, scale. = FALSE)
    df_plot <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])
  } else {
    # Use the first two columns from select.ptx corresponding to PC1 and PC2
    df_plot <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
  }
  
  # Determine axis labels.
  if (run_pca) {
    label_x <- paste0("PC1: ", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 2), " %")
    label_y <- paste0("PC2: ", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 2), " %")
  } else {
    # Use the column names from the DR result if available, otherwise use defaults.
    label_x <- if(!is.null(colnames(df_plot))) colnames(df_plot)[1] else "Dim1"
    label_y <- if(!is.null(colnames(df_plot))) colnames(df_plot)[2] else "Dim2"
  }
  
  # Add the clustering information.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    cluster_factor <- as.factor(clustersX)
    df_plot$Cluster <- cluster_factor
  }
  # Plot using ggplot2
  if (use_cluster) {
    p <<- ggplot(df_plot, aes(x = df_plot[,1], y = df_plot[,2], color = Cluster)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "PCA Plot Colored by Selected Clusters",
           x = label_x,
           y = label_y,
           color = "Clusters") +
      theme_minimal()
  } else {
    p <<- ggplot(df_plot, aes(x = df_plot[,1], y = df_plot[,2])) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "PCA Plot",
           x = label_x,
           y = label_y) +
      theme_minimal()
  }
  return(p)
}

################################################################################
# Sparse PCA visualization function
################################################################################

visualizeSparsePCA <- function() {
  cat("\nPerforming Sparse PCA visualization...\n")
  
  # If unreduced.ptx and saved parameters from SPCA exists, 
  # but the select.ptx is same as unreduced.ptx, rerun SPCA with saved parameters
  run_spca <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_spca <- TRUE
    }
  }
  
  if (run_spca) {
    if (!exists("ncomp_input") || !exists("keepX_vals")) {
      cat("Missing required global parameters (ncomp_input, keepX_vals). Cannot re-run SPCA.\n")
      return(invisible(NULL))
    }
    spca_result <- mixOmics::spca(select.ptx,
                                  ncomp = ncomp_input,
                                  keepX = keepX_vals,
                                  center = FALSE,
                                  scale = FALSE)
    df_plot <- data.frame(SPCA1 = spca_result$x[, 1],
                          SPCA2 = spca_result$x[, 2])
  } else {
    # Use the first two columns from select.ptx (SPCA1 and SPCA2)
    df_plot <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
    colnames(df_plot)[1:2] <- c("SPCA1", "SPCA2")
  }
  
  # Axis labels
  label_x <- colnames(df_plot)[1]
  label_y <- colnames(df_plot)[2]
  
  # Add the clustering information.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    cluster_factor <- as.factor(clustersX)
    df_plot$Cluster <- cluster_factor
  }
  
  if (use_cluster) {
    p <<- ggplot(df_plot, aes(x = df_plot[,1], y = df_plot[,2], color = Cluster)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "Sparse PCA Plot Colored by Selected Clusters",
           x = label_x,
           y = label_y,
           color = "Clusters") +
      theme_minimal()
  } else {
    p <<- ggplot(df_plot, aes(x = df_plot[,1], y = df_plot[,2])) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "Sparse PCA Plot",
           x = label_x,
           y = label_y) +
      theme_minimal()
  }
  
  return(p)
}

################################################################################
# Visualize PLS-DA function:
################################################################################

visualizePLSDA <- function(reproduce = NULL) {
  if (is.na(reproduce$plotDim)) {
    cat("\nPerforming PLS-DA visualization...\n")
    
    # Prompt for selection between 2D or 3D plot
    plotDim <- readline(prompt = "Enter 2 for 2D plot or 3 for 3D plot: ")
    plotDim <- as.numeric(plotDim)
    if (is.na(plotDim) || !(plotDim %in% c(2, 3))) {
      cat("\nInvalid selection, defaulting to 2D plot.\n")
      plotDim <- 2
    }
    reproduce$plotDim <- plotDim
  } else {
    plotDim <- reproduce$plotDim
  }
  
  # Determine whether to re-run PLS-DA using the original data.
  run_plsda <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_plsda <- TRUE
    }
  }
  
  if (run_plsda) {
    # Ensure that necessary global parameters for PLS-DA are available.
    if (!exists("plsda.ncomp") || !exists("temp.plsda.sinfo")) {
      cat("Missing required global parameters (ncomp_input, temp.plsda.sinfo) for PLS-DA.\n")
      return(invisible(NULL))
    }
    ncomp_input <- plsda.ncomp
    # Check if enough components are available
    if (plotDim == 2 && ncomp_input < 2) {
      cat("Error: PLS-DA requires at least two components for a 2D plot.\n")
      return(invisible(NULL))
    } else if (plotDim == 3 && ncomp_input < 3) {
      cat("Error: PLS-DA requires at least three components for a 3D plot.\n")
      return(invisible(NULL))
    }
    
    # Run PLS-DA with options as needed.
    plsda_model <- mixOmics::plsda(select.ptx, temp.plsda.sinfo, ncomp = ncomp_input)
    pls_scores_full <- plsda_model$variates$X[, 1:ncomp_input, drop = FALSE]
    
    if (plotDim == 2) {
      df_plot <- data.frame(Comp1 = pls_scores_full[, 1],
                            Comp2 = pls_scores_full[, 2])
    } else {  # If 3D plot
      df_plot <- data.frame(Comp1 = pls_scores_full[, 1],
                            Comp2 = pls_scores_full[, 2],
                            Comp3 = pls_scores_full[, 3])
    }
    
  } else {
    # Without re-running PLS-DA, use the components from select.ptx.
    if (plotDim == 2) {
      df_plot <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
      colnames(df_plot)[1:2] <- c("Comp1", "Comp2")
    } else {  # If 3D plot
      if (ncol(select.ptx) < 3) {
        cat("Error: 'select.ptx' does not have three columns required for a 3D plot.\n")
        return(invisible(NULL))
      }
      df_plot <- as.data.frame(select.ptx[, 1:3, drop = FALSE])
      colnames(df_plot)[1:3] <- c("Comp1", "Comp2", "Comp3")
    }
  }
  
  # Attach the clustering information.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    df_plot$Cluster <- as.factor(clustersX)
  }
  
  # Branch based on the plot dimension.
  if (plotDim == 2) {
    # 2D plot using ggplot2
    if (use_cluster) {
      p <<- ggplot(df_plot, aes(x = Comp1, y = Comp2, color = Cluster)) +
        geom_point(size = 3) +
        stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
        labs(title = if (exists("clust_method", envir = .GlobalEnv)) {
          paste(get("clust_method", envir = .GlobalEnv), "PLS-DA Plot")
        } else {
          "PLS-DA Plot"
        },
        x = "Component 1",
        y = "Component 2",
        color = "Clusters") +
        theme_minimal()
    } else {
      p <<- ggplot(df_plot, aes(x = Comp1, y = Comp2)) +
        geom_point(size = 3) +
        stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
        labs(title = "PLS-DA Plot",
             x = "Component 1",
             y = "Component 2") +
        theme_minimal()
    }
  } else { # If 3D plot
    if (use_cluster) {
      p <<- plot_ly(data = df_plot, 
                   x = ~Comp1, 
                   y = ~Comp2, 
                   z = ~Comp3,
                   color = ~Cluster,
                   colors = "Set1",
                   type = "scatter3d",
                   mode = "markers") %>%
        layout(title = if (exists("clust_method", envir = .GlobalEnv)) {
          paste(get("clust_method", envir = .GlobalEnv), "PLS-DA 3D Plot")
        } else {
          "PLS-DA 3D Plot"
        },
        scene = list(
          xaxis = list(title = "Component 1"),
          yaxis = list(title = "Component 2"),
          zaxis = list(title = "Component 3")
        ))
    } else {
      p <<- plot_ly(data = df_plot, 
                   x = ~Comp1, 
                   y = ~Comp2, 
                   z = ~Comp3,
                   type = "scatter3d",
                   mode = "markers") %>%
        layout(title = "PLS-DA 3D Plot",
               scene = list(
                 xaxis = list(title = "Component 1"),
                 yaxis = list(title = "Component 2"),
                 zaxis = list(title = "Component 3")
               ))
    }
    
  }
  return(p)
}

################################################################################
# Visualize sPLS-DA function
################################################################################

visualizeSPLSDA <- function() {
  cat("\nPerforming sPLS-DA visualization...\n")
  
  # Determine whether to re-run sPLS-DA by comparing column names.
  run_splsda <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_splsda <- TRUE
    }
  }
  # Reproduce from previously stored parameters if 
  # select.ptx has been changed by another DR method
  if (run_splsda) {
    # Check required global parameters.
    if (!exists("splsda.ncomp") || !exists("splsda.keepX") || !exists("temp.splsda.sinfo")) {
      cat("Missing required global parameters (splsda.ncomp, splsda.keepX, temp.splsda.sinfo) for sPLS-DA.\n")
      return(invisible(NULL))
    }
    # Run sPLS-DA
    ncomp_input <- splsda.ncomp
    keepX <- splsda.keepX
    splsda_model <- mixOmics::splsda(select.ptx, temp.splsda.sinfo, 
                                     ncomp = ncomp_input, keepX = keepX_vals)
    spls_scores_full <- splsda_model$variates$X[, 1:ncomp_input, drop = FALSE]
    if (ncol(spls_scores_full) < 2) {
      cat("Error: sPLS-DA produced less than two components.\n")
      return(invisible(NULL))
    }
    df_plot <- data.frame(Comp1 = spls_scores_full[, 1], 
                          Comp2 = spls_scores_full[, 2])
  } else {
    # Use the first two columns of select.ptx.
    df_plot <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
    colnames(df_plot)[1:2] <- c("Comp1", "Comp2")
  }
  
  # Attach the clustering information.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    df_plot$Cluster <- as.factor(clustersX)
  }
  
  # Create the plot using ggplot2.
  if (use_cluster) {
    p <<- ggplot(df_plot, aes(x = Comp1, y = Comp2, color = Cluster)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "sPLS-DA Scores (Comp1 vs Comp2)",
           x = "Component 1",
           y = "Component 2",
           color = "Clusters") +
      theme_minimal()
  } else {
    p <<- ggplot(df_plot, aes(x = Comp1, y = Comp2)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "sPLS-DA Scores (Comp1 vs Comp2)",
           x = "Component 1",
           y = "Component 2") +
      theme_minimal()
  }
  return(p)
}

################################################################################
# Visualize Kernel PCA
################################################################################

visualizeKernelPCA <- function() {
  cat("\nPerforming Kernel PCA visualization...\n")
  # Determine whether to run kernel PCA based on matching column names.
  run_kpca <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_kpca <- TRUE
    }
  }
  
  if (run_kpca) {
    # Ensure that kernel PCA parameters exist.
    if (!exists("reproduceKernelPCA")) {
      cat("Error: Global parameters 'sigma' and/or 'features' are not set.\n")
      return(invisible(NULL))
    }
    sigmaX <- reproduceKernelPCA$KPCAsigmaX
    features <- reproduceKernelPCA$KPCAnFeatures
    # Run kernel PCA using an RBF kernel.
    kpca_result <- kernlab::kpca(as.matrix(select.ptx),
                                 kernel = "rbfdot",
                                 kpar = list(sigma = sigmaX),
                                 features = features)
    # Extract the rotated (kernel PC) coordinates.
    pc_coords <- kernlab::rotated(kpca_result)
    if (ncol(pc_coords) < 2) {
      cat("Error: Kernel PCA produced less than two components.\n")
      return(invisible(NULL))
    }
    df_dr <- data.frame(PC1 = pc_coords[, 1],
                        PC2 = pc_coords[, 2])
  } else {
    # Fallback: Use the first two columns of select.ptx.
    df_dr <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
    colnames(df_dr)[1:2] <- c("PC1", "PC2")
  }
  
  # Attach the clustering assignment.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    df_dr$Cluster <- as.factor(clustersX)
  }
  
  if (use_cluster) {
    p <<- ggplot(df_dr, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "Kernel PCA of Expression Matrix",
           x = "Kernel PC1", y = "Kernel PC2",
           color = "Clusters") +
      theme_minimal()
  } else {
    p <<- ggplot(df_dr, aes(x = PC1, y = PC2)) +
      geom_point(size = 3, color = "blue") +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "Kernel PCA of Expression Matrix",
           x = "Kernel PC1", y = "Kernel PC2") +
      theme_minimal()
  }
  return(p)
}
################################################################################ 
# UMAP visualization function
################################################################################ 
visualizeUMAP <- function(reproduce = NULL) {
  # Check that the dimensionality reduction (DR) result exists.
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before running UMAP.\n")
    return(invisible(NULL))
  }
  if (is.na(reproduce$plotDim)) {
    cat("\nPerforming UMAP visualization...\n")
    
    # Ask the user if they want a 2D or 3D plot.
    umapDim <- readline(prompt = "Enter 2 for 2D UMAP plot or 3 for 3D UMAP plot: ")
    umapDim <- as.numeric(umapDim)
    if (is.na(umapDim) || !(umapDim %in% c(2, 3))) {
      cat("\nInvalid selection; defaulting to 2D plot.\n")
      umapDim <- 2
    }
    reproduce$plotDim <- umapDim
  } else {
    umapDim <- reproduce$plotDim
  }
 
  # Determine whether to run UMAP based on matching column names.
  run_umap <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_umap <- TRUE
    }
  }
  
  if (run_umap) {
    # Check for required UMAP parameters
    # Use ncomp = 2 for 2D or ncomp = 3 for 3D if no global parameter is available.
    if (!exists("reproduceUMAP")) {
      stop("Global variable 'reproduceUMAP' not found, no saved parameters detected .Please run UMAP first")
    }
    ncomp <- reproduceUMAP$ncomp
    n_neighbors = reproduceUMAP$nNeighbors
    min_dist = reproduceUMAP$minDist
    seed_val = reproduceUMAP$seed
    
    # Run UMAP with the specified parameters.
    umap_result <- uwot::umap(as.matrix(select.ptx),
                              n_components = ncomp,
                              n_neighbors = n_neighbors,
                              min_dist = min_dist,
                              metric = "euclidean",
                              seed = seed_val)
    
    # Create a data frame containing UMAP components.
    if (umapDim == 2) {
      df_umap <- data.frame(UMAP1 = umap_result[, 1],
                            UMAP2 = umap_result[, 2])
    } else {
      if (ncol(umap_result) < 3) {
        cat("Error: UMAP did not produce 3 components, cannot generate a 3D plot.\n")
        return(invisible(NULL))
      }
      df_umap <- data.frame(UMAP1 = umap_result[, 1],
                            UMAP2 = umap_result[, 2],
                            UMAP3 = umap_result[, 3])
    }
    
  } else {
    # Fallback: Use the first 2 or 3 columns of select.ptx.
    if (umapDim == 2) {
      if (ncol(select.ptx) < 2) {
        cat("Error: 'select.ptx' does not have two columns required for a 2D plot.\n")
        return(invisible(NULL))
      }
      df_umap <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
      colnames(df_umap)[1:2] <- c("UMAP1", "UMAP2")
    } else {
      if (ncol(select.ptx) < 3) {
        cat("Error: 'select.ptx' does not have three columns required for a 3D plot.\n")
        return(invisible(NULL))
      }
      df_umap <- as.data.frame(select.ptx[, 1:3, drop = FALSE])
      colnames(df_umap)[1:3] <- c("UMAP1", "UMAP2", "UMAP3")
    }
  }
  
  # Attach the clustering assignment.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_clusters <- FALSE
  } else {
    use_clusters <- TRUE
    df_umap$Cluster <- as.factor(clustersX)
  }
  
  # Plot based on the desired dimensionality.
  if (umapDim == 2) {
    if (use_clusters) {
      p <<- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 3) +
        labs(title = "UMAP Projection",
             x = "UMAP1", y = "UMAP2",
             color = "Clusters") +
        theme_minimal()
    } else {
      p <<- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
        geom_point(size = 3, color = "blue") +
        labs(title = "UMAP Projection",
             x = "UMAP1", y = "UMAP2") +
        theme_minimal()
    }
    print(p)
    final_plot <- p
  } else {
    # 3D plot using plotly.
    if (use_clusters) {
      p <<- plot_ly(data = df_umap, 
                   x = ~UMAP1, 
                   y = ~UMAP2, 
                   z = ~UMAP3, 
                   color = ~Cluster,
                   colors = "Set1",
                   type = "scatter3d", 
                   mode = "markers") %>%
        layout(title = "UMAP 3D Projection",
               scene = list(
                 xaxis = list(title = "UMAP1"),
                 yaxis = list(title = "UMAP2"),
                 zaxis = list(title = "UMAP3")
               ))
    } else {
      p <<- plot_ly(data = df_umap, 
                   x = ~UMAP1, 
                   y = ~UMAP2, 
                   z = ~UMAP3, 
                   type = "scatter3d", 
                   mode = "markers") %>%
        layout(title = "UMAP 3D Projection",
               scene = list(
                 xaxis = list(title = "UMAP1"),
                 yaxis = list(title = "UMAP2"),
                 zaxis = list(title = "UMAP3")
               ))
    }
  }
  return(p)
}

################################################################################
# Visualization function for t-SNE:
################################################################################

visualizeTSNE <- function() {
  # Check that the DR result exists
  if (!exists("select.ptx") || is.null(select.ptx)) {
    cat("\nError: 'select.ptx' is not available. Ensure data is loaded before running t-SNE.\n")
    return(invisible(NULL))
  }
  
  cat("\nPerforming t-SNE visualization...\n")
  
  # If unreduced.ptx exists, but the select.ptx is same as unreduced.ptx, rerun t-SNE
  run_tsne <- FALSE
  if (exists("unreduced.ptx") && !is.null(unreduced.ptx)) {
    if (identical(colnames(select.ptx), colnames(unreduced.ptx))) {
      run_tsne <- TRUE
    }
  }
  
  if (run_tsne) {
    # Ensure required global parameters exist for t-SNE: dims, perplexity, theta, max_iter.
    if (!exists("reproduceTSNE")) {
      stop("Global variable 'reproduceTSNE' not found, no saved parameters detected .Please run t-SNE first")
    }
    # Assign parameters from reproduction variable
    dims2 <- reproduceTSNE$dims
    perplexity <- reproduceTSNE$perplexity
    theta <- reproduceTSNE$theta
    max_iter <- reproduceTSNE$maxIter
    set.seed(reproduceTSNE$seed)
    
    tsne_result <- Rtsne::Rtsne(as.matrix(select.ptx),
                                dims = dims2,
                                perplexity = perplexity,
                                theta = theta,
                                max_iter = max_iter,
                                verbose = TRUE,
                                check_duplicates = FALSE)
    Y <- tsne_result$Y
    df_tsne <- data.frame(TSNE1 = Y[, 1], TSNE2 = Y[, 2])
  } else {
    # Fallback: Use the first two columns of select.ptx.
    df_tsne <- as.data.frame(select.ptx[, 1:2, drop = FALSE])
    colnames(df_tsne)[1:2] <- c("TSNE1", "TSNE2")
  }
  
  # Attach the clustering assignment.
  if (!exists("clustersX") || is.null(clustersX)) {
    cat("\nWarning: Global variable 'clustersX' not found; plotting without cluster colors.\n")
    use_cluster <- FALSE
  } else {
    use_cluster <- TRUE
    df_tsne$Cluster <- as.factor(clustersX)
  }
  
  if (use_cluster) {
    p <<- ggplot(df_tsne, aes(x = TSNE1, y = TSNE2, color = Cluster)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "t-SNE Projection",
           x = "t-SNE 1", y = "t-SNE 2",
           color = "Clusters") +
      theme_minimal()
  } else {
    p <<- ggplot(df_tsne, aes(x = TSNE1, y = TSNE2)) +
      geom_point(size = 3, color = "blue") +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
      labs(title = "t-SNE Projection",
           x = "t-SNE 1", y = "t-SNE 2") +
      theme_minimal()
  }
  return(p)
}
