# ------------------------------------------------------------------------------
# Script:        WrapperFuncs.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Wrapper functions for clusterboot
# ------------------------------------------------------------------------------

################################################################################
# HTKmeans
################################################################################
htkmeansWrapper <- function(x, k) {
  
  # Run HTKmeans
  res <- HTKmeans(x, k = k, standardize = FALSE, lambdas = opt_lambda_global, iter.max = 100)
  
  # Debug: Use the first lambda result and check that it contains the cluster element
  lambdaResult <- res$HTKmeans.out[[1]]
  if (!("cluster" %in% names(lambdaResult))) {
    stop("The HTKmeans result (from the first lambda) does not contain a 'cluster' element.")
  }
  
  # Extract the cluster assignment vector
  assign <- lambdaResult$cluster
  
  # Create a list of logical vectors indicating membership (one per cluster)
  clusterlist <- lapply(1:k, function(i) assign == i)
  
  # Return the result in a list
  list(
    partition   = assign,
    nc          = k,
    clusterlist = clusterlist
  )
}

################################################################################
# FKM
################################################################################

# Wrapper function for FKM
fkmWrapper <- function(x, k) {
  # Run the fuzzy k-means clustering using the default fuzzifier parameter
  res <- FKM(x, k = k)
  
  # Convert the fuzzy membership matrix into hard assignments:
  assign <- apply(res$U, 1, which.max)
  
  # Create clusterlist
  clusterlist <- lapply(1:k, function(i) assign == i)
  
  # Return the list in the format expected by clusterboot.
  return(list(partition = assign, nc = k, clusterlist = clusterlist))
}

################################################################################
# HDBSCAN
################################################################################
hdbscanWrapper <- function(x, k) {
  # Use global minPts
  minPts_val <- hdbscan.minPts
  min_cluster_size_val <- if (exists("hdbscan.min_cluster_size")) hdbscan.min_cluster_size else minPts_val
  
  # Run HDBSCAN
  res <- dbscan::hdbscan(x, minPts = minPts_val, min_cluster_size = min_cluster_size_val)
  
  # Extract cluster assignments
  partition <- as.factor(res$cluster)
  
  # Keep noise as 0, clusters as 1..N
  uniqueClusters <- sort(unique(res$cluster))
  
  # Number of clusters (excluding noise)
  nc <- length(uniqueClusters[uniqueClusters != 0])
  
  # Build cluster list (exclude noise)
  clusterlist <- lapply(uniqueClusters[uniqueClusters != 0], function(i) res$cluster == i)
  
  # Return in format expected by clusterboot
  return(list(partition = partition, nc = nc, clusterlist = clusterlist))
}


################################################################################
# DBSCAN
################################################################################
dbscanWrapper <- function(x, k) {
  # Set DBSCAN parameters 
  eps_val <- dbscan.eps
  minPts_val <- dbscan.minPts
  
  # Run DBSCAN
  res <- dbscan::dbscan(x, eps = eps_val, minPts = minPts_val)
  
  # Extract the cluster partition.
  partition <- res$cluster
  partition <- as.factor(partition)
  levels(partition) <- as.character(as.numeric(levels(partition)) + 1)
  # Determine the number of clusters 
  uniqueClusters <- sort(unique(partition))
  # Exclude noise (0) from the count if present.
  if (0 %in% uniqueClusters) {
    nc <- length(uniqueClusters[uniqueClusters != 0])
  } else {
    nc <- length(uniqueClusters)
  }
  
  # Create a clusterlist
  clusterlist <- lapply(uniqueClusters, function(i) {
    partition == i
  })
  
  # Return the results in a list
  return(list(partition = partition, nc = nc, clusterlist = clusterlist))
}

################################################################################
# Mclust (GMM): 
################################################################################

# If manual model selection: add argument modelName = best_model to Mclust call
gmmWrapper <- function(x, k) {
  mod <<- Mclust(x, G = chosen_G)
  clusters <- mod$classification
  
  return(list(result = mod,
              clusterlist = lapply(1:k, function(i) clusters == i),
              partition = clusters,
              nc = k))
}
################################################################################
# ClusterR (GMM)
################################################################################

gmmWrapper_CR <- function(x, k, full_covariance = FALSE) {
  # full_covariance: if TRUE uses full covariance matrices in the GMM fit
  # Otherwise diagonal covariance
  
  # Fit GMM using ClusterR.
  mod <- ClusterR::GMM(data = x,
                       gaussian_comps = k,
                       dist_mode = "eucl_dist",       
                       seed_mode = "random_subset",
                       km_iter = 100,
                       em_iter = 100,
                       verbose = FALSE,
                       var_floor = 1e-10,
                       seed = global_seed,
                       full_covariance_matrices = full_covariance)
  
  # Number of observations
  n <- nrow(x)
  
  # Prepare a matrix to store the weighted density for each observation and each cluster.
  likelihoods <- matrix(NA, nrow = n, ncol = k)
  
  # For each cluster: compute the probability density weighted by the cluster weight.
  for (j in 1:k) {
    mu_j <- mod$centroids[j, ]
    
    # Extract the covariance matrix 
    if (length(dim(mod$covariance_matrices)) == 2) {
      # Diagonal case: Each row contains the variances.
      sigma_j <- diag(mod$covariance_matrices[j, ])
    } else if (length(dim(mod$covariance_matrices)) == 3) {
      # Full covariance case: A 3D array (the third dimension indexes clusters)
      sigma_j <- mod$covariance_matrices[ , , j]
    } else {
      stop("Unexpected dimensions for covariance_matrices.")
    }
    
    weight_j <- mod$weights[j]
    
    # Compute the density for each observation under the j-th Gaussian component.
    likelihoods[, j] <- weight_j * mvtnorm::dmvnorm(x, mean = mu_j, sigma = sigma_j)
  }
  
  # Assign each observation to the cluster with the highest weighted density.
  clusters <- apply(likelihoods, 1, which.max)
  
  # Return a list 
  return(list(result = mod,
              clusterlist = lapply(1:k, function(i) clusters == i),
              partition = clusters,
              nc = k))
}

################################################################################
# kmeans++
################################################################################

kmeansppWrapper <- function(x, k) {

  # Ensure x (select.npx) is a matrix 
  data_mat <- as.matrix(x)
  
  # Run k-means++
  res <- KMeans_rcpp(data = data_mat,
                     clusters = k,
                     initializer = "kmeans++",
                     num_init = 1,
                     max_iters = 100,
                     verbose = FALSE)
  
  # Extract cluster memberships vector
  clusters <- res$clusters
  
  # Return a list 
  return(list(result = res,
              clusterlist = lapply(1:k, function(i) clusters == i),
              partition = clusters,
              nc = k))
}

################################################################################
# Spectral clustering
################################################################################

speccWrapper <- function(x, k) {
  mod <<- kernlab::specc(x, 
                         centers = k, 
                         kernel = chosen_kernel, 
                         kpar = kpar_list, 
                         nystrom.red = FALSE, 
                         iterations = 200, 
                         mod.sample = 0.75)
  
  clusters <- as.integer(mod)
  
  # Create a list of logical vectors for each cluster membership
  clusterlist <- lapply(1:k, function(i) clusters == i)
  
  # Return the results in a list
  return(list(result = mod,
              clusterlist = clusterlist,
              partition = clusters,
              nc = k))
}



