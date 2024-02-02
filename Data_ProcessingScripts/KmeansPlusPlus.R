### KmeansPlus PLuts 
###SST August -08-04
# This function performs K-means clustering on the input data using the optimal number of clusters based on the distortion_fK criterion.
# Input:
# data: a matrix or data frame of input data
# nclusts: (optional) an integer representing the desired number of clusters. If missing, the function selects the most likely number of clusters from 2 to 20.
# The function first determines the optimal number of clusters using the Optimal_Clusters_KMeans function with the distortion_fK criterion.
# If the nclusts argument is missing, the function selects the most likely number of clusters from 2 to 20 based on the maximum value of the fK criterion.
# Next, the function runs the KMeans_rcpp function with fuzzy set to True and num_init set to 6 to perform fuzzy K-means clustering.
# The function then computes the maximum likelihood for each data point and selects the data points with likelihood below the 25th percentile as noise.
# Finally, the function returns a list containing the trial object from the KMeans_rcpp function, the subset of data points corresponding to the noise, 
# the likelihood below the 25th percentile, the number of clusters used for clustering, and (if nclusts is missing) the indices of the noise data points.
### develop a function to iteratively run through Kmeans
max_likelihood <- function(x) { 
  n_col <- which.max(x)
  p_max <- x[n_col]
  return(p_max)
}
cluster_func <- function(data, prior_clust) { 
  opt_kn = ClusterR::Optimal_Clusters_KMeans(data, max_clusters = 10, plot_clusters = T,
                                             criterion = 'distortion_fK', fK_threshold = 0.85,
                                             initializer = 'optimal_init', tol_optimal_init = 0.2)
  if(!missing(prior_clust)){
    max_clust <- prior_clust /2
    nclusts = which.max(opt_kn[2:max_clust]) + 1 # select the most likely number sless than before 
  } else {
    nclusts = which.max(opt_kn[2:20]) + 1 ## select the most likely number to cluster by 
  }
  trial <- ClusterR::KMeans_rcpp(data, clusters = nclusts, fuzzy = T, num_init = 6,
                                 max_iter = 100, verbose = F, seed = 43)
  likelihood_mat <- as.matrix(trial$fuzzy_clusters)
  max_ <- apply(likelihood_mat, 1, FUN = max_likelihood)
  Q1 <- quantile(max_, 0.25)
  below <- max_ < Q1
  c_below <- which(max_ < Q1)
  sub_data <- data[c_below,]
  likelihood_below <- max_[c_below]
  print('made it here')
  return(list(trial, sub_data, likelihood_below, nclusts))
}

iterative_kmeans <- function(data, nthresh, max_iter) {
  if (missing(data) || missing(nthresh) || missing(max_iter)) {
    #data is a non-scaled dataset
    #nthresh is the number you are ok with being in a noise vector
    #maximum iterations
    stop("missing input")
  }
  #center and scale data 
  dat <- ClusterR::center_scale(data, mean_center = T, sd_scale = T)
  # First iteration of clustering
  it_1 <- cluster_func(data = dat)
  # Second iteration of clustering
  #print('second iteration starts')
  it_2 <- cluster_func(data = it_1[[2]], prior_clust = it_1[[4]])
  
  # Values fixed by the adjustment
  fixed <-
    it_1[[3]] < apply(as.matrix(it_2[[1]]$fuzzy_clusters), 1, max_likelihood)
  print(dim(fixed))
  sub_fixed <- dat[fixed, ]
  print('here')
  # Values not fixed by the adjustment
  unfixed <-
    it_1[[3]] > apply(as.matrix(it_2[[1]]$fuzzy_clusters), 1, max_likelihood)
  #get incides of clusters with a low likelihood membership
  c_below <- which(apply(as.matrix(it_1[[1]]$fuzzy_clusters), 1,
                         max_likelihood) < quantile(apply(
                           as.matrix(it_1[[1]]$fuzzy_clusters), 1, max_likelihood
                         ), 0.25))
  #assign clusters for high likelihood clusters
  c1_clust <- it_1[[1]]$clusters[-c_below]
  c1_dat <-
    dat[-c_below, ] %>% data.table() %>% .[, clusts := c1_clust] %>% .[, noise := F]
  ## assing clusters for low likelihood clusters 
  c2_clusts <- it_2[[1]]$clusters + it_1[[4]]
  c2_dat <-
    dat[c_below,] %>% data.table() %>% .[, clusts := c2_clusts] %>% .[, noise := T]
  #create output dataset
  out_dat <- rbind(c1_dat, c2_dat)
  iterations <- 1
  # Check if number of unfixed values is greater than nthresh
  while (sum(unfixed) > nthresh && iterations <= max_iter) {
    if (iterations == 1) {
      new_dat <- it_2[[2]]
    } else if (iterations > 1) {
      new_dat <- it_3[[2]]
    }
    # If so, print rerunning model
    print('clustering rerunning')
    it_3 <- cluster_func(data = new_dat, it_2[[4]])
    # Values fixed by the adjustment
    fixed <-
      it_2[[3]] < apply(as.matrix(it_3[[1]]$fuzzy_clusters), 1, max_likelihood)
    sub_fixed <- new_dat[fixed, ]
    # Values not fixed by the adjustment
    unfixed <-
      it_2[[3]] > apply(as.matrix(it_3[[1]]$fuzzy_clusters), 1, max_likelihood)
    c_below2 <- which(apply(as.matrix(it_3[[1]]$fuzzy_clusters), 1,
                            max_likelihood) < quantile(apply(
                              as.matrix(it_3[[1]]$fuzzy_clusters), 1, max_likelihood
                            ), 0.25))
    c3_clusts <- it_3[[1]]$clusters + it_1[[4]] + it_2[[4]]
    c3_dat <-
      new_dat[c_below2,] %>% data.table() %>% .[, clusts := c3_clusts] %>% .[, noise := T]
    out_dat <- rbind(c1_dat, c2_dat, c3_dat)
    iterations <- iterations + 1
    if (iterations >= max_iter) {
      print('max iterations exceeded')
    }}
  return(list(out_dat, it_1[[1]], it_2[[1]]))
}