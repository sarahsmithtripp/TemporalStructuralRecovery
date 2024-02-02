# Load required libraries
library(dbscan)
library(cluster)
library(data.table)
library(sgsR)
library(clustra)
library(tidyverse)
library(SamplingBigData)
library(terra)

### here is a place to store this 
# file:///C:/Users/ssmithtr.stu/Downloads/v91i01.pdf
# Generate sample data
set.seed(123)
#read_in_data 
data_ <- arrow::read_feather('D:/Paper2_Clean/Satellite_Data_/Data_Tables/all_indexes.feather')
keep_cols <- grep(pattern = 'x|y|id|NDMI_dist_mag|NDVI_dist_mag|TCA_dist_mag|TCB_dist_mag|TCB_dist_mag_TCG|TCW_dist_mag', x = names(data_), 
invert = T, value = T)
data_fil <- data_ %>% dplyr::select(-c('index','x', 'y', 'NDMI_dist_mag', 'NDVI_dist_mag', 
                                       'TCA_dist_mag', 'TCB_dist_mag', 'TCG_dist_mag', 'TCW_dist_mag')) %>% 
                          filter(!is.na(NDVI_slope)) %>% na.omit()
index_vals <- data_ %>% dplyr::select(-c('NDMI_dist_mag', 'NDVI_dist_mag', 
                                                     'TCA_dist_mag', 'TCB_dist_mag', 'TCG_dist_mag', 'TCW_dist_mag')) %>% 
  filter(!is.na(NDVI_slope)) %>% na.omit()

data_scale <-scale(data_fil) ## scale data for clustering

# iterative kmeans function -----------------------------------------------

source('D:/Paper2_Clean/Data_ProcessingScripts/KmeansPlusPlus.R')

#ctest <- iterative_kmeans(data_fil, nthresh = 50, max_iter = 10)
# save_all_clusts <- c('F:/Quesnel_RPA_Lidar2022/Clusters/all_fire_clusters/kmeans/it_kmeans_highsev.rds')
# saveRDS(ctest, save_all_clusts)
test_sample <- data_fil[sample(nrow(data_fil), 40000), ]
##  test run on  the overall 
#N thresh is set to 1 percent  of the data 
test_c <- iterative_kmeans(data_fil, nthresh = 70000, max_iter = 10)

# Run kmeans ++ from iterative kmeans output -----------------------------------------
## get centroids of iterative kmeans to set the seed for the final clustering 
centroids_1  <- test_c[[2]]$centroids
centroids_2 <- test_c[[3]]$centroids
## CENTROIDS rows of centroids = nclusters, columnds = n_cols 
centroids  <- rbind(centroids_1, centroids_2)
dat <- ClusterR::center_scale(data_fil, mean_center = T, sd_scale = T)
dat_dist <- dist(test) ## calculate a distance matrix
new_kmn <- ClusterR::KMeans_rcpp(as.data.frame(dat), clusters = 13, num_init = 10, max_iters = 100, 
fuzzy = T, CENTROIDS = centroids)
#saveRDS(new_kmn, file ='F:/Quesnel_RPA_Lidar2022/Clusters/all_fire_clusters/kmeans/it_reclust_highsev_Sept.rds' )
#arrow::write_feather(as.data.frame(dat),'F:/Quesnel_RPA_Lidar2022/SatHigh_Sev_NaturalRegn/May22ClusteringData.feather')

# saveRDS(new_kmn, file ='F:/Quesnel_RPA_Lidar2022/Clusters/Exploratory_Clusters/test_reclust_highsev.rds')
# arrow::write_feather(as.data.frame(dat),'F:/Quesnel_RPA_Lidar2022/Clusters/Exploratory_ClustersClusteringData.feather')
# source for combined clusters https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/


# Collapse Like Clusters --------------------------------------------------

#cluster_IDS = arrow::read_feather('F:/Quesnel_RPA_Lidar2022/Clusters/all_fire_clusters/kmeans/high_severite_reclusters/cluster_IDs.feather')
cluster_IDs = readRDS('D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/it_kmeans_highsev.rds')
dat_clust <-as.data.frame(data_fil)
dat_clust$clusters <-  cluster_IDs$clusters
## get a sample to run a test to collapse based on similarities  
library(clues)
in_data <- dat_clust %>% group_by(clusters) %>% 
mutate(count_clust = n()/nrow(.)) %>% sample_frac(., .005, weight = count_clust)
clusts = in_data$clusters
dist_data <- in_data %>% ungroup() %>% dplyr::select(-c('clusters', 'count_clust', 'index'))
dat_dist <- ClusterR::distance_matrix(dist_data, method = 'euclidean')
dunn_unique<- clValid::dunn(Data =dist_data, mem = clusts )
#https://www.rdocumentation.org/packages/fpc/versions/2.2-10/topics/cluster.stats
alt <- fpc::distcritmulti(dist_data, clusts, seed = 123, count = T ) #note sure what this  gives 
alt2 <-  clues::get_Silhouette(y = as.matrix(dist_data), mem = clusts)
#combine based on neighbors 
neighbor_clust <- alt2$neighbor
clust_neighbors <- as.data.frame(cbind(clusts, neighbor_clust)) %>% 
  mutate(average_s = alt2$s) %>% 
  group_by(clusts, neighbor_clust) %>% 
  mutate(nclust= n()) %>% 
  group_by(clusts) %>% 
  mutate(max_n = max(nclust), 
         average_s = mean(average_s, na.rm = T), 
         most_sim = case_when(max_n == nclust ~ neighbor_clust))%>% 
  filter(!is.na(most_sim)) %>% distinct()

##  define sil threshold 
clust_mergers  <- filter(clust_neighbors, average_s < 0.07)
## cluster mergers 
#change clusters 
mapping <- c(1,1, 2, 3, 4, 5, 6, 3, 5, 6, 7, 1, 8)

clusts_new <- as.data.frame(clusts) %>% 
  mutate(new_clusts = mapping[clusts])
alt_clusts <- as.vector(clusts_new$new_clusts)
alt3 <- get_Silhouette(y = as.matrix(dist_data), mem = alt_clusts)
alt3$avg.s

mapping2 <- c(1,2,3,2,4,5,6,4,7,8,1,9)
clusts_new2 <- as.data.frame(clusts) %>% 
  mutate(new_clusts = mapping2[clusts])
alt_clusts2 <- as.vector(clusts_new$new_clusts)
alt4 <- get_Silhouette(y = as.matrix(dist_data), mem = alt_clusts2)
alt4$avg.s ## same as before so going with cluster 3 

new_neighbors <- alt3$neighbor
clust_neighbors2 <- as.data.frame(cbind(alt_clusts, new_neighbors)) %>% 
  mutate(average_s = alt3$s) %>% 
  group_by(alt_clusts, new_neighbors) %>% 
  mutate(nclust= n()) %>% 
  group_by(alt_clusts) %>% 
  mutate(max_n = max(nclust), 
         average_s = mean(average_s, na.rm = T), 
         most_sim = case_when(max_n == nclust ~ new_neighbors))%>% 
  filter(!is.na(most_sim)) %>% distinct()

##add_back to og data 
out_data <- dat_clust %>% 
  mutate(comb_clusters = mapping[clusters], 
         index = index_vals$index, 
         x = index_vals$x, 
         y = index_vals$y)


## test for cluster similarities 
##Get most similiar clusters 
clust_sil <- test$separation.matrix 
clust_sil[clust_sil == 0] <- NA_real_
avg_sil <- test$avg.silwidth
std_separtion <- sd(test$separation)
avg_sil - std_separtion

##Overall, the code aims to identify duplicated clusters based on their range values calculated
# from the minimum and maximum values of centroid distance.
## Clusters where that share another clsuters as their neighbor .
similiar_clusts <-  clust_sil %>% data.frame() %>% 
rowwise() %>% mutate(min_value = min(c_across(everything()), na.rm = T),
 min_column = which.min(c_across(everything()))) %>% rownames_to_column(., var = 'cluster')
cluster_pairs <- dplyr::select(similiar_clusts, c(min_column, cluster)) %>% 
mutate(cluster = as.numeric(cluster))%>% mutate(min_val = min(c_across(everything())), 
max_val =max(c_across(everything())), 
range = paste(min_val, max_val)) #CALCULATE THE MAX 
duplicated_clusters <- cluster_pairs %>% group_by(range) %>% filter(n() > 1)
## 9 = 5


# Create Raster of Cluster Locations  -------------------------------------
data_out <-  merge %>% mutate(id = index) 
## get the index raster to build cluster raster from (this is generated when creating the data folder that is read in )
id_rast <- rast('D:/Paper2_Clean/idx_rast.tif')
id_vals = unique(id_rast[])

matched_IDS = match(id_rast, out_data$index)
filter_data_out <- out_data[matched_IDS[!is.na(matched_IDS)],]
matched_IDS[!is.na(matched_IDS)] <- filter_data_out$comb_clusters
## reproj 
og_clusts <- rast('F:/Quesnel_RPA_Lidar2022/Clusters/SBS_clusters/sbs_clusters.tif')
new_rast <- rast(select(out_data, c(x, y, z = comb_clusters)), type = 'xyz')
crs(new_rast) <- crs(og_clusts)
new_rast <- project(new_rast, og_clusts)
plot(matched_IDS)
# terra::writeRaster(new_rast, filename = 'D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/clusters.tif', 
#                    overwrite = T)
# arrow::write_feather(out_data, sink = 'D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/all_metrics.feather')

