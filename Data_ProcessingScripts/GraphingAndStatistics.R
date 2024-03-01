library(terra)
library(tidyverse)
library(sgsR)
library(data.table)
library(cowplot)
library(gridExtra)
library(ggrepel)#deals with labeling issues 
library(factoextra) #for pca 
library(ggh4x) # extra graphics code 
library(terra)
library(dunn.test)
library(FSA)
library(ggtext)
library(kableExtra)
library(flextable)
library(factoextra)
library(shadowtext)
library(officer)
library(stringr)
set.seed(123)

#location to save outputs
list_save <- list()

# Test for differences in clusters ----------------------------------------
table_function <- function(x,y, data){
  data <- as.data.frame(data)
  sel_dat <- grep(pattern = paste0(x,'|',y), names(data))
  data_use <- data[, sel_dat]
  
  if(nchar(as.character(data_use[1,2])) > 1){
    names(data_use) <- c('y', 'x')
    data_use$y <- as.factor(data_use$y)
    data_use$x <- as.numeric(data_use$x)
  }
  else if(nchar(as.character(data_use[1,2])) == 1){
    names(data_use) <- c('x', 'y')
    data_use$y <- as.factor(data_use$y)
    data_use$x <- as.numeric(data_use$x)
  }
  sum_dat <- data_use %>% group_by(y) %>% dplyr::summarise(mean_dat = mean(x, na.rm = T), 
                                                           sd_dat = sd(x, na.rm = T)) 
  names(sum_dat) <- c('cluster', paste0(x,'_mean'), paste0(x, 'sd'))
  kruskal_obj <- kruskal.test(x ~ y, data= data_use)
  dunn_test_obj <- dunnTest(x ~ y, data= data_use, method = 'sidak', two.sided = T) 
  chi_info<- paste0(expression('$\\chi^2$'),'(' ,kruskal_obj$parameter[1],', N = ', nrow(data_use),
                    ') = ', round(kruskal_obj$statistic, 2),', p = ', round(kruskal_obj$p.value,2))
  sel_pairs <- dunn_test_obj$res %>% dplyr::filter(P.adj >= 0.1) 
  test <- length(unique(sel_pairs$Comparison))
  if(test == 0){
    out_put <- data.frame(x_val = x, 'chi_test' = chi_info, `like_pairs` = 'All Unique', `adj_p` = ' ')
    return(list(out_put, sum_dat))
  }
  else if (test > 0){
    list_vals <- as.character(unique(sel_pairs$Comparison))
    in_out <- c()
    p_out <- c()
    for(i in 1:2){
      if(i == 1){
        in_out <- c(list_vals[i])
        p_out <- round(unique(sel_pairs$P.adj)[i],2)
      }
      else if (i > 1){
        temp <- list_vals[i]
        temp2 <- round(unique(sel_pairs$P.adj)[i],2)
        in_out <- paste(temp,',', in_out ) %>% str_remove(pattern = 'NA ,')
        p_out <- paste(temp2, ',', p_out) %>% str_remove(pattern = 'NA ,')
      }
    }
    in_out
    out_put <-  data.frame(x_val = x, 'chi_test' = chi_info, `like_pairs` = in_out, `adj_p` = p_out)
    return(list(out_put, sum_dat, in_out))
  }
}
 ###

### Set colors for clusters
in_cols = c("#a24936","#ffcf00","#ee6123","#D68C30","#bdb63d","#777949","#006ba6","#7b9d79")
# colors that symbolize years 
year_cols <- c("#22577a","#38a3a5","#57cc99","#80ed99")

## read in raster data 
in_r <- rast('D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/clusters.tif')
values_r <- as.data.frame(in_r, cells = T, na.rm = T)

in_r$strata <- in_r$z

## read in satellite metrics 
data_ <- arrow::read_feather('D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/all_metrics.feather')
#drop values not used in clusiterin 
keep_cols <- grep(pattern = 
                    'x|y|NDMI_dist_mag|NDVI_dist_mag|TCA_dist_mag|TCB_dist_mag|TCB_dist_mag_TCG|TCW_dist_mag', x = names(data_), 
                  invert = T, value = T)
find_index <- data_ %>% base::subset(., select = c(keep_cols, 'index')) %>% 
  filter(!is.na(NDVI_slope)) %>% na.omit()
#filter for empty values 
data_fil <- data_ %>% dplyr::select(-c('x','y')) %>% 
  filter(!is.na(NDVI_slope)) %>% na.omit()
index_PCA <- data_fil$index

# Computer Summary metrics  -----------------------------------------------
data_sum <- data_fil
data_sum$clust <- values_r$z

get_avg_sf <- function(x,metric, name_in) {
  #name_in type of grouping you are looking for 
  x <- as.data.frame(x)
  sel_col <- grep(paste0("^", name_in,"$"), names(x))[1]
  sel_col
  x$clust <- x[,sel_col]
  x_out <- x %>% group_by(clust) %>% pivot_longer(ends_with(metric),
                                                  values_to = 'measure', 
                                                  names_to = 'metric') %>% 
    group_by(metric,clust) %>% 
    summarise(mean_ = mean(measure, na.rm = T), 
              sd_ = sd(measure, na.rm = T)) %>%
    mutate(max_ = max(mean_),
               min_ = min(mean_),
              clust =ifelse(mean_ == max_ | mean_ == min_, clust,NA),
               sd_ = ifelse(mean_ == max_ | mean_ == min_, sd_,"")) %>%
    filter(!is.na(clust))
  return(x_out)
}
slopes <- get_avg_sf(data_sum, 'slope', name_in = 'clust')
regrowth <- get_avg_sf(data_sum, 'regrowth', name_in = 'clust')
measure_5yrs  <- get_avg_sf(data_sum, 'vals_years', name_in = 'clust')

slopes$name <- 'slope'
regrowth$name <- 'regrowth'
measure_5yrs$name <- '5 year measure'


list_save[[1]] <- bind_rows(slopes, regrowth, measure_5yrs)


# PCA of cluster variables ------------------------------------------------

pca_data <- prcomp(data_fil %>%
                     dplyr::select(-c('index', 'comb_clusters', 'clusters')), #remove index and cluster names
                   center = T, scale = T, na.action = na.omit)
pca_ind <- facto_summarize(pca_data, element = 'ind'
                           , result = 'coord', axes = 1:2) #get first two axes 
pca_var_full <- facto_summarize(pca_data, element = 'var', result = 'coord', axes = 1:2) %>% as.data.frame() %>%  
  mutate(scaled_x = summary(pca_data)$importance[2,1] * Dim.1, 
         scaled_y = summary(pca_data)$importance[2,2] * Dim.2, 
         scaled_mag = (scaled_x^2 + scaled_y^2)^0.5, 
         label = paste(name, '=',
                       round(scaled_mag, 2))) %>% arrange(., scaled_mag)
list_save[[2]] <- pca_var_full
#fix labeling for plot
pca_var_full <- pca_var_full %>% mutate(label_short = str_replace_all(name, pattern = "_",
                                                                      replacement = " "),
                                        label_short = str_replace_all(label_short, 
                                                                      pattern = "vals years",
                                                                      replacement = "measure (+5 yrs)"))
#filter to the first five loadings for biplot 
pca_var <- pca_var_full %>% .[(nrow(pca_var_full)-4):nrow(pca_var_full), ]

#add clust values to  the PCA data for coloring 
pca_ind$clust <- values_r$z
##added to look  at OG clustering output (before the clusters are combined)
#pca_ind$alt_clust <- cluster_IDs$clusters
#pca_ind$str_clust <- cluster_IDs$fuzzy_clusters[,1]
pca_sample <- pca_ind %>% group_by(clust) %>% filter(!is.na(clust)) %>% sample_n(1000) 

pca_plot <- ggplot(data = pca_sample, group = as.factor(clust)) +
  geom_point(data = pca_sample,
             aes(x = Dim.1, y = Dim.2, color = as.factor(clust)),
             size = 0.75, alpha = 0.6) +
  stat_ellipse(data = pca_sample, geom = 'polygon',
               aes(x = Dim.1, y = Dim.2, fill = as.factor(clust)),color = 'black', 
               alpha = 0.6, level = 0.4, size = 0.75) +
  coord_equal() + 
  coord_axes_inside(labels_inside = TRUE)  +
  geom_segment(data = pca_var, aes(x = 0, y= 0, xend = 12 * scaled_x, yend = 20 * scaled_y), size = 1.2, 
               arrow = arrow(length = unit(0.5, "cm"))) + 
  geom_label_repel(data = pca_var, aes(x = scaled_x* 15, y= scaled_y * 25, label = label_short)) +
  guides(color = FALSE, 
         fill =guide_legend(ncol = 1, 
                            label.position = 'top',
                            title.position = 'top')) +
  theme(axis.line = element_line(linewidth  = 0.2), axis.text = element_text(size = 15), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_rect(fill="white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3)) + labs(color = 'Spectral Cluster') +
  #scale_color_brewer(palette = 'Paired') + xlim(-6, 6) + ylim(-6, 6) 
  scale_color_manual(values = in_cols) + 
  scale_fill_manual(values = in_cols, 
                    guide = guide_legend(direction = 'horizontal')) + xlim(-8, 8) + ylim(-8, 8) + 
  labs(fill = 'Cluster') + xlab('PC1') + ylab('PC2') 


index_vals <- ggplot(pca_var_full, aes(x = reorder(name, -scaled_mag), y = scaled_mag)) + 
  geom_col(fill = 'magenta4') + shadowtext::geom_shadowtext(aes(y = 0.15, label = label_short), fontface = 'bold' , size = 4.5, angle = 90) +
  theme_bw() + ylab('Scaled PC Loading') +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank())
compound_pc <- cowplot::plot_grid(pca_plot, index_vals, nrow = 2, rel_heights = c(1, 0.4), 
                                  labels = c('A', 'B'))
compound_pc
# save_plot(compound_pc, filename = 'F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/PCA_plot.png',
#           base_width = 8, base_height = 9)

##get the most distinct spectral clsuter 
# clsut <- readRDS('D:/Paper2_Clean/Spectral_Clusters/Kmeans_outputs/it_kmeans_highsev.rds')
# centroids <- clsut[[2]]$centroids
# #change clusters 
# mapping <- c(1,1, 2, 3, 4, 5, 6, 3, 5, 6, 7, 1, 8)
# 
# 
# clusts_new <- as.data.frame(clusts) %>% 
#   mutate(new_clusts = mapping[clusts])
# 
# dist(centroids)
# average_dist <- apply(dist(centroids), 2, FUN = mean)
# names(average_dist) <- mapping[names(average_dist)]
# ordered <- order(average_dist)
# names_out <- order(average_dist)
# real_clust <- mapping[names_out]
# average_dist




# Baseline Graph ----------------------------------------------------------
##Read in full trajectories (compiled in Python)
in_dat <- list.files('D:/Paper2_Clean/Spectral_Clusters/metrics/Trajectories',
                     full.names = T) %>% map(., .f = arrow::read_feather) %>%
  rbindlist()
plot_data <- in_dat %>% filter(Years_indist >=-3 & Years_indist <18) %>% 
  dplyr::select(Diff_Baseline, band, clust, Years_indist) %>% distinct() #%>% 
  mutate(Diff_Baseline = ifelse(clust == 2 & Years_indist == 4, NA, Diff_Baseline))  #clean up for plotting
plot_data
sel_plot <-plot_data %>% subset(clust %in% c(1,3,4,5,8)) %>% 
  subset(band %in%c("TCA",
                    "NBR",
                    "TCW"))
data_sel <- data_sum  %>% 
  subset(clust %in% c(1,3,4,5,8))
library(gg.layers)
# ggplot(data_sel) + 
#   geom_boxplot2(aes(comb_clusters,
#                    TCA_regrowth,
#                    group = comb_clusters),
#                outlier.shape = NA)
### plot data 
base_line_graph <- ggplot(plot_data, aes(Years_indist, group = clust)) +
  #geom_line(aes(y = measure_mean, color = clust_new)) + 
  #geom_ribbon(aes(ymin = lwr_baseline, ymax = upp_baseline, fill = clust), alpha = 0.5) +
  geom_line(aes(y = Diff_Baseline, color = as.factor(clust)), size = 0.8)  + 
  facet_wrap(~band, scales = 'free', ncol = 2)  +
  xlim(-1,10) + 
  geom_hline(yintercept = 0) + xlab(NULL) +
  geom_vline(xintercept = 0, linewidth  = 0.3, linetype = 'longdash' ) + 
  ylab('% Difference Pre-Fire Baseline') + xlab('Years in Disturbance') +
  scale_color_manual(values = in_cols) +
  labs(color = 'Spectral\nCluster') + theme_classic(base_size = 15) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = c(0.7, 0.12), 
        legend.direction = 'horizontal',
        legend.key.size = unit(0.8, "cm"),
        legend.key.width =unit(0.8, "cm"),
        legend.key.height = unit(0.8, "cm"),
        panel.border = 
          element_rect(colour = "black", fill=NA, size=0.4)) +
  xlim(-3, 11)
base_line_graph

# save_plot(base_line_graph, filename = 'F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/Traj_plot.png',
#           base_height = 9, base_width = 10)

# Proportions of Clusters  ------------------------------------------------
cluster_counts <- data.table(values_r) %>% .[, .(Count = .N), by = z] %>% 
  .[, proportion := Count/ sum(Count)]

#Read in cluster counts attributed to their disturbance year 
cluster_counts_attributed <- arrow::read_feather('D:/Paper2_Clean/Spectral_Clusters/metrics/counts_by_year.feather') %>% 
  filter(!is.na(category)) %>% 
  group_by(Dist_Year) %>% 
  mutate(sum_year = sum(years)) %>% 
  ungroup() %>% 
  mutate(prop_year = years/sum_year, 
         clust = category) 
all_years <- ggplot(cluster_counts, 
                    aes(x = 1, y = proportion * 100, fill = as.factor(z))) + 
         geom_bar(stat = 'identity') + scale_fill_manual(values = in_cols) + 
  xlab('All Years') + theme_bw() + ylab('% of Landscape = Cluster') +coord_flip(ylim = c(0,100)) +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x =  element_text(size = 12, face = 'bold'),
        axis.ticks.y = element_line(size = 2), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 17, face = 'bold'),
        panel.border = element_rect(color = "black", 
                                     fill = NA, 
                                     size = 1.5)) #+ 
  scale_y_continuous(breaks = c(0, 25, 50, 100))
all_years
source('D:/Paper2_Clean/Data_ProcessingScripts/Proportions_Graph.R')
#out ## this ithe proportions graph and you can plot it here 
save_plot(out, filename='F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/proportions.png',
          base_width = 6, base_height = 7)

list_save[[5]] <- cluster_counts
list_save[[6]] <- cluster_counts_attributed
# Make Structure Graph  ---------------------------------------------------
#1. read in satellite data 
#2. Read in lidar cluster data 
#3. Compare the two 

sat_clusters <- in_r
sat_df <- as.data.frame(sat_clusters, xy = T, cell = T, na.rm =F)
sat_df <- as.data.table(sat_df)
sat_count <- na.omit(sat_df) %>% 
  group_by(z) %>% summarise(count = n())

##2 Read in lidar data 
#lidar_all <- rast('D:/Paper2_Clean/RPA_data/ABA_Models/Model_Outputs/Lidar_rast2.tif')
lidar_dfs <- lidar_all %>% as.data.frame(xy = T, cell = T, na.rm = F) #%>% 
lidar_dfs <- read_csv2('D:/Paper2_Clean/RPA_data/ABA_Models/Model_Outputs/lidar_models_df.csv') %>% 
  data.table()
lidar_to_base <- sat_df[lidar_dfs, on = c('x','y')] # merge the lidar data and satellite data 


##LidarGraph 
lidar_sat_filtered <- lidar_to_base[complete.cases(lidar_to_base),] 
dim(lidar_sat_filtered)
## stack based on model_type 
stack_clust <- lidar_sat_filtered %>% 
  pivot_longer(., col = c('composition', 'stem_counts','basal_model', 'PBsl', 'bare_ground'),
                            names_to = 'model_type', values_to = 'model_measure') %>% 
  mutate(clust = as.factor(round(z))) 



##clean_lidar 
lidar_sat <- stack_clust %>% 
  filter(!is.na(clust)) # remove points that were not clustered 


# Plotting Lidar  means and deviations  -----------------------------------

lidar_sum <-  lidar_sat %>%
  filter(model_measure < 8000) %>% # filter obvious outliers  
  filter(model_type != 'PBsl') %>%
  mutate(model_measure = case_when(model_type == 'stem_counts' ~ log(model_measure),
                                 TRUE  ~ as.numeric(model_measure)),
         years_indist = as.factor((2022 - dist_year)))

lidar_filter <- lidar_sum  

unique(lidar_filter$dist_year)
lidar_filter <- lidar_sum
lidar_filter$model_type <- factor(lidar_filter$model_type,
                                 labels =  c("Bare~ground~('%')", 
                                    "Basal~Area~(m^2/hectare)",
                                    'Coniferous:Decidious',
                                    "Stem~Counts~(stems/900~m^2)"))
rect_df <- data.frame(dist_year = c(2006, 2010, 2014, 2017),
stem_counts = range(lidar_sat_filtered$stem_counts, na.rm = T),
composition = range(lidar_sat_filtered$composition, na.rm = T),
basal_model= range(lidar_sat_filtered$basal_model, na.rm = T),
bare_ground = range(lidar_sat_filtered$bare_ground, na.rm = T))


structural_sums_year <- lidar_sum %>% group_by(years_indist,  model_type) %>%
  summarise(mean_ = mean(model_measure, na.rm = T), 
            sd_  = sd(model_measure, na.rm = T)) %>% 
  pivot_wider(names_from = years_indist, 
              values_from = c(mean_, sd_))

structural_sums_cluster_year_long <- lidar_sum %>%
  group_by(clust, years_indist, model_type) %>% 
  #antilog stem counts to untransform them 
  mutate(model_measure= case_when(model_type == 'stem_counts' ~ exp(model_measure),
                                TRUE  ~ as.numeric(model_measure)))%>%
  group_by(clust, years_indist, model_type) %>% 
  summarise(mean_ = mean(model_measure, na.rm =T), 
            sd_ = sd(model_measure),
            n_group = n(),
            sem = sd_/(n_group^0.5))
structural_sums_cluster_year_long$model_type <- 
  #add factor labels for plotting 
  factor(structural_sums_cluster_year_long$model_type,
         labels = c("Bare~ground~('%')", 
                    "Basal~Area~(m^2/hectare)",
                    'Coniferous:Decidious',
                    "Stem~Counts~(stems/900~m^2)"))

structural_sums_cluster_year <- structural_sums_cluster_year_long %>% 
  pivot_wider(names_from = years_indist, 
               values_from = c(mean_, sd_))%>% 
  arrange(model_type)



### additional structure graph 
all_dat <- lidar_filter %>%  left_join(mutate(structural_sums_cluster_year_long, years_indist = as.factor(years_indist)), 
                                       by = c('model_type', 'clust', 'years_indist')) %>% 
  mutate(model_measure= case_when(model_type == 'Stem~Counts~(stems/900~m^2)' ~ exp(model_measure),
                           TRUE  ~ as.numeric(model_measure))) %>% filter(!is.na('model_measure')) %>% 
  filter(model_measure < 2000) %>% 
  group_by(years_indist, model_type) %>% 
  mutate(mean_tot = mean(model_measure, na.rm = T), 
         dist_mean = (mean_ - mean_tot)) 
### calculate the total number for each pixel and the proportion for each  year
pie_dat <- dplyr::select(lidar_filter, c('cell', 'years_indist', 'clust')) %>% distinct() %>% 
  group_by(years_indist, clust) %>% 
  summarise(n_year = n()) %>% group_by(clust) %>% 
  mutate(tot_ = sum(n_year),
         prop = n_year/tot_)
## add the total lidar coverage to a file 
list_save[[11]] <- pie_dat

## create a stacked bar chart
#which sums the total number of cluster w/  lidar data for each year acludster
#for each cluster 
lidar_coverage <- ggplot(pie_dat, aes(x = clust, 
                                      y =  prop,
                                      group = clust,
                                      fill = as.factor(years_indist))) + 
  geom_col(position = 'stack') + labs(fill = 'Time Since Disturbance\n(years)')+
  geom_text(aes(y = 0.5, x = clust, label = paste('n =',tot_)),
            size = 6) +
  xlab('Spectral Cluster') + 
  ylab('RPAs Lidar Coverage \nper Spectral Cluster') +
  scale_fill_manual(values = year_cols) + 
  guides(fill= guide_legend(ncol = 4))  +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(fill = NA, color = 'black',size = 1),
        plot.margin = unit(c(1,0.2,0.2,0.2), "cm"),
        legend.position = "top") 
lidar_coverage
# save_plot(lidar_coverage,
#             filename='F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/LidarCoverage.png',
#           base_height = 6)


# PCofA Trajectories  (ADONIS PERMANOVA)-----------------------------------------------------
struct_pc <- dplyr::select(lidar_filter, 
                           c('years_indist', 'clust', 'model_type', 'model_measure', 'cell')) %>%  
  filter(!is.na(model_measure))%>% 
  pivot_wider(id_cols = c('clust', 'years_indist', 'cell'),
              names_from = 'model_type', values_from = 'model_measure',
              values_fn = function(x) {
                out <- median(x, na.rm = T)
              }) %>%  filter(!is.na(`Stem~Counts~(stems/900~m^2)`)) %>% mutate(`Bare~ground~('%')` = 1-`Bare~ground~('%')`)
## save to output dataset 
# Combine Like Trajectories -----------------------------------------------
source('D:/Paper2_Clean/Data_ProcessingScripts/MANOVA_And_PostHoc.R')
#in_cols2 = c("#f6a884","#ffcf00","#7b9d79","#00916e","#E57128")
in_cols2 = c("#a24936","#ffcf00","#777949","#006ba6","#E57128")
new_clustletters <- c('1', '2','C', 'C', 'C', '6','7', 'C')
new_struct <- struct_pc %>% mutate(new_clusts = new_clusters[clust])
# colors that symbolize years 
year_cols5 <- c("#22577a","#38a3a5","#9bc2b2", "#57cc99","#80ed99")


## create a stacked bar chart of lidar data used to combine 
#which sums the total number of cluster w/  lidar data for each year acludster
#for each cluster 
lidar_coverage <- ggplot(pie_dat, aes(x = clust, 
                    y =  prop,
                    group = clust,
                    fill = as.factor(years_indist))) + 
  geom_col(position = 'stack') + 
  labs(fill = 'Time Since Disturbance\n(years)')+
  geom_text(aes(y = 0.5, x = clust, label = paste0('n=',tot_)),
            size = 4) + coord_cartesian(clip = 'off')+
  xlab('Spectral Cluster') + 
  ylab('RPAs Lidar Coverage \nper Spectral Cluster') +
  scale_fill_manual(values = year_cols5) + 
  guides(fill= guide_legend(ncol = 5))  +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(fill = NA, color = 'black',size = 1),
        plot.margin = unit(c(1,0.2,0.2,0.2), "cm"),
                                        legend.position = "top",
        legend.title.align=0.5) 
lidar_coverage
# save_plot(lidar_coverage,
#             filename='F:/Sync/PhD_Writing/Paper2/test/images/March_Clusters/LidarCoverage.png',
#           base_height = 6)

# Clean Up the Lidar Structure Data and Plot Structure Dynamics -----------

clean_struct <- all_dat %>% mutate(new_clusts = new_clustletters[clust]) %>% 
  dplyr::select(c('years_indist', 'model_measure', 'new_clusts', 'model_type')) %>% 
  group_by(model_type, years_indist) %>% 
  mutate(mean_measure = mean(model_measure,na.rm =T ),
         years_indist = as.numeric(ifelse(years_indist==4, 5, paste0(years_indist)))) %>% 
  group_by(model_type, new_clusts, years_indist) %>% 
  mutate(n_group_year = n(),
    mean_clust = mean(model_measure, na.rm = T), 
         sd_clust = sd(model_measure, na.rm = T),
         sem = sd_clust/n_group_year) 
##
struct_mean  <- distinct(dplyr::select(clean_struct,
                       c('mean_measure',
                         'mean_clust', 'sem', 'years_indist','n_group_year',
                         'new_clusts'))) %>% 
  mutate(dist_mean = mean_clust - mean_measure,
         dist_perc = ((mean_clust - mean_measure) / mean_measure )*100) 

head(struct_mean)
struct_mean_n <- mutate(struct_mean, new_groups 
                        = case_when(new_clusts == "1" ~ 1,
                                    new_clusts == "C" ~ 2, 
                                    new_clusts == '2' ~ 3, 
                                    new_clusts == '6' ~ 4, 
                                    new_clusts == '7' ~ 5))
clust_names <- c("Low Density\nConifer",
                 "High Density\nConifer",
                 "Shelterwood to\nLate Growth",
                 "Stem Exclusion to\nStem Loss",
                 "Stem Exclusion to\nIngrowth")
struct_mean_n$new_groups <- 
  #add factor labels for plotting 
  factor(struct_mean_n$new_groups,
         labels = clust_names)

structure_alt2 <- ggplot(data = struct_mean_n) +
  facet_grid(model_type ~ new_groups, scales = 'free_y', 
             labeller = 
               labeller(.rows = label_parsed, .multi_line = T),
             switch = 'y') +
  labs(filler = "Spectra\nCluster") +
  geom_hline(aes(yintercept = 0), linewidth = 0.5)+
  geom_segment( aes(x=as.factor(years_indist),
                    xend=as.factor(years_indist), 
                    y=0, yend=dist_mean), color="grey") +
  geom_point(aes(
    x = as.factor(years_indist),
    y = dist_mean ,
    color = as.factor(new_clusts)
  ), size = 4, 
  position = position_dodge(width = 0.3))  +
  xlab("Time Since Disturbance\n(years)") +
  ylab(NULL) +
  labs(
    color = paste0(
      "Difference From <br> Mean <br>",
      "<span style='font-size: 10pt'>",
      "(+/-SEM)</span>"
    )
  ) +
  scale_color_manual(values = in_cols2) +
  theme_classic(base_size = 15) +
  scale_x_discrete(breaks  =c(5,8,11,12,16),
labels = c('5', '8','11','12','16'))+
  ##add labels that describe the zero line 
  facetted_pos_scales(y = list(
    scale_y_continuous(sec.axis =
                         dup_axis(
                           breaks = c(0.0),
                           labels = c('No Deviance\nFrom Mean'),
                           name = NULL
                         )),
    scale_y_continuous(sec.axis =
                         dup_axis(
                           breaks = c(0.0),
                           labels = c('No Deviance\nFrom Mean'),
                           name = NULL
                         )),
    scale_y_continuous(sec.axis =
                         dup_axis(
                           breaks = c(0),
                           labels = c('No Deviance\nFrom Mean'),
                           name = NULL
                         )),
    scale_y_continuous(sec.axis =
                         dup_axis(
                           breaks = c(0),
                           labels = c('No Deviance\nFrom Mean'),
                           name = NULL
                         ))
  )) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = 'black',
      size = 1
    ), strip.text.y = element_text(size = 10),
    strip.placement = 'outside',
    strip.background.y = element_blank(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(1.15, 0.6),
    plot.margin = unit(c(0.1, 4, 0.1, 0.1), 'cm'),
    axis.text.y.right = element_text(angle = 270, hjust = 0.5, 
                                     size = 10)
  ) + 
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.position = "none")
structure_alt2
# 
save_plot(structure_alt2, 
          filename='F:/Sync/PhD_Writing/Paper2/test/images/March_Clusters/structure_split.png',
          base_width = 12, base_height = 10)


result_data <- struct_mean %>% dplyr::select('model_type', 
                                      'mean_clust', 'sem', 
                                      'new_clusts', 'years_indist') %>% 
  mutate(mean_clusts = ifelse(model_type == 'Stem~Counts~(stems/900~m^2)', round(mean_clust), 
                               mean_clust),
         sem = ifelse(model_type == 'Stem~Counts~(stems/900~m^2)', round(sem), 
                               sem)) %>% 
  rename('mean'= mean_clusts,
         'group' = new_clusts)%>%
  pivot_wider(id_cols = c('group', 'model_type'),
              names_from = 'years_indist', values_from = c('mean', 'sem')) %>% 
  dplyr::select(c('group', 'model_type', 'mean_5', 'sem_5',
                  'mean_8', 'sem_8','mean_11', 'sem_11',
                  'mean_12', 'sem_12',
                  'mean_16', 'sem_16'))
  
new_row <- list( 'Average +/- Standard Deviation By Year')
lidar_sums_table <- result_data %>% group_by(model_type) %>% mutate(group = as.factor(group)) %>% 
  as_grouped_data(groups = 'model_type') %>% 
  flextable() %>% 
  set_header_labels(values = list(model_type = "Model", 
                                  clust = 'Cluster', 
                                  mean_5 = '5 years',
                                  sem_5 = '+/-', 
                                  mean_8 = '8 years',
                                  sem_8 = '+/-',
                                  mean_11 = '11 years', 
                                  sem_11 = '+/-',
                                  mean_12 = '12 years',
                                  sem_12 = '+/-',
                                  mean_16 = '16 years',
                                  sem_16 = '+/-')) %>% 
   colformat_double(digits = 2) %>% 
  add_header_row(values = c('', 'Average +/- SEM By Year'),
                 colwidths = c(2,10), top = TRUE) %>% 
  align(i = 1, j = NULL, part = 'header', align = 'center')

# pivot_wider(names_from = model_type, 
  #              values_from = starts_with("mean_"))
View(structural_sums_cluster_year)

list_save[[7]] <- structural_sums_year
list_save[[8]] <- structural_sums_cluster_year
list_save[[9]] <- result_data


### check for differences in groups by year 
## Interaction Comparisons:
#For interactions, consider pairwise combinations of levels between the two factors.
#Since there are 5 levels in the first factor and 4 in the second, there are 
#5×4 = 32combinations for each level pairing.
#Each of these combinations can interact with each other, 
#so you need to calculate the pairwise combinations of these 32. That's 
#(20 choose 2) = 20 comparisons
filt_aov <- function(dat, filter){
  dat <- filter(lidar_sum, model_type == paste0(filter)) %>% 
    mutate(years_indist = as.factor(years_indist),
           groups = as.factor(groups))
  out <- aov(dat$model_measure ~ dat$groups/factor(dat$years_indist),
             data = dat)
}
lidar_sum <- lidar_sum %>% mutate(groups = new_clustletters[clust])
basal_diff <- filt_aov(lidar_sum, filter = 'basal_model')
basal_diff <-filter(lidar_sum, model_type == 'basal_model') %>% 
  {aov(.$model_measure ~ .$groups/factor(.$years_indist))}

bare_diff <-filter(lidar_sum, model_type == 'bare_ground') %>% 
  {aov(.$model_measure ~ .$groups/factor(.$years_indist))}
comp_diff <-filter(lidar_sum, model_type == 'composition') %>% 
  {aov(.$model_measure ~ .$groups/factor(.$years_indist))}
stem_diff <-filter(lidar_sum, model_type == 'stem_counts') %>% 
  {aov(.$model_measure ~ .$groups/factor(.$years_indist))}


# New Structure Through Time Clusters -------------------------------------

new_struct <- new_struct %>% mutate(groups = new_clustletters[clust])
structure_dat <- dplyr::select(new_struct,-c('cell','clust','new_clusts')) %>% data.table(.) |>
  data.table::melt(id.vars = c('groups', 'years_indist')) |>
  as.data.frame() %>% 
  mutate(years_fct = as.factor(years_indist),
         years_ed = fct_cross(years_fct, groups, 
                              keep_empty = T))
head(structure_dat)
structure_dat$variable <- factor(structure_dat$variable,
                                  labels =  c("Coniferous:Decidious",
                                              "log(Stem~Counts~(stems/900~m^2))",
                                              "Basal~Area~(m^2/hectare)",
                                              "Bare~ground~('%')" ))

Structure_plot_funct <- function(dat, grouping,leg_arg, x_axis = F){ 
  select_ <- filter(dat, variable == grouping)
  ylab <-  paste(grouping)
  if(x_axis == T){ 
    strip = element_text(angle = 45, hjust = 0.5, 
                                    size = 10)
    }
  else{
    strip = element_blank()
  }
  
  ggplot(select_, aes(y= value, group = years_ed)) + 
    geom_boxplot(aes(x= new_clusts , fill = new_clusts, 
                     group = years_ed), width = 1,
                 position = position_dodge2(preserve = 'single'),
                 outlier.shape = NA) +
    facet_wrap(~ years_indist, switch = "x", nrow = 1) +
    #annotate("text", label = ylab, x = 0.5, y = max(select_$value)) +
    theme_classic()+   scale_fill_manual(values = new_colsstruc) +
    guides(size = guide_legend(override.aes = list(size = NULL)),
           color = guide_legend(nrow = 1), 
           linetype = guide_legend(override.aes = list(lty = NULL))) + 
    coord_cartesian()+
    # Customize appearance
    theme(legend.position = leg_arg,
          strip.text = strip,
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x =  element_blank(),# Place the facet label at the bottom
          panel.spacing = unit(2, "lines"),                   # Remove space between facets
          panel.border = element_blank(), # Remove panel borders
          strip.background = element_blank(),                  # Remove strip background
          strip.placement = "outside",
          axis.title.x = element_text(vjust = 0.5),
          plot.margin = unit(c(0.4, 0.0, 0.0, 0.0), "cm"))+
    xlab(NULL) + labs(fill = 'Spectra Clusters Grouped\nby Structural Behaviour')
}
a <- Structure_plot_funct(dat = structure_dat, 
                                grouping = unique(structure_dat$variable)[1],
                          leg_arg = "top")  

b <- Structure_plot_funct(dat = structure_dat, 
                                grouping = unique(structure_dat$variable)[2], 
                                leg_arg = 'none')
c <- Structure_plot_funct(dat = structure_dat, 
                          grouping = unique(structure_dat$variable)[3], 
                          leg_arg = 'none')
d <- Structure_plot_funct(dat = structure_dat, 
                          grouping = unique(structure_dat$variable)[4], 
                          leg_arg = 'none', 
                          x_axis = T) + theme(axis.title.x = element_text(vjust = 6
                                                                          )) + xlab('Years Post Disturbance')
d
plot_labs <- unique(structure_dat$variable) %>% sapply(., FUN = str_replace_all, 
                                                       pattern = "~", 
                                                       replacement = " ")
Structure_plot <- plot_grid(a + theme(legend.position = 'none'),b,c,d, nrow = 4, 
                            align= 'hv', 
                            axis = 'l', 
                            greedy = F, 
                            labels = plot_labs, 
                            label_size = 10, 
                            label_fontface = 'plain', 
                            label_x = 0.35, 
                            label_y = 0.95)
Structure_plot

# save_plot(Structure_plot, filename='F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/new_structures.png',
#           base_width =8, base_height = 10)
#####
#Alternative Structure plot, each cluster group and dynamics through time 

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1] # lower limit
  u <- boxplot.stats(x)$stats[5] # upper limit
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA) # return NA if it's outside those limits
  }
  return(x)
}

dt2 = structure_dat%>% 
  group_by(groups) %>% # filter the outliers by group
  mutate(value = filter_lims(value), # filter the outliers out
         filt = "Filtered",#for knowing that it is filtered when plotting
         value = ifelse(variable == 
                          unique(variable)[2],
                        exp(value), value), ## take the antilog of stems 
         value = ifelse(variable == unique(variable)[4], 
                        1 - value, value),
         variable = fct_relevel(variable, sort)) %>%
  na.omit() # remove outliers

dt2$variable <- factor(dt2$variable,
                       labels = c("Bare~ground~('%')", 
                                  "Bsl~A~(m^2/hct)",
                                  'Conif:Decid',
                                  "Stem~cts~(n/900~m^2)"))

dt2 <- mutate(dt2, new_groups 
                        = case_when(groups == "1" ~ 1,
                                    groups == "C" ~ 2, 
                                    groups == '2' ~ 3, 
                                    groups == '6' ~ 4, 
                                    groups == '7' ~ 5))
clust_names <- c("Low Density\nConifer",
                 "High Density\nConifer",
                 "Shelterwood to\nLate Growth",
                 "Stem Exclusion to\nStem Loss",
                 "Stem Exclusion to\nIngrowth")
dt2$new_groups <- 
  #add factor labels for plotting 
  factor(dt2$new_groups,
         labels = clust_names)

structure_time <- ggplot(dt2 %>% filter(groups %in% c('C', 1))) + 
  geom_boxplot(aes(years_indist, value, group= years_ed,
                   fill = years_fct), 
               position = position_dodge2(preserve = 'single'),
               outlier.colour = NA) + 
  facet_grid(variable~new_groups, scales = 'free', switch = 'both',
             labeller = 
               labeller(.rows = label_parsed, .multi_line = T)) + 
  scale_fill_manual(values= year_cols5) + 
  theme_classic() +
  labs(fill= 'Years Post\nDisturbance') +
  xlab('')+ylab(NULL)+
  scale_x_discrete(drop = F) +
  theme(panel.border = element_rect(fill = NA, color = 'black'),
    legend.position = 'none', 
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 8),
        axis.title.x = element_text(hjust= -0.1, vjust = 10,
                                    size = 10, face = 'bold'),
    axis.text.y = element_text(hjust = 0.2))
structure_time

all_dat <- all_dat %>%
  mutate(year_label = as.factor(if_else(model_type == unique(.$model_type)[2],
                        years_indist,paste(years_indist, 'no_display'))),
         )
# function to suppress labels
delete_no_display <- function(v) {
  if_else(str_detect(v, 'no_display'), '', v)
}
structure_boxplot <- ggplot(all_dat,  
                            aes(x = year_label, y= model_measure)) + 
  geom_boxplot(aes(fill = years_indist), outlier.shape = NA) +
  facet_wrap(~model_type, scales = 'free', labeller = label_parsed,
             nrow = 4,
             switch = 'x') +labs(filler = "Spectra\nCluster") + 
  xlab('All Data') +
  ylab(NULL) + ylab(NULL)+ scale_fill_manual(values = year_cols5) +
  scale_x_discrete(expand = c(0,0), label = delete_no_display) +
  theme_classic() + labs(fill= 'Years in\nDisturbance')+
  theme(panel.border = element_rect(fill = NA, color = 'black'),
        strip.placement = 'outside',
        strip.text = element_blank(),
        axis.title.x = element_text(hjust = -0.1, vjust = 10, size = 10, face = 
                                      'bold'),
        plot.margin = unit(c(0,0,  0.5, 0.5), 'cm'),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.spacing.y = unit(-0.7, 'lines'),
        panel.spacing.x = unit(0.1, 'lines'))
structure_boxplot
library(cowplot)
library(patchwork)
Grouped_clusters <- structure_time +structure_boxplot +
  plot_layout(ncol = 2 ,
            widths = c(2, 0.8)) 

# save_plot(Grouped_clusters, filename='F:/Sync/PhD_Writing/Paper2/test/images/March_Clusters/Cluster_throughtime.png',
#           base_width =8, base_height = 7)


# Structural Grouping Table  ----------------------------------------------
## Here is the data we need from  earlier in this script 
# Most important PC values = PC_var_ful
#Og clusts >- names 'comb_clusters' 
#new groupings -> new_clusters 
#structure data -> new_struct
# spectral metrics  -> data_fil
get_avg_sf <- function(x,metric, name_in) {
  #name_in type of grouping you are looking for 
  x <- as.data.frame(x)
  sel_col <- grep(name_in, names(x))[1]
  x$clust <- x[,sel_col]
  x_out <- x %>% group_by(clust) %>% pivot_longer(ends_with(metric),
                                                  values_to = 'measure', 
                                                  names_to = 'metric') %>% 
    group_by(metric,clust) %>%
    summarise(mean_ = mean(measure, na.rm = T), 
              sd_ = sd(measure, na.rm = T)) %>%
    mutate(max_ = max(mean_),
           min_ = min(mean_)) %>% 
           #clust =ifelse(mean_ == max_ | mean_ == min_, clust,NA),
           #sd_ = ifelse(mean_ == max_ | mean_ == min_, sd_,"")) %>% ## removed because causing errors 
    filter(!is.na(clust))
  return(x_out)
}

## process outline 
spect_dat <- data.table(data_) %>% .[, cell := index]

spect_dat$new_clusts <- as.factor(new_clusters[spect_dat$comb_clusters])
spect_dat$classes <- as.factor(LETTERS[spect_dat$new_clusts])
## reclassify into structure groups 
in_r <- rast()
new_clustnum <-  c(1,2,3,3,3,6,7,3)
new_clustlettersNAN <- c('1','2','C','6','7')
raster_df <- sat_df %>% mutate(z = as.factor(new_clustletters[z])) 
values(sat_clusters) <- raster_df$z
#writeRaster(sat_clusters, filename = 'D:/Paper2_Clean/Spectral_Clusters/Structural_groups.tif')

spect_dat_full <- spect_dat[complete.cases(spect_dat),]

##Get average values 
slopes_in  <- get_avg_sf(spect_dat_full, metric = 'slope', name_in = 'classes')
years_in  <- get_avg_sf(spect_dat, metric = 'vals_years', name_in = 'classes')
regrowth_in  <- get_avg_sf(spect_dat, metric = 'regrowth', name_in = 'classes')


spectral_averages <- bind_rows(slopes_in, years_in, regrowth_in)
spectral_averages<- spectral_averages |> 
  mutate(classes = case_when(clust == 'C' ~ '3-4-5-8', 
                           clust == 'A' ~ '1', 
                           clust == 'B' ~ '2', 
                           clust == 'D' ~ '6', 
                           clust == 'E' ~  '7'))
# Convert the specified columns to numeric
spectral_averages[, c("mean_", "sd_", "max_", "min_")] <- lapply(spectral_averages[, c("mean_", "sd_", "max_", "min_")], as.numeric)

# Divide the specified columns by 100
spectral_averages[, c("mean_", "sd_", "max_", "min_")] <- lapply(spectral_averages[, c("mean_", "sd_", "max_", "min_")],
                                                                 function(x) {round(x/ 100,1)})

spectral_sel <- filter(spectral_averages, metric %in%
                         c(arrange(pca_var_full, desc(scaled_mag))$name)[1:4]) %>% arrange(clust) %>% dplyr::select(-clust)
spect_dat %>% ggplot(.) +
  geom_boxplot(aes(classes, TCA_slope)) + ylim(0, 100)

new_struct$classes <- as.factor(new_struct$new_clusts) 
new_struct <- new_struct |> 
  mutate(classes = ifelse(classes == 3,  '3-4-5-8', classes))
new_struct_sums <- data.table(dplyr::select(new_struct, -c('i.cell')))%>%
  melt(., id.vars = c('classes', 'clust', 'years_indist', 'new_clusts')) %>% 
  .[, .(mean_ = round(mean(value, na.rm = T),1), sd_ = round(sd(value, na.rm = T),1)), by = .(classes, years_indist, variable)] %>% 
  dcast(., classes + variable ~  years_indist, value.var = c('mean_', 'sd_')) 
#reoredr 
col_order <- c(1,2, 3, 7,4,8, 5, 9, 6, 10)
new_struct_sums <- as.data.frame(new_struct_sums) %>% .[, col_order]

spect_struct_big <- cbind(dplyr::select(spectral_sel, -c('max_', 
                                                  'min_')), dplyr::select(new_struct_sums, -classes)) 
dim(spect_struct_big)  
new_order <-  c(4, 1,2, 3, seq(5, 13, by = 1)) 
spect_struct_big <-  spect_struct_big[, new_order]  

# Function to replace underscores with spaces and capitalize words
replace_underscores_and_capitalize <- function(text) {
  if (str_detect(text, pattern = 'sd')) {
    textb <- paste("+/-", str_split(text, pattern = "_")[[1]][3])
  } else {
    texta <- gsub("_+", " ", text)  # Replace one or more underscores with spaces
    textb <- str_to_title(texta)    # Capitalize words
  }
  out <- c(text, textb)
  return(out)
}

text_desc <- read_delim('F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/Cluster_textual_Descriptions.csv', delim = '\\')#[1:5, 1:4]
names(text_desc) <- c('classes', names(text_desc)[2:4])
spect_struct <- full_join(text_desc, spect_struct_big)
names_in <- sapply(names(spect_struct_big), replace_underscores_and_capitalize)
names(spect_struct_big)  <- names_in[seq(2,length(names_in), by =2)]
#names_in[14] <- c("is_last_val_in_group = NA")
names(spect_struct_big) #str(names_in)
spect_struct_big

list_save[[10]] <- spect_struct_big

#save(list_save, file = 'F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/summary_tables.RData')


new_struct_sums_long <- data.table(dplyr::select(new_struct, -c('i.cell')))%>%
  .[, stems := round(exp(.$`Stem~Counts~(stems/900~m^2)`))]%>%
  .[,`Stem~Counts~(stems/900~m^2)`:= NULL] %>% 
  melt(., id.vars = c('clust', 'years_indist', 'new_clusts', 'classes')) %>%
  .[, value := as.numeric(value)] %>% 
  .[,.(mean_ = mean(value, na.rm = T), sd_ = sd(value, na.rm =T)), 
    by = .(new_clusts, years_indist, variable)] %>% 
  .[, .SD[1], by = .(new_clusts, years_indist, variable)] %>% 
  dcast(., new_clusts + variable ~ years_indist, value.var = c('mean_', 'sd_')) %>% 
  arrange(variable) %>% 
  .[, clust := new_clusts] %>% .[, model_type := variable] %>% dplyr::select(., -c('variable', 'new_clusts'))

new_struct_df <-                               
  new_struct_sums_long <- data.table(dplyr::select(new_struct, -c('i.cell')))%>%
  .[, stems := round(exp(.$`Stem~Counts~(stems/900~m^2)`))]%>% .[,`Stem~Counts~(stems/900~m^2)`:= NULL] %>% 
  melt(., id.vars = c('clust', 'years_indist', 'new_clusts', 'classes')) %>% 
  .[,.(mean_ = mean(value, na.rm = T), sd_ = sd(value, na.rm =T)), by = .(new_clusts, years_indist, variable)] %>% 
  .[, .SD[1], by = .(new_clusts, years_indist, variable)] %>% .[variable == 'stems',mean_ := log(mean_)]


