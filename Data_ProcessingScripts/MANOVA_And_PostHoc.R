### This script runs the MANOVA to determine the overlap in our data and the pairwise comparison 
## it also creates a table to import into the final portion of paper 2 

# Read in Libraries and setwd -------------------------------------------------------
 setwd('D:/Paper2_Clean')
# Load the necessary library
# The vegan package contains functions for community ecology analysis, including the adonis function which is used to perform PERMANOVA.
library(vegan)
library(RVAideMemoire)
library(permute)
library(corrplot)
# Assuming 'data' is your data frame and it has columns 'Species', 'Year', 'Variable_A', 'Variable_B', 'Variable_C', and 'Variable_D'
ifelse(exists('struct_pc'),
       struct_pc <- struct_pc,
       struct_pc <- arrow::read_feather(paste0(getwd(),'/Spectral_Clusters/Structural_Data/Structural_Estimates.feather')))

### calculate distance matrix (defaul to read in from file uncomment as per necessary)
# STEP 1: Calculate the distance matrix
# The vegdist function is used to calculate a distance matrix from the data.
# consider columsn 4 to 7 for this calculation, which are where structure values are stored
## UNCOMMENT HERE TO RERUN
#dist_matrix <- vegdist(struct_pc[, 4:7], method = "euclidean") ## calculate euclid

dist_matrix <- readRDS(file = 'D:/Paper2_Clean/Spectral_Clusters/Structural_Data/euc_dist_mat.rds')

# STEP 2: Run the PERMANOVA test using the adonis function
# The adonis function is used to perform the PERMANOVA test.
# The formula argument specifies the model to be used in the analysis, with the distance matrix as the response variable and the species as the explanatory variable.
# The data$clust extracts the 'clust' column from the data frame
# and uses clusters to group in anaylsis use as the grouping variable in the analysis.
## RUN ONCE WILL TAKE ~ 1 day
#permanova_result <- adonis2(dist_matrix ~ struct_pc$clust)

# STEP 3: Create Weighted Sample (for post-hoc analyses)
## create a weighed sample for post-hoc analyses
weighted_sample <- struct_pc %>%
  group_by(clust) %>%
  mutate(n_clust = n(),
         n_prop = n_clust/nrow(.)) %>%
  sample_frac(0.15)

#dist_sample <- vegdist(weighted_sample[, 4:7], method = 'euclidean')
#saveRDS(dist_sample, file = 'D:/Paper2_Clean/Spectral_Clusters/Structural_Data/WeightedSamples/euc_dist_mat.rds')

#STEP 4: Run both a nested MANOVA and the ANOVA
#https://archetypalecology.wordpress.com/2018/02/21/permutational-multivariate-analysis-of-variance-permanova-in-r-preliminary/
dist_sample <- readRDS('D:/Paper2_Clean/Spectral_Clusters/Structural_Data/WeightedSamples/euc_dist_mat.rds')
bird.div<-adonis2(dist_sample~weighted_sample$clust, permutations = 999,
                  perm= how(blocks = weighted_sample$clust,
                            plots = plots(strata = weighted_sample$year))) ## permute within years

bird.div2<-adonis2(dist_sample~weighted_sample$clust * weighted_sample$years_indist, permutations = 999)## permute within years

##create flextable output for manuscript
Test_stat <- rownames(bird.div2[1]) %>% str_remove_all(
  pattern = "weighted_sample"
) %>% str_replace_all(., pattern = "\\$", replacement = " ")
MANOVA <-  data.frame(Df = bird.div2$Df,
                      SumSq = bird.div2$SumOfSqs,
                      `F` = bird.div2$`F`,
                    P_val = bird.div2$`Pr(>F)`)%>%
  round(., digits = 2) %>%
  mutate(P_val = ifelse(P_val ==0, '<< 0.01',P_val))
MANOVA <- cbind(Test_stat, MANOVA) %>%
  data.frame()# add names to the first col

#write.csv(MANOVA,'F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/MANOVA_tab.csv')

# STEP 5: Run Pairwise results for the permanova to combine like clusters
# p  < 0.001, thus we can explore differences by year
pairwise_results <- pairwise.perm.manova(dist_sample, weighted_sample$clust,
                                         test = 'wilks')

pairwise_df <- reshape2::melt(pairwise_results$p.value, na.rm = T) %>%
  mutate(value = round(value, 2))

yint_df <- data.frame(Var1 = seq(2,8, by = 1),
                      labs = seq(1,7, by = 1),
                      Var2 = seq(2, 8, by = 1))

ggplot(data = pairwise_df, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), color = 'grey9') + scale_fill_gradient(low="#ffffff",
                                                     high="#e3221d")+
  geom_text(aes(label = value)) +
  geom_text(data = yint_df, aes(x= Var1, y= Var2, label = labs)) +
  scale_x_discrete(limits = as.character(seq(1,8,by = 1))) +
  scale_y_discrete(limits = as.character(seq(1,8, by = 1))) +
  geom_text(aes(x = 3, y = 6, label = "Post-Hoc Comparisons\nFor Spatiotemporal Differences in\nSpectral Clusters"),
            size = 3)+
  theme(legend.position = "none",
        plot.background = element_rect(fill = NULL),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  xlab('Spectral Cluster')

## Save to a folder
#ggsave(filename ='F:/Sync/PhD_Writing/Paper2/test/images/September_Clusters/cluster_simularities.png', last_plot())
new_clusters <- c(1, 2, 3,3,3,6,7,3)
print("New Clusters Created and added to vector 'new_clusters'")
