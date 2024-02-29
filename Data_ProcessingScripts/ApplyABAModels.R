### build lidar predictions to correlate to satellites
library(tidyverse); library(terra)
library(betareg); library(glmmTMB)

##### SOME FUNCTIONS 
#rename a raster 
renam_funt <- function(r, name){
  names(r) <- paste0(name)
  return(r)
}

## count number of trees in each raster cell function
count_funct <- function(trees, lidar_rast){
  trees$Z <- as.character(trees$Z)
  tree_cn <- terra::rasterize(x = trees,y = lidar_rast,
                              field = "treeID", fun = "length")
  lidar_rast_out <- c(lidar_rast, tree_cn)
  return(lidar_rast_out)
 }

## remove datasets not currently interested in 
remove_dtm <- function(file_list) {
    file_return <- file_list[which(!grepl(pattern = "basal", file_list))]
    return(file_return)
}

## READ In DATA 
#1. Read in models 
file_list <- list.files("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation", 
                        full.names = T, pattern = 'kfol')
file_list
hurdle_models <- readRDS(file_list[[2]])
count_models <- readRDS(file_list[[3]])
composition_model <- readRDS(file_list[[1]])
base_model <- readRDS('F:/Quesnel_RPA_Lidar2022/ABA_validation/bare_ground_model.rds')
# hurdle_models <- readRDS("D:/Quesnel_RPA_Lidar2022/ABA_validation/hurdle_out.rds")
# count_models <- readRDS("D:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_model.rds")
# composition_model <- readRDS("D:/Quesnel_RPA_Lidar2022/ABA_validation/composition_model.rds")


#2. Read in lidar data 
base_dir <- "D:/Paper2_Clean/RPA_data/ABA_Models/Model_Inputs"
blocks_dir_base <- list.dirs(base_dir, recursive = F)  
## isolate lidar data
lidar_dir <- blocks_dir_base %>% 
as.list(strsplit(., ",")) %>%
purrr::map(., .f = list.files,
pattern = '30m.tif$', full.names = T) %>%
Filter(length, .)

### Momentary filter to the area we are interested in applying the models 
# lidar_dir <- lidar_dir[11:12]
# blocks_dir_base <- blocks_dir_base[11:12]
#As a raster
lidar_rast <- lidar_dir %>%purrr::map(., .f = remove_dtm) %>% 
purrr::map(., .f = rast)
## standardize names
names_std <- c("basic_30m_1", "basic_30m_2", "basic_30m_3", "basic_30m_4", 
               "basic_30m_5", "basic_30m_6", "basic_30m_7", "basic_30m_8", 
               "canopy_disp_30m_1", "canopy_disp_30m_2", "canopy_disp_30m_3",
               "canopy_disp_30m_4", "canopy_disp_30m_5", "canopy_disp_30m_6",
               "L_moments_30m_1", "L_moments_30m_2", "L_moments_30m_3", 
               "L_moments_30m_4", "L_moments_30m_5", "L_moments_30m_6", 
               "L_moments_30m_7", "percent_30m_1", "percent_30m_2", 
               "percent_30m_3", "percent_30m_4")

lidar_rast <- purrr::map(lidar_rast, .f = function(x) { 
  set.names(x, names_std)
  return(x)})
lidar_rast[[1]]
#3. Read in ttops shapefiles 
tree_shps <- blocks_dir_base %>%
 as.list(strsplit(., ",")) %>%
purrr::map(., .f = list.files, recursive = TRUE, pattern = "chm_ttops",
full.names = TRUE) %>%
Filter(length, .) %>% purrr::map(., .f = vect) 

##fix CRS of plot  (which is in some weird projection) 
## Set a project crs
crs_proj <- crs(lidar_rast[[1]])

crs_new <- rast(tree_shps[[3]])
crs_old <- crs(rast(tree_shps[[4]]))
crs(tree_shps[[4]]) <- crs_new


## transform any raster that does not meet the necessary standards
## input an index raster 
lidar_rast <- purrr::map(lidar_rast, .f = function(x){
  if (crs(x) == crs_proj) {
    x  <- x
  } else if (crs(x) != crs_proj){
    x <- project(x, crs_proj, res = 30)
  }
  return(x)
})

tree_shps <- purrr::map(tree_shps, .f = function(x){
  if (crs(x) == crs_proj) {
    x  <- x
  } else if (crs(x) != crs_proj){
    x <- project(x, crs_proj)
  }
  return(x)
})




#4. Count number of trees in each raster cell output for ttops function 
lidar_withtrees <-purrr::map2(tree_shps, lidar_rast, .f = count_funct, .progress = TRUE)
names(lidar_withtrees[[4]]) <- names(lidar_withtrees[[1]])

## APPLY MODEL 1 - PROBABILITY OF BASAL AREA 

#1: Subset Lidar layers 
model_1_lidar <- lidar_rast %>% purrr::map(., .f = subset, c("L_moments_30m_6")) %>% purrr::map(., .f = setNames, 
"Lkurt") # %>%purrr::map(., .rf= clamp, lower = -0.1790733, upper = 0.73209, values = TRUE)
#2: predict values  and rename model to something that makes sense 
predict_model_1 <- model_1_lidar %>% purrr::map(., predict, model = hurdle_models[[2]], type = "response")  %>% 
purrr::map(., .f = clamp, lower = 0, upper = 1) %>% purrr::map(., .f = renam_funt, name = "PBsl")

#3. Name model to something PBsl for probability of measured Basal area 

#tests for model 1 (not run) 
# predfun <- function(model, data) {
  # v <- predict(model, data, se.fit=TRUE)
  # cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
# }

# test_model <- values(model_1_lidar[[1]]) * coef(hurdle_models[[1]])[2] + coef(hurdle_models[[1]])[2]
# logit2prob <- function(logit){
  # odds <- exp(logit)
  # prob <- odds / (1 + odds)
  # return(prob)
# }
# test_logit <- logit2prob(test_model)
# range(test_logit, na.rm = TRUE)

#plot(aggregate(predict_model_1[[10]], 17, fun = mean))
#model_1_lidar <- lidar_rast r%% purrr::map(., predict,  model = hurdle_models[[1]])
# save_model_1 <- function(raster, blocks) {
#   writeRaster(raster, filename = paste0(blocks, "/Lidar/processed/05_raster/prob_basal.tif"), 
#   overwrite =T)
# }

## APPLY MODEL 2 - BASAL AREA 
#1. Subset lidar layers 
model_2 <- hurdle_models[[1]]
select_layers <- c("pzabovemean", "zMADmedian", "Lkurt")
model_2_lidar <- lidar_rast %>% purrr::map(., .f = subset, 
c("percent_30m_1",#pzbovemean
 "canopy_disp_30m_3", #zMadMedian
"L_moments_30m_6")) %>% #lKurt
purrr::map(., .f = setNames, select_layers)

#2. Fit model
fit_model_2 <-purrr::map(model_2_lidar, .f = predict, model_2) %>%
  purrr::map(., .f = clamp, lower = 0) %>%
  purrr::map(., .f = renam_funt, name = "basal_model")

# APPLY MODEL 3 - STEM MODELS 
# Define a function to test by hand what the predict function is doing with log link 
model_3 <- count_models[[3]]
tmp <- summary(model_3) #doing this manually
values_mod_3 <- exp(model_3[,'est'])
values_mod_3
model_3_funct <- function(lidar_r) {
  lidar_out <-values_mod_3[1] + (values_mod_3[2] * clamp(lidar_r[select_layers[1]], upper = 44.19 )) +
  (values_mod_3[3] * lidar_r[select_layers[2]]) + 
  (values_mod_3[5] * lidar_r[select_layers[3]]) + 
  (values_mod_3[4] * lidar_r[select_layers[4]])
  class(lidar_out)
return(lidar_out)
}

#1. Select layers for the model 
select_layers <- c('pzabovemean', 'Lkurt', 'ttops', 'zmax', 'CRR')
model_3_lidar <- lidar_withtrees %>% purrr::map(., .f = subset, 
                                                c("percent_30m_1", # perc above mean 
"L_moments_30m_6",  ## Lkurt
"treeID_length",  ## ttops
"basic_30m_2", #zmax 
"canopy_disp_30m_4")) %>%  ## CRR
purrr::map(., .f = setNames, select_layers)

#Check model fit here: 
	# model_3$frame %>% summarize_all(., .funs = c('min', 'max'))

	# model_3_sel <- model_3_lidar[[1]] %>% as.data.frame() %>% summarise_all(., .funs = 
	# c('min', 'max'))
	# model_3_sel
	# compare_fit <-purrr::map(model_3_lidar, .f = model_3_funct)
#2. Fit the model and rename to 'stem counts'
fit_model_3 <- purrr::map(model_3_lidar, .f = predict, model = count_models[[1]]) %>% 
 purrr::map(., .f = renam_funt, name = "stem_counts") %>% purrr::map(., ~  exp(.)) #%>% 
 #purrr::map(., .f = clamp, upper = 1800) ## number of stems that can be in a space 


## APPLY MODEL 4 - COMPOSITION MODEL  
model_4 <- composition_model[[1]]
##testing 
#model_4 <- model_w2
#1. Select lidar layers
select_layers <- c("pzabovemean", "CRR", "Lkurt", "Lcoefvar")
model_4_lidar <- lidar_rast %>% purrr::map(., .f = subset, c("percent_30m_1", # pzabove mean
"canopy_disp_30m_6", # CRR
"L_moments_30m_6",  ## Lkurt
"L_moments_30m_7" # Lcoefvar"
)) %>%purrr::map(., .f = setNames, select_layers)
#2. Predict model
fit_model_4 <-purrr::map(model_4_lidar, .f = predict, model_4) %>%purrr::map(., .f = renam_funt, "composition")#%>%purrr::map(., .f = clamp, lower = 0)

## APPLY MODEL 5 - BARE GROUND MODEL 
model_5 <- base_model[[1]]
#1. Select lidar layers
select_layers <- c("pzabove0.2", "pzabove1.3", "CRR", "zentropy", "Lcoefvar")

model_5_lidar <- lidar_rast %>% purrr::map(., .f = subset, c("percent_30m_2", # pzabove 0.2,
"percent_30m_3", # perc above 1.3
"canopy_disp_30m_4", # CRR
"canopy_disp_30m_5", #zentropy 
"L_moments_30m_7" # Lcoefvar"
)) %>% purrr::map(., .f = setNames, select_layers)
#2. Predict model
fit_model_5 <-purrr::map(model_5_lidar, .f = predict, model_5) %>%
purrr::map(., .f = renam_funt, "bare_ground")#%>%purrr::map(., .f = clamp, lower = 0)

######  join modeling data 
# Get year of disturbance based on site number
site_list <- unlist(blocks_dir_base) #%>% str_extract(., pattern = "Site[[:digit:]][[:digit:]].*") ## Get the site ID 
dist_year_list <- 
  c('2003','2006','2006','2006','2014','2010','2011', '2014','2018','2018','2018', 
    '2011',## for the west fire, this was flow in 2021, so shift by 1 year to make sure that the correct 
    '2011',## 'time' since disturbance is applied 
    '2011') ## these are the disturbance years
attach_site_df <- data.frame(Plot_ID = site_list,  dist_year = dist_year_list) ## create a dataframe

model_data_frame <-purrr::map2(fit_model_2, fit_model_3,.f = function(x,y) c(x,y)) %>%
 purrr::map2(., .y = fit_model_4, .f = function(x,y) c(x,y)) %>%
  purrr::map2(., .y = predict_model_1, .f = function(x,y) c(x,y)) %>% 
  purrr::map2(., .y = fit_model_5, .f = function(x, y) c(x,y))

### recreate as a raster
for(i in 1:length(model_data_frame)){
  model_data_frame[[i]]$dist_year <- rep(dist_year_list[i], ncell(model_data_frame[[i]]))
}
out <- sprc(model_data_frame)
m <- merge(out)
writeRaster(m, filename = 'D:/Paper2_Clean/RPA_data/ABA_Models/Model_Outputs/lidar_rast2.tif',
            overwrite = T)
plot(m)

#### SAVE RASTERS OF MODELS 
# for(i in 1:10){
#   writeRaster(model_data_frame[[i]], filename = paste0(blocks_dir_base[[i]], '/Lidar/four_applied_models.tif'),
#   overwrite = TRUE)
# }
model_data_frame <- model_data_frame %>% purrr::map(., .f = as.data.frame, xy = TRUE)


## Add Site # and Disturbance year to TABLE of estimated values 
for(i in 1:length(model_data_frame)) {  
  model_data_frame[[i]]$Plot_ID = attach_site_df[i, 1]
  model_data_frame[[i]]$dist_year = attach_site_df[i,2]
}
lidar_models <- bind_rows(model_data_frame) %>% 
mutate(meas_bsl = ifelse(.$PBsl >= 0.5, 'Bsl', 'No_Bsl'))
dim(lidar_models)
write.csv2(lidar_models, 
           "D:/Paper2_Clean/RPA_data/ABA_Models/Model_Outputs/lidar_models_df.csv")



