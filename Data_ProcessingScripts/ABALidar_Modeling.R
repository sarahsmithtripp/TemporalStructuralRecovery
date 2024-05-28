library(tidyverse); library(corrplot); library(data.table); library(R.utils); library(caret); library(Hmisc)
set.seed(123)

## For proportional models 
library(AER)
library(betareg)
library(leaps)


# Global Functions --------------------------------------------------------

#calculate the mean for pixel 
myfunc <- function(x){
  x_out <- mean(x/100,na.rm = T)
}

# Read in FieldData -------------------------------------------------------

field_data <- read.csv('F:/Quesnel_RPA_Lidar2022/ABA_validation/plot_2022_2024_liam_averagedplots.csv')

lidar_plots2022 <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/lidar_plot_metrics_ttops.csv") %>% filter(!grepl('22 | 53_1', Plot_ID)) %>% 
  dplyr::select(-'X.1')
lidar_plots2024 <-read.csv('F:/Quesnel_RPA_Lidar2022/ABA_validation/Field2024Data/lidar_plot_metrics.csv') %>% 
  dplyr::select(-c('X.1', 'zvar'))
setdiff(names(lidar_plots2024), names(lidar_plots2022))
## Additional Data Gathered from RPA lidar in same forest type, but not fire disturbed 
## added to insure that model captured the full range of forest  types 

##read liam data 
liam_lidar <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/FieldData/liam_lidar_metrics.csv")
liam_field <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/FieldData/liam_bsl.csv") %>% 
  dplyr::rename(Plot_ID = PlotID) %>% dplyr::select(-c(X))

lidar_plots <- rbind(lidar_plots2022, lidar_plots2024)

metrics <- (dplyr::select(lidar_plots, -c(Plot_ID, zmin, n, X)))
M <- cor(metrics, method = 'spearman')
M_mod <- M 
M_mod[upper.tri(M)] <- 0
diag(M_mod) <- 0
corrplot::corrplot(M_mod)

drop_metrics <- apply(M_mod, 2, function(x) any(x > 0.7))
drop_metrics

metrics_dropped <- metrics[, !drop_metrics] %>% cbind(lidar_plots['Plot_ID'])

# Composition and Bare Ground Models --------------------------------------
source('D:/Paper2_Clean/Data_ProcessingScripts/CompositionModel.R')

## Graphs to Pull out 
## Update Predictions
decid_model$prediction <- predict(in_model, decid_model, type = "response")

decid_model <- decid_model%>% mutate(disturbance_year = case_when(startsWith(.$Plot_ID,'22_') ~ 2003,
                                                                  startsWith(.$Plot_ID,'26_') ~ 2006,
                                                                  startsWith(.$Plot_ID,'27_') ~ 2006,
                                                                  startsWith(.$Plot_ID,'72_') ~ 2011,
                                                                  startsWith(.$Plot_ID,'35_') ~ 2006,
                                                                  startsWith(.$Plot_ID,'35_') ~ 2014,
                                                                  startsWith(.$Plot_ID,'74_') ~ 2014,
                                                                  startsWith(.$Plot_ID,'53_') ~ 2010,
                                                                  startsWith(.$Plot_ID,'92_') ~ 2018,
                                                                  startsWith(.$Plot_ID,'94_') ~ 2018,
                                                                  startsWith(.$Plot_ID,'88_') ~ 2018))

decid_model$residual <- decid_model$y - decid_model$prediction
predictions_comp <- ggplot(decid_model, aes(group= as.factor(disturbance_year))) +
  geom_point(aes(log10(y+1), log10(prediction+1)))+  scale_y_continuous(limits = c(0, 0.3)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab('% measured coniferous') +
  theme_classic(base_size = 12) + ylab('% modeled coniferous') + labs(color = 'Disturbance Year')
predictions_comp

comp_graph <- ggplot(decid_model) +
  geom_point(aes(log10(y+1), log10(prediction+1)), size = 3)+  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0,0.3)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab('Measured (Coniferous:Decidious)') +
  egg::theme_article() + ylab('Predicted (Coniferous:Decidious)')

# Save Deciduous Model Outputs  
# save_list <- list(in_model, params_df_new, comp_graph)
# saveRDS(final_model, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/composition_model.rds")
# saveRDS(save_list, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/composition_kfolmodel.rds")


# Bare Ground Model ---------------------------------------------------------
source('D:/Paper2_Clean/Data_ProcessingScripts/Bare_Groud_Model.R')

#Graphs to pull out 
## Update Predictions
bare_ground_model$prediction <- predict(in_model, bare_ground_model, type = "response")

predictions_bare_ground_model <- ggplot(bare_ground_model)  +
  geom_point(aes(log10(Prop.Bare+1), log10(prediction+1)))+  scale_y_continuous(limits = c(0, 0.3)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab('% measured bare ground') +
  theme_classic(base_size = 12) + ylab('% modeled bare ground')
predictions_bare_ground_model


# final_model <- in_model
# save_list <- list(in_model, params_df_new_bar, predictions_bare_ground_model)
# saveRDS(save_list, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/bare_ground_model.rds")


# Basal Area Model  -------------------------------------------------------
## read in specific basal area field plot data (in own excel sheet because it has DBH and tree heights etcs.)
BA_plots <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/FieldData/average_plot_metrics.csv") %>% filter(!grepl('22', Plot_ID))
lidar_plots <- rbind(dplyr::select(lidar_plots, -'ttops'), liam_lidar) ## Remove ttops which are not necessary for BA model 
bsl_df <- rbind(liam_field, field_plots[,c('Plot_ID', 'Basal_perhect')])

## These models are gamma hurdle models because the data is zero-inflated 
## original script set a new seed 
set.seed(999)
source('D:/Paper2_Clean/Data_ProcessingScripts/BA_model.R')

#graphs 
xl <-  expression("Measured Basal Area (m"^2*"/ hectare)")
yl <- expression("Predicted Basal Area (m"^2*"/ hectare)")
predictions_bsl <- ggplot(basal_df, aes(color = as.factor(disturbance_year))) + geom_point(aes(Basal_perhect, predicted), color = 'black', size = 3) +
  #geom_point(aes(Basal_perhect, predicted_kfold, shape = 'Kfold Model Estimate'), color = 'blue', size = 3 )
geom_abline(slope = 1, intercept = 0) + labs(color = "Disturbance Year") +
xlab(xl) + ylab(yl) +
xlim(0,30) + ylim(0,30) +egg::theme_article()

residuals <- ggplot(basal_df %>% mutate(id = rownames(.)),
 aes(color = as.factor(disturbance_year))) +
geom_point(aes(x = Plot_ID, y= residual)) + theme_bw()

hurdle_models <- list(in_model, in_model_zero, params_df_new, params_df, predictions_bsl)

# saveRDS(hurdle_models, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/hurdle_out.rds")
# saveRDS(hurdle_models, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/hurdle_out.rds")

# Stem Count Model  -------------------------------------------------------
set.seed(123) ## reset seed to stay consistent
# read in correct stem count data 
## Read data into R (drop site 22 from model)
field_plots <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/FieldData/stems_per_plot.csv") %>% filter(!grepl('22', Plot_ID))
lidar_plots <- read.csv("D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/lidar_plot_metrics_ttops.csv") %>% filter(!grepl('22', Plot_ID))

field_plots <- field_plots %>% mutate(stems_perm2 = round(.$con_stems/900))
#source('D:/Paper2_Clean/Data_ProcessingScripts/CountModels.R') # Will take some time to run 
graph_data  <- input_data

predictions_stems <- ggplot(graph_data) + 
  geom_point(aes(con_stems, prediction), size = 3) +
  geom_abline(slope = 1, intercept = 0) + theme_minimal() + 
  xlab('Measured (# Stems / Plot)') + ylab('Predicted (# Stems / Plot)')  + egg::theme_article() 


residuals <- ggplot(stems_model %>% mutate(id = rownames(.)),
                    aes(color = as.factor(disturbance_year))) +
  geom_point(aes(x = Plot_ID, y= residual)) + theme_bw()

in_model <- fit_zibinom2
in_model$fit$par <- params_df_named$est
final_model <- fit_zibinom2
#saveRDS(final_model, file = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/stem_count_model.rds")
#future::save_rds(list_save, pathname = "D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/stem_count_kfold.rds")
#thing <- readRDS("F:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_kfold.rds")



