# this script builds stem count models and is called in ABALidar_Modeling.R
# models use a poisson dataset because the dataset is coun t 

## set time to get full time for this script to run 
library(tictoc)
tic()

## Create input data 
stems_model <-  field_data %>% dplyr::select(con_stems, Plot_ID) %>% mutate(Plot_ID = str_to_sentence(Plot_ID)) %>%right_join(rbind(lidar_plots2022, lidar_plots2024), by = 'Plot_ID')


## run stepwise glm 
input_data <- stems_model %>% 
mutate(con_stems = replace_na(.$con_stems, 0), disturbance_year = case_when(startsWith(.$Plot_ID,'22_') ~ 2003,
                                      startsWith(.$Plot_ID,'26_') ~ 2006,
                                      startsWith(.$Plot_ID,'27_') ~ 2006,
                                      startsWith(.$Plot_ID,'72_') ~ 2011,
                                      startsWith(.$Plot_ID,'35_') ~ 2006, 
                                      startsWith(.$Plot_ID,'35_') ~ 2014,
                                      startsWith(.$Plot_ID,'74_') ~ 2014,
                                      startsWith(.$Plot_ID,'53_') ~ 2010,
                                      startsWith(.$Plot_ID,'92_') ~ 2018,
                                      startsWith(.$Plot_ID,'94_') ~ 2018,
                                      startsWith(.$Plot_ID,'88_') ~ 2018, endsWith(.$Plot_ID, '2003') ~ 2003,
                                      endsWith(.$Plot_ID, '2009') ~ 2009, 
                                      endsWith(.$Plot_ID, '2015') ~ 2015, 
                                      endsWith(.$Plot_ID, '2017') ~ 2017)) 


## poisson because it is a count dataset 
library(glmmTMB)
library(bbmle)
library(buildmer)
fit_ <- buildglmmTMB(con_stems ~ pzabovemean + Lkurt + L2 + L3 + ttops + CRR +zmean + VCI + pzabove1.3 + Lcoefvar+ ziqr, data=input_data)
## model is super zero inflacted so fit a zero inflated model -- 0 are both part of the data and also oversaturated 
fit_zibinom <- glmmTMB(con_stems ~ ziqr + CRR + ttops + VCI , 
data= input_data, ziformula = ~ Lkurt, family = nbinom2)
summary(fit_zibinom)
wt1 <- 1/(predict(fit_zibinom, type = "response")^2)
fit_zibinom2w1 <- update(fit_zibinom, weights = wt1)


AICtab(fit_zibinom, fit_zibinom2w1)
## model checkign
# https://aosmith.rbind.io/2017/12/21/using-dharma-for-residual-checks-of-unsupported-models/ 
library(DHARMa)
sim_dat <- simulate(fit_zibinom2w1 , nsim = 1000) %>% collapse()
sim_mod2 <- simulateResiduals(fit_zibinom, n= 100, plot = F)
#compare the residuals for fit_zibinom2
predict_resp <- predict(fit_zibinom2w1, type = 'response')
resid_resp <- resid(fit_zibinom2w1)

## fit_zibion2 is the best, whic includes the pzabove the mean lkurtosis and ttops
pchisq(2 * (logLik(fit_zibinom2w1) - logLik(fit_zibinom)), df = 4, lower.tail = FALSE)
# Dispersion statistic
E2 <- resid(fit_zibinom2w1, type = "pearson")
N  <- nrow(input_data)
p  <- length(coef(fit_zibinom2w1)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

logLik(fit_zibinom2w1)
logLik(fit_zibinom)
##Create null model to get pseudo R squared
null_mod <- glmmTMB(con_stems ~ 1, data = input_data)
McFaddensR <-1-logLik(fit_zibinom2w1)/logLik(null_mod)
MuMIn::r.squaredGLMM(fit_zibinom)
# Fit K Folds to Zipp Binomial Model  -------------------------------------
library(caret)
library(parallel)
folds <-createMultiFolds(y=input_data$con_stems,k=5,times=10)
r_square_val <- list()
new_mod_params <- list()

kfold_fit_function <- function(fold, in_data) {
  train <- in_data[fold,]
  test <- in_data[-fold,]
  count_zero <- length(which(train$con_stems == 0))
  
  if (count_zero >= 9) {
    out <- c('too many zeros')
    print(out)
  } else if (count_zero < 9) {
    skip_to_next <- FALSE
    out <- c('')
    
    tryCatch({
      fit_zibinom2_k <- glmmTMB(con_stems ~ ziqr + CRR + ttops + VCI, 
                                data = train, ziformula = ~ Lkurt, family = nbinom2)
      wt1 <- 1 / (predict(fit_zibinom2_k, type = "response")^2)
      fit_zibinom2w1 <- update(fit_zibinom2_k, weights = wt1)
      test$wt1 <- rep(1, nrow(test))
      
      r_square_val <- cor(test$con_stems, predict(fit_zibinom2w1, newdata = test, type = "response", allow.new.levels = TRUE))^2
      rmse_val <- sqrt(mean((predict(fit_zibinom2w1, newdata = test, type = "response", allow.new.levels = TRUE) - test$con_stems)^2, na.rm = TRUE))
      
      if (!is.na(r_square_val)) {
        new_mod_params <- t(fit_zibinom2w1$fit$parfull)
        out <- list(r_square_val, new_mod_params, rmse_val)
      } else {
        out <- c('mod_not_fit')
      }
    }, error = function(e) {
      skip_to_next <- TRUE
    })
    
    return(out)
  }
}

## Set up a cluster with one worker running on another machine
library(furrr)
library(parallel)
cl <- detectCores()
plan(multicore, workers = cl)
t <-  kfold_fit_function(fold = folds[[1]], in_data = input_data)
kfold_mods <-  map(folds, .f = kfold_fit_function, in_data = input_data) ## Should take about 4 min 
fitted <- which(!grepl('too many', kfold_mods))
selected_kfolds <- kfold_mods[fitted]
rmse <- map(selected_kfolds, .f =function(x) as.data.frame(x[3][1][[1]])) %>% 
  map(., .f = function(x) x %>% mutate(new = paste0(x[1][1][1]),
                                       new = as.numeric(ifelse(new == "", NA, new))) %>% dplyr::select(new)) %>% bind_rows()%>% 
  mutate(new = as.numeric(.$new))
mean_rmse <- mean(rmse$new, na.rm = T)

r_square_vals <- map(selected_kfolds, .f =function(x) as.data.frame(x[1][1][[1]])) %>% 
  map(., .f = function(x) x %>% mutate(new = paste0(x[1][1][1]),
                                       new = as.numeric(ifelse(new == "", NA, new))) %>% dplyr::select(new)) %>% bind_rows()%>% 
  mutate(new = as.numeric(.$new))
mean_r_square <- mean(r_square_vals$new, na.rm = T)

fixed_var <- attributes(unlist(fit_zibinom2_k$modelInfo$terms$cond$fixed))$term.labels
zi_var <- attributes(unlist(fit_zibinom2_k$modelInfo$terms$zi$fixed))$term.labels 

params_df_new <- map(selected_kfolds, .f = function(x) x[[2]] %>% as.data.frame()) %>% purrr::list_rbind()
params_df_named <- params_df_new %>%summarise_all(., .funs = c(mean, sd), na.rm = T) %>% 
  pivot_longer(cols = everything(.), names_to = c(".value",'type'), 
               names_sep = "_.") %>% t()  %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>% 
  data.frame(.) %>%  dplyr::rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = c('(Intercept)', fixed_var, '(Intercept)', zi_var, 'Beta_adj'),
         r_square = mean_r_square, 
         df = fit_zibinom2_k$modelInfo$nobs - (length(fixed_var) + length(zi_var)),
         model_family = fit_zibinom2_k$modelInfo$family$family)

# Create New Model Inputs  ------------------------------------------------
in_model <- fit_zibinom2_k
in_model$fit$par <- params_df_named$est
input_data$wt1 <- rep(1, nrow(input_data))
input_data$prediction <- predict(fit_zibinom2_k, newdata = input_data, type = "response", allow.new.levels = TRUE)
toc()

# final_model <- fit_zibinom2
# list_save <- list(in_model, predictions_stems, params_df_named)
# #saveRDS(final_model, file = "D:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_model.rds")
# future::save_rds(list_save, pathname = "F:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_kfold.rds")
# thing <- readRDS("F:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_kfold.rds")
# 

