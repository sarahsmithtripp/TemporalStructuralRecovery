# this script builds stem count models and is called in ABALidar_Modeling.R
# models use a poisson dataset because the dataset is coun t 

## set time to get full time for this script to run 
library(tictoc)
tic()

## Drop hyper correlated variables 
metrics <- (dplyr::select(lidar_plots, -c(Plot_ID, zmin, n, X)))
M <- cor(metrics, method = 'spearman')
M_mod <- M 
M_mod[upper.tri(M)] <- 0
diag(M_mod) <- 0
#corrplot(M_mod)

drop_metrics <- apply(M_mod, 2, function(x) any(x > 0.7))

metrics_dropped <- metrics[, !drop_metrics]
## Create input data 
stems_model <- left_join(metrics_dropped, dplyr::select(lidar_plots, c('zmax', 'Plot_ID'))) %>% 
left_join(., dplyr::select(field_plots, c('Plot_ID', 'con_stems'))) 

## run stepwise glm 
input_data <- dplyr::select(stems_model, - c("X.1")) %>% 
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
                                      startsWith(.$Plot_ID,'88_') ~ 2018))


## poisson because it is a count dataset 
library(glmmTMB)
library(bbmle)

## model is super zero inflacted so fit a zero inflated model -- 0 are both part of the data and also oversaturated 
fit_zibinom<- glmmTMB(con_stems ~ pzabovemean + Lkurt + zmax + CRR + ttops + (1|disturbance_year), 
data= input_data, ziformula = ~ CRR, family = nbinom2)
fit_zibinom1 <- update(fit_zibinom, . ~ pzabovemean + Lkurt + zmax + ttops + (1|disturbance_year))
fit_zibinom2 <- update(fit_zibinom, . ~ pzabovemean + Lkurt + zmax + ttops, se = TRUE)
wt1 <- 1/(predict(fit_zibinom2, type = "response")^2)
fit_zibinom2w1 <- update(fit_zibinom2, weights = wt1)
fit_zibinom2a <- update(fit_zibinom, . ~ pzabovemean + Lkurt + ttops) 
fit_zibinom2b <- update(fit_zibinom, . ~ zmax + Lkurt + ttops)
fit_zibinom4 <- update(fit_zibinom, . ~ pzabovemean + Lkurt + ttops + (1|disturbance_year))

fit_zibinom3 <- update(fit_zibinom, . ~ pzabovemean + Lkurt + CRR + ttops + (1|disturbance_year))

fit_zibinom5 <- update(fit_zibinom, . ~ pzabovemean + Lkurt + (ttops| disturbance_year))
fit_zipois <- update(fit_zibinom, family = poisson)

AICtab(fit_zibinom, fit_zibinom1, fit_zibinom2, fit_zipois, fit_zibinom3, fit_zibinom4, fit_zibinom2a, fit_zibinom2b)
## model checkign
# https://aosmith.rbind.io/2017/12/21/using-dharma-for-residual-checks-of-unsupported-models/ 
library(DHARMa)
sim_dat <- simulate(fit_zibinom2 , nsim = 1000) %>% collapse()
sim_mod2 <- simulateResiduals(fit_zibinom2, n= 100, plot = F)
sim_mod4 <- simulateResiduals(fit_zibinom4)
#compare the residuals for fit_zibinom2
predict_resp <- predict(fit_zibinom2, type = 'response')
resid_resp <- resid(fit_zibinom2)

## fit_zibion2 is the best, whic includes the pzabove the mean lkurtosis and ttops
pchisq(2 * (logLik(fit_zibinom2) - logLik(fit_zibinom4)), df = 4, lower.tail = FALSE)
# Dispersion statistic
E2 <- resid(fit_zibinom2, type = "pearson")
N  <- nrow(input_data)
p  <- length(coef(fit_zibinom2)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

logLik(fit_zibinom2)
logLik(fit_zibinom4)
##Create null model to get pseudo R squared
null_mod <- glmmTMB(con_stems ~ 1, data = input_data)
McFaddensR <-1-logLik(fit_zibinom2)/logLik(null_mod)
MuMIn::r.squaredGLMM(fit_zibinom2)
# Fit K Folds to Zipp Binomial Model  -------------------------------------
library(caret)
library(parallel)
folds <-createMultiFolds(y=input_data$con_stems,k=5,times=1000)
r_square_val <- list()
new_mod_params <- list()

kfold_fit_function <- function(fold, in_data)
{
  train <- in_data[fold,]
  test <- in_data[-fold,]
  count_zero <- length(which(train$con_stems == 0))
  if (count_zero >= 9) {
    out <- c('too many zeros')
    print(out)
  }
  else if (count_zero < 9) {
    skip_to_next = FALSE
    out <- c('')
    tryCatch({
      fit_zibinom2_k <-
        glmmTMB::glmmTMB(con_stems ~ pzabovemean + Lkurt + zmax + ttops, se = TRUE, ziformula = ~ CRR, 
                  family = nbinom2, data = train)
      r_square_val <-
        cor(test$con_stems, predict(fit_zibinom2_k, test)) ^ 2
      rmse_val <- ((mean((predict(fit_zibinom2_k, test) - test$con_stems), na.rm = T)^2) ^ 0.5)
      ##ad to df
      if (is.na(r_square_val) == FALSE) {
        new_mod_params <- fit_zibinom2_k$fit$parfull %>%
          t()
        out <- list(r_square_val, new_mod_params, rmse_val)
      }
      else if (is.na(r_square_val) == TRUE) {
        out <- c('mod_not_fit')
      }},
      error = function(e) {
        skip_to_next = T
      })
    return(out)
  }
}
## Set up a cluster with one worker running on another machine
library(furrr)
cl <- detectCores()
plan(multicore, workers = cl)
kfold_mods <- future_map(folds, .f = kfold_fit_function, in_data = input_data) ## Should take about 4 min 
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

fixed_var <- attributes(unlist(fit_zibinom2$modelInfo$terms$cond$fixed))$term.labels
zi_var <- attributes(unlist(fit_zibinom2$modelInfo$terms$zi$fixed))$term.labels 

params_df_new <- map(selected_kfolds, .f = function(x) x[[2]] %>% as.data.frame()) %>% purrr::list_rbind()
params_df_named <- params_df_new %>%summarise_all(., .funs = c(mean, sd), na.rm = T) %>% 
  pivot_longer(cols = everything(.), names_to = c(".value",'type'), 
               names_sep = "_.") %>% t()  %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>% 
  data.frame(.) %>%  dplyr::rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = c('(Intercept)', fixed_var, '(Intercept)', zi_var, 'Beta_adj'),
         r_square = mean_r_square, 
         df = fit_zibinom2$modelInfo$nobs - (length(fixed_var) + length(zi_var)),
         model_family = fit_zibinom2$modelInfo$family$family)

# Create New Model Inputs  ------------------------------------------------
in_model <- fit_zibinom2
in_model$fit$par <- params_df_named$est

toc()

# final_model <- fit_zibinom2
# list_save <- list(in_model, predictions_stems, params_df_named)
# #saveRDS(final_model, file = "D:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_model.rds")
# future::save_rds(list_save, pathname = "F:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_kfold.rds")
# thing <- readRDS("F:/Quesnel_RPA_Lidar2022/ABA_validation/stem_count_kfold.rds")
# 

