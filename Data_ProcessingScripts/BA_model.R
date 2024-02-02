#This script uses a gamma hurdel model to fit a basal area model to the dataset because the data is zero inflated 
# metrics to drop is rerun w/in the scrip becasue lidar plots is augmented with data from collaborator 
#Augmented data details:
  #collected same year in  similiar forest type and location 
  #added to make sure that the data was not limited to the low basal areas we sampled in the field 

## Drop hyper correlated variables 
## plot correlations to see which variables should be dropped 
metrics <- (dplyr::select(lidar_plots, -c(Plot_ID, zmin, n, X)))
M <- cor(metrics, method = 'spearman')
M_mod <- M 
M_mod[upper.tri(M)] <- 0
diag(M_mod) <- 0
#corrplot(M_mod)

drop_metrics <- apply(M_mod, 2, function(x) any(x > 0.5))
drop_metrics
metrics_dropped <- metrics[, !drop_metrics]

## this gives us 6 variables to develop models from which is GREAT!!! 
##Create Model Dataframe 
basal_df <- lidar_plots %>% #left_join(metrics_dropped, dplyr::select(lidar_plots, c('zkurt', 'Plot_ID'))) %>% 
  left_join(., dplyr::select(bsl_df, c('Plot_ID', 'Basal_perhect'))) %>% 
  mutate(y = ifelse(is.na(.$Basal_perhect), 0, 1)) 

# Originally Attemped to a Fit a Lasso Model  -----------------------------
##Build models with these 6 variables LASSO model 
library(glmnet)
x <- as.matrix(dplyr::select(basal_df, names(metrics_dropped)))
y <- basal_df$y

# Fitting the model (Ridge: Alpha = 0)
## sourced from: https://ricardocarvalho.ca/post/ridge/
## and https://www.science.smith.edu/~jcrouser/SDS293/labs/lab10-r.html

## select best nu,ber of lamda for validation
# cv.lasso<- cv.glmnet(x, y, family='binomial', nfolds=5, alpha=0, parallel=TRUE, standardize=TRUE, 
# type.measure='deviance', 
# exact = T)
# 
# # Results
# plot(cv.lasso)
# plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
# bestlambda <- cv.lasso$lambda.min
# cv.lasso$lambda.1se
# lasso_coef = predict(cv.lasso, type = "coefficients", s = bestlambda) # Display coefficients using lambda chosen by CV
# lasso_coef
# coef(cv.lasso, s=cv.lasso$lambda.min)
# 
# lasso_model <- c(lasso_coef)

# Abandoned and fit normal regression model -------------------------------
##Abandon and fit to normal regression model 
#load packages
library(caret)
library(MASS)
library(rms)
library(egg)
library(pROC)
library(bestglm)
library(MASS)
library(MuMIn)

## mostly sourced from https://stats.oarc.ucla.edu/r/dae/logit-regression/
input_data <- dplyr::select(basal_df, !c('Basal_perhect', 'Plot_ID', 'X', 'n'))
train = data.frame(input_data)
dt = sort(sample(nrow(train), nrow(train)*.8))
train_df = train[dt,]
test = train[-dt,]
model <- glm(y ~  ., data = input_data, family = binomial) %>% stepAIC(direction = 'both' ,
trace = T, k = 2)
coef(model)
summary(model)
null_mod_m1 <- glm(y ~ 1, data = input_data, family = binomial)
AIC(model, null_mod_m1)
BIC(model, null_mod_m1)

###UPDATE TEst the prediction outputs
test_mod_m1 <- glm(y ~ Lkurt, data = train, family = binomial) %>% stepAIC(direction = 'both', trace = T, k = 2)
prediction_test <- predict(test_mod_m1, newdata = test, type = 'link', se = T)
## look at the probability output of the data
prediction <- predict(model, newdata = input_data, type = 'link', se = T)
L_kurt_r <- range(input_data$Lkurt)

input_data$prediction <- predict(model, newdata = input_data, type = 'response')

standard_errors_ofModel <- function(fit) {
    out_df <- within(fit, {
        PredictedProb <- plogis(fit)
    LL <- plogis(fit - (1.96 * se.fit))
    UL <- plogis(fit + (1.96 * se.fit))
    })
return(out_df)
}

out_prediction <- cbind(input_data, standard_errors_ofModel(prediction))
test_out <- cbind(test, standard_errors_ofModel(prediction_test))
with(model, null.deviance - deviance)
with(model, pchisq(null.deviance-deviance, df.null - df.residual, lower.tail = FALSE))
logLik(model)


##assess ## of times right 
test_out <- mutate(test_out, predicted_outcome = ifelse(PredictedProb >= 0.5, 1, 0))
test_out$right_wrong <- ifelse(test_out$predicted_outcome == test_out$y, 'Model Right', 'Model Wrong')
z <- (test_out$right_wrong == "Model Right")
sum(z)
#interacccuracy 
22/28

new_fit <- lrm(y ~ Lkurt, data = train)
AUC <- 0.79 ## this is the probability that a random positive is positioned to the right of a random neagit
logistic_plot <- ggplot(test_out, aes(x = Lkurt, y = PredictedProb)) + 
geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.2) + geom_line() + 
geom_point(aes(x = Lkurt, y = y, color = right_wrong), size = 3) + scale_color_manual(values = c("#413629", "#f14a0e")) +
labs(colour = "") + ylab("Probability of Basal Area > 0") + xlab("Estimate of L-moment Kurtosis") +
theme_article(base_size =  12, base_family = "sans") +theme(legend.position = "top")
logistic_plot
## k fol validation  on the logistic regresssion model

folds <-createMultiFolds(y=input_data$y,k=5,times=20)
auc_value<-as.numeric()
auc_old <- as.numeric()
params <- list()
for(i in 1:length(folds)){
  train <- input_data[folds[[i]],]
  test <- input_data[-folds[[i]],]
  model_kFold<-glm(y~ Lkurt,
                   family = binomial(link=logit), 
                   data=train )
  model_pre<-predict(model_kFold,type='response', newdata=test)
  bin <- ifelse(model_pre <= 0.5, 0, 1)
  model_right <- ifelse(test[,'y'] == bin, 1, 0)
  auc_value <- append(auc_value, sum(model_right)/nrow(test))
  ##ad to df
  params[[i]] <- data.frame(int = model_kFold$coefficients[1], 
                            lkrt = model_kFold$coefficients[2])
  ## test og model 
  old_model_pre <- predict(test_mod_m1, newdata=test)
  bin <- ifelse(auc_old <= 0.5, 0, 1)
  old_model_right <- ifelse(test[,'y'] == bin, 1, 0)
  auc_old <- append(auc_old, sum(model_right)/nrow(test))
}

mean(auc_value)
mean(auc_old)
params_df <- bind_rows(params) %>% summarise_all(., .funs = c(mean, sd)) %>% 
  pivot_longer(cols = everything(.), names_to = c(".value",'type'), 
               names_sep = "_.") %>% t() %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>% data.frame(.) %>% 
  dplyr::rename('est' = 'X1', 'stdev'= 'X2') %>% mutate(metric = c(names(model_kFold$coefficients)), 
                                                        auc = mean(auc_value), 
                                                        df = c(model_kFold$df.null))



## input model for 1 is better using lasso to fit glm with lasso model
## sourced from: https://seananderson.ca/2014/05/18/gamma-hurdle/
## second model --- non zero data (the shit that matters )
m2_data <- dplyr::select(basal_df, -c(y, Plot_ID)) %>% 
filter(., Basal_perhect > 0)
# garbage model
m2 <- glm(Basal_perhect ~  ., data = m2_data, family = Gamma(link = log)) 
summary(m2)
#plot(m2)
###weighted to deal with the dispersion 
n <- 18 
null_mod_m2 <- glm(Basal_perhect ~ 1, data = m2_data, family =Gamma(link = log))
m2_w1 <- glm(Basal_perhect ~  pzabovemean + zMADmedian + L3 + Lkurt, data = m2_data, family = Gamma(link = log), 
weights = runif(n)) 
m2_w1_drop <- glm(Basal_perhect ~ pzabovemean + zMADmedian + Lkurt, data = m2_data, family = Gamma(link = log), 
weights = runif(n))
logLik(m2_w1, m2_w1_drop) 
AIC(m2_w1, m2_w1_drop,  m2)

m2_extra <- glm(Basal_perhect ~ pzabovemean + Lkurt, data = m2_data, 
family = Gamma(link = log))

r.squaredGLMM(m2_extra)


### McFaddens R-squared for model 
McFaddensR <-1-logLik(m2_w1_drop)/logLik(null_mod_m2)
MuMIn::r.squaredGLMM(m2_w1_drop)
label_R <- expression('R*2')
with(summary(m2_w1_drop), 1-deviance/null.deviance)

## stay with simple model with significant terms 
shape_gamma <- gamma.shape(m2_w1_drop)

# Model Cross Validation --------------------------------------------------
## k fol validation  on the logistic regresssion model
folds <-createMultiFolds(y=m2_data$Basal_perhect,k=5,times=20)
r_square_val <- as.numeric()
rmse_val <- as.numeric()
new_mod_params <- list()
for(i in 1:length(folds)){
  train <- m2_data[folds[[i]],]
  test <- m2_data[-folds[[i]],]
  n <- nrow(train)
  m2_w1_drop_test <- glm(Basal_perhect ~ pzabovemean + zMADmedian + Lkurt, 
                         data = train, family = Gamma(link = log), 
                  weights = runif(n))
  model_pre <-predict(m2_w1_drop_test,type='response', newdata=test)
  r_square_val <- append(r_square_val, MuMIn::r.squaredGLMM(m2_w1_drop)[1,1])
  rmse_val <- append(rmse_val, 
                     ((mean((predict(m2_w1_drop_test, test) - test$Basal_perhect),
                            na.rm = T)^2) ^ 0.5))
  ##ad to df
  new_mod_params[[i]] <- data.frame(int = m2_w1_drop_test$coefficients[1], 
                                    pszabov = m2_w1_drop_test$coefficients[2],
                                    zMADmed = m2_w1_drop_test$coefficients[3],
                                    Lkurt = m2_w1_drop_test$coefficients[4])

}
mean(r_square_val)
mean(rmse_val)
##Create Export Data Frame 
params_df_new <- bind_rows(new_mod_params) %>% summarise_all(., .funs = c(mean, sd), na.rm = T) %>% 
pivot_longer(cols = everything(.), names_to = c(".value",'type'), 
             names_sep = "_.") %>% t() %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>% 
  data.frame(.) %>%  dplyr::rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = names(m2_w1_drop_test$coefficients),
         r_square = mean(r_square_val), 
         df = m2_w1_drop_test$df.null,
         model_family = paste0(m2_w1_drop_test$family$family, m2_w1_drop_test$family$link))
params_long <- params_df_new %>% mutate(int_upp = est + stdev,
                                        int_lwr = est - stdev) %>% 
  pivot_longer(cols = starts_with('int'), names_to = 'direction', values_to = 'confidence')

# Create New Kfold_Model  -------------------------------------------------
## New non zero gamma model 
in_model <- m2_w1_drop_test
in_model$coefficients <- params_df_new$est

#new bernoulli model 
in_model_zero <- model_kFold
in_model_zero$coefficients <- params_df$est

# Non-zero Gamma model:
##add disturbance_year
basal_df <- basal_df%>% mutate(disturbance_year = case_when(startsWith(.$Plot_ID,'22_') ~ 2003,
                                      startsWith(.$Plot_ID,'26_') ~ 2006,
                                      startsWith(.$Plot_ID,'27_') ~ 2006,
                                      startsWith(.$Plot_ID,'72_') ~ 2011,
                                      startsWith(.$Plot_ID,'35_') ~ 2006, 
                                      startsWith(.$Plot_ID,'35_') ~ 2014,
                                      startsWith(.$Plot_ID,'74_') ~ 2014,
                                      startsWith(.$Plot_ID,'53_') ~ 2010,
                                      startsWith(.$Plot_ID,'92_') ~ 2018,
                                      startsWith(.$Plot_ID,'94_') ~ 2018,
                                      startsWith(.$Plot_ID,'88_') ~ 2018)) %>% 
                                      mutate(disturbance_year = ifelse(is.na(.$disturbance_year), "Undisturbed",.$disturbance_year) )
basal_df$predicted <- predict(in_model, basal_df,type = 'response')
basal_df$residual <- basal_df$Basal_perhect - basal_df$predicted


