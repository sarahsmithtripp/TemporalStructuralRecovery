# this script runs the proportional models for composition based on field data 
#it is sourced in the script ABALidar_Modeling.R

field_met <- data.frame(Plot_ID = field_data$Plot_ID, 
                        rounded_out = paste(field_data$rounded_out)) %>% 
  mutate(Plot_ID = str_to_sentence(Plot_ID))
decid_model <- filter(field_met, Plot_ID %in% lidar_plots$Plot_ID) %>%
  mutate(y= case_when(rounded_out == 0.0 ~ 0.0001,   
                      rounded_out == 1.0 ~ 0.9999, 
                      is.na(rounded_out)  ~ 0.0001, 
                      TRUE ~ as.numeric(as.character(rounded_out)))) %>% 
  left_join(lidar_plots) 

## run stepwise glm
input_data <- dplyr::select(decid_model, - c('Plot_ID', "rounded_out", 'X')) %>% na.omit()
test <- leaps::regsubsets(y ~ . ,
               data =input_data,    # 1 best model for each number of predictors
               nvmax = )# NULL for no limit on number of variables)
summary(test)$rsq
m1 <- input_data[(names(input_data)[summary(test)$which[9,]])] %>% na.omit()
cor(input_data$Lskew, input_data$Lkurt)
summary(test)
model <- betareg(y ~ ., data = na.omit(m1))
in_w1 <- abs(1/resid(model))
model_2 <- betareg(y ~ pzabovemean + CRR + Lkurt + Lcoefvar, data = na.omit(input_data))
in_w2 <- abs(1/resid(model))
in_w2a <- 1/input_data$CRR
in_w2a  <- 1/fitted(model)^2
in_w2b <- 1/resid(model)^2
model_w2 <- betareg(y ~zkurt + pzabovemean+ CRR + Lkurt + Lcoefvar, data = na.omit(input_data), weights = in_w2, link = c('logit'))
model_w2a <- betareg(y ~zkurt + pzabovemean+ CRR + Lkurt + Lcoefvar, data = na.omit(input_data), weights = in_w2a)
model_w2b <- betareg(y ~ zkurt + pzabovemean+ CRR + Lkurt + Lcoefvar, data = na.omit(input_data), weights = in_w2b)
BIC(model_w2,model_w2a, model_w2b)
#plot(model_w2a)
summary(model_w2)

model_3 <- betareg(y ~ pzabove1.3 + CRR, data = input_data)
in_w3 <- abs(1/resid(model_3))
model_w3 <- betareg(y ~ pzabove1.3 + CRR, data = input_data, weights = in_w3)
model_w4 <- betareg(y ~ pzabove1.3 + zMADmedian + VCI, data = input_data)
summary(model_w4)
## best model is model_w2


# fit k-folds to percent decidoous  ----------------------------------------
library(caret)
folds <-createMultiFolds(y=m1$y,k=5,times=100)
r_square_val <- as.numeric()
rmse_val <- as.numeric()
new_mod_params <- list()
for(i in 1:length(folds)){
  train <- m1[folds[[i]],]
  test <- m1[-folds[[i]],]
  n_low <-  length(which(train$y < 0.001))
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(
    {kfold_model_2 <- betareg(y ~ ., data = train)
            r_square_val <- append(r_square_val, kfold_model_2$pseudo.r.squared)
            rmse_val <- append(rmse_val,
                               (mean((predict(kfold_model_2, test) - test$y), na.rm = T)^2)^0.5)

             ##ad to df
            new_mod_params[[i]] <- kfold_model_2$coefficients$mean %>%
              t() %>%
              data.frame(., phi = kfold_model_2$coefficients$precision)}, error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }
}
mean(r_square_val[1:100], na.rm =T)
mean(rmse_val[1:100], na.rm = T)
params_df_new <- bind_rows(new_mod_params[1:100]) %>% summarise_all(., .funs = c(mean, sd), na.rm = T) %>%
  pivot_longer(cols = everything(.), names_to = c(".value",'type'),
               names_sep = "_.") %>% t() %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>%
  data.frame(.) %>%  rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = c(names(kfold_model_2$coefficients$mean), 'phi'),
         r_square = mean(r_square_val),
         df = kfold_model_2$df.null,
         model_family = paste0('betaregression log link'))

# Create New Updated Model  -----------------------------------------------
in_model <- model
##Update Coefficients and R Square
in_model$coefficients$mean <- params_df_new$est[1:10]
in_model$coefficients$precision <- params_df_new$est[6]
in_model$pseudo.r.squared <- params_df_new$r_square[1]
summary(in_model)

m1$fit <- fitted(model)
plot(m1$fit, m1$y)

decid_model <- right_join(decid_model, input_data['n'])
decid_model$prediction <- predict(in_model, decid_model, type = "response")
