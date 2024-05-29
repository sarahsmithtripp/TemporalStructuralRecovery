# this script uses field data to calculate the bare ground proportion for ABA models 

# Bare ground model  -----------------------------------------------
#Big bare ground model

field_met <- data.frame(Plot_ID = field_data$Plot_ID, 
           Prop.Bare = paste(field_data$Prop.Bare)) %>% 
  mutate(Plot_ID = str_to_sentence(Plot_ID))
bare_ground_model <- left_join(select(field_met, Prop.Bare, Plot_ID), metrics_dropped, by = c('Plot_ID' = 'Plot_ID')) %>% 
  mutate(y = round((as.numeric(Prop.Bare)/100), 4))

## run stepwise betaregression 
input_data <- dplyr::select(bare_ground_model, - c('Plot_ID', 'Prop.Bare')) %>% na.omit()
test <- leaps::regsubsets(y ~ . ,
                          data =input_data,    # 1 best model for each number of presdictors
                          nvmax = 11)
summary(test)$rsq
m1 <- input_data[c('y',(names(input_data)[summary(test)$which[8,]]))] %>% na.omit()
m1
mods_1 <- betareg(y ~  . , data =m1)
w1 <- 1/resid(mods_1)^2
mods_1w <- update(mods_1, weights = w1)
#plot(mods_1)
summary(mods_1)
plot(mods_1w)
summary(mods_1w)
mods_2 <- betareg(y ~ pzabovemean+ zcv + pzabove2 + Lcoefvar + ttops, data = m1)
summary(mods_2)
w1 <- 1/resid(mods_2)^2
mods_2w <- update(mods_2, weights = w1)
#plot(mods_2w)
summary(mods_2w)

### Fit in Kfolds to get model accuracy
folds <-createMultiFolds(y=input_data$y,k=5,times=100)
bare_r_square_val <- as.numeric()
bare_new_mod_params <- list()
for(i in 1:length(folds)){
  train <- m1[folds[[i]],]
  test <- m1[-folds[[i]],]

  # Note that print(b) fails since b doesn't exist

  tryCatch(
    {kfold_mod1 <- betareg(y ~., data = train)
        w1 <- 1/resid(kfold_mod1)^2
      kfold_mod1w <- update(kfold_mod1, weights=w1, control = betareg.control(start = mods_1w$coefficients,
                                                                              fstol =1e-3))
            bare_r_square_val <- append(bare_r_square_val, kfold_mod1w$pseudo.r.squared)
            ##ad to df
            bare_new_mod_params[[i]] <- kfold_mod1w$coefficients$mean %>%
              t() %>%
              data.frame(., phi = kfold_mod1w$coefficients$precision)}, error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }
}
mean(bare_r_square_val[1:100], na.rm =T)
sd(bare_r_square_val, na.rm = T)
params_df_new_bar <- bind_rows(bare_new_mod_params[1:100]) %>% summarise_all(., .funs = c(mean, sd), na.rm = T) %>%
  pivot_longer(cols = everything(.), names_to = c(".value",'type'),
               names_sep = "_.") %>% t() %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>%
  data.frame(.) %>%  rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = c(names(mods_1w$coefficients$mean), 'phi'),
         r_square = mean(bare_r_square_val, na.rm = T),
         df = mods_1w$df.null,
         model_family = paste0('betaregression log link'))
##

in_model <- mods_1w
##Update Coefficients and R Square
in_model$coefficients$mean <- params_df_new_bar$est[1:10]
in_model$coefficients$precision <- params_df_new_bar$est[7]
in_model$pseudo.r.squared <- mean(bare_r_square_val[1:100], na.rm =T)
fitted(in_model)
dim(input_data)
bare_ground_model
