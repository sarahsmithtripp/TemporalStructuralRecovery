# this script uses field data to calculate the bare ground proportion for ABA models 

# Bare ground model  -----------------------------------------------
#Big bare ground model
bare_ground_model <- left_join(field_plots, lidar_plots, by = c('Plot_ID' = 'Plot_ID')) %>%
select(-c('Prop.Con', 'Prop.Decidious', 'decid_con_rat', 'rounded_out', 'X', 'n', 'Plot_ID'))

mods_1 <- betareg(Prop.Bare ~ zsd + zcv + zskew + VCI, data = bare_ground_model)
mods_2 <- betareg(Prop.Bare ~ zsd + zskew + VCI + L2, data = bare_ground_model)
mods_3 <- betareg(Prop.Bare ~ pzabove0.2 + pzabove1.3 + CRR + zentropy + Lcoefvar, data = bare_ground_model)
#(mods_3)
mods_4 <- betareg(Prop.Bare ~ zmax + CRR + zentropy + Lcoefvar, data = bare_ground_model)
mods_5 <- betareg(Prop.Bare ~ pzabove1.3+ zentropy + VCI, data = bare_ground_model)
### best model is model 3

### Fit in Kfolds to get model accuracy
folds <-createMultiFolds(y=bare_ground_model$Prop.Bare,k=5,times=100)
bare_r_square_val <- as.numeric()
bare_new_mod_params <- list()
for(i in 1:length(folds)){
  train <- bare_ground_model[folds[[i]],] %>% select(-'X.1')
  test <- bare_ground_model[-folds[[i]],] %>% select(-'X.1')

  # Note that print(b) fails since b doesn't exist

  tryCatch(
    {kfold_mod1 <- betareg(Prop.Bare ~pzabove0.2 + pzabove1.3 + CRR + zentropy + Lcoefvar, data = train)
            bare_r_square_val <- append(bare_r_square_val, kfold_mod1$pseudo.r.squared)
            ##ad to df
            bare_new_mod_params[[i]] <- kfold_mod1$coefficients$mean %>%
              t() %>%
              data.frame(., phi = kfold_mod1$coefficients$precision)}, error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }
}
mean(bare_r_square_val[1:100], na.rm =T)
sd(bare_r_square_val, na.rm = T)
params_df_new_bar <- bind_rows(bare_new_mod_params[1:100]) %>% summarise_all(., .funs = c(mean, sd), na.rm = T) %>%
  pivot_longer(cols = everything(.), names_to = c(".value",'type'),
               names_sep = "_.") %>% t() %>% .[2:nrow(.), ] %>% apply(., 2, as.numeric) %>%
  data.frame(.) %>%  rename('est' = 'X1', 'stdev' = 'X2') %>%
  mutate(param_names = c(names(mods_3$coefficients$mean), 'phi'),
         r_square = mean(bare_r_square_val, na.rm = T),
         df = mods_3$df.null,
         model_family = paste0('betaregression log link'))
##

in_model <- mods_3
##Update Coefficients and R Square
in_model$coefficients$mean <- params_df_new_bar$est[1:6]
in_model$coefficients$precision <- params_df_new_bar$est[7]
in_model$pseudo.r.squared <- mean(bare_r_square_val[1:100], na.rm =T)
