library(tidyverse)
library(kableExtra)
library(flextable)
### run the other scripting models 
### Read in model info 
list_rds <- list.files(path= 'D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation', pattern  = 'kfol', full.names = T) %>% 
  map(., .f = readRDS)
bare_ground <- readRDS('D:/Paper2_Clean/RPA_data/ABA_Models/ABA_validation/bare_ground_model.rds')
### plot the graphs 
graphs <- map(list_rds, .f = function(x) x[[3]])


a <- graphs[[1]] # composition
b <- list_rds[[2]][[5]] # Basal area 
c <- list_rds[[3]][[2]] + xlim(0, 1500) ## stems 
d <- bare_ground[[3]] + egg::theme_article()
library(cowplot)
metrics_3 <- plot_grid(b,a,c, nrow= 1)

#reduce margins
margins_left <- c(1,-1,-1,1)
margins_right <- c(1,1,-1,1)

a <- a + theme(plot.margin=unit(margins_left, "cm"), text = element_text(size = 14))
b <- b + theme(plot.margin=unit(margins_left, "cm"), text = element_text(size = 14))
c <- c + theme(plot.margin=unit(margins_right, "cm"), text = element_text(size = 14))

d <- d + theme(plot.margin=unit(margins_right, "cm"), text = element_text(size = 14)) + 
  ylab('Predicted % Bare Ground') + xlab('Measured % Bare Ground')

#plot
cp <- plot_grid(a, b, c, d, labels = LETTERS[1:4], label_size = 12, align="hv")
#save_plot('C:/Users/ssmithtr.stu/Desktop/test.png', cp)
##Update may 17th, add bare ground models 
metrics_4 <- plot_grid(d, b, c, a, nrow = 2, labels = c('A', 'B', 'C', 'D'))

### make the model kable 
#small 
dfs <- bind_rows(list_rds[[1]][[2]], list_rds[[2]][[3]], list_rds[[3]][[3]])
#big
### add names 
bare_ground[[2]]$name <- '% Bare Ground'
list_rds[[2]][[3]]$name <- 'Basal Area'
list_rds[[3]][[3]]$name <- 'Stem Count'
list_rds[[1]][[2]]$name <- 'Proportion Coniferous:Decidous'
dfs <- bind_rows(list_rds[[2]][[3]], list_rds[[3]][[3]],
                 list_rds[[1]][[2]], bare_ground[[2]])
#dfs[[3]] <- list_rds[[3]][[3]]
dfs_bind <- bind_rows(dfs) %>% relocate(model_family, .before = est) %>% 
  relocate(name, .before = model_family) %>% 
  relocate(df, .before = est) %>% relocate(param_names, .before = est) %>% 
  group_by(model_family, df)
dfs_bind[,5:7] <- apply(dfs_bind[,5:7], 2, FUN = round, digits = 2)
#dfs[[4]] <- list_rds[[2]][[1]]

set_flextable_defaults(font.family = 'sans',
  font.size = 11, font.color = "#303030",
  padding = 3, table.layout = "autofit")
names(dfs_bind) <- c('Model', 'Family', 'df', 'parameters', 
                     'est', 'stdev', 'R sq')

cluster_tests <- dfs_bind %>% flextable() %>%  merge_v(j = c(1:3, 7))
  # kable(.) %>% collapse_rows(1:3, row_group_label_position = 'stack') %>% 
  # kable_classic(full_width = F, html_font = 'calibri')