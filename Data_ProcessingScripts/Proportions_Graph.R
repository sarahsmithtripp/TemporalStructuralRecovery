n_ <- ggplot(cluster_counts, 
                    aes(x = 1, y = proportion * 100, fill = as.factor(z))) + 
  geom_bar(stat = 'identity') + scale_fill_manual(values = in_cols) + 
  xlab('All\nYears') + theme_bw() +
  coord_flip(ylim = c(0,100)) +
  theme(legend.position = 'none', 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y =  element_text(size = 12, face = 'bold', angle = 0,
                                     vjust = 0.4),
        axis.ticks.x = element_line(size = 2),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1)) 
 n_
selected_data <- cluster_counts_attributed %>%
  #subset(Dist_Year %in% c(1992, 2006, 2010, 2014, 2017, 2018)) %>% 
  subset(Dist_Year %in% c(1988, 1992, 2006, 2007, 2010, 2014, 2017)) %>% 
  group_by(Dist_Year) %>% 
  mutate(x_val = cur_group_id()) %>% ungroup()
selected_data <- cluster_counts_attributed %>%
  #subset(Dist_Year %in% c(1992, 2006, 2010, 2014, 2017, 2018)) %>% 
  subset(Dist_Year %in% c(1988, 1992, 2006, 2007, 2010, 2014, 2017)) %>% 
  group_by(Dist_Year) %>% 
  mutate(x_val = cur_group_id()) %>% ungroup()
some_years <- ggplot(selected_data,
                     aes(x =x_val, y = prop_year, fill = as.factor(clust))) + 
  scale_x_continuous(breaks = sort(unique(selected_data$x_val)), 
                     labels = sort(unique(selected_data$Dist_Year))) + 
  geom_col(position = 'fill') + scale_fill_manual(values = in_cols) + 
  labs(fill = 'Cluster')+ 
  theme_bw() + coord_flip(ylim = c(0,1)) + 
  theme(axis.title.y = element_text(size = 8), 
        axis.ticks.x = element_line(size = 2), 
        #axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, face = 'bold'),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) +
  scale_y_continuous(breaks = c( 0, 0.25, 0.5, 0.75, 1
                              ), 
                     labels = seq(0,100, by = 25)) + xlab(' ') + 
  ylab('% of Landscape in Cluster')

some_years
plotb <- some_years + theme(legend.position = 'none')
legend = cowplot::get_legend(some_years + theme(legend.position = 'right') + guides(colour = guide_legend(ncol = 2)))
proportions <- plot_grid(n_, plotb, rel_heights = c(0.4, 0.9), labels = c('A', 'B'), vjust = 1,
                         label_size = 16, ncol =1 )
out <- plot_grid(proportions, legend, rel_widths = c(1, 0.1))                           
out


