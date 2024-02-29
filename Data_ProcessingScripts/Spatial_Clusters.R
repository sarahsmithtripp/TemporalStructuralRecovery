library(terra)

###
#This creates a nice table to input into the paper that describes the spatial aggregation of the structural groups on the landscape 
###
groups <- rast('D:/Paper2_Clean/Spectral_Clusters/Structural_groups.tif')
levels(groups)
metrics_calc <- c("lsm_c_contig_cv",
                  "lsm_c_contig_mn",
                  "lsm_c_iji",
                  "lsm_c_cpland",
                  "lsm_c_ai",
                  "lsm_c_pland",
                  "lsm_c_cohesion",
                  "lsm_l_mutinf",
                  "lsm_l_relmutinf")
library(landscapemetrics)
class_metrics <- calculate_lsm(groups, 
                               what = metrics_calc)

write_csv(class_metrics, 
          file = paste0(getwd(),'/images/class_structure.csv'))

## read in details csv
class_descript <-read_xlsx("images/September_Clusters/landscape_metrics.xlsx")

joined<- left_join(class_metrics, class_descript) %>% 
  mutate(value= round(value, 2), 
         mimimum = round(minimum), 
         maximum = round(as.numeric(str_replace(maximum,
                                     pattern = "\u221E",
                                     replacement = 'Inf')))) %>%
  filter(!is.na(label)) %>% 
  dplyr::select(-c('level', 'layer', 'id', 'names'))
joined_pivot <- pivot_wider(joined, id_cols = c('metric',
                                        'Description', 'type',
                                        'minimum','maximum', 'unit'),
                            names_from = class, values_from =value) %>% 
  group_by(type)

## test flextable here 
library(magrittr)
t1 <- flextable(as_grouped_data(joined_pivot, groups = "type")) %>%
  set_header_labels(., values = list(
  label = "Metric",
  Description = "Description*",
  type = "Level",
  '1' = '1',
  '2' = '2',
  '3' = '6',
  '4' = '7',
  '5' = 'C',
  'NA' = "Landscape",
  minimum = "min",
  maximum = "max",
  unit = "unit")) %>%
  add_footer_lines(.,
                   values =  c('*Hesselbarth, M. H. K., Sciaini, M.,
                               With, K. A., Wiegand, K., & Nowosad, J. (2019).
                               landscapemetrics: An open-source R tool to calculate landscape metrics.
                               Ecography, 42(10), 1648–1657. https://doi.org/10.1111/ecog.04617')) %>%
  autofit() %>%
  #colformat_double(digits = 2) %>%
  padding(padding = 0.5, part = 'all') %>%
  add_header_row(values = c("", "Range"," ", "Structural Group", ''),
                 colwidths = c(3,2, 1, 5, 1), top = TRUE)  %>%
  flextable::set_caption(caption = "Spatial Metrics for
                         Structural Trends Across the landscape") #%>%
t1 
write_csv(joined_pivot, file = "F:/Sync/PhD_Writing/Paper2/test/images/class_structure.csv")

