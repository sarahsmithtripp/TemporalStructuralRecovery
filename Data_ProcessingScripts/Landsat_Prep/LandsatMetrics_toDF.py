"""
Script takes the metrics calculated and masked in prior files and then turns them into a dataframe 
It drops metrics that are not used in clustering analyses 
It alos filters out areas with bad data 
Finally it combines all of this data into one dataframe to send to R 

"""
import rasterio as rio
import pandas as pd
import numpy as np
import glob
from osgeo import gdal
## Raster data 
masked_dir = 'F:/NewBap/Stacks/cropped/SBS_metrics/'
masked_r = glob.glob(masked_dir + "/masked/*.tif")
masked_rl = [] # make empty list to fill
names_indices = list(["dist_mag", "regrowth", "vals_years", "slope"])
bands = ['NBR', 'NDMI', 'NDVI', 'TCA', 'TCB', 'TCG', 'TCW']
band_names = [f"{band}_{name}" for band in bands for name in names_indices]

#Open the multiband raster
in_file = 'F:/NewBap/Stacks/cropped/SBS_metrics/masked/Full_Rast/indices.tif'
##These files are combined and stacked in R and then read in here to turn into a DF 
## See Masked_PythonMetrics_ToFires for sources
#in_file = 'F:/NewBap/Stacks/cropped/SBS_metrics/masked/TCA.tif'
with rio.open(in_file) as src: 
        arr = src.read()
        shp = src.shape
        transform = src.transform
        #get the numer of rows and columns in the raster 
        nrows, ncols = arr.shape[1:]
        # use the transform object to get the arras of x any y coordinates 
        x = np.arange(ncols) * transform[0] + transform[2]
        y = np.arange(nrows) * transform[4] + transform[5]
        # Get the x,y coordinates of each pixel
        x, y = np.meshgrid(x, y)
        # Stack the x,y coordinates and the pixel values from each band into a single array
        stacked = np.column_stack((x.flatten(), y.flatten(), arr.reshape((arr.shape[0], -1)).T))
        df = pd.DataFrame(stacked)
        df.columns = ['x', 'y'] + band_names
        print(stacked.shape)
        df_test = df.dropna(subset = ['NBR_dist_mag'])
        df_reset =df_test.reset_index()
        df_reset.to_feather(f'F:/NewBap/Stacks/cropped/SBS_metrics/masked/Full_Rast/all_indexes2.feather')

## comment out the rest of this code because it has been replaced by the above
## This code is left in case I need to go back and look at it
## This code is also used to create the dataframes for the individual bands
# ## This code is also used to create the dataframes for the individual bands
#     ### older code that I dropped 
#     #for in_df, band_name in zip(masked_rl, bands): 
#     for masked_file, band_name in zip(masked_r[4:7], bands[4:7]): 
#         with rio.open(masked_file) as src:
#             arr = src.read()
#             shp = src.shape
#             transform = src.transform
#             #get the numer of rows and columns in the raster 
#             nrows, ncols = arr.shape[1:]
#             # use the transform object to get the arras of x any y coordinates 
#             x = np.arange(ncols) * transform[0] + transform[2]
#             y = np.arange(nrows) * transform[4] + transform[5]
#             # Get the x,y coordinates of each pixel
#             x, y = np.meshgrid(x, y)
#             # Stack the x,y coordinates and the pixel values from each band into a single array
#             stacked = np.column_stack((x.flatten(), y.flatten(), arr.reshape((arr.shape[0], -1)).T))
#             print(stacked.shape)
#             # Create a pandas dataframe from the stacked array with columns for the x and y coordinates
#             # Create a list of column names using names_indices and slice it to match the number of bands
#             band_names = names_indices
#             print(band_names)
#             df = pd.DataFrame(stacked)

#             df.columns = ['x', 'y'] + band_names
#             #masked_rl.append(df)
#             df.to_feather(f'F:/NewBap/Stacks/cropped/SBS_metrics/masked_df/{band_name}.feather')


#     ### pivot these data frames and select only columns of interest 
#     import pandas as pd 
#     #read in one dataframe to know what we are looking at 
#     bands = ['NDMI', 'NDVI', 'TCA', 'TCB', 'TCG', 'TCW']
#     store = []
#     names_indices = list(["pre_disturbance", "post_disturbance", "dist_mag", 
#                                     "recov_80", "regrowth", "vals_years", "slope"])
#     for band_name in (bands): 
#         df = pd.read_feather(f'F:/NewBap/Stacks/cropped/SBS_metrics/masked_df/{band_name}.feather')
#         df_test = df.dropna(subset = ['dist_mag'])
#         # Replace 'suffix' with the desired suffix you want to add
#         df_filtered = df_test.drop(columns = ["pre_disturbance", "post_disturbance", "recov_80", ])
#         columns_to_rename = df_filtered.columns[2:6]  # Selects the columns from index 2 to 7 (columns are zero-indexed)
#         new_column_names = [col + '_' + band_name for col in columns_to_rename]  # Creates new column names with the suffix
#         df_filtered.rename(columns=dict(zip(columns_to_rename, new_column_names)), inplace=True)
#         df_reset =df_test.reset_index()
#         df_reset.to_feather(f'F:/NewBap/Stacks/cropped/SBS_metrics/masked_df/filtered/{band_name}.feather')

#     ##Read in dataframes and bind to send to R 
#     import modin.pandas as pd
#     import glob 
#     import os
#     in_df = 'F:/NewBap/Stacks/cropped/SBS_metrics/masked_df/filtered'
#     dataframes= []
#     bands = ['NBR','NDMI', 'NDVI', 'TCA', 'TCB', 'TCG', 'TCW']
#     for filename, band_name in zip(os.listdir(in_df), bands):
#         print(filename)
#         if filename.endswith('.feather'):  # Adjust the file extension if needed
#             file_path = os.path.join(in_df, filename)
#             df = pd.read_feather(file_path)
#             df_filtered = df.drop(columns = ["pre_disturbance", "post_disturbance", "recov_80", ])
#             columns_to_rename = df_filtered.columns[2:6]  # Selects the columns from index 2 to 7 (columns are zero-indexed)
#             new_column_names = [col + '_' + band_name for col in columns_to_rename]  # Creates new column names with the suffix
#             df_filtered.rename(columns=dict(zip(columns_to_rename, new_column_names)), inplace=True)
#             dataframes.append(df_filtered)
#     combined_df = pd.concat(dataframes, ignore_index=True)
