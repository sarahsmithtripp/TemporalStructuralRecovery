"""
This script executes the following tasks:
it executes a script stored in the IRSS github organization that then calculates a suite of post-disturbance metrics 
It takes masked rasters with the greatest year of change and then runs the metrics on those
In this script the last year of raster band is the year of disturbance
"""
#load librarys 
import rasterio
import numpy as np 
import os 
import timeit 
from scipy.stats import mstats
# get functions 
exec(open('D:/Paper2_Clean/Data_ProcessingScripts/Landsat_Prep/Metrics_Calculations/prepost_metrics.py').read())


def metrics_calc(file_name): 
    print('starting file' + file_name)
    startt = timeit.default_timer()
    ### all satellite data can be downloaded at: https://github.com/saveriofrancini/bap
    path_in ='F:/BapSummer/_2022_stacked_metrics/west/masked_GYC/Masked_WGYC/' + file_name + '.tif'
    rast = rasterio.open(path_in) 
    rast_r = rast.read()
    obs = rast_r[:-1, :, :]
    dist = rast_r[-1, :, :]
    np_pre_disturbance = np_metrics(dist, obs, start = 1984, max_range = 5, slope = False)
   

    #save metrics 
    # write tif
    save = rast.profile
    save.update(count = np_pre_disturbance.shape[0])
    path_out = "F:/BapSummer/_2022_stacked_metrics/west/_metrics/" + file_name + ".tif"
    with rasterio.open(path_out, 'w', **save) as dst:
        for band_idx in range(np_pre_disturbance.shape[0]):
            dst.write(np_pre_disturbance[band_idx], indexes = band_idx + 1)
    stopt = timeit.default_timer()
    print(f"np runtime: {stopt - startt}") 
    print(file_name + " done!")
    return path_out

##update June 6, move onto new bap measures 
def metrics_calc_slope(file_name): 
    print('starting file' + file_name)
    startt = timeit.default_timer()
    path_in ='F:/NewBap/Stacks/cropped/' + file_name + 'cleaned.tif'
    print(path_in)
    rast = rasterio.open(path_in) 
    rast_r = rast.read()
    obs = rast_r[:-1, :, :]
    dist = rast_r[-1, :, :]
    np_pre_disturbance = np_metrics(dist, obs, start = 1984, max_range = 5, slope = True)

    #save metrics 
    # write tif
    save = rast.profile
    save.update(count = np_pre_disturbance.shape[0])
    path_out = "F:/NewBap/Stacks/cropped/SBS_metrics/" + file_name + "wslope.tif"
    with rasterio.open(path_out, 'w', **save) as dst:
        for band_idx in range(np_pre_disturbance.shape[0]):
            dst.write(np_pre_disturbance[band_idx], indexes = band_idx + 1)
    stopt = timeit.default_timer()
    print(f"np runtime: {stopt - startt}") 
    print(file_name + " done!")
    return path_out

bands = ("NBR", "NDMI", "NDVI", "TCA", "TCB", "TCG", "TCW")

metrics_calc_slope('TCG')
metrics_calc_slope('TCW')
metrics_calc_slope('NDMI')
metrics_calc_slope('NBR')
metrics_calc_slope('TCA')

merged_df = pd.concat([df.set_index(['x', 'y']) for df in dataframes], axis=1).reset_index()
merged_df.head()


in_r =rio.open('F:/BapSummer/_2022_stacked_metrics/west/nbr_west.tif').read(1)

unique = np.arange(1, in_r.shape[1] * in_r.shape[0] + 1)
unique_cell = unique.reshape(in_r.shape[0], in_r.shape[1])
unique_cell.shape


save = rio.open('F:/BapSummer/_2022_stacked_metrics/west/nbr_west.tif').profile
save.update(count =1)
path_out = "F:/BapSummer/_2022_stacked_metrics/west/index_raster.tif"
with rio.open(path_out, 'w', **save) as dst:
    dst.write(unique_cell, indexes = 1)


merged_df['ID'] = unique
## ID represents the values of the raster we want to evantually write back to 
merged_filterd = merged_df.dropna(subset=['dist_mag_NBR'])


merged_filterd.reset_index().to_feather('F:/BapSummer/_2022_stacked_metrics/west/_metrics/df/merged_metrics/metrics_may21.feather')