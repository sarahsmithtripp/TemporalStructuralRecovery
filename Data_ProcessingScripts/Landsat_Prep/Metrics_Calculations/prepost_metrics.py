"""
This module contains various approaches for computing pre/post disturbance values.

"""
import rasterio
import rioxarray
import timeit

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import numpy.ma as ma 
from scipy.stats import mstats


def xarray_metrics(disturbance_array, observations, start, end, pre_range):
    """
    xarray implementation of pre/post max values

    NOTE: this is the layman's algo. It loops over years and gets the 
    max, then combines the values at the end. There is likely a ~faster~ way
    to do this. Possibly through `apply_ufunc` or split-apply-combine workflow.

    
    Parameters
    -----------
    disturbance_array: xr.DataArray
        2D (x, y) array of disturbance years. A value X at some pixel 
        (i,j) represents a disturbance during year X at spatial loc (i,j).
        x/y dimensions must match x/y dimensions of `observations`.
    
    observations: xr.DataArray
        3D (z (time), x, y) array of spectral/index values. x/y 
        dimensions must match `disturbance_array`
    
    start: int
        The start year of the timeseries. pre-values will only
        be computed for years `start`+1 because otherwise, there
        is no previous value to look to.
    
    end: int
        The final year to compute metrics for.
    
    pre_range: int
        The # of years before the disturbance year to look back at 
        for pre- values. 
    """
    dist_max = xr.zeros_like(disturbance_array)

    for dist_year in range(start, end+1):
        disturbance_mask = xr.where(disturbance_array == float(dist_year), True, False)
        dist_year_index = dist_year - start
        if dist_year_index != 0:
            if dist_year_index >= pre_range:
                unmasked_window = observations[dist_year_index-pre_range:dist_year_index,:,:]
            else:
                reduce_pre_range = dist_year_index
                unmasked_window = observations[dist_year_index-reduce_pre_range:dist_year_index,:,:]

            max_previous = unmasked_window.max(dim="band", skipna=True) # calculate max here 
            max_previous_sliced = xr.where(disturbance_mask, max_previous, 0)
            dist_max = dist_max + max_previous_sliced

    return dist_max



def xarray_metrics_faster(disturbance_array, observations, start, end):
    """
    xarray implementation of pre/post max values but attempting to
    be faster than `xr_metrics`
    """
    # NOTE: tried to expanding into new xarray that holds indexes (similar to np)
    # and implementation was much slower. Doesn't rule it out but good to note.

    # TODO: try to get an `apply_ufunc` version going.
    return 0


def np_metrics(disturbance_array, observations, start, max_range, slope):
    """
    numpy implementation of pre/post-disturbance max values

    NOTE: fast, fast, fast. Taken directly from Lukas.

    Parameters
    -----------
    disturbance_array: np.array
        2D (x, y) array of disturbance years. A value X at some pixel 
        (i,j) represents a disturbance during year X at spatial loc (i,j).
        x/y dimensions must match x/y dimensions of `observations`.
    
    observations: np.array
        3D (z (time), x, y) array of spectral/index values. x/y 
        dimensions must match `disturbance_array`
    
    start: int
        The start year of the timeseries. pre-values will only
        be computed for years `start`+1 because otherwise, there
        is no previous value to look to.

    max_range: int
        The max # of years to calculate indices for.  Set for the 
        last year of the data if interested in the current year  
        
    slope: bool
        calculate thiel-sen slope on the dataframe for the range of years. Operation is very slow
        Boolean to avoid time if use does not want 
    """
    #Set up data 
    pre_range = 2
    post_range = 2
    disturbance_year_index = (disturbance_array-start).astype(int)


    #pre-disturbance
    ##always want two years before distubance b/c indexing is poor otherwise 
    invalid_years = np.logical_or(disturbance_year_index < -2, disturbance_year_index > observations.shape[0]-pre_range)  # make sure that we are not trying to access anything before the time series starts or after it ends
    disturbance_year_index[invalid_years] = 1 # for now just select the second year for these locations, will overwrite later
    sel_years = np.stack([disturbance_year_index + i for i in range(-pre_range, 0)])
    years_to_summarize = np.take_along_axis(observations, sel_years, axis=0)
    pre_disturbance = np.nanmax(years_to_summarize, axis=0)  # if all are nan in a certain pixel, the output will also be nan
    pre_disturbance[invalid_years] = np.nan
    
    #post disturbance 
    invalid_years = np.logical_or(disturbance_year_index < -2, disturbance_year_index > observations.shape[0]-post_range)  # make sure that we are not trying to access anything before the time series starts or after it ends
    disturbance_year_index[invalid_years] = 1 # for now just select the second year for these locations, will overwrite later
    sel_years = np.stack([disturbance_year_index + i for i in range(0, post_range)])
    years_to_summarize =  np.take_along_axis(observations, sel_years, axis=0)
    post_disturbance = np.nanmin(years_to_summarize, axis = 0)
    post_disturbance[invalid_years] = np.nan
    
    #disturbance magnitude
    dist_mag = pre_disturbance - post_disturbance

    #80% recovery 
    recov_80 = pre_disturbance * 0.8 

    #value for the maximum year of interest
    invalid_years = np.logical_or(disturbance_year_index < -2, disturbance_year_index > observations.shape[0]-max_range)  # make sure that we are not trying to access anything before the time series starts or after it ends
    disturbance_year_index[invalid_years] = 1
    max_years = np.stack([disturbance_year_index + i for i in range(max_range-1, max_range)])
    years_to_summarize = np.take_along_axis(observations, max_years, axis = 0)
    vals_years = np.nanmean(years_to_summarize, axis = 0)
    vals_years[invalid_years] = np.nan
    regrowth = vals_years-post_disturbance
    
    #Theil Sen Slope ### SLOW SLOW SLOW B/C for Loop Operation 
    if slope == True:
        #option to leave out slope b/c slow 
        sel_years = np.stack([disturbance_year_index + i for i in range(1, max_range)])
        years_to_summarize = np.take_along_axis(observations, sel_years, axis = 0)
        slope = np.round(theil_sen_slope(years_to_summarize))
        stacked_metrics =  np.stack([pre_disturbance, post_disturbance, dist_mag, 
                                recov_80, regrowth, vals_years, slope])
    else: 
        stacked_metrics =  np.stack([pre_disturbance, post_disturbance, dist_mag, 
                                recov_80, regrowth, vals_years])
    return stacked_metrics


def theil_sen_slope(arr):
    # Get the number of bands (i), rows (k), and columns (j)
    i, k, j = arr.shape

    # Reshape the array into a 2D array with k*j rows and i columns
    flat_arr = arr.transpose(1, 2, 0).reshape(k*j, i)

    # Initialize a list to store the slopes
    slopes = []

    # Compute the slope for each kj element
    for idx in range(k*j):
        # Get the x and y values for the current kj element
        x = np.arange(i)
        y = flat_arr[idx, :]

        # Compute the pairwise differences between all x and y values
        dx = np.subtract.outer(x, x)
        dy = np.subtract.outer(y, y)

        # Mask the differences where x and y have the same value or either x or y is NaN
        mask = (dx != 0) & ~np.isnan(dx) & ~np.isnan(dy)
        dx = dx[mask]
        dy = dy[mask]

        # If no valid pairs are left after masking, append a NaN slope
        if dx.size == 0 or dy.size == 0:
            slopes.append(np.nan)
            continue

        # Compute the slope as the median of the pairwise differences
        slope = np.nanmedian(dy / dx)
        slopes.append(slope)

    # Reshape the slopes back into a 3D array with shape (k, j)
    slope_arr = np.array(slopes).reshape(k, j)

    return slope_arr 

def theil_sen_faster(arr):
    """
    faster implementation of theil sen slopes that doesn't pass through a for loop 
    This is really struggling, might have to think of an alternative slope calculation 
    e.g. ransac https://vitalflux.com/ransac-regression-explained-with-python-examples/
    https://stackoverflow.com/questions/20343500/efficient-1d-linear-regression-for-each-element-of-3d-numpy-array
    """
    n_min = 3
    mask_arr = ma.masked_invalid(arr)
    valid_idx = mask_arr.count(axis=0) >= n_min

    new_r = mask_arr[:, valid_idx]
    arr_in = ~new_r.mask
    #get non NA 
    arr_clean = arr_in.data
    x = arr.data
    startt = timeit.default_timer()
    slopes = np.zeros(arr_in.shape[1])
    count = []
    # Theilslopes for each valid pixel
    for i in range(arr_in.shape[1]):
        if new_r[:, i].sum() >= n_min:
            slopes[i] = mstats.theilslopes(new_r[:, i])[0]
        else: 
            count.append(1)
    stopt = timeit.default_timer()
    print(f"np runtime: {stopt - startt}") # I get 52 seconds on the references raster  
    # # Create output grid with original dimensions
    # out = np.ma.masked_all_like(mask_arr.ma.stack[0]) ## unsure what this code is doing here 
    # # Fill in the valid indices
    # out[valid_idx] = slopes

    return 0
if __name__ == '__main__':
    
    start_year = 1984
    end_year = 2022
    path_to_test_tif = "test_real.tif"

    # Read in the TIF as Dask array and run through xarray metric func
    observations_and_map = rioxarray.open_rasterio(path_to_test_tif, chunks="auto")
    xr_disturbance_map = observations_and_map.sel(band=40)
    observations = observations_and_map.drop(band=40)
    startt = timeit.default_timer()
    xr_pre_disturbance = xarray_metrics(xr_disturbance_map, observations, start_year, end_year, pre_range)
    stopt = timeit.default_timer()
    print(f"xr runtime: {stopt - startt}")
    # write tif
    xr_pre_disturbance.rio.to_raster("xr_test.tif")

    # Read in TIF as np array and run through np metric func
    data = rasterio.open(path_to_test_tif)
    arr = data.read()
    dist_year = arr[-1, ...]
    obs = arr[:-1, ...]
    startt = timeit.default_timer()
    np_pre_disturbance = np_metrics(dist_year, obs, start_year, max_range = 5)
    stopt = timeit.default_timer()
    print(f"np runtime: {stopt - startt}")  

    # write tif
    save = data.profile
    save.update(count = np_pre_disturbance.shape[0])
    with rasterio.open("np_test.tif", 'w', **save) as dst:
        for band_idx in range(np_pre_disturbance.shape[0]):
            dst.write(np_pre_disturbance[band_idx], indexes = band_idx + 1)
    