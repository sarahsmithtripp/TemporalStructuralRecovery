"""
This script reads in a  collection of files, masks  the values that are zero, and the areas  that are not disturbed, 
It then stacks them by year
"""
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from shapely.geometry import box
import os
import glob
import numpy as np
import numpy.ma  as ma
from rasterio.warp import calculate_default_transform, reproject, Resampling

"""
Example of how to use this script 
## run on NBR
# Directory where your raster files are stored
directory = 'F:/NewBap/NBR_new/'

# Use glob to get all raster files in the directory
raster_files = glob.glob(os.path.join(directory, '*.tif'))
dist = rasterio.open(f'//FRST-FRM-2232B/bc/10S/change_metrics/SRef_10S_Greatest_Change_Year.dat', 'r')

# Open the first file to get the metadata
with rasterio.open(raster_files[0]) as src0:
    meta = src0.meta

bounds  = src0.bounds
bbox = gpd.GeoDataFrame({"id":1, "geometry": [box(*bounds)]})
coords = getFeatures(bbox)
out_img, out_transform = mask(dist, shapes=coords, crop=True)


# Update metadata to reflect the number of layers
meta.update(count = len(raster_files))

# Write it to a new file
with rasterio.open('F:/NewBap/Stacks/NBR_stack.tif', 'w', **meta) as dst:
    for id, layer in enumerate(raster_files, start=1):
        print(id)
        with rasterio.open(layer) as src1:
            r_read = src1.read()
            mask1 = (r_read[0] == 0) | (out_img[0] == 0)
            r_mask = ma.masked_array(r_read, mask1)
            r_masked_filled = r_mask.filled(-99999)
            dst.write_band(id, r_masked_filled[0])
"""

"""
functions to read in 
"""
def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def mask_stack(file_name):# Directory where your raster files are stored
    directory = f'F:/NewBap/{file_name}_new/'

    # Use glob to get all raster files in the directory
    raster_files = glob.glob(os.path.join(directory, '*.tif'))
    ## disturbance data also available at https://opendata.nfis.org/mapserver/nfis-change_eng.html
    dist = rasterio.open(f'//FRST-FRM-2232B/bc/10S/change_metrics/SRef_10S_Greatest_Change_Year.dat', 'r')

    # Open the first file to get the metadata
    with rasterio.open(raster_files[0]) as src0:
        meta = src0.meta

    bounds  = src0.bounds
    bbox = gpd.GeoDataFrame({"id":1, "geometry": [box(*bounds)]})
    coords = getFeatures(bbox)
    out_img, out_transform = mask(dist, shapes=coords, crop=True)


    # Update metadata to reflect the number of layers
    meta.update(count = len(raster_files))

    # Write it to a new file
    with rasterio.open(f'F:/NewBap/Stacks/{file_name}_stack.tif', 'w', **meta) as dst:
        for id, layer in enumerate(raster_files, start=1):
            print(id)
            with rasterio.open(layer) as src1:
                r_read = src1.read()
                # if out_img.shape != r_read.shape:
                #     base_array = np.full((out_img.shape[1], out_img.shape[2]), np.nan) 
                #     with rasterio.open(f'F:/NewBap/{file_name}_new/adj/mosaic_{id}.tif', 'w', **meta) as dst:
                #         reproject( 
                #         source=r_read,
                #         destination=base_array,
                #         src_transform=src1.transform,
                #         src_crs=src1.crs,
                #         dst_transform=out_transform,
                #         dst_crs=src0.crs,
                #         resampling=Resampling.nearest)
                #         dst.write(base_array, 1)
                #     r_read = rasterio.open(f'F:/NewBap/{file_name}_new/adj/mosaic_{id}.tif').read()
                mask1 = (r_read[0] == 0) | (out_img[0] == 0)
                r_mask = ma.masked_array(r_read, mask1)
                r_masked_filled = r_mask.filled(-99999)
                dst.write_band(id, r_masked_filled[0])


