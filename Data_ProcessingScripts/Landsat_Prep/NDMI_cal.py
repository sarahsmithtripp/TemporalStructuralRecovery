"""
This script reads NBR and SWIR1 as a raster from the cleaned metrics. I use cleaned because 
I don't care about the null values in either 
"""
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from shapely.geometry import box
import os
import glob
import numpy as np
import numpy.ma  as ma

