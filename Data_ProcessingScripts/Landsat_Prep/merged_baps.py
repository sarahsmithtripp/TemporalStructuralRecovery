"""
This script searches the directory to merge and mosaic a collection of BAP rasters downloaded from GEE. It merges everything
by year. It rights this all to a single year for the processing and stores this to a new raster representing a single year. 
"""
import json
import rasterio
from rasterio.merge import merge
import shutil
with open("F:/NewBap/call_list.json") as f:
    raster_files = json.load(f)


import os
import glob

# Define the directory where your files are
directory = "F:/NewBap/ndmi_new/raw_NDMI"

# Use glob to recursively search the directory and return all files that contain 'TCB' and end with '.tif'
tcb_files = glob.glob(os.path.join(directory, '**/*TCB*.tif'), recursive=True)

years = range(1984, 2022)

# Loop over the years
def mosaic_rasters(in_files, years, file_name):
    for year in years:
        ## get all open dataset so you can close them at the end 
        opened_ = []
        for raster in in_files:
            src = rasterio.open(raster, 'r')
            opened_.append(src)
            lyrs = src.descriptions
            subdataset = f'{year}' 
            if subdataset in lyrs:
                src_subdataset = src.read(lyrs.index(subdataset)+1)
                # Assuming the band you want to merge is the first one, i.e., index 1
                    # Write the result to a new file
                out_ = src.meta.copy()
                out_.update({'height':src_subdataset.shape[0], 
                            'width':src_subdataset.shape[1], 
                            'count':1}) 
                name_save  = src.name.split("/", 5)[-1].split('.tif')[0]
                with rasterio.open(f"F:/NewBap/{file_name}_new/temp/{name_save}_mosaic{year}.tif", "w", **out_ ) as tmp : 
                    tmp.write(src_subdataset,1)
            else:
                print(f"Subdataset does not exist: {subdataset}")
        read_in =  glob.glob(os.path.join(f'F:/NewBap/{file_name}_new/temp', '*.tif'))
        src_files_to_mosaic = []
        for r in read_in:
            src = rasterio.open(r)
            opened_.append(src)
            src_files_to_mosaic.append(src)
        # Merge the rasters
        mosaic, out_trans = rasterio.merge.merge(src_files_to_mosaic)

        # Write the result to a new file
        out_meta = src.meta.copy()
        out_meta.update({"driver": "GTiff",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": out_trans
                        }
                        )
        with rasterio.open(f"F:/NewBap/{file_name}_new/mosaic_{year}.tif", "w", **out_meta) as dest:
            dest.write(mosaic)
            opened_.append(dest)
        for i in opened_:
            i.close()
        for r in read_in:
            os.remove(r)

mosaic_rasters(in_files = tcb_files, years = range(1988, 2022), file_name= 'TCB')
# Use glob to recursively search the directory and return all files that contain 'TCB' and end with '.tif'
nbr = glob.glob(os.path.join(directory, '**/*NBR*.tif'), recursive=True)
mosaic_rasters(in_files = nbr, years = years, file_name= 'NBR')
ndvi = glob.glob(os.path.join(directory, '**/*NDVI*.tif'), recursive=True)
mosaic_rasters(in_files = ndvi, years = range(2004, 2023), file_name= 'NDVI')

ndmi = glob.glob(os.path.join(directory, '**/*ndmi*.tif'), recursive=True)
mosaic_rasters(in_files = ndmi, years = years, file_name= 'ndmi')

tcg = glob.glob(os.path.join(directory, '**/*TCG*.tif'), recursive=True)
mosaic_rasters(in_files = tcg, years = years, file_name= 'TCG')

tcw = glob.glob(os.path.join(directory, '**/*TCW*.tif'), recursive=True)
mosaic_rasters(in_files = tcw, years = years, file_name= 'TCW')

tca = glob.glob(os.path.join(directory, '**/*TCA*.tif'), recursive=True)
mosaic_rasters(in_files = tca, years = years, file_name= 'TCA')

