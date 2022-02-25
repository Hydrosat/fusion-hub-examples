'''

author: Joe McGlinchy (Hydrosat, Inc.)
'''

import rasterio as rio
from rasterio.env import Session
#from fiona.crs import from_epsg # deprecated in favor of pyproj
from pyproj.crs import CRS
from rasterio.mask import mask
from shapely.geometry import mapping, box, Point, Polygon
import requests
import multiprocessing as mp
import rioxarray as rxr
import xarray as xr
import pandas as pd
import os
import urllib

class FH_Hydrosat(object):
    
    def __init__(self, items, geometry=None, crs=None):
        
        self.items = items
        self.item_href = [i.to_dict()['assets']['lst']['href'] for i in items]
        self.item_desc = [i.to_dict()['links'][0]['href'] for i in items]
        self.datetime = [i.to_dict()['properties']['datetime'] for i in items]
        self.geometry = geometry
        self.crs = crs
        
                         
    def _update_lists(self, idx_ls):
        ''' update class attributes by index'''
        
        self.item_href = [_item for i,_item in enumerate(self.item_href) if i in idx_ls]
        self.item_desc = [_item for i,_item in enumerate(self.item_desc) if i in idx_ls]
        self.datetime = [_item for i,_item in enumerate(self.datetime) if i in idx_ls]
        
        
    def download_single_asset(self, idx, local_folder=None):
        ''' this function downloads a single asset referenced by index'''
        
        dl_asset = self.item_href[idx]
        dl_url = self.item_desc[idx]
        
        if local_folder is not None:
            outfile = os.path.join(local_folder, os.path.basename(dl_url))
        else:
            outfile = os.path.basename(dl_url)
            
        # download file
        return urllib.request.urlretrieve(dl_asset, outfile)
        
        
    def download_multiple_assets(self):
        # TODO
        pass
    
    def create_animation(self):
        # TODO
        pass
    
    def remove_below_data_perc(self):
        # TODO
        pass
    
    def stack(self, chunks=2048):
        ''' this function stacks the data files and adds a time dimension '''
            
        ds = xr.concat([rxr.open_rasterio(f) for f in self.item_href], dim='time')
         
        #datetimes2 = pd.to_datetime(self.datetime, format='%Y-%m-%dT%H:%M:%S.%fZ', utc=True)
        datetimes2 = pd.to_datetime(self.datetime, infer_datetime_format=True, utc=True) # more general construction
        datetimes2 = [d.to_pydatetime() for d in datetimes2]

        ds2 = ds.assign_coords(time=('time', datetimes2))
        
        # assign crs
        with rio.open(self.item_href[0]) as src:
            raster_crs = src.profile['crs']
        
        ds2.rio.write_crs(raster_crs)
        
        return ds2.chunk({'x':chunks, 'y':chunks})
        
        