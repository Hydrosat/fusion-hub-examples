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
from functools import partial
import numpy as np

from matplotlib import pyplot as plt, animation

# This is needed to display graphics calculated outside of jupyter notebook
from IPython.display import HTML, display

# needed on Linux 
from multiprocessing import set_start_method
set_start_method("spawn", force=True)


class FH_StackedDataset(object):
    
    def __init__(self, data_array):
        # super().__init__() # until correctly subclassing the xarray.core.DataArray class
        
        self.ds = data_array
        
    def remove_below_data_perc(self, ds, thresh=0.1):
        ''' This function removes items from the xarray DataArray which have data coverage below a specified percentage threshold.'''
        
        # check params
        if type(ds) is not xr.DataArray:
            raise TypeError('parameter ds must be of type xarray.DataArray')
            
        if len(ds.shape) != 4:
            raise TypeError('parameter ds must be of shape (time, band, x, y)')
            
        if (thresh > 1.0) or (thresh < 0.0000001):
            raise ValueError('parameter thresh must be in range (0, 1.0]')
            
        #########################################################################
        ########### remove the entries meeting the threshold criteria ###########
        data_count = ds.count(dim=('x', 'y', 'band')).values # number of valid pixels
        max_data = ds.shape[2] * ds.shape[3]               # max number of valid pixels  (all)
        data_perc = data_count / max_data                    # valid pixels as percentage
        valid_idx = np.where(data_perc > thresh)[0]          # apply the threshold
        valid_ds = ds.isel(time=valid_idx)                   # select valid entries
        
        return valid_ds
    
    def create_animation(self, figsize=(12,6), cmap='inferno', interval=200, vmin=280, vmax=320, save_ani=False, anipath='animation.gif'):
        # TODO
        # Get a handle on the figure and the axes
        fig, ax = plt.subplots(figsize=figsize)

        # Plot the initial frame. 
        cax = self.ds[0,0,:,:].plot(
            add_colorbar=True,
            cmap=cmap,
            vmin=vmin, vmax=vmax,
            cbar_kwargs={
                'extend':'neither'
            }
        )

        # Next we need to create a function that updates the values for the colormesh, as well as the title.
        def animate(frame):
            cax.set_array(self.ds[frame,0,:,:].values.flatten())
            ax.set_title("Time = " + str(self.ds.coords['time'].values[frame])) # original had [:13] at the end, to shorten title

        # Finally, we use the animation module to create the animation.
        ani = animation.FuncAnimation(
            fig,             # figure
            animate,         # name of the function above
            frames=self.ds.shape[0],       # Could also be iterable or list
            interval=interval     # ms between frames
        )
        
            
        if save_ani:
            ani.save(anipath)
        
        # close the figure
        plt.close()
        
        return ani
        
    
    
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
            
        # download file if it doesn't exist
        if not os.path.exists(outfile):
            return urllib.request.urlretrieve(dl_asset, outfile)
        else:
            return
        
        
    def download_multiple_assets(self, idx_list, nproc=2, local_folder=None):
        # TODO
        ''' this function downloads multiple assets by index.'''
        if nproc > 6:
            raise ValueError('specify nproc <=6')
            
        if type(idx_list) is not list:
            raise TypeError('idx_list should be list of indexes. If only a single value, use a list or call download_single_asset().')
            
        if len(idx_list) == 1:
            download_single_asset(idx_list[0], local_folder)
            
        if (local_folder is not None) and (not os.path.exists(local_folder)):
            os.makedirs(local_folder)
            
        else:
            # TODO 
            # call multiprocess with functools.partial
            dl_func = partial(self.download_single_asset,local_folder=local_folder)
            
            with mp.get_context("spawn").Pool(nproc) as pool:
                print(f'using {nproc} processes to download {len(idx_list)} assets')
                pool.map(dl_func, idx_list)
            
        return    
    
    
    def remove_below_data_perc(self, ds, thresh=0.1):
        # TODO
        ''' This function removes items from the xarray DataArray which have data coverage below a specified percentage threshold.'''
        
        # check params
        # is ds a data array?
        if type(ds) is not xr.DataArray:
            raise TypeError('parameter ds must be of type xarray.DataArray')
            
        if len(ds.shape) != 4:
            raise TypeError('parameter ds must be of shape (time, band, x, y)')
            
        if (thresh > 1.0) or (thresh < 0.0000001):
            raise ValueError('parameter thresh must be in range (0, 1.0]')
            
        #########################################################################
        ########### remove the entries meeting the threshold criteria ###########
        data_count = ds.count(dim=('x', 'y', 'band')).values # number of valid pixels
        max_data = ds.shape[2] * data.shape[3]               # max number of valid pixels  (all)
        data_perc = data_count / max_data                    # valid pixels as percentage
        valid_idx = np.where(data_perc > thresh)[0]          # apply the threshold
        valid_ds = ds.isel(time=valid_idx)                   # select valid entries
        
        return valid_ds
        
    
    def stack(self, chunks=2048, cache=False):
        ''' this function stacks the data files and adds a time dimension '''
            
        ds = xr.concat([rxr.open_rasterio(f, cache=cache, chunks=chunks) for f in self.item_href], dim='time')
         
        #datetimes2 = pd.to_datetime(self.datetime, format='%Y-%m-%dT%H:%M:%S.%fZ', utc=True)
        datetimes2 = pd.to_datetime(self.datetime, infer_datetime_format=True, utc=True) #more general conversion
        datetimes2 = [d.to_pydatetime() for d in datetimes2]

        ds2 = ds.assign_coords(time=('time', datetimes2))
        
        # assign crs
        with rio.open(self.item_href[0]) as src:
            raster_crs = src.profile['crs']
        
        ds2.rio.write_crs(raster_crs)
        
        return FH_StackedDataset(ds2.chunk({'x':chunks, 'y':chunks})) # return the class above, which has some added functionality
        
        
