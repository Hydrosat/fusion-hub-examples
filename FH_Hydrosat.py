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
import geopandas as gpd
import datetime
from odc import stac

from matplotlib import pyplot as plt, animation

# This is needed to display graphics calculated outside of jupyter notebook
from IPython.display import HTML, display

# needed on Linux 
from multiprocessing import set_start_method
set_start_method("spawn", force=True)

#classless function to unpack combined_qa asset from v19+
# write a value dictionary for the combinations of qa
qa_dict = {
    'Prepared High-Resolution t0': {'position': 0, 'value': 1},
    'Prepared High-Resolution t1': {'position': 1, 'value': 1},
    'Prepared Low-Resolution t0': {'position': 2, 'value': 1},
    'Prepared Low-Resolution t1': {'position': 3, 'value': 1},
    'Sharpened High-Resolution t0': {'position': 4, 'value': 1},
    'Sharpened High-Resolution t1': {'position': 5, 'value': 1},
    'Sharpened Low-Resolution t0': {'position': 6, 'value': 1},
    'Sharpened Low-Resolution t1': {'position': 7, 'value': 1},
    'Fused t1': {'position': 8, 'value': 1}
}

def plot_lst_and_qa(hdst_item, mask_val=None, keep_val=None):
    
    mask_href = hdst_item.to_dict()['assets']['combined_qa']['href']
    with rio.open(mask_href) as src:
        qa = src.read(1)

    lst_href = hdst_item.to_dict()['assets']['lst']['href']
    with rio.open(lst_href) as src:
        lst = src.read(1)
     
    ## apply masking or keeping of values
    # if a mask_val was provided, use it to set those pixels to np.nan so they won't be used 
    # for analysis or visualizaton. This can be done by numpy's `where` function, which operates as
    # np.where(condition, value_if_true, value_if_false)
    # if there are multiple values in a list to mask, we will mask each one sequentially.
    if mask_val is not None:
        
        # if a list of mask values was provided, sequentially set the lst data where the QA data 
        # matches those values to np.nan
        if type(mask_val) == list:
            for mv in mask_val:
                lst = np.where(qa == mask_val, np.nan, lst)
        elif type(mask_val) == float:
            
            # if only one value supplied to mask, us np.where() to set the lst data
            # to np.nan where that condition is true
            lst = np.where(qa == mask_val, np.nan, lst)
        else:
            raise ValueError("mask_val must be either float or list of floats")
    
    # if a keep_val was provided, use it to set pixels not equal to those values to np.nan so 
    # they won't be used for analysis or visualizaton. This can be done with boolean logic and indexing
    # directly into the array. if there are multiple values, we will sequentially OR the mask as more values
    # are used.
    if keep_val is not None:
        if type(keep_val) == list:
            
            # if it's a list, there should be multiple values. iterate through them.
            for i,kv in enumerate(keep_val):
                # build keep mask
                if i == 0:
                    keep_mask = qa == kv
                else:
                    keep_mask = keep_mask | (qa == kv)
            
            # take the complement of the keep_mask to retain only lst data where the QA values match
            # the provided list.
            lst[~keep_mask] = np.nan
                
        elif type(keep_val) == float:
            
            # if only one value provided to keep, set the lst data not equal to that value to np.nan
            lst[qa != keep_val] = np.nan
        else:
            raise ValueError("keep_val must be either float or list of floats")
        

    ## display the data side-by-side
    fig, ax = plt.subplots(1,2, figsize=(10,5))
    im = ax[0].imshow(lst, cmap='inferno')
    ax[0].set_title('HDST')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    #### make a unique value colormap for the QA data
    cmap = plt.cm.cividis  # define the colormap
    # extract all colors from the cmap
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    num_vals = int(np.unique(qa).shape[0])
    bounds = np.linspace(0, num_vals, num_vals)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


    im = ax[1].imshow(qa, cmap=cmap, norm=norm)
    ax[1].set_title('COMBINED_QA')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional', 
                              ticks=bounds+0.5, boundaries=bounds, format='%1i')

    cbar.ax.set_yticklabels(np.unique(qa))

    [a.axis('off') for a in ax]

    plt.show()
    
    return


def unpack_qa_value(qa_val=0):

    bit_desc_full = []
    # if it is valued 0, it is a clear pixel
    if qa_val == 0:
        bit_desc_full.append('valid pixel in all inputs')
        
    else:        
        # get the bit string for the current qa value
        bit_str = format(int(qa_val), '#010b')
        bit_vals = bit_str.split('b')[1]

        bit_desc = []
        for pos,b in enumerate(bit_vals[::-1]):

            if b=='1':
                for k,v in qa_dict.items():
                    if v['position'] == pos:
                        bit_desc.append(k)
                        
        bit_desc_full = bit_desc
    
    return bit_desc_full

# classless function to generate lst_mask values from a set of STAC items
def get_point_lst_mask_vals(item_hr, item_cr_t0, item_cr_tn, pt, tol=1500):
    ''' this function samples the lst_mask assets for a point location
    
    item_hr: STAC item for high res observation
    item_cr_t0: STAC item for coarse res observation at t0
    item_cr_tn: STAC item for coarse res observation at t1
    '''
    # check pt param type
    if type(pt) is not Point:
        raise(TypeError, "input pt must be of type shapely.geometry.Point")

    point_df = gpd.GeoDataFrame({'geometry':[pt]}, crs=CRS.from_epsg(4326))

    # project the point to raster CRS
    ds = rxr.open_rasterio(item_hr.to_dict()['assets']['lst']['href'])
    raster_crs = CRS.from_wkt(ds.spatial_ref.crs_wkt)
    point_df_utm = point_df.to_crs(raster_crs)
    set_x, set_y = point_df_utm['geometry'][0].x, point_df_utm['geometry'][0].y
    
    vals = []
    for _i, i in enumerate((item_hr, item_cr_t0, item_cr_tn)):
        ds = rxr.open_rasterio(i.to_dict()['assets']['lst_mask']['href'])
        
        # swap qa band values
        ds = xr.where(ds == 1, 0, 1)
        
        if _i == 1:
            ds = xr.where(ds == 1, 2, ds)
        if _i == 2:
            ds = xr.where(ds == 1, 4, ds)
        
        val = ds.isel(band=0).sel(x=set_x, y=set_y, method='nearest', tolerance=tol).values
        vals.append(val)
    
    return np.array(vals)
                               
    
def get_rel_items_for_pred(item_hr, prepped_hr_items, prepped_cr_items):
    
    pred_doy = item_hr.to_dict()['id'].split('_')[2]
    pred_str = datetime.datetime.strptime(pred_doy, '%Y%j').strftime('%Y-%m-%d')
    from_doy = item_hr.to_dict()['id'].split('_')[4]
    from_str = datetime.datetime.strptime(from_doy, '%Y%j').strftime('%Y%m%d')

    high_res_item = [i for i in prepped_hr_items if from_str in i.to_dict()['id']][0]
    coarse_res_item_t0 = [i for i in prepped_cr_items if from_doy in i.to_dict()['assets']['lst']['href']][0]
    coarse_res_item_tn = [i for i in prepped_cr_items if pred_doy in i.to_dict()['assets']['lst']['href']][0]


    return (high_res_item, coarse_res_item_t0, coarse_res_item_tn, pred_str)

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
    
    def __init__(self, items, geometry=None, crs=None, asset='lst'):
        
        self.items = items
        self.item_href = [i.to_dict()['assets'][asset]['href'] for i in items]
        self.item_desc = [os.path.basename(href.split('?')[0]) for href in self.item_href]
        self.datetime = [i.datetime for i in items]
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
        
        # illegal character set
        sets = [':','*','?','"','<','>','|']
        
        if local_folder is not None:
            outfile = os.path.join(local_folder, os.path.basename(dl_url)) #+ '.tif')
            for char in outfile:
                if char in sets:
                    outfile = outfile.replace(char,'_')
        else:
            outfile = os.path.basename(dl_url) #+ '.tif'
            for char in outfile:
                if char in sets:
                    outfile = outfile.replace(char,'_')
            
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
            self.download_single_asset(idx_list[0], local_folder)
            
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
        ''' this function stacks the data files and adds a time dimension 
        updated to use odc stac!'''
            
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
    
        
    def stack2(self, chunks=2048, cache=False):
        ''' this function stacks the data files and adds a time dimension. Updated to ensure matching x-y dims. '''
        
        base_ds = rxr.open_rasterio(self.item_href[0], cache=cache, chunks=chunks)
        restof_ds = [rxr.open_rasterio(f, cache=cache, chunks=chunks).rio.reproject_match(base_ds) for f in self.item_href[1:]]
        ds = xr.concat([base_ds] + restof_ds, dim='time')
         
        #datetimes2 = pd.to_datetime(self.datetime, format='%Y-%m-%dT%H:%M:%S.%fZ', utc=True)
        datetimes2 = pd.to_datetime(self.datetime, infer_datetime_format=True, utc=True) #more general conversion
        datetimes2 = [d.to_pydatetime() for d in datetimes2]

        ds2 = ds.assign_coords(time=('time', datetimes2))
        
        # assign crs
        with rio.open(self.item_href[0]) as src:
            raster_crs = src.profile['crs']
        
        ds2.rio.write_crs(raster_crs)
        
        return FH_StackedDataset(ds2.chunk({'x':chunks, 'y':chunks})) # return the class above, which has some added functionality
    
    
    def stack_odc(self, bands=['lst', 'combined_qa'], cache=False):
        ''' warning, this is currently much slower than either of stack or stack2'''
        
        return stac.load(self.items, bands=bands, cache=cache)
    
        
    def _extract_point_val(self, href, set_x=None, set_y=None, tol=20, band=0):
        """ construct a pandas DataFrame which is a time series for a single pixel across the search result.
            pt: shapely.geometry.Point
            tol: float. specified to retrieve nearest pixel
        """
        
        try:
            # open the raster dataset to sample
            ds = rxr.open_rasterio(href, chunks=2048, cache=False)
            val = ds.isel(band=band).sel(x=set_x, y=set_y, method='nearest', tolerance=tol).values
        except Exception as e:
            val = np.nan
        
        return val
            
    
    def _extract_point_val_window(self, href, set_x=None, set_y=None, tol=20, pad=1, band=0):
        """ construct a pandas DataFrame which is a time series for a single pixel across the search result. The value
            is computer from the area mean from +/- `pad` value around the center pixel.
            
            pt: shapely.geometry.Point
            tol: float. specified to retrieve nearest pixel
        """
        
        try:
            # open the raster dataset to sample
            ds = rxr.open_rasterio(href, chunks=2048, cache=False)
            xslice = slice(set_x-pad*tol, set_x+pad*tol)
            yslice = slice(set_y+pad*tol, set_y-pad*tol)
            val = np.nanmean(ds.isel(band=band).sel(x=xslice, y=yslice).mean(dim=('x', 'y')).values)

        except Exception as e:
            val = np.nan

        return val
    
    # TODO
    def _extract_area_val(self, href, poly_df, stat='mean', clip_dict={'all_touched':True, 'crop':True}):
        """ construct a pandas DataFrame which is a time series for pixels in the geometry across the search result."""
        try:
            # open the raster dataset to sample
            ds = rxr.open_rasterio(href, chunks=2048, cache=False)
            clipped_ds = ds.rio.clip(poly_df.geometry, *clip_dict)
            
            if stat=='mean':
                val = np.nanmean(clipped_ds.values)
            elif stat=='max':
                val = np.nanmax(clipped_ds.values)
            elif stat=='min':
                val = np.nanmin(clipped_ds.values)
            elif stat=='std':
                val = np.nanstd(clipped_ds.values)
            elif stat=='var':
                val = np.nanvar(clipped_ds.values)
            elif stat=='median':
                val = np.nanmedian(clipped_ds.values)
            else:
                val = np.nan

        except Exception as e:
            val = np.nan

        return val
        
    
    
    def point_time_series_from_items(self, pt, tol=20, nproc=2, pad=None, band=0):
        """ construct a pandas DataFrame which is a time series for a single pixel across the search result."""
        
        # check pt param type
        if type(pt) is not Point:
            raise(TypeError, "input pt must be of type shapely.geometry.Point")
            
        point_df = gpd.GeoDataFrame({'geometry':[pt]}, crs=CRS.from_epsg(4326))
        
        # project the point to raster CRS
        ds = rxr.open_rasterio(self.item_href[0])
        raster_crs = CRS.from_wkt(ds.spatial_ref.crs_wkt)
        point_df_utm = point_df.to_crs(raster_crs)
        set_x, set_y = point_df_utm['geometry'][0].x, point_df_utm['geometry'][0].y
        
        # return a time series using multiprocessing and helper function
        # call multiprocess with functools.partial
        if pad is None:
            sample_func = partial(self._extract_point_val, set_x=set_x, set_y=set_y, tol=tol, band=band)
        else:
            if type(pad) is not int:
                raise TypeError("parameter pad should be of type int")
            sample_func = partial(self._extract_point_val_window, set_x=set_x, set_y=set_y, tol=tol, pad=pad, band=band)

        with mp.get_context("spawn").Pool(nproc) as pool:
            print(f'using {nproc} processes to sample {len(self.item_href)} assets')
            vals = pool.map(sample_func, self.item_href)
            
        return list(vals)
            
            
    # TODO
    def area_time_series_from_items(self, poly, nproc=2, stat='mean', clip_dict={'all_touched':True, 'crop':True}):
        """ construct a pandas DataFrame which is a time series for a single pixel across the search result."""
        
        # check pt param type
        if type(pt) is not Polygon:
            raise(TypeError, "input pt must be of type shapely.geometry.Polygon")
            
        valid_stats = ('mean', 'std', 'var', 'median', 'min', 'max')
        if stat not in valid_stats:
            raise(ValueError, f"parameter 'stat' must be in {valid_stats}")
            
        poly_df = gpd.GeoDataFrame({'geometry':[poly]}, crs=CRS.from_epsg(4326))
        
        # reproject the polygon to raster CRS
        ds = rxr.open_rasterio(self.item_href[0])
        raster_crs = CRS.from_wkt(ds.spatial_ref.crs_wkt)
        poly_df_utm = poly_df.to_crs(raster_crs)
        
        sample_func = partial(self._extract_area_val_window, poly=poly_df, stat=stat, clip_dict=clip_dict)

        with mp.get_context("spawn").Pool(nproc) as pool:
            print(f'using {nproc} processes to sample {len(self.item_href)} assets')
            vals = pool.map(sample_func, self.item_href)
            
        return list(vals)
            

                               
                               
    
