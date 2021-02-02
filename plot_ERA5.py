import argparse as ap

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
from cartopy.mpl.gridliner import (LATITUDE_FORMATTER, LONGITUDE_FORMATTER,
                                   Gridliner)
def read_data(path, var=None):
    dset = xr.open_dataset(path)
    dset = fix_longitude(dset)
    datavar=list(dset.keys())

    if len(datavar) > 1 and var==None:
        raise(IndexError('netCDF file contains more than one variable and varialbe was not provided'))
    elif var!=None:
        datavar=var
    else:
        datavar=''.join(datavar)
    dset.attrs['varName']=datavar
    return dset


def fix_longitude(dset):
    """
    DESCRIPTION:
    ============
        Center longitude around 180 degree in stead of 0

    """

    lon_attrs=dset.longitude.attrs
    dset = dset.assign_coords(longitude=(((dset.longitude + 180) % 360) - 180))
    dset=dset.sortby('longitude')
    dset.longitude.attrs=lon_attrs
    return dset

def draw_map(ax=None, top_labels=False, right_labels=False):
    if ax==None:
        ax=plt.gca()
    gl = ax.gridlines(transform = ccrs.PlateCarree(), draw_labels = True, linestyle ='--')
    gl.top_labels = top_labels
    gl.right_labels = right_labels
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.coastlines('110m', color='gray', alpha=0.8)
    return ax 

def plot_contour(da,ax=None, colors='black'):
    if ax == None:
        ax = plt.gca()
    da = da.isel(time=0)
    CS = ax.contour(da['longitude'],da['latitude'],da,transform=ccrs.PlateCarree(), colors=colors)
    ax.clabel(CS, CS.levels, inline=True, fmt='%d')
    return ax

def plot_contourf(da, ax=None, cmap='viridis'):
    if ax == None:
        ax = plt.gca()
    da=da.isel(time=0)
    ax.contourf(da['longitude'],da['latitude'],da ,transform=ccrs.PlateCarree(), cmap=cmap)

    return ax 

def plot_pcolormesh(da, ax=None, cmap='viridis'):
    """
    Create pcolormesh plot

    """

    if ax==None:
        ax=plt.gca()
 
    da.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(), cmap=cmap)

    return ax
 

if __name__=='__main__':
    parser=ap.ArgumentParser(description='Command line interface for plotting ERA5 data')
    parser.add_argument('path', help='path to ERA5 dataset')
    parser.add_argument('--method', '--m', default='pcolormesh', help='which ploting function to use')
    parser.add_argument('--title', '--t', default=None, help='title of plot')
    parser.add_argument('--filename', '--fn', default=None, help='Filename to which plot should be stored')
    parser.add_argument('--variable', '--v', default=None, 
            help='If there are more than one variable in the dataset variable name has to be specified')
    args = parser.parse_args()
    path=args.path
    title=args.title
    file_name=args.filename
    var=args.variable
    method=args.method
    dset=read_data(path,var)
    da=dset[dset.varName]
    fig,ax =plt.subplots(subplot_kw={'projection':ccrs.Robinson(central_longitude=90.0)})
    print(method)
    ax = draw_map(ax)
    if method=='pcolormesh':
        ax=plot_pcolormesh(da,ax=ax)
    elif method=='contour':
        ax=plot_contour(da, ax=ax)
    elif method=='contourf':
        ax=plot_contourf(da,ax=ax)
    else:
        raise(ValueError('{} is not a valid plotting method choose either pcolormesh, contour, contourf'.format(method)))
    
    plt.show()
    if file_name.endswith('.pdf'):
        plt.savefig(file_name, bbox_inches='tight')
    elif file_name.endswith('png'):
        plt.savefig(file_name, dpi=300,bbox_inches='tight')
    else: 
        raise(ValueError("invalid file format {}".format(file_name)))
