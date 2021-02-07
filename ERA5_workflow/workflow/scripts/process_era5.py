import xarray as xr
import time

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

def calculate_monthly_clm(path):
    ds=read_data(path)
    ds=ds.groupby('time.month').mean(keep_attrs=True)
    ds.attrs['history']=time.ctime() +' xarray.Dataset.groupby(\'time.month\').mean(keep_attrs=True) ' + ds.attrs['history']
    return ds

def resample_seasonal(path,season):
    ds=read_data(path)
    if season=='DJF':
        ds=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        ds=ds.sel(time=(ds.time.dt.season=='DJF'))
    elif season=='MAM':
        ds=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        ds=ds.sel(time=(ds.time.dt.season=='MAM'))
    elif season=='JJA':
        ds=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        ds=ds.sel(time=(ds.time.dt.season=='JJA'))
    elif season=='SON':
        ds=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        ds=ds.sel(time=(ds.time.dt.season=='SON'))
    return ds
    