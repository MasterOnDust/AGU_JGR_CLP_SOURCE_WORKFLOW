import xarray as xr
import DUST

def create_timeseries(ds, lon0=None,lon1=None,lat0=None,lat1=None):
    with xr.set_options(keep_attrs=True):
        ds = DUST.utils.utils.region_slice(ds, lon0, lon1, lat0,lat1)
    ts = ds[ds.varName].sum(dim=['lon','lat'], keep_attrs=True)
    
    ts = ts.to_dataset(name=ds.varName)
    ts.attrs=ds.attrs
    if lon0:
        ts.attrs['lon0']=lon0
    if lat0:
        ts.attrs['lat0']=lat0
    if lon1:
        ts.attrs['lon1']=lon1
    if lat1:
        ts.attrs['lat1']=lat1
    return ts