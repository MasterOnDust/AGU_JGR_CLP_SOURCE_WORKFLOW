import xarray as xr
import DUST
from netCDF4 import num2date, date2num
import numpy as np
from IPython import embed

def create_timeseries(path, lon0=None,lon1=None,lat0=None,lat1=None):
    ds = xr.open_dataset(path)
    with xr.set_options(keep_attrs=True):
        ds = DUST.utils.utils.region_slice(ds, lon0, lon1, lat0,lat1)
    ems=ds.Emission*ds.area #multiply by area to get units of kg
    ems = ems.sum(dim=['lon','lat'])
    ems.attrs['units']='kg'
    ems.attrs['long_name']=ds.Emission.attrs['long_name']
    ems = ems.to_dataset(name='Emission')
    ems.attrs=ds.attrs
    ems.attrs['lon0']=lon0
    ems.attrs['lat0']=lat0
    ems.attrs['lon1']=lon1
    ems.attrs['lat1']=lat1
    ems.time.attrs=ds.time.attrs
    return ems

def resample_emission_flux(paths, lon0=None,lon1=None,
                           lat0=None,lat1=None, frequency='m'):
    """
    DESCRIPTION:
    ===========
        Resample data to monthly or yearly frequency, if frequency is different from 
        monthly ('m') cum_emission variable will be set as the Emission variable.
        
    
    """
    dsets=[]
    d0= DUST.read_flexdust_output(paths[0])
    start_date=d0.startdate
    starthour=d0.starthour
    for path in paths:
        temp_ds=DUST.read_flexdust_output(path)
        temp_ds=DUST.utils.utils.region_slice(temp_ds, lon0,lon1,lat0,lat1)

        if frequency=='m':
            resampled_ds=temp_ds.Emission.resample(time='m').sum().to_dataset(name='Emission')
            resampled_ds['Emission'].attrs['cell_methods']='time: summed'
            resampled_ds['Emission'].attrs['units']='kg/m2' 
            resampled_ds['Emission'].attrs['long_name']='monlthy simulated accumulated dust emissions'
            
            daily = temp_ds.time.resample(time='d').asfreq()
            month_start = list(daily.where(daily.dt.is_month_start,drop=True))
            month_start[0] = temp_ds.time[0]

            month_start= [date2num(month.dt,units='days since 1979-1-1', calendar='standard') for month in month_start]
            time_bnds=[date2num(t.time.dt, units='days since 1979-1-1', calendar='standard') for t in resampled_ds.time]
       
            time=xr.DataArray(data=month_start,
                    dims=['time'],
                    coords=dict(time=month_start),
                    attrs=dict(
                        units='days since 1979-1-1',
                        bounds='time_bnds')
            )
                
            resampled_ds=resampled_ds.assign(time=time)
            bnds=np.array((month_start,time_bnds)).transpose()
            time_bnds=xr.DataArray(data=bnds,
                    dims=['time','nbnd'],
                    coords=dict(time=resampled_ds.time,
                    nbnd=np.arange(bnds.shape[1])) ,
                    attrs=dict(
                    units='days since 1979-1-1')
                )
            resampled_ds=resampled_ds.assign(time_bnds=time_bnds)
            dsets.append(resampled_ds)
        else:
            resampled_ds=temp_ds.cum_emission.to_dataset(name='Emission')
            resampled_ds['Emission'].attrs['cell_methods']='time: summed'
            resampled_ds['Emission'].attrs['units']='kg/m2'   
            resampled_ds['Emission'].attrs['long_name']='yearly accumulated simulated dust emissions'
            t0 = date2num(temp_ds.time[0].dt, units='days since 1979-1-1', calendar='standard')
            time=xr.DataArray(data=t0,
                    dims=['time'],
                    coords=dict(time=[t0]),
                    attrs=dict(
                        units='days since 1979-1-1',
                        bounds='time_bnds'
                    )
            )
            
            resampled_ds=resampled_ds.assign(time=time)
            t1 = date2num(temp_ds.time[-1].dt, units='days since 1979-1-1', calendar='standard')
            bnds=np.array([t0, t1]).reshape(1,2)
            time_bnds=xr.DataArray(data=bnds,
                    dims=['time','nbnd'],
                    coords=dict(time=[t0],
                    nbnd=np.arange(bnds.shape[1])) ,
                    attrs=dict(
                        units='days since 1979-1-1'
                    )
            )
            resampled_ds=resampled_ds.assign(time_bnds=time_bnds)
            resampled_ds.time.attrs=time.attrs
            dsets.append(resampled_ds)
    ds=xr.concat(dsets, dim='time',data_vars=['Emission','time_bnds'],combine_attrs='override')
    ds=ds.assign(cum_emission=ds.Emission.sum(dim='time', keep_attrs=True))
    ds['cum_emission'].attrs['units']='kg/m2'
    ds['cum_emission'].attrs['long_name']='total simulated emissions'
    ds['area']=d0['area']
    ds['soil']=d0['soil']
    edate = str(temp_ds.time[0].dt.strftime('%Y').values)
    sdate = str(d0.time[0].dt.strftime('%Y').values)
    ds.attrs['title']='FLEXDUST simulated dust emissions {}-{}'.format(sdate, edate)
    ds.attrs['source']='FLEXDUST'
    ds.attrs['reference']='10.1002/2016JD025482'
    ds.attrs['soil_map'] = 'ISRIC'
    ds.attrs['lon0']=lon0
    ds.attrs['lat0']=lat0
    ds.attrs['lon1']=lon1
    ds.attrs['lat1']=lat1
    ds.attrs['varName']='Emission'
    
    return ds
    
    
                                
        
            