import xarray as xr
import numpy as np

def calculate_temperature_gradient(temp_dataset, n_latitude_bins=17):
    """
    DESCRIPTION:
    ===========
        Calculates the gradient in zonally averaged temperature between 30 and 
        120 degress East. 
    
    
    """
    ds = temp_dataset[temp_dataset.varName].sel(longitude=slice(30,120), 
                                                latitude=slice(80,0)).mean(dim='longitude',keep_attrs=True)
    
    ds_grouped=ds.groupby_bins('latitude',17).mean(keep_attrs=True)
    ds_grouped['latitude_bins'] = [intervall.mid for intervall in ds_grouped.latitude_bins.values]
    temp_gradient = ds_grouped.differentiate('latitude_bins')
    temp_gradient['latitude_bins'].attrs['long_name']='Midpoint of latitude bin'
    temp_gradient = temp_gradient.to_dataset(name='t2m_gradient')
    temp_gradient['anomalies'] = temp_gradient['t2m_gradient'] - temp_gradient['t2m_gradient'].mean(dim='time')
    temp_gradient['t2m_gradient'].attrs['units']='K/deg(latitude)'
    temp_gradient['t2m_gradient'].attrs['long_name']='Temperature gradient'
    temp_gradient['anomalies'].attrs['units']='K/deg(latitude)'
    temp_gradient['anomalies'].attrs['long_name']='Temperature gradient anomalies'
    temp_gradient['t2m'] = ds_grouped
    return temp_gradient
    
    
    