################################################################################
# Created by Ove Haugvaldstad                                                  #
# 23.02.2021                                                                   #
#                                                                              #
################################################################################

import pandas as pd
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
 
def select_years_to_composite(timeseries,criterion='1-std'):
    """
    DESCRIPTION:
    ===========
        Chooses which years to be included in the composites, the selection criterion can 
        either be 2, 1 or 0.5 standard devivation around the mean.
    
    USAGE:
    ======
        weak_years, strong_years = select_years_to_composite(timeseires, criterion='1-std')
    
    PARAMETERS:
    ==========
        timeseries : either xarray DataArray or pandas.series
        criterion : valid values '0.5-std', '1-std' and '2-std'
        
    """
    
    if criterion == '05-std':
        c = 0.5
    elif criterion == '1-std':
        c = 1
    elif criterion == '2-std':
        c = 2
    else:
        raise(ValueError('Invalid criterion provided: {}'.format(criterion)))
    
    std = c*timeseries.std()
    mean = timeseries.mean()
    
    if isinstance(timeseries, xr.core.dataarray.DataArray):
        if 'year' in timeseries.dims:
            strong_years =timeseries.where(timeseries > (mean+std), drop=True).year.values
            weak_years=timeseries.where(timeseries < (mean-std), drop=True).year.values
        else:
            strong_years = timeseries.where(timeseries > (mean+std), drop=True).time.dt.year.values
            weak_years = timeseries.where(timeseries < (mean-std), drop=True).time.dt.year.values
    elif isinstance(timeseries,pd.core.series.Series):
        strong_years = timeseries.where(timeseries > (mean+std)).dropna.index.year.values
        weak_years = timesereis.where(timeseries < (mean-std)).dropna.index.year.values
    else:
        raise(ValueError('Invalid datatype provided'))
    return weak_years,strong_years

    
def create_composite(years_to_composite,contour_f=None,
                               contour=None, 
                               quiver=None, 
                               pcolormesh=None,
                               composite_anomailies=False,
                                composite_anomalies_kw = {},
                                calc_std=True):
    """
    DESCRIPTION:
    ===========
        Create a composite based on an array of years. Return dataset containing composite data
        , depending on wether the data is provided as contour_f, contour or quiver, determine how
        the comopsite data will interface with the plotting functions. If composite_anomalies is set 
        to True provided the anomalies of the composite compared to climatetology is also calculated. 
    
    USAGE:
    ======
        composite_data = create_composite([2005,2007,2009], contour_f=temperature_data, quiver=winds)
    
    PARAMETERS:
    ==========
        years_to_composite : list/array of years make composite of 
        contour_f : 2D xarray DataArray, will be interfaced as contourf plot 
        contour : 2D xarray DataArray, will be interfaced as contour plot
        pcolormesh : 2D xarray DataArray will be interfaced as pcolormesh plot 
                        (cannot be combined with contourf)
        quiver : xarray Dataset containing both u and v winds.
        climatology_anomalies : Boolean default (False) if true the anomaly of the composi is also computed
       
    RETURNS:
    =======
        xarray Dataset containing the calculated composit, with interface ready for plotting
    """

            
    
    
    if all(v is None for v in [contour,contour_f,pcolormesh, quiver]):
        raise(ValueError('No data provided'))
    if all(v is not None for v in [contour_f, pcolormesh]):
        raise(ValueError('Both contour_f and pcolormesh cannot be provided at the same time'))
    if len(years_to_composite)==0:
        raise(ValueError('Years to composite cannot be of 0 size'))
    variables = {}
    composite_data = xr.Dataset()
    composite_data.attrs['title'] = 'Composite mean data'
    if isinstance(contour,xr.core.dataarray.DataArray):
        contour_varName = contour.name
        variables[contour_varName] = contour
        contour['time'] = contour.time.dt.year
        composite = _average_composite(years_to_composite, contour)
        composite_data[contour_varName] = composite
        composite_data[contour_varName+'_clim'] = _calculate_climatology(contour, **composite_anomalies_kw)
        composite_data.attrs['contour'] = contour_varName
    
    if quiver:
        variables['quiver'] = quiver
        quiver['time'] = quiver.time.dt.year
        
        composite = _average_composite(years_to_composite,quiver)
        composite_data['u'] = composite['u']
        composite_data['v'] = composite['v']
        composite_data['u_clim'] = _calculate_climatology(quiver['u'], **composite_anomalies_kw)
        composite_data['v_clim'] = _calculate_climatology(quiver['v'], **composite_anomalies_kw)
        composite_data.attrs['quiver'] = ['u','v']
    
    if isinstance(contour_f, xr.core.dataarray.DataArray):
        contourf_varName = contour_f.name
        variables[contourf_varName] = contour_f
        contour_f['time'] = contour_f.time.dt.year
        composite=_average_composite(years_to_composite, contour_f)
        composite_data[contourf_varName] = composite
        composite_data[contourf_varName+'_clim'] = _calculate_climatology(contour_f, **composite_anomalies_kw)
        composite_data.attrs['contourf'] = contourf_varName
    
    if isinstance(pcolormesh, xr.core.dataarray.DataArray):
        pcolormesh_varName = pcolormesh.name
        variables[pcolormesh_varName] = pcolormesh
        pcolormesh['time'] = pcolormesh.time.dt.year
        composite = _average_composite(years_to_composite,pcolormesh)
        composite_data[pcolormesh_varName] = composite
        composite_data[pcolormesh_varName+'_clim'] = _calculate_climatology(pcolormesh, **composite_anomalies_kw)
        composite_data.attrs['pcolormesh'] = pcolormesh_varName
    if composite_anomailies:
        for key, item in variables.items():
            if key=='quiver':
                climatology = _calculate_climatology(item['u'], **composite_anomalies_kw)
                composite_data['u_anomalies']=composite_data['u'] - climatology
                climatology = _calculate_climatology(item['v'], **composite_anomalies_kw)
                composite_data['v_anomalies']=composite_data['v'] - climatology
            else:
                climatology = _calculate_climatology(item, **composite_anomalies_kw)
                composite_data[key+'_anomalies'] =composite_data[key] - climatology

        composite_data.attrs['anomalies'] = 'True' 
        composite_data.attrs['climatology'] = [climatology.attrs['start_year'], climatology.attrs['end_year']]
    else:
        composite_data.attrs['anomalies'] = 'False' 
    
    composite_data.attrs['years_composited']=years_to_composite
    
    if calc_std:
        for key, item in variables.items():
            if key=='quiver':
                composite_data['u_std'] = _calculate_composite_std(years_to_composite,item['u'])
                composite_data['v_std'] = _calculate_composite_std(years_to_composite,item['v'])
            else:
                 composite_data[key+'_std'] = _calculate_composite_std(years_to_composite,item)
#     from IPython import embed; embed() 
    return composite_data    
        
def _average_composite(years_to_composite,data):
    """averages the composite years"""
    if len(years_to_composite) ==1:
        data = data.sel(time=years_to_composite)
    else:
        data = data.sel(time=years_to_composite).mean(dim='time',keep_attrs=True)
    return data

def _calculate_composite_std(years_to_composite, data):
    """calculate the standard diviation of the composite years"""
    data = data.sel(time=years_to_composite).std(dim='time', keep_attrs=True)
    return data
    

def _calculate_climatology(data, start_year=None, end_year=None):
    """calclulates climatology for calculating """
    if start_year == None:
        if isinstance(data.time[0].values,np.datetime64):
            start_year = str(data.time[0].dt.year.values)
        else:
            start_year = int(data.time[0].values)
    if end_year == None:
        if isinstance(data.time[0].values,np.datetime64):
            end_year = str(data.time[-1].dt.year.values)
        else:
            end_year = int(data.time[-1].values)
    data=data.sel(time=slice(start_year,end_year)).mean(dim='time', keep_attrs=True)
    data.attrs['start_year'] = start_year
    data.attrs['end_year'] = end_year
    return data        
    
    
        