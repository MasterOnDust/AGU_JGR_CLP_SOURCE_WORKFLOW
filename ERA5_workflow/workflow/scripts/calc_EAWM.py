#!/usr/bin/env python

import xarray as xr
from .process_era5 import read_data

def calc_EAWM_zonal_300hPa(u_300):


    siberian_high_domain=((u_300.longitude >= 80) & (u_300.longitude <= 140) & 
                (u_300.latitude <= 60) & (u_300.latitude >= 50))
    alutian_low_domain=((u_300.longitude >= 110) & (u_300.longitude <= 170) & 
                (u_300.latitude <= 37.5) & (u_300.latitude >= 27.5))
    alutian_low = u_300.where(alutian_low_domain,drop=True).mean(dim=['longitude','latitude'])
    siberian_high = u_300.where(siberian_high_domain,drop=True).mean(dim=['longitude','latitude'])

    eawmi=alutian_low['u']-siberian_high['u'] 
    eawmi = eawmi.to_dataset(name='EAWMI')
    eawmi['EAWMI'].attrs['long_name']='East Asian Winter Monsoon index'
    eawmi.attrs['reference'] = " Jhun, J., & Lee, E. (2004)"
    eawmi.attrs['doi'] = "https://doi.org/10.1175/1520-0442(2004)017%3C0711:ANEAWM%3E2.0.CO;2"
    return eawmi

def calc_MO_index(mslp_data):
    mslp_Nemuro=mslp_data.msl.sel(longitude=145,latitude=43.75)
    mslp_Irkutsk=mslp_data.msl.sel(longitude=105,latitude=52.5) 
    standardized_mo =(((mslp_Irkutsk-mslp_Irkutsk.mean(dim='time'))/mslp_Irkutsk.std()) - 
                    ((mslp_Nemuro-mslp_Nemuro.mean(dim='time'))/mslp_Nemuro.std()))
    mo = standardized_mo.to_dataset(name='MO')
    mo.MO.attrs['long_name']='Standardized MO index'
    mo.attrs['reference'] = "https://doi.org/10.1029/2008JD010824"
    mo.attrs['title'] = 'Standardized MO index, based the difference in mslp between Irkutsk and Nemuro'
    
    return mo

