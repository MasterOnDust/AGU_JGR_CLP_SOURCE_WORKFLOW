import xarray as xr
from .process_era5 import read_data

def calc_EAWM_zonal_300hPa(zonal_300hPa_path, frequency):
    ds =read_data(zonal_300hPa_path)
    if frequency=='monthly':
        u_300=ds
    elif frequency =='DJF':
        u_300=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        u_300=u_300.sel(time=(u_300.time.dt.season=='DJF'))
    elif frequency=='MAM':
        # 'Q-NOV' year ends in november 
        u_300=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        u_300=u_300.sel(time=(u_300.time.dt.season=='MAM'))
    elif frequency=='JJA':
        # 'Q-NOV' year ends in november 
        u_300=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        u_300=u_300.sel(time=(u_300.time.dt.season=='JJA'))
    elif frequency=='SON':
        # 'Q-NOV' year ends in november 
        u_300=ds.resample(time='Q-NOV').mean(keep_attrs=True)
        u_300=u_300.sel(time=(u_300.time.dt.season=='JJA'))
    else:
        raise(ValueError("Invalid sampling freqency {}".format(frequency)))

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



    

