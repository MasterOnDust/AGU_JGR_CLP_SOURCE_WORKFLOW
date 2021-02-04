import numpy as np
import xarray as xr

from eofs.xarray import Eof
from .process_era5 import read_data

def calc_NAO_station(mslp_anomalies):
    ds_reykjavik=mslp_anomalies.sel(longitude=-21.93, latitude=64.13, method="nearest")
    ds_Lisbon=mslp_anomalies.sel(longitude=-9.14, latitude=37.71, method="nearest")
    Lisbon_NAO=ds_Lisbon.groupby('time.month')/ds_Lisbon.groupby('time.month').std()
    reykjavik_NAO=ds_reykjavik.groupby('time.month')/ds_reykjavik.groupby('time.month').std()
    return Lisbon_NAO-reykjavik_NAO

def calc_NAO_EOF(monthly_mslp_path, freqency):
    """
    DESCRIPTION
    ===========
        Calculates the EOF based NAO using the EOFS python package. Returning 
        the scaled 1th principal component 
    """
    monthly_mslp=read_data(monthly_mslp_path)
    if freqency=='monthly':
        mlsp_data=monthly_mslp
    elif freqency =='DJF':
        mlsp_data=monthly_mslp.resample(time='Q-NOV').mean(keep_attrs=True)
        mlsp_data=mlsp_data.sel(time=(mlsp_data.time.dt.season=='DJF'))
    elif freqency=='MAM':
        # 'Q-NOV' year ends in november 
        mlsp_data=monthly_mslp.resample(time='Q-NOV').mean(keep_attrs=True)
        mlsp_data=mlsp_data.sel(time=(mlsp_data.time.dt.season=='MAM'))
    elif freqency=='JJA':
        # 'Q-NOV' year ends in november 
        mlsp_data=monthly_mslp.resample(time='Q-NOV').mean(keep_attrs=True)
        mlsp_data=mlsp_data.sel(time=(mlsp_data.time.dt.season=='JJA'))
    elif freqency=='SON':
        # 'Q-NOV' year ends in november 
        mlsp_data=monthly_mslp.resample(time='Q-NOV').mean(keep_attrs=True)
        mlsp_data=mlsp_data.sel(time=(mlsp_data.time.dt.season=='JJA'))
    else:
        raise(ValueError("Invalid sampling freqency {}".format(freqency)))
    
    mslp_data_ano= mlsp_data['msl']-mlsp_data['msl'].mean(dim='time')

    bolean=((mslp_data_ano.longitude >= -90) & (mslp_data_ano.longitude <= 40) & 
                (mslp_data_ano.latitude <= 80) & (mslp_data_ano.latitude >= 20))
    mslp_data_ano=mslp_data_ano.where(bolean, drop=True)

    coslat = np.cos(np.deg2rad(mslp_data_ano.coords['latitude'].values)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(mslp_data_ano, weights=wgts)

    # Retrieve the leading EOF, expressed as the covariance between the leading PC
    # time series and the input SLP anomalies at each grid point.
    eof1_corr = solver.eofsAsCorrelation(neofs=1)
    eof1_cov = solver.eofsAsCovariance(neofs=1)
    out_data=solver.pcs(1,1).copy()

    out_data['eof1_corr']=eof1_corr
    out_data['eof1_cov']=eof1_cov

    return out_data



