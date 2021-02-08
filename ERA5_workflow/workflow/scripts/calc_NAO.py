#!/usr/bin/env python

import numpy as np
from eofs.xarray import Eof

def calc_NAO_station(mslp_anomalies):
    """
    DESCRIPTION
    ===========
        Calculates the Station based NAO index, which is defined as the 
        normalized difference in sea level pressure between Reykjavik and 
        Lisbon.

        Input:
            mslp_anaomalies: xarray.dataset containing mean sea level pressure anomailies. 

        Returns:
            NAO station index : norm(Lisbon mean sea level pressure anomalies) - 
                                norm(Reykjavik mean sea level pressure anomalies) 


    """

    ds_reykjavik=mslp_anomalies.sel(longitude=-21.93, latitude=64.13, method="nearest")
    ds_Lisbon=mslp_anomalies.sel(longitude=-9.14, latitude=37.71, method="nearest")
    Lisbon_NAO=ds_Lisbon.groupby('time.month')/ds_Lisbon.groupby('time.month').std()
    reykjavik_NAO=ds_reykjavik.groupby('time.month')/ds_reykjavik.groupby('time.month').std()
    return Lisbon_NAO-reykjavik_NAO

def calc_NAO_EOF(mslp_data_ano):
    """
    DESCRIPTION
    ===========
        Calculates the EOF based NAO using the EOFS python package.
        Input takes mean sea level pressure anomalies.
        
        Input:
            mslp_anaomalies: xarray.dataset containing mean sea level pressure anomailies. 

        Returns: 
            out_data: the scaled 1th principal component time series 
    """

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
    out_data=out_data.to_dataset(name='NAO_EOF')
    out_data['NAO_EOF'].attrs['standard_name']='EOF_North_atlantic_oscillation_index'
    out_data['NAO_EOF'].attrs['long_name']='Normalized 1th principal component of MSLP'
    out_data['eof1_corr']=eof1_corr
    out_data['eof1_cov']=eof1_cov
    return out_data



