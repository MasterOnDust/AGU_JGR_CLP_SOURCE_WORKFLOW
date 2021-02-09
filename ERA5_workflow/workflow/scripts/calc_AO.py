import numpy as np
from eofs.xarray import Eof

def calc_AO_EOF(geopot_data_ano):
    """
    DESCRIPTION
    ===========
        Calculates the AO using the EOFS python package.
        Input takes 1000hPa geopotential height anomalies.
        
        Input:
            mslp_anaomalies: xarray.dataset containing mean sea level pressure anomailies. 

        Returns: 
            out_data: the scaled 1th principal component time series 
    """

    bolean=((geopot_data_ano.latitude <= 90) & (geopot_data_ano.latitude >= 20))
    geopot_data_ano=geopot_data_ano.where(bolean, drop=True)

    coslat = np.cos(np.deg2rad(geopot_data_ano.coords['latitude'].values)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(geopot_data_ano, weights=wgts)

    # Retrieve the leading EOF, expressed as the covariance between the leading PC
    # time series and the input SLP anomalies at each grid point.
    eof1_corr = solver.eofsAsCorrelation(neofs=1)
    eof1_cov = solver.eofsAsCovariance(neofs=1)
    out_data=solver.pcs(1,1).copy()
    out_data=out_data.to_dataset(name='AO_EOF')
    out_data['AO_EOF'].attrs['standard_name']='arctic_oscillation_index'
    out_data['AO_EOF'].attrs['long_name']='Normalized 1th principal component of 1000hPa geopotential height'
    out_data['eof1_corr']=eof1_corr
    out_data['eof1_cov']=eof1_cov
    return out_data
