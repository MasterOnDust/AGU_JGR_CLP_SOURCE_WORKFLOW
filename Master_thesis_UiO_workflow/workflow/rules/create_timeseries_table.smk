
SDATE=config['sdate']
EDATE=config['edate']

SDATE_m = config['m_sdate']
EDATE_m = config['m_edate']

rule calc_receptor_correlations:
    input:
        NAO_EOF_path = 'results/nao/era5.single_level.NAO_EOF.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        AO_EOF_path = 'results/ao/era5.1000hPa.AO_EOF.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        APVI_path = 'results/apv/era5.500hPa.APVI_index.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        MO_path = 'results/eawmi/era5.single_level.EAWM_MO.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        EAWMI_path = 'results/eawmi/era5.single_level.U300hPa_EAWM.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Precipitation_path = 'results/precip/era5.single_level.total_precipitation.{location}.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.csv',
        NAO_station_path = 'results/nao/era5.single_level.NAO_station.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Temp_gradient = 'results/2mt_gradient/era5.single_level.2m_temperature_gradient.East_asia.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        receptor_data = expand('results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               region=['taklamakan','mongolia','north_west', 'total'] ,allow_missing=True)
    
    output:
        outpath='results/correlations_receptor/{kind}/{kind}.{location}.{psize}.correlations.{season}.csv'
    run:
        import pandas as pd
        import xarray as xr
        def read_data(path, sdate_m, edate_m):
            ds = xr.open_dataset(path)
            ds['time'] = ds.time.dt.year 
            ds = ds.sel(time=slice(sdate_m, edate_m))
            return ds
        season=wildcards.season
        NAO_EOF = read_data(input.NAO_EOF_path, config['m_sdate'], config['m_edate'])
        AO_EOF = read_data(input.AO_EOF_path, config['m_sdate'], config['m_edate'])
        APVI = read_data(input.APVI_path, config['m_sdate'], config['m_edate'])
        MO = read_data(input.MO_path, config['m_sdate'], config['m_edate'])
        NAO_station = read_data(input.NAO_station_path, config['m_sdate'], config['m_edate'])
        Temp_gradient = read_data(input.Temp_gradient, config['m_sdate'], config['m_edate'])
        temp_gradient_anomalies = Temp_gradient.anomalies.sel(latitude_bins=slice(40,50)).mean(dim='latitude_bins')
        EAWMI = read_data(input.EAWMI_path, config['m_sdate'], config['m_edate'])
        precip = pd.read_csv(input.Precipitation_path, index_col=0)
        precip.index = pd.to_datetime(precip.index).year
        precip = precip['tp'].loc[config['m_sdate'] : config['m_edate']]
        
        indicies = {
            'NAO_EOF_'+season : NAO_EOF.NAO_EOF.sel(mode=0).values,
            'AO_EOF_'+season : AO_EOF.AO_EOF.sel(mode=0).values,
            'APVI_'+season : APVI.APVI.values,
            'MO_'+season : MO.MO.values,
            'EAWMI_'+season : EAWMI.EAWMI.values,
            'Temp_gradient_anomalies_'+season : temp_gradient_anomalies.values,
            'NAO_station_'+season : NAO_station.NAO,
            
        }
        df_indices = pd.DataFrame(indicies, index=AO_EOF.time.values)
        df_indices = df_indices.join(precip)
        taklamakan =  xr.open_dataset(input.receptor_data[0])
        mongolia = xr.open_dataset(input.receptor_data[1])
        north_west = xr.open_dataset(input.receptor_data[2])
        total = xr.open_dataset(input.receptor_data[3])
        
        receptor_data = {
            'taklamakan' : taklamakan[taklamakan.varName].values,
            'mongolia' : mongolia[mongolia.varName].values,
            'north_west' : north_west[north_west.varName].values,
            'total' : total[total.varName].values
        }
        df_receptor_data = pd.DataFrame(receptor_data, index=taklamakan.year.values)
        data = df_receptor_data.join(df_indices)
        corrs = data.corr().loc[['taklamakan','mongolia','north_west', 'total'], df_indices.keys()]
        corrs.to_csv(output.outpath)