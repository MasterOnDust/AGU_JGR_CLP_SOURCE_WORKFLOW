
SDATE=config['sdate']
EDATE=config['edate']

SDATE_m = config['m_sdate']
EDATE_m = config['m_edate']

rule calc_receptor_correlations:
    input:
        NAO_EOF_path_djf = 'results/nao/era5.single_level.NAO_EOF.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        AO_EOF_path_djf = 'results/ao/era5.1000hPa.AO_EOF.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        APVI_path_djf = 'results/apv/era5.500hPa.APVI_index.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        MO_path_djf = 'results/eawmi/era5.single_level.EAWM_MO.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        SHI_path_djf = 'results/eawmi/era5.single_level.SHI.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        EAWMI_path_djf = 'results/eawmi/era5.single_level.U300hPa_EAWM.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Precipitation_path_djf = expand('results/precip/era5.single_level.total_precipitation.{location}.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.csv',
                                    location=config['receptors'].keys() ,allow_missing=True),
        NAO_station_path_djf = 'results/nao/era5.single_level.NAO_station.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Temp_gradient_djf = 'results/2mt_gradient/era5.single_level.2m_temperature_gradient.East_asia.DJF.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        NAO_EOF_path = 'results/nao/era5.single_level.NAO_EOF.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        AO_EOF_path = 'results/ao/era5.1000hPa.AO_EOF.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        APVI_path = 'results/apv/era5.500hPa.APVI_index.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        MO_path = 'results/eawmi/era5.single_level.EAWM_MO.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        SHI_path = 'results/eawmi/era5.single_level.SHI.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        EAWMI_path = 'results/eawmi/era5.single_level.U300hPa_EAWM.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Precipitation_path = expand('results/precip/era5.single_level.total_precipitation.{location}.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.csv',
                            location=config['receptors'].keys() ,allow_missing=True),
        NAO_station_path = 'results/nao/era5.single_level.NAO_station.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        Temp_gradient = 'results/2mt_gradient/era5.single_level.2m_temperature_gradient.East_asia.MAM.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        emission_data = expand(config['flexdust_results']+'/emission_flux.time_series.{region}.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                        region=['taklamakan','mongolia','north_west', 'total'],allow_missing=True),
        receptor_data_wetdep_2micron = expand('results/model_results/time_series/wetdep/wetdep.{location}.total.2micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys() ,allow_missing=True),
        receptor_data_drydep_2micron = expand('results/model_results/time_series/drydep/drydep.{location}.total.2micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys() ,allow_missing=True),
        receptor_data_wetdep_20micron = expand('results/model_results/time_series/wetdep/wetdep.{location}.total.20micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys() ,allow_missing=True),
        receptor_data_drydep_20micron = expand('results/model_results/time_series/drydep/drydep.{location}.total.20micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys() ,allow_missing=True)
        
    output:
        outpath='results/timeseries_table.csv'
    run:
        import pandas as pd
        import xarray as xr
        from thesis_toolbox.composites.create_composites import detrend_timeseries
        def read_data(path, sdate_m, edate_m, detrend=False):
            ds = xr.open_dataset(path)
            ds['time'] = ds.time.dt.year 
            if detrend:
                ds[ds.varName] = detrend_timeseries(ds[ds.varName])
            else:
                ds = ds.sel(time=slice(sdate_m, edate_m))
            
            return ds
        def prepare_precip_data(paths, sdate_m,edate_m, season, detrend=False):
            df = pd.DataFrame()    
            for precip_path in paths:
                receptor_name = precip_path.split('/')[-1].split('.')[3] 
                precip = pd.read_csv(precip_path, index_col=0)
                precip.index = pd.to_datetime(precip.index).year
                if detrend:
                    df['Precip '+receptor_name+' '+season] = detrend_timeseries(precip['tp'].loc[sdate_m : edate_m])
                else:
                    df['Precip '+receptor_name+' '+season] = precip['tp'].loc[sdate_m : edate_m]
            return df
        def prepare_depostion_data(paths,psize, detrend=False):
            df = pd.DataFrame()
            for path in paths:
                receptor_name = path.split('/')[-1].split('.')[1]
                kind = path.split('/')[-1].split('.')[0]
                ds = xr.open_dataset(path)
                if detrend:
                    df[f'{receptor_name} {kind} {psize}'] = detrend_timeseries(ds[ds.varName].values)
                else:
                    df[f'{receptor_name} {kind} {psize}'] = ds[ds.varName].values
            df.index = ds.year.values
            return df
        # Read MAM indices 
        NAO_EOF_MAM = read_data(input.NAO_EOF_path, config['m_sdate'], config['m_edate'])
        AO_EOF_MAM = read_data(input.AO_EOF_path, config['m_sdate'], config['m_edate'])
        APVI_MAM = read_data(input.APVI_path, config['m_sdate'], config['m_edate'])
        MO_MAM = read_data(input.MO_path, config['m_sdate'], config['m_edate'])
        SHI_MAM = read_data(input.SHI_path, config['m_sdate'], config['m_edate'])
        NAO_station_MAM = read_data(input.NAO_station_path, config['m_sdate'], config['m_edate'])
        Temp_gradient_MAM = read_data(input.Temp_gradient, config['m_sdate'], config['m_edate'])
        temp_gradient_anomalies_MAM = Temp_gradient_MAM.anomalies.sel(latitude_bins=slice(40,50)).mean(dim='latitude_bins')
        EAWMI_MAM = read_data(input.EAWMI_path, config['m_sdate'], config['m_edate'])

        #Read DJF indicies
        NAO_EOF_DJF = read_data(input.NAO_EOF_path_djf, config['m_sdate'], config['m_edate'])
        AO_EOF_DJF = read_data(input.AO_EOF_path_djf, config['m_sdate'], config['m_edate'])
        APVI_DJF = read_data(input.APVI_path_djf, config['m_sdate'], config['m_edate'])
        MO_DJF = read_data(input.MO_path_djf, config['m_sdate'], config['m_edate'])
        SHI_DJF = read_data(input.SHI_path_djf, config['m_sdate'], config['m_edate'])
        NAO_station_DJF = read_data(input.NAO_station_path_djf, config['m_sdate'], config['m_edate'])
        Temp_gradient_DJF = read_data(input.Temp_gradient_djf, config['m_sdate'], config['m_edate'])
        temp_gradient_anomalies_DJF = Temp_gradient_DJF.anomalies.sel(latitude_bins=slice(40,50)).mean(dim='latitude_bins')
        EAWMI_DJF = read_data(input.EAWMI_path_djf, config['m_sdate'], config['m_edate'])

        precip_djf = prepare_precip_data(input.Precipitation_path_djf,config['m_sdate'], config['m_edate'],'DJF')
        precip_mam = prepare_precip_data(input.Precipitation_path,config['m_sdate'],config['m_edate'],'MAM')

        
        indicies = {
            'NAO EOF MAM' : NAO_EOF_MAM.NAO_EOF.sel(mode=0).values,
            'AO EOF MAM' : AO_EOF_MAM.AO_EOF.sel(mode=0).values,
            'APVI MAM' : APVI_MAM.APVI.values,
            'MO MAM' : MO_MAM.MO.values,
            'SHI MAM' : SHI_MAM.SHI.values,
            'EAWMI MAM' : EAWMI_MAM.EAWMI.values,
            'Temp gradient anomalies MAM' : temp_gradient_anomalies_MAM.values,
            'NAO station MAM' : NAO_station_MAM.NAO,
            'NAO EOF DJF' : NAO_EOF_DJF.NAO_EOF.sel(mode=0).values,
            'AO EOF DJF' : AO_EOF_DJF.AO_EOF.sel(mode=0).values,
            'APVI DJF' : APVI_DJF.APVI.values,
            'MO DJF' : MO_DJF.MO.values,
            'SHI DJF' : SHI_DJF.SHI.values,
            'EAWMI DJF' : EAWMI_DJF.EAWMI.values,
            'Temp gradient anomalies_DJF' : temp_gradient_anomalies_DJF.values,
            'NAO station DJF' : NAO_station_DJF.NAO
            
        }
        df_indices = pd.DataFrame(indicies, index=AO_EOF_DJF.time.values)
        df_indices = df_indices.join([precip_mam,precip_djf])
        taklamakan =  xr.open_dataset(input.emission_data[0])
        mongolia = xr.open_dataset(input.emission_data[1])
        north_west = xr.open_dataset(input.emission_data[2])
        total = xr.open_dataset(input.emission_data[3])
        
        emission_data = {
            'Emissions taklamakan' : taklamakan[taklamakan.varName].values,
            'Emissions mongolia' : mongolia[mongolia.varName].values,
            'Emissions north_west' : north_west[north_west.varName].values,
            'Emissions total' : total[total.varName].values
        }

        df_emission_data = pd.DataFrame(emission_data, index=taklamakan.time.dt.year.values)
        data = df_emission_data.join(df_indices)
        depo_wet_2micron = prepare_depostion_data(input.receptor_data_wetdep_2micron,'2micron')
        depo_wet_20micron = prepare_depostion_data(input.receptor_data_wetdep_20micron,'20micron')
        depo_dry_2micron = prepare_depostion_data(input.receptor_data_drydep_2micron,'2micron')
        depo_dry_20micron = prepare_depostion_data(input.receptor_data_drydep_20micron,'20micron')
        data = data.join([depo_wet_2micron, depo_wet_20micron,depo_dry_2micron,depo_dry_20micron])
        data.to_csv(output.outpath)

