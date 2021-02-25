SDATE = config['sdate']
EDATE = config['edate']

rule windspeed_geopot_wind_composite:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        geopot=config['intermediate_files'] + '/era5.{plevel}.GeopotHeight.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/windspeed_geopot/era5.wind_geopot.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        if wildcards.weak_strong == 'strong':
            years_to_compsite = strong_years
        else:
            years_to_compsite = weak_years
        comp_type = wildcards.comp_type
        if comp_type == 'std_anom':
            composite_anomailies=True
            calc_std=True
        elif comp_type=='anom':
            composite_anomailies=True
            calc_std=False
        elif comp_type=='std':
            composite_anomailies=False
            calc_std=True
        else:
            composite_anomailies=False
            calc_std=False
        wind_u = xr.open_dataset(input.u_wind)
        wind_v = xr.open_dataset(input.v_wind)
        windspeed = xr.Dataset()
        winds=xr.Dataset()
        winds['u'] = wind_u['u']
        winds['v'] = wind_v['v']
        windspeed['hws']=np.sqrt(wind_u['u']**2 + wind_v['v']**2)
        windspeed.attrs['varName']='hws'
        windspeed['hws'].attrs['units'] = 'm/s'
        windspeed['hws'].attrs['long_name'] = 'wind speed'
        geopot = xr.open_dataset(input.geopot)
        geopot= geopot[geopot.varName]
        # This is kind of a bad hack, I'll make a new rule that remove empty datasets
        if len(years_to_compsite) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:     
            comp = create_composite(years_to_compsite, contour_f=windspeed['hws'], contour=geopot,
                                            quiver=winds, composite_anomailies=composite_anomailies,
                                            calc_std=calc_std, composite_anomalies_kw={'start_year':params.clim_sdate, 
                                                                                       'end_year':params.clim_edate})
        comp.to_netcdf(output.outpath)

rule mean_sea_level_pressure_and_wind:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        msl=config['intermediate_files'] + '/era5.single_level.mean_sea_level_pressure.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/msl_pressure_wind/era5.wind_msl.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        if wildcards.weak_strong == 'strong':
            years_to_compsite = strong_years
        else:
            years_to_compsite = weak_years
        comp_type = wildcards.comp_type
        if comp_type == 'std_anom':
            composite_anomailies=True
            calc_std=True
        elif comp_type=='anom':
            composite_anomailies=True
            calc_std=False
        elif comp_type=='std':
            composite_anomailies=False
            calc_std=True
        else:
            composite_anomailies=False
            calc_std=False
        wind_u = xr.open_dataset(input.u_wind)
        wind_v = xr.open_dataset(input.v_wind)
        windspeed = xr.Dataset()
        winds=xr.Dataset()
        winds['u'] = wind_u['u']
        winds['v'] = wind_v['v']
        windspeed['hws']=np.sqrt(wind_u['u']**2 + wind_v['v']**2)
        windspeed.attrs['varName']='hws'
        windspeed['hws'].attrs['units'] = 'm/s'
        windspeed['hws'].attrs['long_name'] = 'wind speed'
        msl= xr.open_dataset(input.msl)
        msl = msl[msl.varName]
        if len(years_to_compsite) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:
            comp = create_composite(years_to_compsite, contour_f=windspeed['hws'], contour=msl,
                                            quiver=winds, composite_anomailies=composite_anomailies,
                                            calc_std=calc_std, composite_anomalies_kw={'start_year':params.clim_sdate, 
                                                                                       'end_year':params.clim_edate})
        comp.to_netcdf(output.outpath)

rule snow_cover_and_temperature:
    input:
        snow_depth=config['intermediate_files'] + '/era5.single_level.snow_depth.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        temperature_2m=config['intermediate_files'] + '/era5.single_level.2m_temperature.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/snow_temperature/era5.temp_2m_snow.single_level.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']   
    threads: 1
    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        if wildcards.weak_strong == 'strong':
            years_to_compsite = strong_years
        else:
            years_to_compsite = weak_years
        
        snow_depth = xr.open_dataset(input.snow_depth)
#         from IPython import embed; embed()
        snow_depth = snow_depth[snow_depth.varName]
        temperature_2m = xr.open_dataset(input.temperature_2m)
        temperature_2m = temperature_2m[temperature_2m.varName]
        comp_type = wildcards.comp_type
        if comp_type == 'std_anom':
            composite_anomailies=True
            calc_std=True
        elif comp_type=='anom':
            composite_anomailies=True
            calc_std=False
        elif comp_type=='std':
            composite_anomailies=False
            calc_std=True
        else:
            composite_anomailies=False
            calc_std=False
        if len(years_to_compsite) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:
            comp = create_composite(years_to_compsite, pcolormesh=snow_depth, contour=temperature_2m,
                                                composite_anomailies=composite_anomailies,
                                                calc_std=calc_std, composite_anomalies_kw={'start_year':params.clim_sdate, 
                                                                                       'end_year':params.clim_edate})
        comp.to_netcdf(output.outpath)

rule precipitation_and_mslp:
    input:
        mslp=config['intermediate_files'] + '/era5.single_level.mean_sea_level_pressure.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        total_precipitation=config['intermediate_files'] + '/era5.single_level.total_precipitation.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/snow_temperature/composites/precip_mslp/era5.precip_mslp.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        
        if wildcards.weak_strong == 'strong':
            years_to_compsite = strong_years
        else:
            years_to_compsite = weak_years
        precip = xr.open_dataset(input.total_precipitation)
        precip = precip[precip.varName]
        mslp = xr.open_dataset(input.mslp)
        mslp = mslp[mslp.varName]
        comp_type = wildcards.comp_type
        if comp_type == 'std_anom':
            composite_anomailies=True
            calc_std=True
        elif comp_type=='anom':
            composite_anomailies=True
            calc_std=False
        elif comp_type=='std':
            composite_anomailies=False
            calc_std=True
        else:
            composite_anomailies=False
            calc_std=False
        
        if len(years_to_compsite) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:
            comp =  create_composite(years_to_compsite, pcolormesh=total_precipitation, contour=mslp,
                                            composite_anomailies=composite_anomailies,
                                            calc_std=calc_std, composite_anomalies_kw={'start_year':params.clim_sdate, 
                                                                                       'end_year':params.clim_edate})
        comp.to_netcdf(output.outpath)