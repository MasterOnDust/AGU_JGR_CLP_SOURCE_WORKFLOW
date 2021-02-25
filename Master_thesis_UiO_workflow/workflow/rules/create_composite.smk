SDATE = config['sdate']
EDATE = config['edate']

rule windspeed_geopot_wind_composite:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        geopot=config['intermediate_files'] + '/era5.{plevel}.GeopotHeight.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/windspeed_geopot/era5.wind_geopot.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'

    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
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
        if wildcards.weak_strong == 'strong':
            comp = create_composite(strong_years, contour_f=windspeed['hws'], contour=geopot,
                                            quiver=winds, composite_anomailies=composite_anomailies,
                                            calc_std=calc_std)
        else:     
            comp = create_composite(weak_years, contour_f=windspeed['hws'], contour=geopot,
                                            quiver=winds, composite_anomailies=composite_anomailies,
                                            calc_std=calc_std)
        comp.to_netcdf(output.outpath_low)

rule mean_sea_level_pressure_and_wind:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        msl=config['intermediate_files'] + '/era5.single_level.mean_sea_level_pressure.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/msl_pressure_wind/era5.wind_msl.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'

    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
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
        if wildcards.weak_strong == 'strong':
            comp = create_composite(strong_years, contour_f=windspeed['hws'], contour=msl,
                                                quiver=winds, composite_anomailies=composite_anomailies,
                                                calc_std=calc_std)
        else:
            comp = create_composite(weak_years, contour_f=windspeed['hws'], contour=msl,
                                            quiver=winds, composite_anomailies=composite_anomailies,
                                            calc_std=calc_std)
        comp.to_netcdf(output.outpath)

rule snow_cover_and_temperature:
    input:
        snow_depth=config['intermediate_files'] + '/era5.single_level.snow_depth.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        temperature_2m=config['intermediate_files'] + '/era5.single_level.2m_temperature.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/snow_temperature/composites/era5.temp_2m_snow.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{comp_type}.{threshold}.{weak_strong}.{sdate}-{edate}.nc'
    run:
        from scripts.create_composites import select_years_to_composite, create_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        
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
        
        if wildcards.weak_strong == 'strong':
            comp = create_composite(strong_years, pcolormesh=snow_depth, contour=temperature_2m,
                                                composite_anomailies=composite_anomailies,
                                                calc_std=calc_std)
        else:
            comp =  create_composite(weak_years, pcolormesh=snow_depth, contour=temperature_2m,
                                            composite_anomailies=composite_anomailies,
                                            calc_std=calc_std)
        comp.to_netcdf(output.outpath)