SDATE = config['sdate']
EDATE = config['edate']

rule windspeed_geopot_wind_composite:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        geopot=config['intermediate_files'] + '/era5.{plevel}.GeopotHeight.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/windspeed_geopot_{plevel}/era5.wind_geopot.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{c}_{criterion}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    run:
        from thesis_toolbox.composites.setup_thesis_data import geopot_wind_composite
        from thesis_toolbox.composites.create_composites import select_years_to_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],int(wildcards.c),wildcards.criterion)
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
        # This is kind of a bad hack, I'll make a new rule that remove empty datasets
        if len(weak_years) == 0 or len(strong_years) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:     
            comp = geopot_wind_composite(geopot, wind_u, wind_v, weak_years,strong_years,wildcards.plevel, 
                                    wildcards.season, wildcards.location,wildcards.kind)
        comp.to_netcdf(output.outpath)

rule mean_sea_level_pressure_and_wind:
    input:
        u_wind=config['intermediate_files'] + '/era5.{plevel}.u_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        v_wind=config['intermediate_files'] + '/era5.{plevel}.v_component_of_wind.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        msl=config['intermediate_files'] + '/era5.single_level.mean_sea_level_pressure.{season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc',
        timeseries='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/msl_pressure_wind_{plevel}/era5.wind_msl.{plevel}.composite.{kind}.{psize}.{season}.{location}.{region}.{c}_{criterion}.{sdate}-{edate}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    run:
        from thesis_toolbox.composites.setup_thesis_data import mslp_wind_composite
        from thesis_toolbox.composites.create_composites import select_years_to_composite
        import numpy as np
        timeseries = xr.open_dataset(input.timeseries)
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],int(wildcards.c),wildcards.criterion)
        wind_u = xr.open_dataset(input.u_wind)
        wind_v = xr.open_dataset(input.v_wind)
        windspeed = xr.Dataset()
        winds=xr.Dataset()
        winds['u'] = wind_u['u']
        winds['v'] = wind_v['v']
        # windspeed['hws']=np.sqrt(wind_u['u']**2 + wind_v['v']**2)
        # windspeed.attrs['varName']='hws'
        # windspeed['hws'].attrs['units'] = 'm/s'
        # windspeed['hws'].attrs['long_name'] = 'wind speed'
        msl= xr.open_dataset(input.msl)

        if len(weak_years) == 0 or len(strong_years) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        else:
            comp = mslp_wind_composite(msl, wind_u, wind_v, weak_years, strong_years, wildcards.plevel, 
            wildcards.season, wildcards.location, wildcards.kind)
        comp.to_netcdf(output.outpath)