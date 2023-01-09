SDATE = config['sdate']
EDATE = config['edate']
import glob
def get_MERRA2_paths(wildcards):

    paths = []
    for i in range(SDATE, EDATE+1):
        paths.append(glob.glob(config['MERRA2_path']+'/dust/*{}*.nc4'.format(i)))
    return [item for sublist in paths for item in sublist]

def get_forcing_paths(wildcards):
    outpaths = {}
    if wildcards.forcing=='merra2':
        paths = []
        for i in range(SDATE, EDATE+1):
            paths.append(glob.glob(config['MERRA2_path']+'/meteo/*{}*.nc4'.format(i)))
        outpaths['paths'] = [item for sublist in paths for item in sublist]
    else:
        outpaths['u_wind']=config['intermediate_files'] + f'/era5.{wildcards.plevel}.u_component_of_wind.{wildcards.season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc'
        outpaths['v_wind']=config['intermediate_files'] + f'/era5.{wildcards.plevel}.v_component_of_wind.{wildcards.season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc'
        outpaths['msl']=config['intermediate_files'] + f'/era5.single_level.mean_sea_level_pressure.{wildcards.season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc'
        outpaths['geopot']=config['intermediate_files'] + f'/era5.{wildcards.plevel}.GeopotHeight.{wildcards.season}.'+ str(SDATE)+'-'+str(EDATE)+'.nc'
    
    return outpaths
rule windspeed_geopot_wind_composite_merra:
    input:
        unpack(get_forcing_paths),
        MERRA2_path=get_MERRA2_paths
    output:
        outpath='results/composites_MERRA2/windspeed_geopot_{plevel}/{forcing}.wind_geopot.{plevel}.composite.{kind}.{size_bin}.{season}.{location}.{c}_{criterion}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    wildcard_constraints:
        size_bin='clay|fine-silt',
        forcing='era5|merra2'
    run:
        from thesis_toolbox.composites.setup_thesis_data import geopot_wind_composite
        from thesis_toolbox.composites.setup_MERRA2_composite import merra2_geopot_wind_composite 
        from thesis_toolbox.composites.create_composites import select_years_to_composite
        import numpy as np
        from thesis_toolbox.utils import get_locations_CLP, create_deposition_timeseries_MERRA2
        import pandas as pd
        import xarray as xr
        locs = get_locations_CLP()
        if wildcards.location == 'BAODE':
            lon,lat = locs.loc['BADOE']
        else:
            lon,lat = locs.loc[wildcards.location]
        paths = sorted(input.MERRA2_path)
        ds = xr.open_mfdataset(paths,concat_dim=['time'], combine='nested')
        if wildcards.size_bin=='clay':
            sizebin = ['002']
        elif wildcards.size_bin=='fine-silt':
            sizebin = ['005']
        else:
            sizebin = None

        timeseries = create_deposition_timeseries_MERRA2(ds, lon,lat,kind=wildcards.kind, sizebins=sizebin)
        timeseries=timeseries.resample(time='Q-NOV').sum(keep_attrs=True)
        timeseries=timeseries.sel(time=(timeseries.time.dt.season == 'MAM'))
        timeseries = timeseries.assign_coords(time=timeseries.time.dt.year)
        timeseries = timeseries.compute()
        weak_years, strong_years = select_years_to_composite(timeseries,int(wildcards.c),wildcards.criterion)

        # This is kind of a bad hack, I'll make a new rule that remove empty datasets
        if len(weak_years) == 0 or len(strong_years) == 0:
            
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        elif wildcards.forcing == 'era5':     
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
            comp = geopot_wind_composite(geopot, wind_u, wind_v, weak_years,strong_years,wildcards.plevel, 
                                    wildcards.season, wildcards.location,wildcards.kind)
        else:
            paths = sorted(input.paths)
            ds = xr.open_mfdataset(paths,
                       concat_dim=['time'],parallel=True, combine='nested')
                
            ds=ds.resample(time='Q-NOV').mean()
            ds=ds.sel(time=(ds.time.dt.season == wildcards.season))
            plevel = float(wildcards.plevel[:-3])
            comp = merra2_geopot_wind_composite(ds, weak_years, strong_years,plevel, wildcards.season, wildcards.location, wildcards.kind)
        comp.to_netcdf(output.outpath)

rule mean_sea_level_pressure_and_wind_merra:
    input:
        unpack(get_forcing_paths),
        MERRA2_path=get_MERRA2_paths
    output:
        outpath='results/composites_MERRA2/msl_pressure_wind_{plevel}/{forcing}.wind_msl.{plevel}.composite.{kind}.{size_bin}.{season}.{location}.{c}_{criterion}.nc'
    params:
        clim_sdate = config['clim_date0'],
        clim_edate = config['clim_date1']
    threads: 1
    wildcard_constraints:
        size_bin='clay|fine-silt',
        forcing='era5|merra2'
    run:
        from thesis_toolbox.composites.setup_MERRA2_composite import merra2_mslp_wind_composite 
        from thesis_toolbox.composites.setup_thesis_data import mslp_wind_composite
        from thesis_toolbox.composites.create_composites import select_years_to_composite
        import numpy as np
        from thesis_toolbox.utils import get_locations_CLP,create_deposition_timeseries_MERRA2
        import pandas as pd
        locs = get_locations_CLP()
        if wildcards.location == 'BAODE':
            lon,lat = locs.loc['BADOE']
        else:
            lon,lat = locs.loc[wildcards.location]
        paths = sorted(input.MERRA2_path)
        ds = xr.open_mfdataset(paths,concat_dim=['time'], combine='nested')
        if wildcards.size_bin=='clay':
            sizebin = ['002']
        elif wildcards.size_bin=='fine-silt':
            sizebin = ['005']
        else:
            sizebin = None

        timeseries = create_deposition_timeseries_MERRA2(ds, lon,lat,kind=wildcards.kind, sizebins=sizebin)
        timeseries=timeseries.resample(time='Q-NOV').sum(keep_attrs=True)
        timeseries=timeseries.sel(time=(timeseries.time.dt.season == 'MAM'))
        timeseries = timeseries.assign_coords(time=timeseries.time.dt.year)
        timeseries = timeseries.compute()
        weak_years, strong_years = select_years_to_composite(timeseries,int(wildcards.c),wildcards.criterion)


        if len(weak_years) == 0 or len(strong_years) == 0:
            comp = xr.Dataset()
            comp.attrs['years_composited']=0
        
        elif wildcards.forcing == 'era5': 
            wind_u = xr.open_dataset(input.u_wind)
            wind_v = xr.open_dataset(input.v_wind)
            windspeed = xr.Dataset()
            winds=xr.Dataset()
            winds['u'] = wind_u['u']
            winds['v'] = wind_v['v']
            msl= xr.open_dataset(input.msl)
            comp = mslp_wind_composite(msl, wind_u, wind_v, weak_years, strong_years, wildcards.plevel, 
            wildcards.season, wildcards.location, wildcards.kind)
        else:
            paths = sorted(input.paths)
            ds = xr.open_mfdataset(paths,
                       concat_dim=['time'],parallel=True, combine='nested')
            ds=ds.resample(time='Q-NOV').mean()
            ds=ds.sel(time=(ds.time.dt.season == wildcards.season))
            plevel = float(wildcards.plevel[:-3])
            
            comp=merra2_mslp_wind_composite(ds, weak_years, strong_years, plevel, wildcards.season, wildcards.location, wildcards.kind)
        comp.to_netcdf(output.outpath)

rule calc_composites_MERRA2:
    input:
        expand('results/composites_MERRA2/msl_pressure_wind_{plevel}/{forcing}.wind_msl.{plevel}.composite.{kind}.{size_bin}.{season}.{location}.{c}_{criterion}.nc'
                ,plevel=['850hPa'], kind=['total','drydep'],size_bin=['clay'], season=['DJF', 'MAM'], location=['SACOL','LINGTAI','BAODE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],
                forcing=['merra2'],c=['4'], criterion=['rank']),
        expand('results/composites_MERRA2/windspeed_geopot_{plevel}/{forcing}.wind_geopot.{plevel}.composite.{kind}.{size_bin}.{season}.{location}.{c}_{criterion}.nc',
                plevel=['500hPa'], kind=['total','drydep'], size_bin=['clay'], season=['DJF','MAM'], location=['SACOL','LINGTAI','BAODE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],
                forcing=['merra2'],c=['4'], criterion=['rank'])
    threads: 1