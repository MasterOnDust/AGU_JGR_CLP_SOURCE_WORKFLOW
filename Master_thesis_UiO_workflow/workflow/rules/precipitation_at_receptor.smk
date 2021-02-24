from IPython import embed

def determine_path(frequency):
    if frequency =='monthly':
        path = config['download_path']
    else:
        path = config['intermediate_files']
    return path

rule precipitation_at_receptor:
    input:
        precip_data=config['intermediate_files']+"/era5.single_level.total_precipitation.{frequency}.{sdate}-{edate}.nc"
    output:
        outpath=config['Precipitation'] + "/era5.single_level.total_precipitation.{receptor}.{frequency}.{sdate}-{edate}.csv" 
    wildcard_constraints:
        frequency='DJF|MAM|JJA|SON',
        receptor='|'.join(config['receptors'].keys())
    run:
        ds = xr.open_dataset(input.precip_data)
        location = config['receptors'][wildcards.receptor]
        ds = ds.sel(longitude=location['lon'], latitude=location['lat'], method='nearest')
        ds.attrs['title'] = 'Total precipitation at ' + wildcards.receptor
        ds.attrs['institution']='ECWMF'
        ds.attrs['source']='ERA5 reanalysis'
        ds.to_dataframe().to_csv(output.outpath)

rule precipitation_at_receptor_monthly:
    input:
        precip_data=config['download_path'] +"/era5.single_level.total_precipitation.{frequency}.{sdate}-{edate}.nc"
    output:
        outpath=config['Precipitation'] + "/era5.single_level.total_precipitation.{receptor}.{frequency}.{sdate}-{edate}.csv" 
    wildcard_constraints:
        frequency='monthly',
        receptor='|'.join(config['receptors'].keys())
    run:
        ds = xr.open_dataset(input.precip_data)
        location = config['receptors'][wildcards.receptor]
        ds = ds.sel(longitude=location['lon'], latitude=location['lat'], method='nearest')
        ds.attrs['title'] = 'Total precipitation at ' + wildcards.receptor
        ds.attrs['institution']='ECWMF'
        ds.attrs['source']='ERA5 reanalysis'
        ds.to_dataframe().to_csv(output.outpath)
        
rule precipitation_in_source_region_series:
    input:
        precip_data=lambda wildcards: determine_path(wildcards.frequency) +"/era5.single_level.total_precipitation.{frequency}.{sdate}-{edate}.nc"
    output:
        outpath=config['Precipitation'] + "/era5.single_level.total_precipitation.{source_region}.{frequency}.{sdate}-{edate}.csv" 
    
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON',
        source_region='|'.join(config['source_regions'].keys())
    
    run:
        ds=xr.open_dataset(input.precip_data)
        source_region=config['source_regions'][wildcards.source_region]
        ds=ds.sel(longitude=slice(source_region['lon0'], source_region['lon1']),
                 latitude=slice(source_region['lat1'], source_region['lat0']))
        ds=ds.mean(dim=['longitude','latitude'], keep_attrs=True)
        ds.attrs['title'] = 'Averaged precipitation over ' + wildcards.source_region
        ds.attrs['institution']='ECWMF'
        ds.attrs['source']='ERA5 reanalysis'
        ds.to_dataframe().to_csv(output.outpath)
 
rule calculate_precipitation_anomalies:
    input: 
        precip_data=config['Precipitation'] + "/era5.single_level.total.precipitation.{area}.{frequency}.{sdate}-{edate}.csv"
    output:
        outpath=config['Precipitation']+"/era5.single_level.total.precipitation.{area}.{frequency}_anomalies.{sdate}-{edate}.csv"
    wildcard_constraints:
        frequency='DJF|MAM|JJA|SON'
    run:
        import pandas as pd
        df = pd.read_csv(input.precip_data,index_col=0)
        mean_precip=df['tp'].mean()
        anomalies = df['tp']-mean_precip
        anomalies.to_csv(output.outpath)
        