from scripts.process_era5 import calculate_monthly_clm,read_data,resample_seasonal

rule calculate_seasonal_anomaly:
    input:
        season = config['intermediate_files']+"/era5.{level}.{varName}.{season}.{sdate}-{edate}.nc",
        season_avg = config['intermediate_files']+"/era5.{level}.{varName}.{season}_avg.{sdate}-{edate}.nc"
    output:
        outpath=config['intermediate_files']+"/era5.{level}.{varName}.{season}_anomaly.{sdate}-{edate}.nc"
    run:
        ds_clim=read_data(input.season_avg)
        ds_season=read_data(input.season)
        monthly_anomaly=ds_season[ds_season.varName] - ds_clim[ds_clim.varName]
        monthly_anomaly=monthly_anomaly.to_dataset(name=ds_season.varName)
        monthly_anomaly.to_netcdf(output.outpath)

rule calculate_anomaly_monthly:
    input:
        monthly = config['download_path']+"/era5.{level}.{varName}.monthly.{sdate}-{edate}.nc",
        monthly_clim =config['intermediate_files']+"/era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"

    output:
        outpath=config['intermediate_files']+"/era5.{level}.{varName}.monthly_anomaly.{sdate}-{edate}.nc"
    threads: 1    
    run:
        ds_clim=xr.open_dataset(input.monthly_clim)
        ds_monthly=xr.open_dataset(input.monthly)
        monthly_anomaly=ds_monthly.groupby('time.month') - ds_clim
        monthly_anomaly.to_netcdf(output.outpath)
    
        
        
rule calculate_monthly_climatology:
    input:
        src='workflow/scripts/process_era5.py',
        path=config['download_path']+"/era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        outpath=config['intermediate_files']+"/era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"
    run:
        ds= calculate_monthly_clm(input.path)
        ds.to_netcdf(output.outpath)


rule resample_seasonal:
    input:
        path=config['download_path']+"/era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        outpath=config['intermediate_files']+"/era5.{level}.{varName}.{season}.{sdate}-{edate}.nc"
    wildcard_constraints:
        season='DJF|MAM|JJA|SON'
    run:
        ds=resample_seasonal(input.path, wildcards.season)
        ds.to_netcdf(output.outpath)