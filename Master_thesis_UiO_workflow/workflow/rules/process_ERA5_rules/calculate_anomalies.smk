from scripts.process_era5 import read_data

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