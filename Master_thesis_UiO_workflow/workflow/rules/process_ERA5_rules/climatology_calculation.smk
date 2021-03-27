from thesis_toolbox.process_era5 import calculate_monthly_clm,read_data,resample_seasonal


        
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
        import xarray as xr
        ds = xr.open_dataset(input.path)
        ds=resample_seasonal(ds, wildcards.season)
        ds.to_netcdf(output.outpath)