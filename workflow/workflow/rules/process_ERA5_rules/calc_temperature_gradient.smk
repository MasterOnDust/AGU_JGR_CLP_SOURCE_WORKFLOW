
rule temperature_gradient:
    input:
        temp_path = config['intermediate_files']+'/era5.single_level.2m_temperature.{season}.{sdate}-{edate}.nc'
    output:
        outpath = 'results/2mt_gradient/era5.single_level.2m_temperature_gradient.East_asia.{season}.{sdate}-{edate}.nc'
    params:
        n_bins=config['n_bins_temperature_gradient']
    threads: 1
    run:
        from thesis_toolbox.calc_indicies.calculate_temperature_gradient import calculate_temperature_gradient
        from thesis_toolbox.process_era5 import read_data
        temp_ds = read_data(input.temp_path)
        temp_ds = calculate_temperature_gradient(temp_ds,n_latitude_bins=params.n_bins)
        temp_ds.attrs['title']='Temperature gradient East Asia'
        temp_ds.attrs['source']='ERA5 reanalysis'
        temp_ds.to_netcdf(output.outpath)
    
