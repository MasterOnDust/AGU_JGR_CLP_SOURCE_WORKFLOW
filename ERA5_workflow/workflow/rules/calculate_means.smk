
rule annual_average:
    input:
        config['download_path']+"/era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        config['intermediate_files']+"/era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    threads: 1
    shell:
        """
        ncra --mro -d time,,,12,12 {input} {output}
        """

rule seasonal_average:
    input:
        path=config['intermediate_files']+"/era5.{level}.{varName}.{season}.{sdate}-{edate}.nc"
    output:
        outpath=config['intermediate_files']+"/era5.{level}.{varName}.{season}_avg.{sdate}-{edate}.nc"
    threads:1
    run:
        ds=xr.open_dataset(input.path)
        ds=ds.mean(dim='time',keep_attrs=True)
        ds.to_netcdf(output.outpath)


rule whole_period_average:
    input:
        config['intermediate_files']+"/era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    output:
        config['intermediate_files']+"/era5.{level}.{varName}.avg.{sdate}-{edate}.nc"