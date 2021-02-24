
rule APVI_index:
    input:
        geo_pot_anomaly=config['intermediate_files']+"/era5.500hPa.GeopotHeight.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=config['APVI_index_outpath']+"/era5.500hPa.APVI_index.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    
    run:
        geo_pot_anomaly=xr.open_dataset(input.geo_pot_anomaly)
        bolean=((geo_pot_anomaly.longitude >= 60) & (geo_pot_anomaly.longitude <= 150) & 
                (geo_pot_anomaly.latitude <= 90) & (geo_pot_anomaly.latitude >= 45))
        
        APVI_index=geo_pot_anomaly.where(bolean,drop=True).mean(dim=['longitude','latitude'])
        APVI_index=APVI_index.rename(Z='APVI')
        APVI_index.attrs['title']='Asian polar vortex intensity index'
        APVI_index['APVI'].attrs['long_name']='Asian polar vortex intensity index'
        APVI_index.attrs['institution']='ECWMF'
        APVI_index.attrs['source']='ERA5 reanalysis'
        APVI_index.attrs['reference'] = "10.1007/s00382-020-05507-9"
        APVI_index.to_netcdf(output.outpath)
           