
##################################################################################################
# Rules for calculating the NAO index, both using the station based definition and the EOF based #
# definition.                                                                                    #
##################################################################################################

# Import scripts

from thesis_toolbox.calc_indicies.calc_NAO import calc_NAO_EOF, calc_NAO_station

DOWNLOADS=config["download_path"]
DATA_FOLDER=config["intermediate_files"]
NAO=config['nao_output_path']
rule eof_based_NAO:
    input:
        mslp_data_anomaly=DATA_FOLDER+"/era5.single_level.mean_sea_level_pressure.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=NAO+"/era5.single_level.NAO_EOF.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        mslp_data_anomaly = xr.open_dataset(input.mslp_data_anomaly)
        NAO_eof = calc_NAO_EOF(mslp_data_anomaly['msl'])

        NAO_eof.attrs['institution']='ECWMF'
        NAO_eof.attrs['source']='ERA5 reanalysis'
        NAO_eof.attrs['title']='EOF based NAO index ' + wildcards.frequency
        NAO_eof.attrs['reference'] = "https://ajdawson.github.io/eofs/latest/examples/nao_xarray.html, https://doi.org/10.1016/j.jmarsys.2008.11.026"
        NAO_eof.to_netcdf(output.outpath)
        
rule calculate_NAO_station:
    input:
        mslp_data_anomaly=DATA_FOLDER+"/era5.single_level.mean_sea_level_pressure.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=NAO+"/era5.single_level.NAO_station.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        ds_anomaly=xr.open_dataset(input.mslp_data_anomaly)
        ds_NAO=calc_NAO_station(ds_anomaly)
        ds_NAO=ds_NAO.rename(msl='NAO')
        ds_NAO.NAO.attrs['long_name']='Normalized sea level pressure (Lisbon - Reykjavik)' 
        ds_NAO.NAO.attrs['standard_name']='Station based NAO index'
        ds_NAO.attrs['title']='Normalized sea level pressure Lisbon - Reykjavik ' + wildcards.frequency
        ds_NAO.attrs['institution']='ECWMF'
        ds_NAO.attrs['source']='ERA5 reanalysis'
        ds_NAO.attrs['comment']='NAO index calculated according to the station based definition'
        ds_NAO.attrs['reference']="10.1126/science.269.5224.676"
        ds_NAO.to_netcdf(output.outpath)