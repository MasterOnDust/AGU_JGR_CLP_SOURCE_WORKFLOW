
##################################################################################################
# Rules for calculating the NAO index, both using the station based definition and the EOF based #
# definition.                                                                                    #
##################################################################################################

# Import scripts

from scripts.calc_NAO import calc_NAO_EOF, calc_NAO_station

DOWNLOADS=config["download_path"]
DATA_FOLDER=config["intermediate_files"]
rule eof_based_NAO:
    input:
        src='src/calc_NAO.py',
        mslp_data=DOWNLOADS+"/era5.single_level.mean_sea_level_pressure.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"/era5.single_level.NAO_EOF.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        NAO_eof = calc_NAO_EOF(input.mslp_data, wildcards.frequency)
        NAO_eof.to_netcdf(output.outpath)
        
rule calculate_NAO_station:
    input:
        src='src/calc_NAO.py',
        mslp_data_anomaly=DATA_FOLDER+"/era5.single_level.mean_sea_level_pressure.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"/era5.single_level.NAO.station.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        ds_anomaly=xr.open_dataset(input.mslp_data_anomaly)
        ds_NAO=calc_NAO_station(ds_anomaly)
        ds_NAO=ds_NAO.rename(msl='NAO')
        ds_NAO.NAO.attrs['long_name']='North atlantic ocilation index'
        ds_NAO.to_netcdf(output.outpath)