
from scripts.calc_EAWM import calc_EAWM_zonal_300hPa, calc_MO_index

rule calc_EAWMI_zonal_u300hPa:
    input:
        u_300_data=config['intermediate_files']+"/era5.300hPa.u_component_of_wind.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=config['EAWMI_path']+"/era5.single_level.U300hPa_EAWM.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    threads: 1
    run:
        u_300_ds = xr.open_dataset(input.u_300_data) 
        eawmi=calc_EAWM_zonal_300hPa(u_300_ds)
        eawmi.to_netcdf(output.outpath)

rule calc_MO_mslp_index:
    input:
        mslp_data = config['intermediate_files']+"/era5.single_level.mean_sea_level_pressure.{frequency}.{sdate}-{edate}.nc"
    output:
        outpath = config['EAWMI_path']+"/era5.single_level.EAWM_MO.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        djf_mslp = xr.open_dataset(input.mslp_data)
        MO= calc_MO_index(djf_mslp)
        MO.to_netcdf(output.outpath)