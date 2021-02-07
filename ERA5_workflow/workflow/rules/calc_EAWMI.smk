
from scripts.calc_EAWM import calc_EAWM_zonal_300hPa

rule calc_EAWMI_zonal_u300hPa:
    input:
        u_300_data=config['intermediate_files']+"/era5.300hPa.u_component_of_wind.monthly.{sdate}-{edate}.nc"
    output:
        outpath=config['intermediate_files']+"/era5.single_level.U300hPa_EAWM.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    threads: 1
    run: 
        eawmi=calc_EAWM_zonal_300hPa({input.u_300_data})
        eawmi.to_netcdf(output.outpath)