from scripts.calc_AO import calc_AO_EOF


DOWNLOADS=config["download_path"]
DATA_FOLDER=config["intermediate_files"]
AO=config['ao_output_path']
rule eof_based_AO:
    input:
        src='workflow/scripts/calc_AO.py',
        geoPot_data_anomaly=DATA_FOLDER+"/era5.1000hPa.GeopotHeight.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=AO+"/era5.1000hPa.AO_EOF.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        geoPot_data_anomaly = xr.open_dataset(input.geoPot_data_anomaly)
        AO_eof = calc_AO_EOF(geoPot_data_anomaly['Z'])

        AO_eof.attrs['institution']='ECWMF'
        AO_eof.attrs['source']='ERA5 reanalysis'
        AO_eof.attrs['title']='EOF based AO index ' + wildcards.frequency
        AO_eof.attrs['reference'] = "https://ajdawson.github.io/eofs/latest/examples/nao_xarray.html, https://doi.org/10.1175/1520-0442(2000)013<1000:AMITEC>2.0.CO;2"
        AO_eof.to_netcdf(output.outpath)