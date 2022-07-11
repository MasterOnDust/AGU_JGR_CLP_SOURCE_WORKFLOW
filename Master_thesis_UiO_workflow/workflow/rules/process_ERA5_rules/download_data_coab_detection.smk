
rule download_EA_3hourly:
    output:
        outpath=DOWNLOADS+"/CAOB_data/era5.single_level.{ERA5_vname}.3hourly.{sdate}-{edate}.nc"
    
    message: "downloading {output}"

    notebook: 
        "../../notebooks/download_ERA5/EA_singlelevels_download.ipynb.py"