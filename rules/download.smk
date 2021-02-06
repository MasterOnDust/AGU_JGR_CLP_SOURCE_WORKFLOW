rule download_monthly_era5_pressure:
    input: "download_era5_monthly.py"
    output: DATA_FOLDER + "era5.{plevel}hPa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DATA_FOLDER
    shell: 
        """
        [[ -d {params.data_dir} ]] || mkdir {params.data_dir}
        python {input} {wildcards.plevel} {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

rule download_monthly_era5_single_level:
    input: "download_era5_monthly.py"
    output: DATA_FOLDER + "era5.single_level.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DATA_FOLDER
    shell:
        """
        [[ -d {params.data_dir} ]] || mkdir {params.data_dir}
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """