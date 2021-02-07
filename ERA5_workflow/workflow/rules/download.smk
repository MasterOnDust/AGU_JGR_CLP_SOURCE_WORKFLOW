
# Load config file:

DOWNLOADS=config["download_path"]
VARNAME_PRESSURE_LEVEL="|".join(config['variables_pressure_levels'].keys())
VARNME_SURFACE="|".join(config['variables_single_level'])



rule download_monthly_era5_pressure:
    input: "workflow/scripts/download_era5_monthly.py"
    output: DOWNLOADS + "/era5.{plevel}hPa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DOWNLOADS
    wildcard_constraints:
        varName_ERA5=VARNAME_PRESSURE_LEVEL + '|' + VARNME_SURFACE
    shell: 
        """
        python {input} {wildcards.plevel} {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

rule download_monthly_era5_single_level:
    input: "workflow/scripts/download_era5_monthly.py"
    output: DOWNLOADS + "/era5.single_level.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DOWNLOADS
    wildcard_constraints:
        varName_ERA5=VARNAME_PRESSURE_LEVEL + '|' + VARNME_SURFACE
    shell:
        """
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

