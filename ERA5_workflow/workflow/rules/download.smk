
# Load config file:

DOWNLOADS=config["download_path"]
VARNAMES="|".join(config['variables'])

rule download_monthly_era5_pressure:
    input: "workflow/scripts/download_era5_monthly.py"
    output: DOWNLOADS + "/era5.{plevel}hPa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DOWNLOADS
    wildcard_constraints:
        varName_ERA5=VARNAMES
    shell: 
        """
        python {input} {wildcards.plevel} {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        ncks -O --msa -d longitude,181.,360. -d longitude,0.,180.0 {output} {output}
        ncap2 -O -s 'where(longitude > 180) longitude=longitude-360' {output} {output}
        """

rule download_monthly_era5_single_level:
    input: "workflow/scripts/download_era5_monthly.py"
    output: DOWNLOADS + "/era5.single_level.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DOWNLOADS
    wildcard_constraints:
        varName_ERA5=VARNAMES
    shell:
        """
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        ncks -O --msa -d longitude,181.,360. -d longitude,0.,180.0 {output} {output}
        ncap2 -O -s 'where(longitude > 180) longitude=longitude-360' {output} {output}
        """

