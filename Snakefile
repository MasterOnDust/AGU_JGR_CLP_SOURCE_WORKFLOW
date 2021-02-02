
#####################################################################################
# Snakefile for downloading and processing monthly ERA5 data                        #
#                                                                                   #
# Created by Ove Haugvaldstad                                                       #
#                                                                                   #
#####################################################################################

#-------------------------------CONFIGURATION---------------------------------------#

PLEVELS=[1000,850, 500, 250]
VARNAME_PRESSURE_LEVEL=['temperature','Geopotential'] # Name of variables to download pressure levels 
#VARNME_SURFACE=['2m_temperature', 'surface_pressure', 'snow_depth']
VARNME_SURFACE=['surface_pressure']
TEMPORALRES='monthly'
SDATE=2017
EDATE=2019

SEASONS=['DJF','MAM','JJA','SON']

DATA_FOLDER='/opt/uio/flexpart/ERA5/'

FIGIS=DATA_FOLDER+'/figs/'
#----------------------------------CONSTANTS----------------------------------------#

g0=9.80665 #m/s


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
from plot_ERA5 import plot_contour, plot_contourf,read_data, draw_map



wildcard_constraints:
    season='DJF|MAM|JJA|SON',
    varName_ERA5='|'.join(VARNAME_PRESSURE_LEVEL+VARNME_SURFACE)

rule all:
    input: 
        expand(DATA_FOLDER+"era5.{plevel}hpa.{varName}.monthly.{sdate}-{edate}.nc", plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL,sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.monthly.{sdate}-{edate}.nc", varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.annual.{sdate}-{edate}.nc", varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.{season}_avg.{sdate}-{edate}.nc", varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hpa.{varName}.{season}_avg.{sdate}-{edate}.nc", plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL,sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hpa.GeopotHeight.monthly.{sdate}-{edate}.nc", plevel=PLEVELS, sdate=SDATE, edate=EDATE),  
        expand(DATA_FOLDER+"era5.{plevel}hpa.GeopotHeight.{season}_avg.{sdate}-{edate}.nc", plevel=PLEVELS, sdate=SDATE, edate=EDATE, season=SEASONS)           
rule plot:
    input:
        temperature_files=DATA_FOLDER+"era5.{level}.temperature.{season}_avg.{sdate}-{edate}.nc",
        geopot_files=DATA_FOLDER+"era5.{level}.GeopotHeight.{season}_avg.{sdate}-{edate}.nc"
    output:
        FIGIS+"era5.{level}.temp_geopot.{season}_avg.{sdate}-{edate}.png",

    run:
        """
        ds_t=read_data({input.temperature_files})
        ds_geopot=read_data({input.geopot_files})
        
        """
        

rule create_anomaly:
    input:
        monthly = DATA_FOLDER+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc",
        annual_avg =DATA_FOLDER+"era5.{level}.{varName}.avg.{sdate}-{edate}.nc"

    output:
        DATA_FOLDER+"era5.{level}.{varName}.anomaly.{sdate}-{edate}.nc"
    threads: 1    
    shell:
        """
        ncbo -d time,,11 {input.monthly} {input.annual_avg} {output}
        """


rule annual_average:
    input:
        DATA_FOLDER+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        DATA_FOLDER+"era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    threads: 1
    shell:
        """
        ncra --mro -d time,,,12,12 {input} {output}
        """

rule seasonal_average:
    input:
        DATA_FOLDER+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        DATA_FOLDER+"era5.{level}.{varName}.{season}_avg.{sdate}-{edate}.nc"
    threads:1
    shell:
        """
        if [[ {wildcards.season} == DJF ]]; then
            ncra -d time,11,,12,3 {input} {output}
        elif [[ {wildcards.season} == MAM ]]; then
            ncra -d time,2,,12,3 {input} {output}
        elif [[ {wildcards.season} == JJA ]]; then
            ncra -d time,5,,12,3 {input} {output}
        else
            ncra -d time,8,,12,3 {input} {output}
        fi
        """


rule whole_period_average:
    input:
        DATA_FOLDER+"era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    output:
        DATA_FOLDER+"era5.{level}.{varName}.avg.{sdate}-{edate}.nc"

rule geopotential_to_geopotential_height:
    input: DATA_FOLDER + "era5.{plevel}hpa.Geopotential.monthly.{sdate}-{edate}.nc"
    output: DATA_FOLDER + "era5.{plevel}hpa.GeopotHeight.monthly.{sdate}-{edate}.nc"
    params: g=g0
    shell:
        """
        ncap2 -v -s'z=z/{params.g}f*0.1f;' {input} {output} 
        ncatted -a units,z,o,c,'dam' {output}
        ncatted -a long_name,z,o,c,'Geopotential height' {output}  
        ncatted -a standard_name,z,o,c,'geopotential height' {output}
        ncrename -v z,Z {output}
        """




rule download_monthly_era5_pressure:
    input: "download_era5_monthly.py"
    output: DATA_FOLDER + "era5.{plevel}hpa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DATA_FOLDER
    shell: 
        """
        [[ -d dir ]] || mkdir dir {params.data_dir}
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
        [[ -d dir ]] || mkdir dir {params.data_dir}
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

rule clean:
    params: data_folder=DATA_FOLDER
    shell:
        """
        find {params.data_folder} -type f \! -name "*monthly*.nc" -exec rm {{}} \;
        rm {params.data_folder}*GeopotHeight*.nc
        """