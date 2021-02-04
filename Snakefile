
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
VARNME_SURFACE=['surface_pressure','2m_temperature', 'mean_sea_level_pressure']
TEMPORALRES='monthly'
SDATE=1979
EDATE=2019

SEASONS=['DJF','MAM','JJA','SON']

DATA_FOLDER='/opt/uio/flexpart/ERA5/'

FIGIS=DATA_FOLDER+'figs/'
#----------------------------------CONSTANTS----------------------------------------#

g0=9.80665 #m/s


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
from src.plot_ERA5 import plot_geopot_temperature
from src.process_era5 import calculate_monthly_clm,read_data
from src.calc_NAO import calc_NAO_station,calc_NAO_EOF
from src.calc_EAWMI import calc_EAWM_zonal_300hPa
from IPython import embed


wildcard_constraints:
    season='DJF|MAM|JJA|SON',
    varName_ERA5='|'.join(VARNAME_PRESSURE_LEVEL+VARNME_SURFACE) + '|u_component_of_wind',
    fig_file_extention='png|pdf'

rule all:
    input:
        expand(DATA_FOLDER+"era5.{plevel}hPa.{varName}.monthly.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL,sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.monthly.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.annual.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.{season}_avg.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hPa.{varName}.{season}_avg.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hPa.GeopotHeight.monthly.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, sdate=SDATE, edate=EDATE),  
        expand(DATA_FOLDER+"era5.{plevel}hPa.GeopotHeight.{season}_avg.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.300hPa.u_component_of_wind.monthly.{sdate}-{edate}.nc",
                        sdate=SDATE, edate=EDATE),
        expand(FIGIS+"era5.{plevel}hPa.temp_geopot.DJF_JJA_comp_avg.{sdate}-{edate}.png",
                        sdate=SDATE, edate=EDATE,plevel=PLEVELS),
        expand(FIGIS+"era5.single_level.temp2m_sp.DJF_JJA_comp_avg.{sdate}-{edate}.png",
                        sdate=SDATE, edate=EDATE,plevel=PLEVELS),

        DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.monthly_anomaly.1979-2019.nc",
        DATA_FOLDER+"era5.single_level.NAO.monthly.1979-2019.nc"
rule eof_based_NAO:
    input:
        src='src/calc_NAO.py',
        mslp_data=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.single_level.NAO_EOF.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        NAO_eof = calc_NAO_EOF(input.mslp_data, wildcards.frequency)
        NAO_eof.to_netcdf(output.outpath)
        
rule calculate_NAO_station:
    input:
        src='src/calc_NAO.py',
        mslp_data_anomaly=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.monthly_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.single_level.NAO_station.{frequency}.{sdate}-{edate}.nc"
    run:
        ds_anomaly=xr.open_dataset(input.mslp_data_anomaly)
        ds_NAO=calc_NAO_station(ds_anomaly)
        ds_NAO.to_netcdf(output.outpath)

rule calc_EAWMI_zonal_u300hPa:
    input:
        u_300_data=DATA_FOLDER+"era5.300hPa.u_component_of_wind.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.single_level.U300hPa_EAWM.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    threads: 1
    run: 
        eawmi=calc_EAWM_zonal_300hPa({input.u_300_data})
        eawmi.to_netcdf(output.outpath)
        

    

rule calculate_monthly_anomaly:
    input:
        monthly = DATA_FOLDER+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc",
        monthly_clim =DATA_FOLDER+"era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"

    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.monthly_anomaly.{sdate}-{edate}.nc"
    threads: 1    
    run:
        ds_clim=xr.open_dataset(input.monthly_clim)
        ds_monthly=xr.open_dataset(input.monthly)
        monthly_anomaly=ds_monthly.groupby('time.month') - ds_clim
        monthly_anomaly.attrs['varName']='msl'
        monthly_anomaly.to_netcdf(output.outpath)
        

rule plot_suface_pressure_and_temperature:
    input:
        script='plot_ERA5.py',
        surface_pressure_winter=DATA_FOLDER+"era5.single_level.surface_pressure.DJF_avg.{sdate}-{edate}.nc",
        surface_pressure_summer=DATA_FOLDER+"era5.single_level.surface_pressure.JJA_avg.{sdate}-{edate}.nc",
        surface_temperature_winter=DATA_FOLDER+"era5.single_level.2m_temperature.DJF_avg.{sdate}-{edate}.nc",
        surface_temperature_summer=DATA_FOLDER+"era5.single_level.2m_temperature.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGIS+"era5.single_level.temp2m_sp.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}"
    message: "Creating figure at {output}"
    threads: 1
    run:
        fig,ax = plt.subplots(figsize=(14,6),nrows=2, subplot_kw=dict(projection=ccrs.PlateCarree()))

        ax,im = plot_geopot_temperature(ax,input.surface_pressure_winter, input.surface_pressure_summer,
                            input.surface_temperature_winter, input.surface_temperature_summer,extent=[-120,150,10,90])
        
        fig.suptitle('Averaged surface pressure and temperature {}-{}'.format(
            wildcards.sdate,wildcards.edate),y=0.94, fontsize=18)
        fig.subplots_adjust(bottom=0.15)
        cbar_ax = fig.add_axes([0.17, 0.05, 0.67, 0.05])
        cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', format='%d')
        cbar.set_label('K')

        plt.savefig(output.out_path, bbox_inches='tight')


rule plot_geopotential_and_temperature:
    input:
        script='plot_ERA5.py',
        temperature_files_winter=DATA_FOLDER+"era5.{level}.temperature.DJF_avg.{sdate}-{edate}.nc",
        temperature_files_summer=DATA_FOLDER+"era5.{level}.temperature.JJA_avg.{sdate}-{edate}.nc",
        geopot_files_winter=DATA_FOLDER+"era5.{level}.GeopotHeight.DJF_avg.{sdate}-{edate}.nc",
        geopot_files_summer=DATA_FOLDER+"era5.{level}.GeopotHeight.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGIS+"era5.{level}.temp_geopot.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}",
    message: "Creating figure at {output}"
    threads: 1
    run:
        
        fig,ax = plt.subplots(figsize=(14,6),nrows=2, subplot_kw=dict(projection=ccrs.PlateCarree()))

        ax,im = plot_geopot_temperature(ax, input.temperature_files_winter, 
                            input.temperature_files_summer,
                            input.geopot_files_winter, input.geopot_files_summer,extent=[-120,150,10,90])
        
        fig.suptitle('Averaged geopotential height and temperature {}-{} at {}'.format(
            wildcards.sdate,wildcards.edate,wildcards.level),y=0.94, fontsize=18)
        fig.subplots_adjust(bottom=0.15)
        cbar_ax = fig.add_axes([0.17, 0.05, 0.67, 0.05])
        cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', format='%d')
        cbar.set_label('K')

        plt.savefig(output.out_path, bbox_inches='tight')
        
        
rule calculate_monthly_climatology:
    input:
        src='src/process_era5.py',
        path=DATA_FOLDER+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"
    run:
        ds= calculate_monthly_clm(input.path)
        ds.to_netcdf(output.outpath)

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
    input: DATA_FOLDER + "era5.{plevel}hPa.Geopotential.monthly.{sdate}-{edate}.nc"
    output: DATA_FOLDER + "era5.{plevel}hPa.GeopotHeight.monthly.{sdate}-{edate}.nc"
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
    input: "src/download_era5_monthly.py"
    output: DATA_FOLDER + "era5.{plevel}hPa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DATA_FOLDER
    shell: 
        """
        python {input} {wildcards.plevel} {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

rule download_monthly_era5_single_level:
    input: "src/download_era5_monthly.py"
    output: DATA_FOLDER + "era5.single_level.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=DATA_FOLDER
    shell:
        """
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """
rule setup_folder_structure:
    output:
        data_folder=DATA_FOLDER,
        figs=FIGIS
    shell:
        """
        [[ -d {output.data_folder} ]] || mkdir {output.data_folder}
        [[ -d {output.figs} ]] || mkdir {output.figs}
        """
rule clean:
    params: 
        data_folder=DATA_FOLDER,
        figs=FIGIS
    shell:
        """
        find {params.data_folder} -type f \! -name "*monthly*.nc" -exec rm {{}} \;
        rm {params.data_folder}*GeopotHeight*.nc
        rm -r {params.figs}
        """