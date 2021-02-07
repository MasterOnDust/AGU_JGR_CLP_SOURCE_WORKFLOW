
#####################################################################################
# Snakefile for downloading and processing monthly ERA5 data                        #
#                                                                                   #
# Created by Ove Haugvaldstad                                                       #
#                                                                                   #
#####################################################################################

# TODO:
#   MAKE NEW folders so I don't get trouble with Ambiguous rules 
#-------------------------------CONFIGURATION---------------------------------------#

PLEVELS=[850, 500, 250]
VARNAME_PRESSURE_LEVEL=['temperature','Geopotential'] # Name of variables to download pressure levels 
#VARNME_SURFACE=['2m_temperature', 'surface_pressure', 'snow_depth']
VARNME_SURFACE=['2m_temperature', 'mean_sea_level_pressure', 'snow_depth']                                                                                                           
TEMPORALRES='monthly'
SDATE=1979
EDATE=2019

SEASONS=['DJF','MAM','JJA','SON']

DATA_FOLDER='/opt/uio/flexpart/ERA5/'

MONTHLY_DATA=DATA_FOLDER+'monthly_data/'

FIGS=DATA_FOLDER+'figs/'
#----------------------------------CONSTANTS----------------------------------------#

g0=9.80665 #m/s


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
from src.plot_ERA5 import (plot_geopot_temperature,draw_map_orogaphic_east_asia,
                plot_contour,plot_contourf)
from src.process_era5 import calculate_monthly_clm,read_data,resample_seasonal

from src.calc_NAO import calc_NAO_station,calc_NAO_EOF
from src.calc_EAWM import calc_EAWM_zonal_300hPa
from IPython import embed


wildcard_constraints:
    season='DJF|MAM|JJA|SON',
    varName_ERA5='|'.join(VARNAME_PRESSURE_LEVEL+VARNME_SURFACE) + '|u_component_of_wind',
    fig_file_extention='png|pdf'

rule all:
    input:
        expand(MONTHLY_DATA+"era5.{plevel}hPa.{varName}.monthly.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL,sdate=SDATE, edate=EDATE),
        expand(MONTHLY_DATA+"era5.single_level.{varName}.monthly.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.annual.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE),
        expand(DATA_FOLDER+"era5.single_level.{varName}.{season}_avg.{sdate}-{edate}.nc", 
                        varName=VARNME_SURFACE, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hPa.{varName}.{season}_avg.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, varName=VARNAME_PRESSURE_LEVEL, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(DATA_FOLDER+"era5.{plevel}hPa.GeopotHeight.{season}_avg.{sdate}-{edate}.nc", 
                        plevel=PLEVELS, sdate=SDATE, edate=EDATE, season=SEASONS),
        expand(MONTHLY_DATA+"era5.300hPa.u_component_of_wind.monthly.{sdate}-{edate}.nc",
                        sdate=SDATE, edate=EDATE),
        expand(FIGS+"era5.{plevel}hPa.temp_geopot.DJF_JJA_comp_avg.{sdate}-{edate}.png",
                        sdate=SDATE, edate=EDATE,plevel=PLEVELS),
        expand(FIGS+"era5.single_level.temp2m_mslp.DJF_JJA_comp_avg.{sdate}-{edate}.png",
                        sdate=SDATE, edate=EDATE,plevel=PLEVELS),

        DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.monthly_anomaly.1979-2019.nc",
        DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.DJF.1979-2019.nc",
        DATA_FOLDER+"era5.single_level.NAO.station.DJF.1979-2019.nc"
        # FIGS+"era5.single_level.mslp.DJF_orographic_avg.1979-2019.png"
        #DATA_FOLDER+"era5.single_level.NAO.monthly.1979-2019.nc"
rule eof_based_NAO:
    input:
        src='src/calc_NAO.py',
        mslp_data=MONTHLY_DATA+"era5.single_level.mean_sea_level_pressure.monthly.{sdate}-{edate}.nc"
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
        mslp_data_anomaly=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.{frequency}_anomaly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.single_level.NAO.station.{frequency}.{sdate}-{edate}.nc"
    wildcard_constraints:
        frequency='monthly|DJF|MAM|JJA|SON'
    run:
        ds_anomaly=xr.open_dataset(input.mslp_data_anomaly)
        ds_NAO=calc_NAO_station(ds_anomaly)
        ds_NAO=ds_NAO.rename(msl='NAO')
        ds_NAO.NAO.attrs['long_name']='North atlantic ocilation index'
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
        

rule calculate_seasonal_anomaly:
    input:
        season = DATA_FOLDER+"era5.{level}.{varName}.{season}.{sdate}-{edate}.nc",
        season_avg = DATA_FOLDER+"era5.{level}.{varName}.{season}_avg.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.{season}_anomaly.{sdate}-{edate}.nc"
    run:
        ds_clim=read_data(input.season_avg)
        ds_season=read_data(input.season)
        monthly_anomaly=ds_season[ds_season.varName] - ds_clim[ds_clim.varName]
        monthly_anomaly=monthly_anomaly.to_dataset(name=ds_season.varName)
        monthly_anomaly.to_netcdf(output.outpath)

rule calculate_anomaly_monthly:
    input:
        monthly = MONTHLY_DATA+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc",
        monthly_clim =DATA_FOLDER+"era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"

    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.monthly_anomaly.{sdate}-{edate}.nc"
    threads: 1    
    run:
        ds_clim=xr.open_dataset(input.monthly_clim)
        ds_monthly=xr.open_dataset(input.monthly)
        monthly_anomaly=ds_monthly.groupby('time.month') - ds_clim
        monthly_anomaly.to_netcdf(output.outpath)
        
rule plot_orographic_map_mslp:
    input:
        mean_sea_level_pressure=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.{season}_avg.{sdate}-{edate}.nc",
    output:
        out_path=FIGS+"era5.single_level.mslp.{season}_orographic_avg.{sdate}-{edate}.{fig_file_extention}"
    threads: 1
    run:
        fig,ax =plt.subplots(figsize=(10,10),
             subplot_kw={'projection':ccrs.Orthographic(central_latitude=45,central_longitude=70)})
        ax.coastlines('110m')
        gl=ax.gridlines(transform = ccrs.PlateCarree(), linestyle ='--', draw_labels=False,linewidth=2)
        ds_mslp=read_data(input.mean_sea_level_pressure)
        da_mslp=ds_mslp[ds_mslp.varName].isel(time=0)/100
        da_mslp.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),
                levels=[990,995,1000,1005,1010,1015,1020,1025, 1030],extend='both')
        fig.suptitle('Averaged {} mean sea level pressure {}-{}'.format(
            wildcards.season,wildcards.sdate,wildcards.edate),y=0.94, fontsize=18)
        if wildcards.fig_file_extention=='png':
            plt.savefig(output.out_path, bbox_inches='tight')
        else:
            plt.savefig(output.out_path, bbox_inches='tight')

rule plot_mslp_and_temperature:
    input:
        script='src/plot_ERA5.py',
        mean_sea_level_pressure_winter=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.DJF_avg.{sdate}-{edate}.nc",
        mean_sea_level_pressure_summer=DATA_FOLDER+"era5.single_level.mean_sea_level_pressure.JJA_avg.{sdate}-{edate}.nc",
        surface_temperature_winter=DATA_FOLDER+"era5.single_level.2m_temperature.DJF_avg.{sdate}-{edate}.nc",
        surface_temperature_summer=DATA_FOLDER+"era5.single_level.2m_temperature.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGS+"era5.single_level.temp2m_mslp.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}"
    message: "Creating figure at {output}"
    threads: 1
    run:
        fig,ax = plt.subplots(figsize=(14,6),nrows=2, subplot_kw=dict(projection=ccrs.PlateCarree()))

        ax,im = plot_geopot_temperature(ax,input.mean_sea_level_pressure_winter, input.mean_sea_level_pressure_summer,
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
        script='src/plot_ERA5.py',
        temperature_files_winter=DATA_FOLDER+"era5.{level}.temperature.DJF_avg.{sdate}-{edate}.nc",
        temperature_files_summer=DATA_FOLDER+"era5.{level}.temperature.JJA_avg.{sdate}-{edate}.nc",
        geopot_files_winter=DATA_FOLDER+"era5.{level}.GeopotHeight.DJF_avg.{sdate}-{edate}.nc",
        geopot_files_summer=DATA_FOLDER+"era5.{level}.GeopotHeight.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGS+"era5.{level}.temp_geopot.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}",
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
        path=MONTHLY_DATA+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.monthly_clim.{sdate}-{edate}.nc"
    run:
        ds= calculate_monthly_clm(input.path)
        ds.to_netcdf(output.outpath)

rule annual_average:
    input:
        MONTHLY_DATA+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        DATA_FOLDER+"era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    threads: 1
    shell:
        """
        ncra --mro -d time,,,12,12 {input} {output}
        """

rule seasonal_average:
    input:
        path=DATA_FOLDER+"era5.{level}.{varName}.{season}.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.{season}_avg.{sdate}-{edate}.nc"
    threads:1
    run:
        ds=xr.open_dataset(input.path)
        ds=ds.mean(dim='time',keep_attrs=True)
        ds.to_netcdf(output.outpath)


rule whole_period_average:
    input:
        DATA_FOLDER+"era5.{level}.{varName}.annual.{sdate}-{edate}.nc"
    output:
        DATA_FOLDER+"era5.{level}.{varName}.avg.{sdate}-{edate}.nc"

rule resample_seasonal:
    input:
        path=MONTHLY_DATA+"era5.{level}.{varName}.monthly.{sdate}-{edate}.nc"
    output:
        outpath=DATA_FOLDER+"era5.{level}.{varName}.{season}.{sdate}-{edate}.nc"
    wildcard_constraints:
        season='DJF|MAM|JJA|SON'
    run:
        ds=resample_seasonal(input.path, wildcards.season)
        ds.to_netcdf(output.outpath)
        


rule geopotential_to_geopotential_height:
    input: MONTHLY_DATA + "era5.{plevel}hPa.Geopotential.monthly.{sdate}-{edate}.nc"
    output: MONTHLY_DATA + "era5.{plevel}hPa.GeopotHeight.monthly.{sdate}-{edate}.nc"
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
    output: MONTHLY_DATA + "era5.{plevel}hPa.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=MONTHLY_DATA
    shell: 
        """
        python {input} {wildcards.plevel} {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """

rule download_monthly_era5_single_level:
    input: "src/download_era5_monthly.py"
    output: MONTHLY_DATA + "era5.single_level.{varName_ERA5}.monthly.{sdate}-{edate}.nc"
    message: "downloading {output}"
    threads: 1
    params: data_dir=MONTHLY_DATA
    shell:
        """
        python {input} -1 {wildcards.varName_ERA5} {wildcards.sdate} {wildcards.edate} --fn {output}
        ncks -O --mk_rec_dim time {output} {output}
        """
rule setup_folder_structure:
    output:
        data_folder=DATA_FOLDER,
        figs=FIGS,
        monthly_data=MONTHLY_DATA
    shell:
        """
        [[ -d {output.data_folder} ]] || mkdir {output.data_folder}
        [[ -d {output.monthly_data} ]] || mkdir {output.monthly_data}
        [[ -d {output.figs} ]] || mkdir {output.figs}
        """
rule clean:
    params: 
        data_folder=DATA_FOLDER,
        figs=FIGS
    shell:
        """
        find {params.data_folder} -type f \! -name "*monthly*.nc" -exec rm {{}} \;
        rm {params.data_folder}*GeopotHeight*.nc
        rm -r {params.figs}
        """