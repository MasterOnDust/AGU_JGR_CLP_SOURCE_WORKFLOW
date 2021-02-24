from scripts.plot_ERA5 import (plot_geopot_temperature,draw_map_orogaphic_east_asia,
                plot_contour,plot_contourf)

FIGS=config['figs_path']

rule plot_orographic_map_mslp:
    input:
        mean_sea_level_pressure=config['intermediate_files']+"/era5.single_level.mean_sea_level_pressure.{season}_avg.{sdate}-{edate}.nc",
    output:
        out_path=FIGS+"/era5.single_level.mslp.{season}_orographic_avg.{sdate}-{edate}.{fig_file_extention}"
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
        script='workflow/scripts/plot_ERA5.py',
        mean_sea_level_pressure_winter=config['intermediate_files']+"/era5.single_level.mean_sea_level_pressure.DJF_avg.{sdate}-{edate}.nc",
        mean_sea_level_pressure_summer=config['intermediate_files']+"/era5.single_level.mean_sea_level_pressure.JJA_avg.{sdate}-{edate}.nc",
        surface_temperature_winter=config['intermediate_files']+"/era5.single_level.2m_temperature.DJF_avg.{sdate}-{edate}.nc",
        surface_temperature_summer=config['intermediate_files']+"/era5.single_level.2m_temperature.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGS+"/era5.single_level.temp2m_mslp.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}"
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
        script='workflow/scripts/plot_ERA5.py',
        temperature_files_winter=config['intermediate_files']+"/era5.{level}.temperature.DJF_avg.{sdate}-{edate}.nc",
        temperature_files_summer=config['intermediate_files']+"/era5.{level}.temperature.JJA_avg.{sdate}-{edate}.nc",
        geopot_files_winter=config['intermediate_files']+"/era5.{level}.GeopotHeight.DJF_avg.{sdate}-{edate}.nc",
        geopot_files_summer=config['intermediate_files']+"/era5.{level}.GeopotHeight.JJA_avg.{sdate}-{edate}.nc"
    output:
        out_path=FIGS+"/era5.{level}.temp_geopot.DJF_JJA_comp_avg.{sdate}-{edate}.{fig_file_extention}",
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