
rule plot_geopot_200hPa_windspeed_composites:
    input:
        data_path = 'results/composites/windspeed_geopot/era5.wind_geopot.200hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/figs/windspeed_geopot_200hPa/era5.wind_geopot.200hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.{extension}'
    wildcard_constraints:
        extension = 'png|pdf'
    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_large_scale
        import xarray as xr
        import numpy as np
        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            ds = ds.sel(longitude=slice(0,180), latitude=slice(90,0))
            ds = ds.squeeze()
            fig = plt.figure(figsize=(12,12.3))
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.03, hspace=0.05)
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
    
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_large_scale(ax1)
            im = ds.hws_clim.where(ds.hws_clim >= 20, drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=np.arange(20,70,5),
                                                                    cmap='viridis_r', extend='max', add_colorbar=False,ax =ax1)
            CS = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(0), ax=ax1,colors='black', linewidths=.7, levels=26, 
                                        add_labels=False, alpha=.7, vmax=1245,vmin=1120)
            CS1 = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(0), ax=ax1,colors='black', 
                                         linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax1.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1, zorder=1030)
            ax1.set_xlabel('')
            ax1.set_ylabel('')

            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Wind speed [m/s]',shrink=0.9)
            # # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            im = ds.hws.plot.contourf(transform=ccrs.PlateCarree(),
                                                                cmap='bwr', add_colorbar=False, ax=ax2, levels=np.linspace(-9,9, 16))
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)


            # # # Anomalies
            # ax3 = fig.add_subplot(spec[2,0], projection=ccrs.PlateCarree(central_longitude))
            cax2 = fig.add_subplot(spec[1,1])
            fig.colorbar(im, cax=cax2, label='Composite difference 200hPa winds [m/s] \n  (strong years - weak years)',
                         shrink=0.9, ticks=im.levels)

            CS = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', linewidths=.7, 
                                   add_labels=False, alpha=.7, vmin=-6, vmax=6, levels=26)
            CS1 = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black',
                                    linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax2.clabel(CS1, fmt='%.2f', colors='black', fontsize=12, inline=1)
            if config['add_fig_title']:
                fig.suptitle('200hPa geopotential height and windspeed at 200hPa'.format(wildcards.season), 
                             y=.915,size=16)
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()



rule plot_geopot_500hPa_windspeed_composites:
    input:
        data_path = 'results/composites/windspeed_geopot/era5.wind_geopot.500hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/figs/windspeed_geopot/era5.wind_geopot.500hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.{extension}'
    wildcard_constraints:
        extension = 'png|pdf'
    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_large_scale
        import xarray as xr
        import numpy as np
        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            lon0 = config['receptors'][wildcards.location]['lon'] 
            lat0 = config['receptors'][wildcards.location]['lat']
            ds = ds.sel(longitude=slice(0,180), latitude=slice(90,0))
            ds = ds.squeeze()
            fig = plt.figure(figsize=(12,12.3))
            if wildcards.season == 'DJF':
                clevs=np.arange(12,37,3)
                extend='max'
            else:
                clevs = np.arange(12,33,3)
                extend='neither'
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.03, hspace=0.05)
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_large_scale(ax1)
            im = ds.hws_clim.where(ds.hws_clim > 10, drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=clevs, 
                                                                    cmap='viridis_r', extend=extend, add_colorbar=False,ax =ax1)
            CS = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black', 
                                        linewidths=.7, levels=np.arange(510,592,2), add_labels=False, alpha=.7)
            CS1 = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black', 
                                         linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax1.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1, zorder=1030)
            ax1.set_xlabel('')
            ax1.set_ylabel('')
            ax1.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Wind speed [m/s]',shrink=0.9)
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            cax2 = fig.add_subplot(spec[1,1])
            im = ds.hws.plot.contourf(transform=ccrs.PlateCarree(),  levels= np.linspace(-7.5,7.5,16),
                                                                cmap='bwr', add_colorbar=False, ax=ax2)

            fig.colorbar(im, cax=cax2, 
                         label='Composite difference 500hPa winds [m/s] \n  (strong years - weak years)',
                        ticks=im.levels)
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            CS = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', linewidths=.7, 
                                   add_labels=False, alpha=.7,  vmin=-6, vmax=6, levels=20)
            CS1 = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', linewidths=2, levels=CS.levels[::4], 
                                    add_labels=False)
            ax2.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1,zorder=1030)
            Q = ax2.quiver(ds.longitude[::22], ds.latitude[::22], ds.u[::22,::22], 
                           ds.v[::22,::22],transform=ccrs.PlateCarree(),color='saddlebrown', 
                      scale_units='xy', zorder=1002, minlength=2, pivot='middle')
            qk=ax2.quiverkey(Q, 0.93,0.96, U=2, label='2 m/s', labelpos='E', coordinates='axes', color='black')
            qk.text.set_backgroundcolor('w')
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            ax2.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            if config['add_fig_title']:
                fig.suptitle('500hPa geopotential height and windspeed at 500hPa {}'.format(wildcards.season)
                             , y=.915,size=16)
            
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()

rule plot_snow_depth_composite:
    input:
        data_path = 'results/composites/snow_temperature/era5.temp_2m_snow.single_level.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/figs/snow_temperature/era5.temp_2m_snow.single_level.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.{extension}'
    wildcard_constraints:
        extension = 'png|pdf'
    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_large_scale
        import xarray as xr
        import numpy as np
        from matplotlib import cm
        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            lon0 = config['receptors'][wildcards.location]['lon'] 
            lat0 = config['receptors'][wildcards.location]['lat']
            fig = plt.figure(figsize=(12,12.3))
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.03, hspace=0.05)
            # Climatology
            ds = ds.sel(longitude=slice(0,180), latitude=slice(90,0))

            ds = ds.assign(t2m=ds.t2m)
            ds = ds.assign(t2m_clim=ds.t2m_clim-273)
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            ax1.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            map_large_scale(ax1)
            ds.t2m_clim.plot.contour(transform=ccrs.PlateCarree(), vmin=-10, vmax=30, levels=17, colors='black', 
                         linewidths=.4, ax=ax1, add_labels=False)
            im = ds.sd_clim.where(ds.sd_clim > 0.005,drop=True).plot.pcolormesh(transform=ccrs.PlateCarree(), 
                                                            add_colorbar=False, add_labels=False, ax =ax1, 
                                                              vmax=1, cmap='viridis_r', extend='max')
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Snow depth [m]',shrink=0.9, extend='max')
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            bwr_dicreate = cm.get_cmap('bwr', 17)
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            CS=ds.t2m.plot.contour(transform=ccrs.PlateCarree(),  colors='black', 
                                     linewidths=.4, ax=ax2, add_labels=False, vmin=-3, vmax=3, levels=np.arange(-3,3.5,0.5))
            im = ds.sd.plot.pcolormesh(transform=ccrs.PlateCarree(), add_colorbar=False, add_labels=False, ax =ax2, 
                                           cmap=bwr_dicreate, vmin=-0.2, vmax=0.2)
            cax2 = fig.add_subplot(spec[1,1])
            fig.colorbar(im,cax=cax2, label='Snow depth anomalies [m]', extend='both')
            ax2.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            plt.savefig(output.outpath, bbox_inches='tight')
            if config['add_fig_title']:
                fig.suptitle('Snowdepth and 2m temeperature {}'.format(wildcards.season), y=.915,size=16)
        else:
            open(output.outpath, 'a').close()
            
rule plot_composite_mslp_wind_vector850hPa:
    input:
        data_path='results/composites/msl_pressure_wind_850hPa/era5.wind_msl.850hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc',
        orography='downloads/ERA5_orography.nc'
    output:
        outpath = 'results/composites/figs/msl_pressure_wind_850hPa/era5.windvectors_msl.850hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.png'
    threads: 1
    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_large_scale
        import xarray as xr
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np

        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            oro = xr.open_dataset(input.orography)
            oro = oro.sel(longitude=slice(69,105), latitude=slice(40,27)).isel(time=0)
            lon0 = config['receptors'][wildcards.location]['lon'] 
            lat0 = config['receptors'][wildcards.location]['lat']
            ds = ds.sel(longitude=slice(0,180), latitude=slice(90,0))
            ds = ds.squeeze()
            ds = ds.assign(msl_clim=ds.msl_clim/100, msl = ds.msl/100)
            
            if wildcards.season == 'DJF':
                contour_levels = np.linspace(998, 1032, 16) 
                clevs_comp = np.linspace(-8, 8,17)
                extend='max'
            else:
                contour_levels = np.linspace(1004, 1024, 16)
                clevs_comp = np.linspace(-6,6,17) 
                extend='neither'
            
            fig = plt.figure(figsize=(12,12.3))
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.03, hspace=0.05)
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            ax1.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            map_large_scale(ax1)
            im = ds.msl_clim.plot.contourf(transform=ccrs.PlateCarree(),levels=contour_levels , cmap='viridis_r'
                                                        , extend=extend, add_colorbar=False,ax =ax1)

            Q = ax1.quiver(ds.longitude[::22], ds.latitude[::22], ds.u_clim[::22,::22], 
                           ds.v_clim[::22,::22],transform=ccrs.PlateCarree(),color='saddlebrown', 
                      scale_units='xy', zorder=1002, minlength=2, pivot='middle')
            qk = ax1.quiverkey(Q, 0.93,0.96, U=4, label='4 m/s', labelpos='E',color='black',coordinates='axes')
            oro.z.where(oro.z > 33000,drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=2, colors='lightgrey'
                                                                    , add_colorbar=False,ax =ax1, zorder=1100, hatches='/')

            qk.text.set_backgroundcolor('w')
            ax1.set_xlabel('')
            ax1.set_ylabel('')
            ax1.set_title('')
            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Mean Sea Level Pressure [hPa]',shrink=0.9, format='%d', ticks=im.levels)
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            im = ds.msl.plot.contourf(transform=ccrs.PlateCarree(),levels=16, vmin=-6, vmax=6,
                                                                cmap='bwr', add_colorbar=False, ax=ax2)
            oro.z.where(oro.z > 33000,drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=2, colors='lightgrey'
                                                                    , add_colorbar=False,ax =ax2, zorder=1100, hatches='/')
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            ax2.set_title('')
            cax2 = fig.add_subplot(spec[1,1])
            fig.colorbar(im, cax=cax2, 
                 label='Mean Sea Level Pressure \n  (strong years - weak years)',
                ticks=im.levels)

            Q = ax2.quiver(ds.longitude[::22], ds.latitude[::22], ds.u[::22,::22], 
                           ds.v[::22,::22],transform=ccrs.PlateCarree(),color='saddlebrown', 
                      scale_units='xy', zorder=1002, minlength=2, pivot='middle')
            qk=ax2.quiverkey(Q, 0.93,0.96, U=2, label='2 m/s', labelpos='E', coordinates='axes', color='black')
            qk.text.set_backgroundcolor('w')
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            if config['add_fig_title']:
                fig.suptitle('Mean sea level Pressure and 850 hPa windvectors {}'.format(wildcards.season), y=.915,size=16)
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()
            
rule plot_composite_mslp_precipitation:
    input:
        data_path = 'results/composites/precip_mslp/era5.precip_mslp.single_level.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc'
    threads: 1
    output:
        outpath = 'results/composites/figs/precip_mslp/era5.precip_mslp.single_level.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.png'
    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_close_up
        import xarray as xr
        from matplotlib import cm
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np
        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            lon0 = config['receptors'][wildcards.location]['lon'] 
            lat0 = config['receptors'][wildcards.location]['lat']
            ds = ds.sel(longitude=slice(70,140), latitude=slice(60,25))
            ds = ds.squeeze()
            ds = ds.assign(msl_clim=ds.msl_clim/100, msl = ds.msl/100)
            ds = ds.assign(tp_clim=ds.tp_clim*1000, tp=ds.tp*1000)

            fig = plt.figure(figsize=(11.5,11.2))
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.01, hspace=0.07)
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_close_up(ax1)
            im = ds.tp_clim.where(ds.tp_clim > 1, drop=True).plot.contourf(transform=ccrs.PlateCarree(),  vmin=1,
                                            cmap='viridis_r', extend='neither', add_colorbar=False,ax =ax1, vmax=15, levels=15)
            CS = ds.msl_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black', linewidths=.7, 
                                          levels=np.arange(1008,1028,1), add_labels=False, alpha=.7)
            CS1 = ds.msl_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black', linewidths=2, 
                                           levels=[1008,1010,1015,1020,1025], add_labels=False)
            ax1.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1, zorder=1030)
            ax1.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            ax1.set_xlabel('')
            ax1.set_ylabel('')

            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Total Precipitation [mm]',shrink=0.9)
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_close_up(ax2)
            im = ds.tp.plot.contourf(transform=ccrs.PlateCarree(),  vmin=-5, vmax=5, levels=16,
                                            cmap='bwr', extend='both', add_colorbar=False,ax =ax2)
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            cax2 = fig.add_subplot(spec[1,1])
            fig.colorbar(im, cax=cax2, label='Composite difference precipitation [mm] \n strong years - weak years',shrink=0.9, ticks=im.levels)
            CS = ds.msl.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', linewidths=.7,
                                     add_labels=False, alpha=.7, vmin=-4, vmax=4, levels=20)
            CS1 = ds.msl.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', linewidths=2, 
                                      levels=CS.levels[::4], add_labels=False)
            ax2.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1,zorder=1030)
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            ax2.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            if config['add_fig_title']:
                fig.suptitle('Mean sea level Pressure and mean total precipitation {}'.format(wildcards.season), y=.915,size=16)
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()

rule plot_composite_mslp_200hPa_wind:
    input:
        data_path='results/composites/msl_pressure_wind_200hPa/era5.wind_msl.200hPa.composite.{kind}.{psize}.{season}.{location}.{region}.{std}.{sdate}-{edate}.nc'
    output:
        outpath = 'results/composites/figs/msl_pressure_wind_200hPa/era5.wind_msl.200hPa.composite.{kind}.{psize}.{seasons}.{location}.{region}.{std}.{sdate}-{edate}.png'
    threads: 1    

    run:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from scripts.plot_ERA5 import map_large_scale
        import xarray as xr
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np

        ds = xr.open_dataset(input.data_path)
        if 'longitude' in ds.dims: # Check that input data isn't empty
            lon0 = config['receptors'][wildcards.location]['lon'] 
            lat0 = config['receptors'][wildcards.location]['lat']
            ds = ds.sel(longitude=slice(0,180), latitude=slice(90,0))
            ds = ds.assign(msl_clim=ds.msl_clim/100, msl = ds.msl/100)
            ds = ds.squeeze()

            fig = plt.figure(figsize=(12,12.3))
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[10,0.3], wspace=0.03, hspace=0.05)
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_large_scale(ax1)
            im = ds.hws_clim.where(ds.hws_clim >= 20, drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=np.arange(20,60,5),
                                                                    cmap='viridis_r', extend='neither', add_colorbar=False,ax =ax1)
            CS = ds.msl_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black',
                                          linewidths=.7, levels=np.arange(1004,1028,1), add_labels=False, alpha=.7)
            CS1 = ds.msl_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black', 
                                           linewidths=2, levels=[1005,1010,1015,1020,1025], add_labels=False)
            ax1.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1, zorder=1030)
            ax1.set_xlabel('')
            ax1.set_ylabel('')

            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            ax1.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            cax1 = fig.add_subplot(spec[0,1])
            fig.colorbar(im,cax=cax1,label='Wind speed [m/s]',shrink=0.9)
            # # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            im = ds.hws.plot.contourf(transform=ccrs.PlateCarree(),
                                                                cmap='bwr', add_colorbar=False, ax=ax2, levels=np.linspace(-9,9, 16))
            ax2.text(x=0.01, y=0.03,
                     s='Composite {} {} {} \n weak year = {} \n strong years = {}'.format(wildcards.location,wildcards.psize, 
                                                                                       wildcards.kind ,ds.weak_years, ds.strong_years), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)


            # # # Anomalies
            cax2 = fig.add_subplot(spec[1,1])
            fig.colorbar(im, cax=cax2, label='Composite difference (strong years - weak years) [m/s]',shrink=0.9, ticks=im.levels)

            CS = ds.msl.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', 
                                     linewidths=.7, add_labels=False, alpha=.7, vmin=-6, vmax=6, levels=20)
            CS1 = ds.msl.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', 
                                      linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax2.clabel(CS1, fmt='%.2f', colors='black', fontsize=12, inline=1)
            ax2.scatter(lon0,lat0,transform=ccrs.PlateCarree(), marker='*', color='black')
            if config['add_fig_title']:
                fig.suptitle('Mean sea level Pressure and windspeed at 200hPa {}'.format(wildcards.season), y=.915,size=16)
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()