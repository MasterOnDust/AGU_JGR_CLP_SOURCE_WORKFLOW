
rule plot_geopot_500hPa_windspeed_composites:
    input:
        data_path = 'results/composites/windspeed_geopot/era5.wind_geopot.500hPa.composite.{kind}.{psize}.MAM.{location}.{region}.anom.05-std.{weak_strong}.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/figs/windspeed_geopot/era5.wind_geopot.500hPa.composite.{kind}.{psize}.MAM.{location}.{region}.anom.05-std.{weak_strong}.{sdate}-{edate}.png'
    
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
            fig = plt.figure(constrained_layout=True, figsize=(10,13.5))
            spec = fig.add_gridspec(ncols=2, nrows=3, width_ratios=[16,0.5])
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_large_scale(ax1)
            im = ds.hws_clim.where(ds.hws_clim > 10, drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=np.arange(10,33,3), 
                                                                    cmap='viridis_r', extend='neither', add_colorbar=False,ax =ax1)
            CS = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black',
                                        linewidths=.7, levels=np.arange(510,592,2), add_labels=False, alpha=.7)
            CS1 = ds.Z_clim.plot.contour(transform=ccrs.PlateCarree(), ax=ax1,colors='black',
                                         linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax1.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1, zorder=1030)
            Q = ax1.quiver(ds.longitude[::22], ds.latitude[::22], ds.u_clim[::22,::22],
                           ds.v_clim[::22,::22],transform=ccrs.PlateCarree(),color='saddlebrown', 
                      scale_units='xy', zorder=1002, minlength=2, pivot='middle')
            ax1.quiverkey(Q, 0.9,0.8, U=20, label='20 m/s', labelpos='E')
            ax1.set_xlabel('')
            ax1.set_ylabel('')
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)

            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)
            im = ds.hws.where(ds.hws >= 10, drop=True).plot.contourf(transform=ccrs.PlateCarree(), levels=np.arange(10,33,3), 
                                                                cmap='viridis_r', add_colorbar=False, ax=ax2)

            CS = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black',
                                   linewidths=.7, levels=np.arange(510,592,2), add_labels=False, alpha=.7)
            CS1 = ds.Z.plot.contour(transform=ccrs.PlateCarree(), ax=ax2,colors='black', 
                                    linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax2.clabel(CS1, fmt='%d', colors='black', fontsize=12, inline=1,zorder=1030)
            Q = ax2.quiver(ds.longitude[::22], ds.latitude[::22], ds.u[::22,::22], 
                           ds.v[::22,::22],transform=ccrs.PlateCarree(),color='saddlebrown', 
                      scale_units='xy', zorder=1002, minlength=2, pivot='middle')
            ax2.quiverkey(Q, 0.9,0.8, U=20, label='20 m/s', labelpos='E')
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            ax2.text(x=0.01, y=0.03, s='years composited = {}'.format(ds.years_composited), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            # Colorbar windspeed and composites
            cax1 = fig.add_subplot(spec[0:2,1])
            fig.colorbar(im,cax=cax1,label='Wind speed [m/s]',shrink=0.9)
            # Anomalies
            ax3 = fig.add_subplot(spec[2,0], projection=ccrs.PlateCarree())
            cax2 = fig.add_subplot(spec[2,1])
            map_large_scale(ax3)
            im = ds.hws_anomalies.where((ds.hws_anomalies > 0.4) | (ds.hws_anomalies < -0.4), drop=True).plot.contourf(transform=ccrs.PlateCarree(), 
                                                                    cmap='bwr', extend='neither', levels=10, add_colorbar=False, ax=ax3)
            fig.colorbar(im, cax=cax2, label='Wind speed anomalies [m/s]',shrink=0.9)
            CS = ds.Z_anomalies.plot.contour(transform=ccrs.PlateCarree(), ax=ax3,colors='black', 
                                             linewidths=.7, levels=np.arange(-3,3.3,0.3), add_labels=False, alpha=.7)
            CS1 = ds.Z_anomalies.plot.contour(transform=ccrs.PlateCarree(), ax=ax3,colors='black', 
                                              linewidths=2, levels=CS.levels[::4], add_labels=False)
            ax3.clabel(CS1, fmt='%.2f', colors='black', fontsize=12, inline=1)
            ax3.text( x=0.03,y=0.94, s='c)', fontsize=16, transform=ax3.transAxes)

            ax3.set_xlabel('')
            ax3.set_ylabel('')
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()

rule plot_snow_depth_anomalies:
    input:
        data_path = 'results/composites/snow_temperature/era5.temp_2m_snow.single_level.composite.{kind}.{psize}.MAM.{location}.{region}.anom.05-std.{weak_strong}.{sdate}-{edate}.nc'
    output:
        outpath='results/composites/figs/snow_temperature/era5.temp_2m_snow.single_level.composite.{kind}.{psize}.MAM.{location}.{region}.anom.05-std.{weak_strong}.{sdate}-{edate}.png'
    
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
            ds = ds.assign(t2m=ds.t2m - 273)
            ds = ds.assign(t2m_clim=ds.t2m_clim-273)
            fig = plt.figure(constrained_layout=True, figsize=(10,13.5))
            spec = fig.add_gridspec(ncols=2, nrows=3, width_ratios=[16,0.5])
            # Climatology
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree())
            map_large_scale(ax1)
            ds.t2m_clim.plot.contour(transform=ccrs.PlateCarree(), levels=np.arange(-10, 30, 4), colors='black', 
                                     linewidths=.4, ax=ax1, add_labels=False)
            im = ds.sd_clim.where(ds.sd_clim > 0.005,drop=True).plot.pcolormesh(transform=ccrs.PlateCarree(), add_colorbar=False, add_labels=False, ax =ax1, 
                                          vmax=1, cmap='viridis_r', extend='max')
            ax1.text(x=0.01, y=0.03, s='Climatology = {}-{}'.format(ds.climatology[0],ds.climatology[1]), 
                    bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax1.transAxes, zorder=1050)
            ax1.text( x=0.03,y=0.94, s='a)', fontsize=16, transform=ax1.transAxes)
            # Composite difference
            ax2 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree())
            map_large_scale(ax2)

            ax2.text( x=0.03,y=0.94, s='b)', fontsize=16, transform=ax2.transAxes)
            ax2.text(x=0.01, y=0.03, s='years composited = {}'.format(ds.years_composited), 
                                bbox={'facecolor': 'lightpink',  'pad': 3}, transform=ax2.transAxes, zorder=1050)
            ds.t2m.plot.contour(transform=ccrs.PlateCarree(), levels=np.arange(-10, 30, 4), colors='black', 
                                     linewidths=.4, ax=ax2, add_labels=False)
            im = ds.sd.where(ds.sd > 0.005,drop=True).plot.pcolormesh(transform=ccrs.PlateCarree(), add_colorbar=False, add_labels=False, ax =ax2, 
                                          vmax=1, cmap='viridis_r',extend='max')

            # Colorbar windspeed and composites
            cax1 = fig.add_subplot(spec[0:2,1])
            fig.colorbar(im,cax=cax1,label='Snow depth [m]',shrink=0.9)
            # Anomalies
            ax3 = fig.add_subplot(spec[2,0], projection=ccrs.PlateCarree())
            cax2 = fig.add_subplot(spec[2,1])
            map_large_scale(ax3)
            im = ds.sd_anomalies.plot.pcolormesh(transform=ccrs.PlateCarree(),add_colorbar=False, add_labels=False,ax=ax3, vmin=-0.25, vmax=0.25,
                                            cmap='bwr')
            ds.t2m_anomalies.plot.contour(transform=ccrs.PlateCarree(),add_labels=False, colors='black', ax= ax3, levels=np.arange(-2,2.1,0.5), linewidths=0.8)
            fig.colorbar(im,cax=cax2, label='Snow depth anomalies [m]')
            plt.savefig(output.outpath, bbox_inches='tight')
        else:
            open(output.outpath, 'a').close()