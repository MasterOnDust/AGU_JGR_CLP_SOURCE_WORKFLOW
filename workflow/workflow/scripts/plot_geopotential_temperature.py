import argparse as ap

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr


from plot_ERA5 import plot_contour, plot_contourf,read_data, draw_map



if __name__ == "__main__":
    parser=ap.ArgumentParser(description='Command line interface for plotting ERA5 data')
    parser.add_argument('--geopotential_height','--gh', help='path geopotential height data', nargs=(1,2))
    parser.add_argument('--temperature_data', '--td', help='path to temperature data', nargs=(1,2))
    parser.add_argument('--title', '--t', default=None, help='title of plot')
    parser.add_argument('--filename', '--fn', default=None, help='Filename to which plot should be stored')
    parser.add_argument('--NH_polar_plot', action='store_true')
    args = parser.parse_args()
    temperature=args.temperature_data
    geopotential_height=args.geopotential_height
    title=args.title
    file_name=args.filename
    pole_plot=args.NH_polar_plot
    t_dset=read_data(temperature)
    geop_deset=read_data(geopotential_height)
    
    da_temp=t_dset[t_dset.varName]
    da_geopot=geop_deset[geop_deset.varName]
    if pole_plot:
        fig,ax =plt.subplots(figsize=(10,12),subplot_kw={'projection':ccrs.Orthographic(0, 90)})
    # fig,ax =plt.subplots(figsize=(10,7),subplot_kw={'projection':ccrs.Mercator(central_longitude=90.0)})
    # fig,ax =plt.subplots(figsize=(10,6),subplot_kw={'projection':ccrs.NorthPolarStereo(central_longitude=90.0)})
    # fig,ax =plt.subplots(figsize=(10,6),subplot_kw={'projection':ccrs.Mollweide(central_longitude=90.0)})
    ax = draw_map(ax, top_labels=True, right_labels=True)
    # ax.set_extent([-120,150,10,90])
    ax=plot_contourf(da_temp)
    ax=plot_contour(da_geopot)
    
    plt.show()

    