"""
Run all the notebooks for generating figures and such

"""


rule plot_ao_mo_composite:
    input:
        expand(config['intermediate_files']+'/era5.850hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind'],
                season=['DJF','MAM']),
        expand(config['intermediate_files']+'/era5.500hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind','GeopotHeight'],
                season=['DJF','MAM']),
        config['old_base']+'/results/ao/era5.1000hPa.AO_EOF.DJF.1979-2019.nc',
        config['old_base']+'/results/eawmi/era5.single_level.EAWM_MO.DJF.1979-2019.nc',
        config['old_base']+'/downloads/ERA5_orography.nc'
    output:
        path_500hpa='figs/winter_MO_AO_composite_850hPa.pdf',
        path_850hpa='figs/winter_MO_AO_composite_500hPa.pdf'
    notebook:
        'notebooks/AO_MO_composite.ipynb'


rule plot_circulation_clim:
    input:
        expand(config['intermediate_files']+'/era5.850hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind'],
                season=['DJF','MAM']),
        expand(config['intermediate_files']+'/era5.single_level.{variable}.{season}.1979-2019.nc',
                variable = ['mean_sea_level_pressure'],
                season=['DJF','MAM']),
        expand(config['intermediate_files']+'/era5.500hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind','GeopotHeight'],
                season=['DJF','MAM']),
        expand(config['intermediate_files']+'/era5.200hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind','GeopotHeight'],
                season=['DJF','MAM']),
        oro = 'downloads/ERA5_orography.nc'
    output:
        outpath = 'figs/climatology_1999-2019.pdf'
    notebook:
        'notebooks/Circulation_climatology.py.ipynb'

rule plot_emission_eof_analysis:
    input:
        config['old_base']+'/results/model_results/intermediate_results/emission_flux.china.MAM.1999-2019.nc',
        data850hpa= expand(config['intermediate_files']+'/era5.850hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind'],
                season=['DJF','MAM']),
        data_single_level=expand(config['intermediate_files']+'/era5.single_level.{variable}.{season}.1979-2019.nc',
                variable = ['mean_sea_level_pressure'],
                season=['DJF','MAM']),
        oro = 'downloads/ERA5_orography.nc'
    output:
        outpath='figs/eofs/emissions_eof_analysis.pdf',
        composite_path_djf='figs/eofs/emissions_composite_djf.pdf',
        composite_path_mam='figs/eofs/emissions_composite_mam.pdf',
        elow_plot = 'figs/eofs/eof_pca_elow_plot.png'

    notebook:
        'notebooks/EOF_analysis.py.ipynb'

