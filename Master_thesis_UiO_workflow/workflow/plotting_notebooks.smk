"""
Run all the notebooks for generating figures and such

"""


rule plot_source_contrib_composite:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['total_deposition','drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        micron2='figs/composite_source_contrib_2micron_tot_dep.pdf',
        micron20='figs/composite_source_contrib_20micron_tot_dep.pdf'
    notebook:
        'notebooks/Deposition_composite_difference.py.ipynb'

rule plot_ao_mo_composite:
    input:
        expand(config['intermediate_files']+'/era5.850hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind'],
                season=['DJF','MAM']),
        expand(config['intermediate_files']+'/era5.500hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind','GeopotHeight'],
                season=['DJF','MAM']),
        'results/ao/era5.1000hPa.AO_EOF.DJF.1979-2019.nc',
        'results/eawmi/era5.single_level.EAWM_MO.DJF.1979-2019.nc'
    output:
        path_500hpa ='figs/winter_MO_AO_composite.pdf',
        path_850hpa = 'figs/winter_MO_AO_composite_500h.pdf'
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


rule plot_850hPa_composite:
    input:
        data850hpa= expand(config['intermediate_files']+'/era5.850hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind'],
                season=['DJF','MAM']),
        data_single_level=expand(config['intermediate_files']+'/era5.single_level.{variable}.{season}.1979-2019.nc',
                variable = ['mean_sea_level_pressure'],
                season=['DJF','MAM']),
        total_depo=expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['total_deposition','drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate']),
        oro = 'downloads/ERA5_orography.nc'
    output:
        silt_composite_path_djf='figs/Composite_maps/mslp_850hPa_20micron_DJF.pdf',
        clay_composite_path_djf='figs/Composite_maps/mslp_850hPa_2micron_DJF.pdf', 
        silt_composite_path_mam='figs/Composite_maps/mslp_850hPa_20micron_MAM.pdf',
        clay_composite_path_mam='figs/Composite_maps/mslp_850hPa_2micron_MAM.pdf'   
    notebook:
        'notebooks/Depostion_circulation_composite_850hPa_and_mslp.ipynb'

rule plot_dust_loading_trajectories:
    input:
        trajec_files = expand('results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc',
                         size=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        path_map = 'figs/average_dust_transport_trajectories.pdf',
        path_vertical_profile = 'figs/average_dust_transport_height.pdf'
    notebook:
        'notebooks/Dust_loading_trajectory_analysis.py.ipynb'

rule plot_dust_emissions:
    input:
        'results/timeseries_table.csv',
        'results/model_results/intermediate_results/emission_flux.china.MAM.1999-2019.nc'
    output:
        ems_map='figs/emission_map_1999_2019.pdf',
        time_series='figs/emission_timeseries_1999_2019.pdf'
    notebook:
        'notebooks/DustEmissions.py.ipynb'

rule plot_source_contribution:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['total_deposition','drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        clay_plot = 'figs/2micron_total_depositon_source_contribution.pdf',
        silt_plot ='figs/20micron_total_depositon_source_contribution.pdf'
    notebook:
        'notebooks/Dust_depostion_source_contribution.py.ipynb'

rule plot_weak_strong_dust_trajecs:
    input:
        expand('results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc',
                         size=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate']),
        'results/timeseries_table.csv'

    output:
        weak_strong_clay_drydep='figs/2_micron_drydep_weak_strong_trajecs.pdf',
        weak_strong_silt_drydep='figs/20_micron_drydep_weak_strong_trajecs.pdf',
        weak_strong_silt_wetdep='figs/20_micron_wetdep_weak_strong_trajecs.pdf',
        weak_strong_clay_wetdep='figs/2_micron_wetdep_weak_strong_trajecs.pdf'
    notebook:
        'notebooks/Weak_strong_dust_loading_trajectories.py.ipynb'

rule plot_500hPa_composite:
    input:
        expand(config['intermediate_files']+'/era5.500hPa.{variable}.{season}.1979-2019.nc',
                variable = ['u_component_of_wind', 'v_component_of_wind','GeopotHeight'],
                season=['DJF','MAM']),

        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','YINCHUAN', 'LANTIAN', 'LUOCHUAN'],
              kind=['total_deposition','drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        silt_composite_path_djf='figs/Composite_maps/geopot_ws_500hPa_20micron_DJF.pdf',
        clay_composite_path_djf='figs/Composite_maps/geopot_ws_500hPa_2micron_DJF.pdf', 
        silt_composite_path_mam='figs/Composite_maps/geopot_ws_500hPa_20micron_MAM.pdf',
        clay_composite_path_mam='figs/Composite_maps/geopot_ws_500hPa_2micron_MAM.pdf'
    
    notebook:
        'notebooks/Depostion_500hPa_composte.py.ipynb'