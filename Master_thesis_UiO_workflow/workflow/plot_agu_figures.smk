LOCS_AGU_PAPER = ['SACOL','LINGTAI','BADOE','SHAPOTOU','LANTIAN', 'LUOCHUAN']

rule plot_source_contribution_agu:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              kind='total_deposition',sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        clay_plot = 'figs/agu/2micron_total_deposition_source_contribution.pdf',
        silt_plot ='figs/agu/20micron_total_deposition_source_contribution.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution.py.ipynb'

rule plot_source_contribution_wetdep_agu:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              kind='wetdep',sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        clay_plot = 'figs/agu/2micron_wet_deposition_source_contribution.pdf',
        silt_plot ='figs/agu/20micron_wet_deposition_source_contribution.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_wetdep.py.ipynb'

rule plot_source_contribution_drydep_agu:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              kind='drydep',sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        clay_plot = 'figs/agu/2micron_dry_deposition_source_contribution.pdf',
        silt_plot ='figs/agu/20micron_dry_deposition_source_contribution.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_drydep.py.ipynb'

rule plot_source_contribution_drydep_wet_combo_agu:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        clay_plot = 'figs/agu/2micron_dry_deposition_wetdep_combo_source_contribution.pdf',
        # silt_plot ='figs/agu/20micron_dry_deposition_source_contribution.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_drydep_wetdep_combo.py.ipynb'


rule plot_500hPa_composite_agu:
    input:
        expand('results/composites/windspeed_geopot_500hPa/era5.wind_geopot.500hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True)
    output:
        silt_composite_path_djf='figs/agu/Composite_maps/geopot_ws_500hPa_{kind}_20micron_DJF.pdf',
        clay_composite_path_djf='figs/agu/Composite_maps/geopot_ws_500hPa_{kind}_2micron_DJF.pdf', 
        silt_composite_path_mam='figs/agu/Composite_maps/geopot_ws_500hPa_{kind}_20micron_MAM.pdf',
        clay_composite_path_mam='figs/agu/Composite_maps/geopot_ws_500hPa_{kind}_2micron_MAM.pdf'
    
    notebook:
        'notebooks/AGU_paper_figures/Deposition_500hPa_composite.py.ipynb'

rule plot_correlation_matrix_agu:
    input:
        'results/timeseries_table.csv'
    output:
        outpath_clay='figs/agu/correlation_matrix_clay.pdf',
        outpath_silt='figs/agu/correlation_matrix_silt.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Plot_correlation_matrix.py.ipynb'

rule plot_850hPa_composite_agu:
    input:
        outpath=expand('results/composites/msl_pressure_wind_850hPa/era5.wind_msl.850hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
        psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True),
        oro = 'downloads/ERA5_orography.nc'
    output:
        silt_composite_path_djf='figs/agu/Composite_maps/mslp_850hPa_{kind}_20micron_DJF.pdf',
        clay_composite_path_djf='figs/agu/Composite_maps/mslp_850hPa_{kind}_2micron_DJF.pdf', 
        silt_composite_path_mam='figs/agu/Composite_maps/mslp_850hPa_{kind}_20micron_MAM.pdf',
        clay_composite_path_mam='figs/agu/Composite_maps/mslp_850hPa_{kind}_2micron_MAM.pdf'   
    notebook:
        'notebooks/AGU_paper_figures/Depostion_circulation_composite_850hPa_and_mslp.py.ipynb'

rule plot_dust_loading_trajectories_agu:
    input:
        trajec_files = expand('results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc',
                         size=['2micron','20micron'],   loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],
              kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    output:
        path_map = 'figs/agu/average_dust_transport_trajectories.pdf',
        path_vertical_profile = 'figs/agu/average_dust_transport_height.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Dust_loading_trajectory_analysis.py.ipynb'

rule create_deposition_histogram:
    input:
        expand('results/model_results/intermediate_results/timeseries/{kind}/{kind}.{loc}.{psize}.Day.{year}.csv',
        kind=['drydep','wetdep'],loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],psize=['2micron','20micron'], 
        year=[str(y) for y in range(1999,2020)])

rule plot_all_agu:
    input:
        'figs/agu/Composite_maps/mslp_850hPa_drydep_20micron_DJF.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_drydep_2micron_DJF.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_drydep_20micron_MAM.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_drydep_2micron_MAM.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_total_deposition_20micron_DJF.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_total_deposition_2micron_DJF.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_total_deposition_20micron_MAM.pdf',
        'figs/agu/Composite_maps/mslp_850hPa_total_deposition_2micron_MAM.pdf',
        rules.plot_source_contribution_agu.output,
        rules.plot_source_contribution_drydep_agu.output,
        rules.plot_source_contribution_wetdep_agu.output,
        'figs/agu/Composite_maps/geopot_ws_500hPa_total_deposition_20micron_DJF.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_total_deposition_2micron_DJF.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_total_deposition_20micron_MAM.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_total_deposition_2micron_MAM.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_drydep_20micron_DJF.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_drydep_2micron_DJF.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_drydep_20micron_MAM.pdf',
        'figs/agu/Composite_maps/geopot_ws_500hPa_drydep_2micron_MAM.pdf',
        rules.plot_correlation_matrix_agu.output,
        rules.plot_dust_loading_trajectories_agu.output
        