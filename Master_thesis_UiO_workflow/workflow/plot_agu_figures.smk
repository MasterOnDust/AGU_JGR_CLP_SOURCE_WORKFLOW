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
        'notebooks/AGU_paper_figures/Plot_correlation_matrix_v2.py.ipynb'

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
        loc=['SACOL','LINGTAI','BADOE','SHAPOTOU','LANTIAN', 'LUOCHUAN'], 
        year=[str(y) for y in range(1999,2020)], allow_missing=True)
    output:
        outpath = 'figs/agu/{kind}_{psize}_deposition_histogram.pdf'

    notebook:
        'notebooks/AGU_paper_figures/Create_deposition_histogram.ipynb'

rule plot_combo_plotv2:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
            loc=LOCS_AGU_PAPER,
            kind=['total_deposition','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'],
            allow_missing=True)
    output:
        combopath='figs/agu/fraction_source_contrib_combo_v2_{psize}.png'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_drydep_wetdep_combo_v2.py.ipynb'

rule plot_combo_plotv4:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
            loc=LOCS_AGU_PAPER,
            kind=['total_deposition','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'],
            psize=['2micron','20micron'])
    output:
        combopath='figs/agu/fraction_source_contrib_combo_v4.png'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_drydep_wetdep_combo_v4.py.ipynb'



rule plot_composite_combo:
    input:
        expand('results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
            loc=LOCS_AGU_PAPER,
            kind=['total_deposition','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'],
            psize=['2micron','20micron'])
    output:
        source_contrib_diff_path='figs/agu/fraction_source_contrib_combo_composite.png'
    notebook:
        'notebooks/AGU_paper_figures/Deposition_composite_difference.py.ipynb'


rule plot_composite_combo_850hPa:
    input:
        outpath=expand('results/composites/msl_pressure_wind_850hPa/era5.wind_msl.850hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
        psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True),
        oro = 'downloads/ERA5_orography.nc'
    output:
        composite_facet_plot='figs/agu/combo_850hPa_composite_figure_{kind}.png'
    notebook:
        'notebooks/AGU_paper_figures/Depostion_circulation_composite_850hPa_and_mslp_v2.py.ipynb'

rule plot_psize_composite_combo_850hPa:
    input:
        outpath=expand('results/composites/msl_pressure_wind_850hPa/era5.wind_msl.850hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
        kind=['drydep','wetdep'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True),
        oro = 'downloads/ERA5_orography.nc'
    output:
        composite_facet_plot='figs/agu/combo_psize_850hPa_composite_figure_{psize}.png'
    notebook:
        'notebooks/AGU_paper_figures/Depostion_circulation_composite_850hPa_and_mslp_v3.py.ipynb'

rule plot_composite_combo_500hPa:
    input:
        expand('results/composites/windspeed_geopot_500hPa/era5.wind_geopot.500hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
             psize=['2micron','20micron'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True)
    output:
        composite_facet_plot='figs/agu/combo_500hPa_composite_figure_{kind}.png'
    notebook:
        'notebooks/AGU_paper_figures/Deposition_500hPa_composite_v2.py.ipynb'


rule plot_all_agu:
    input:
        rules.plot_source_contribution_agu.output,
        rules.plot_source_contribution_drydep_agu.output,
        rules.plot_source_contribution_wetdep_agu.output,
        rules.plot_correlation_matrix_agu.output,
        expand(rules.create_deposition_histogram.output, kind=['wetdep','drydep'],psize=['20micron','2micron']),
        expand(rules.plot_850hPa_composite_agu.output, kind=['drydep','wetdep','total_deposition']),
        expand(rules.plot_500hPa_composite_agu.output, kind=['drydep','wetdep','total_deposition']),
 

        