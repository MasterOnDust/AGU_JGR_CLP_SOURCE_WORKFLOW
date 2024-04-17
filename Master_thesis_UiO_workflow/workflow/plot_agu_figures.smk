LOCS_AGU_PAPER = ['SHAPOTOU','SACOL','BAODE','LUOCHUAN','LINGTAI','LANTIAN']



rule plot_correlation_matrix_agu:
    input:
        'results/timeseries_table.csv'
    output:
        outpath_clay='figs/agu/correlation_matrix_clay.pdf',
        outpath_silt='figs/agu/correlation_matrix_silt.pdf'
    notebook:
        'notebooks/AGU_paper_figures/Plot_correlation_matrix_v2.py.ipynb'

rule plot_dust_loading_trajectories_agu:
    input:
        trajec_files = expand(config['old_base']+'/results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc',
                         size=['2micron','20micron'],   loc=['SACOL','LINGTAI','BAODE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],
              kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate']),
        emi_dust = config['intermediate_results_models']+'/emission_flux.china.MAM.1999-2019.nc'
    output:
        path_map = 'figs/agu/average_dust_transport_trajectories.pdf',
        path_vertical_profile = 'figs/agu/average_dust_transport_height.pdf'
    
    notebook:
        'notebooks/AGU_paper_figures/Dust_loading_trajectory_analysis.py.ipynb'

rule create_deposition_histogram:
    input:
        wetdep_2micron = expand(config['intermediate_results_models']+'/timeseries/wetdep/wetdep.{loc}.2micron.Day.{year}.csv',
        loc=LOCS_AGU_PAPER,
        year=[str(y) for y in range(1999,2020)], allow_missing=True),
        wetdep_20micron = expand(config['intermediate_results_models']+'/timeseries/wetdep/wetdep.{loc}.20micron.Day.{year}.csv',
        loc=LOCS_AGU_PAPER,
        year=[str(y) for y in range(1999,2020)], allow_missing=True),
        drypdep_2micron = expand(config['intermediate_results_models']+'/timeseries/drydep/drydep.{loc}.2micron.Day.{year}.csv',
        loc=LOCS_AGU_PAPER, 
        year=[str(y) for y in range(1999,2020)], allow_missing=True),
        drydep_20micron = expand(config['intermediate_results_models']+'/timeseries/drydep/drydep.{loc}.20micron.Day.{year}.csv',
        loc=LOCS_AGU_PAPER,
        year=[str(y) for y in range(1999,2020)], allow_missing=True),

    params:
        locs =  ['SHAPOTOU','SACOL','BAODE','LUOCHUAN','LINGTAI','LANTIAN']
    output:
        outpath = 'figs/agu/deposition_histograms.pdf'

    notebook:
        'notebooks/AGU_paper_figures/Create_deposition_histogram_v2.ipynb'


rule plot_source_contribution_agu:
    input:
        expand(config['old_base']+'/results/model_results/time_series/{kind}/{kind}.{loc}.total.{psize}.MAM.{sdate}-{edate}.nc',
            loc=LOCS_AGU_PAPER,
            kind=['drydep','wetdep'],sdate=config['m_sdate'], edate=config['m_edate'],
            psize=['2micron','20micron'])
    output:
        combopath='figs/agu/fraction_source_contrib_combo_v4.png'
    notebook:
        'notebooks/AGU_paper_figures/Dust_deposition_source_contribution_drydep_wetdep_combo_v4.py.ipynb'




rule plot_psize_composite_combo_850hPa:
    input:
        outpath=expand('results/composites/msl_pressure_wind_850hPa/era5.wind_msl.850hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
        kind=['drydep','total_deposition'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True),
        oro = config['orography_file']
    output:
        composite_facet_plot='figs/agu/combo_psize_850hPa_composite_figure_{psize}.png'
    notebook:
        'notebooks/AGU_paper_figures/Depostion_circulation_composite_850hPa_and_mslp_v4.py.ipynb'


rule plot_composite_combo_500hPa_v2:
    input:
        expand('results/composites/windspeed_geopot_500hPa/era5.wind_geopot.500hPa.composite.{kind}.{psize}.{season}.{loc}.total.4_rank.{sdate}-{edate}.nc',
             kind=['drydep','total_deposition'],   loc=LOCS_AGU_PAPER,
              sdate=config['m_sdate'], edate=config['m_edate'],season=['MAM','DJF'],allow_missing=True),
    output:
        composite_facet_plot='figs/agu/combo_500hPa_composite_figure_{psize}.png'
    notebook:
        'notebooks/AGU_paper_figures/Deposition_500hPa_composite_v3.py.ipynb'


rule plot_wetdep_height_at_arrival:
    input:
        trajec_files = expand(config['old_base']+'/results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc',
                         size=['2micron','20micron'],   loc=['SACOL','LINGTAI','BAODE','SHAPOTOU','LANTIAN', 'LUOCHUAN'],
              kind=['wetdep'],sdate=config['m_sdate'], edate=config['m_edate'])
    params:
        assumed_cloud_base=2500
    
    output:
        facet_plot_wetdep_at_site = 'figs/agu/wetdep_height_at_arrival.pdf',
    notebook:
        'notebooks/AGU_paper_figures/Wetdep_incloud_vs_below_cloud.py.ipynb'


rule plot_all_agu:
    input:
        rules.plot_ao_mo_composite.output,
        rules.plot_correlation_matrix_agu.output,
        rules.plot_dust_loading_trajectories_agu.output,
        rules.create_deposition_histogram.output,
        rules.plot_source_contribution_agu.output,
        expand(rules.plot_composite_combo_500hPa_v2.output, psize=['2micron','20micron']),
        expand(rules.plot_psize_composite_combo_850hPa.output,psize=['2micron','20micron'])
 

        