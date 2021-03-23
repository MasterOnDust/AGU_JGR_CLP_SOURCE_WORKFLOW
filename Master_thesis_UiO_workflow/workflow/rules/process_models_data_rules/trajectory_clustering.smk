
def get_trajectory_paths(path0,years, kind, location, size):
    paths=[]
    for year in years:
        temp_path=path0+'/'+kind+'/'+size +'/'+str(year)
        for date_tag in ['0331_00','0430_00','0531_00']:
            paths.append(temp_path+'/'+ str(year)+date_tag+location+'/output/trajectories.txt')
    return paths


rule cluster_trajectories:
    input:
        paths=lambda wildcards: get_trajectory_paths(config['flexpart_output'],
                                                     [str(year) for year in range(int(wildcards.sdate), int(wildcards.edate)+1)], 
                                                     wildcards.kind, 
                                                     wildcards.location, wildcards.size)
    output:
        outpath = 'results/model_results/cluster_trajectories/trajectories.{kind}.{location}.{size}.{niters}.{sdate}-{edate}.nc'
    run:
        from scripts.process_model_data_scripts.process_trajectories import (read_trajectories,
                                                                    cluster_trajectories)
        import pandas as pd
        lons, lats, heights = read_trajectories(input.paths)
        ds = cluster_trajectories(lons, lats, heights, niters=int(wildcards.niters))
        ds.to_netcdf(output.outpath)
        
        
rule cluster_trajectories_composites:
    input:
        timeseries = 'results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.' + str(config['m_sdate'])+'-'+str(config['m_edate'])+'.nc',
        paths_trajectories = lambda wildcards: get_trajectory_paths(config['flexpart_output'],
                            [str(year) for year in range(int(config['m_sdate']), int(config['m_sdate'])+1)], 
                            wildcards.kind, 
                            wildcards.location, 
                            wildcards.psize)
    output:
        outpath = 'results/model_results/figs/{composite_kind}.{kind}.{location}.{region}.{psize}.{threshold}.{fig_file_extention}'
    threads: 1
    run:
        from scripts.create_composites import select_years_to_composite, detrend_timeseries
        from scripts.process_model_data_scripts.process_trajectories import (read_trajectories,
                                                                    cluster_trajectories, plot_trajectories)
        import matplotlib.pyplot as plt 
        timeseries = xr.open_dataset(input.timeseries)
        timeseries[timeseries.varName] = detrend_timeseries(timeseries[timeseries.varName]) 
        weak_years, strong_years = select_years_to_composite(timeseries[timeseries.varName],wildcards.threshold)
        
        
        if wildcards.composite_kind =='strong':
            strong_years = get_trajectory_paths(config['flexpart_output'], strong_years, wildcards.kind, wildcards.location, wildcards.psize)
            lons,lats, height = read_trajectories(strong_years, wildcards.kind)
            trajec = cluster_trajectories(lons,lats,height)
        elif wildcards.composite_kind == 'weak': 
            weak_years = get_trajectory_paths(config['flexpart_output'],weak_years, wildcards.kind, wildcards.location, wildcards.psize)
            lons,lats, height = read_trajectories(weak_years, wildcards.kind) 
            trajec = cluster_trajectories(lons,lats,height)
        else:
            lons, lats,height = read_trajectories(input.paths_trajectories, wildcards.kind)
            trajec = cluster_trajectories(lons,lats,height)
        plot_trajectories(trajec)
        ax = plt.gca()
        ax.set_title('{} years {} {} {}'.format(wildcards.composite_kind,wildcards.location, wildcards.psize, wildcards.kind))
        plt.savefig(output.outpath, bbox_inches='tight')

        

    