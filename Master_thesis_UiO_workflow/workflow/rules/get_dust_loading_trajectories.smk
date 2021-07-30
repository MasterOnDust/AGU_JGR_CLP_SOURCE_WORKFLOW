
def get_trajectory_paths(w):
    paths=[]
    years = [str(year) for year in range(int(w.sdate), int(w.edate)+1)]
    for year in years:
        temp_path=config['flexpart_path']+'/'+w.kind+'/'+w.size +'/'+year
        for date_tag in ['0331_00','0430_00','0531_00']:
            paths.append(temp_path+'/'+ year+date_tag+w.loc+'/output/trajectories.txt')
    return paths



rule build_dust_loading_trajectories:
    input:
        'results/timeseries_table.csv',
        source_contrib_paths= lambda wildcards: get_paths(config['source_contrib_path'], wildcards.sdate, wildcards.edate,
                                           wildcards.kind, wildcards.loc, wildcards.size),
        paths_trajecs = get_trajectory_paths
    output:
        outpath='results/model_results/trajectories/dust_loading_traj_{kind}_{size}_{loc}_{sdate}-{edate}.nc'
    wildcard_constraints:
        size='2micron|20micron',
        kind='drydep|wetdep'

    params:
        threshold = 1e-9,
        use_dask = True
    run:
        import xarray as xr
        import numpy as np
        from fpcluster.read_trajectories import load_trajectories, read_trajectories
        paths_trajecs = input.paths_trajecs
        paths_source_contrib = input.source_contrib_paths
        time_stamps = []
        ddep = []
        threshold = params.threshold
        use_dask = params.use_dask
        for path in paths_source_contrib:
            if use_dask:
                ds = xr.open_dataset(path,chunks={'time':30})
            else: 
                ds = xr.open_dataset(path)
            ds = ds[ds.varName].sum(dim=['btime','lon','lat'])
            ddep.append(ds.values)
            time_stamps.append(ds.time.values)
        ddep = np.concatenate(ddep)
        time_stamps = np.concatenate(time_stamps)
        time_stamps = time_stamps[ddep >= threshold]
        ddep = ddep[ddep >= threshold]
        ds = read_trajectories(paths_trajecs,kind=wildcards.kind,timestamps=time_stamps)
        da = xr.DataArray(data=ddep, dims=['time'], coords={'time':time_stamps})
        da = da[~da.get_index("time").duplicated()] #avoid duplicate indices?
        ds = ds.sel(time=ds.time[~ds.get_index("time").duplicated()])
        ds = ds.assign({wildcards.kind:da})
        ds.attrs['varName'] = wildcards.kind
        ds.to_netcdf(output.outpath)