

def get_paths(path0,sdate, edate, kind, location, size):
    years = [str(year) for year in range(int(sdate), int(edate)+1)]
    paths=[]
    if size=='2micron':
        tag='FINE-SILT'
    else:
        tag='COARSE-SILT'
    for year in years:
        temp_path=path0+'/'+kind+'/'+size+'/'+year+'/'+kind
        for dstart,dend in zip(['0306','0331','0430'],['0331', '0430','0531']):
            paths.append(temp_path+'_'+location+'_'+tag+'_'+year+dstart+'-'+year+dend+'.nc')

    return paths

rule resample_source_contrib:
    input:
        paths= lambda wildcards: get_paths(config['source_contrib_path'], wildcards.sdate, wildcards.edate,
                                           wildcards.kind, wildcards.location, wildcards.size)
        
    output:
        outpath=config['intermediate_results_models']+'/{kind}/{kind}.{location}.{size}.monthly.{sdate}-{edate}.nc'
    wildcard_constraints:
        location='|'.join(config['receptors'].keys()),
        size='2micron|20micron',
        kind='drydep|wetdep'
    threads: 1
    run: 
        from dust.utils.resample import resample_monthly, concatenate_monthly
        import time
        paths=input.paths
        dsets_list=[resample_monthly(xr.open_dataset(path)) for path in paths]

        dsets=concatenate_monthly(dsets_list)


        dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),'resample_source_contrib') + dsets.attrs['history']

        dsets.to_netcdf(output.outpath)
def get_flexpart_input_paths(w):
    first_piece = f'/{w.kind}/{w.size}/{w.year}/{w.year}' 
    sec_piece = f'_00{w.location}/output'
    path_folder = expand(config['flexpart_path']+first_piece+'{edate}'+sec_piece, edate=['0331', '0430','0531'])
    paths = []
    for folder in path_folder:
        ncfile = glob_wildcards(folder+'/{ncfile}.nc')
        paths.append(folder+'/{}.nc'.format(''.join(ncfile.ncfile))) 
    return paths


rule process_flexpart_output:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_paths = lambda wildcards: get_flexpart_input_paths(wildcards)
    output:
        config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0306-{year}0331.nc',
        config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0331-{year}0430.nc',
        config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0430-{year}0531.nc'
    wildcard_constraints:
        tag='FINE-SILT|COARSE-SILT'
    threads:1
    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1']
    run:
        import dust
        import xarray as xr 
        import os
        from DUST.process_data_dust import process_per_pointspec, process_per_timestep, create_output
        # print(input.flexdust_path)
        flexdust_ds = DUST.read_flexdust_output(input.flexdust_path+'/')['dset']
        flexdust_ds = flexdust_ds.sel(lon=slice(params.x0,params.x1), lat=slice(params.y0,params.y1))
        outpath = config['source_contrib_path']+f'/{wildcards.kind}/{wildcards.size}/{wildcards.year}'
        for fp_path in input.flexpart_paths:
            ds = xr.open_dataset(fp_path, chunks={'time':40, 'pointspec':15})
            ds, out_data, surface_sensitivity = process_per_pointspec(ds, flexdust_ds, params.x0, params.x1, params.y0, 
                                            params.y1, height=None)
            relcom_str=str(ds.RELCOM[0].values.astype('U35')).strip().split(' ',2)[1:]
            ds.attrs['relcom']=[s.strip() for s in relcom_str]
            out_ds = create_output(out_data,surface_sensitivity,ds)
            ds.close()
            spec_com = ds.spec001_mr.attrs['long_name']
            f_name = out_ds.attrs['varName']
            shape_dset = out_ds[f_name].shape
            encoding = {'zlib':True, 'complevel':9, 'chunksizes' : (1,10, shape_dset[2], shape_dset[3]),
            'fletcher32' : False,'contiguous': False, 'shuffle' : False}
            outFile_name = os.path.join(outpath,out_ds.attrs['filename'])
            print('writing to {}'.format(outFile_name))
            out_ds.to_netcdf(outFile_name, encoding={f_name:encoding, 'surface_sensitivity':encoding})


rule create_deposition_timeseries:
    input:
        paths= lambda wildcards: get_paths(config['source_contrib_path'], wildcards.year, wildcards.year,
                                           wildcards.kind, wildcards.location, wildcards.size)
        
    output:
        outpath=config['intermediate_results_models']+'/timeseries/{kind}/{kind}.{location}.{size}.{freq}.{year}.csv'
    wildcard_constraints:
        location='|'.join(config['receptors'].keys()),
        size='2micron|20micron',
        kind='drydep|wetdep',
        freq='3H|Day|Weekly'
    run:
        import xarray as xr
        paths = input.paths
        if wildcards.freq=='Day':
            freq='D'
        elif wildcards.freq=='Weekly':
            freq='W'
        else:
            freq=wildcards.freq
        def pre_process(dset):
            dset = dset.drop_vars('surface_sensitivity')
            dset = dset.drop_dims('numpoint')
            with xr.set_options(keep_attrs=True):
                if dset.ind_receptor==4 or dset.ind_receptor==3:
                    # Depostion is accumulative
                    dset = dset.sum(dim=['btime','lon','lat'])
                else:
                    # Concentration is not accumulative
                    dset = dset.sum(dim=['btime','lon','lat'])
            return dset
        
        ds = xr.open_mfdataset(paths,concat_dim=['time'], parallel=True, 
                chunks={'time':40}, preprocess=pre_process)
        if ds.ind_receptor==4 or ds.ind_receptor==3:
            da_dep = ds[wildcards.kind].resample(time=freq).sum()
        else:
            da_dep = ds[wildcards.kind].resample(time=freq).mean()
        df = da_dep.to_pandas()
        df.to_csv(output.outpath)