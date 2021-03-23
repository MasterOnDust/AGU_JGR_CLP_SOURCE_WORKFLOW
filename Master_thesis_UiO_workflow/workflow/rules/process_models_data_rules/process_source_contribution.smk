

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
        from DUST.utils.resample import resample_monthly, concatenate_monthly
        import time
        paths=input.paths
        dsets_list=[resample_monthly(xr.open_dataset(path)) for path in paths]

        dsets=concatenate_monthly(dsets_list)


        dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),'resample_source_contrib') + dsets.attrs['history']

        dsets.to_netcdf(output.outpath)