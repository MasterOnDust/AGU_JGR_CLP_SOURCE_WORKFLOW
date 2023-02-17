

rule calc_total_deposition:
    input:
        wetdep_path=config['intermediate_results_models']+'/wetdep/wetdep.{location}.{psize}.monthly.{sdate}-{edate}.nc',
        drydep_path=config['intermediate_results_models']+'/drydep/drydep.{location}.{psize}.monthly.{sdate}-{edate}.nc'
    output:
        outpath=config['intermediate_results_models']+'/total_deposition/total_deposition.{location}.{psize}.monthly.{sdate}-{edate}.nc'
    
    run:
        from thesis_toolbox.process_model_output.calc_totaldeposition import get_total_depostion
        import time
        
        wetdep=xr.open_dataset(input.wetdep_path)
        drydep=xr.open_dataset(input.drydep_path)
        totdep = get_total_depostion(drydep,wetdep)
        totdep.attrs['history'] = time.ctime() + 'total deposition calculated from ' + input.drydep_path + ' ' + input.wetdep_path
        totdep.attrs['filename']=output.outpath 
        totdep.to_netcdf(output.outpath)
        
rule sum_spring_deposition:
    input:
        depo_data=config['intermediate_results_models']+'/{kind}/{kind}.{location}.{psize}.monthly.{sdate}-{edate}.nc'
    output:
        outpath='results/model_results/{kind}/{kind}.{location}.{psize}.MAM.{sdate}-{edate}.nc'
    wildcard_constraints:
        location='|'.join(config['receptors'].keys())
    run:
        ds = xr.open_dataset(input.depo_data)
        spring_depo=ds[ds.varName].groupby('time.year').sum(dim='time',keep_attrs=True)
        spring_depo=spring_depo.to_dataset(name=ds.varName)
        spring_depo.attrs=ds.attrs
        spring_depo.attrs['filename']=output.outpath
        spring_depo['RELLAT']=ds['RELLAT1']
        spring_depo['RELLNG']=ds['RELLNG1']
        spring_depo.to_netcdf(output.outpath)
        
rule source_contribution_source_region_timeseries:
    input:
        depo_data='results/model_results/{kind}/{kind}.{location}.{psize}.MAM.{sdate}-{edate}.nc'
    output:
        outpath='results/model_results/time_series/{kind}/{kind}.{location}.{region}.{psize}.MAM.{sdate}-{edate}.nc'
    wildcard_constraints:
        region='|'.join(config['source_regions'].keys()) + '|total'
        
    run:
        from thesis_toolbox.process_model_output.process_source_contribution import create_timeseries
        ds=xr.open_dataset(input.depo_data)
        if wildcards.region=='total':
            ds = create_timeseries(ds)
        else:
            region=config['source_regions'][wildcards.region]
            ds=create_timeseries(ds,region['lon0'], region['lon1'], region['lat0'], region['lat1'])
        
        ds.attrs['title']='FLEXDUST/FLEXPART simulated dust deposition' + wildcards.region
        
        ds.to_netcdf(output.outpath)