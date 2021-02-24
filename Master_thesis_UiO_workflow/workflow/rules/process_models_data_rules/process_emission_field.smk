
############################################################################
#DUST source regions:                                                      #
#      lat          long                                                   #
#;A: 37-42  ; 100-110    North west                                        #
#;B: 43-50  ; 95-110     Mongolia                                          #
#;C: 43-47  ; 80-90     Jungger Basin                                      #
#;D: 36-42  ; 75-90    Taklamakan                                          #
#;E: 35-40 ; 90-100     Qaidam basin                                       #
#;F: 42-45 ; 112-124    North east                                         #
############################################################################

from scripts.process_model_data_scripts.emission_flux_series import resample_emission_flux, create_timeseries

rule resample:
    input:
        script='workflow/scripts/emission_flux_series.py',
        paths=expand(config['flexdust_path']+'/FLEXDUST1999_2019/{year}/FLEXDUST_{year}0301_{year}0531.nc', 
                    year=[year for year in range(config['sdate'], config['edate']+1)])
    
    output:
        outpath=config['intermediate_results_models']+'/emission_flux.china.{frequency}.{sdate}-{edate}.nc'
    
    wildcard_constraints:
        frequency='MAM|monthly'
    run:
        if wildcards.frequency=='monthly':
            freq='m'
        else:
            freq=wildcards.frequency
        ds = resample_emission_flux(input.paths, config['domain']['lon0'],config['domain']['lon1'],
                                   config['domain']['lat0'], config['domain']['lat1'],
                                   frequency=freq)
        ds.to_netcdf(output.outpath)
        
rule source_region_timeseries:
    input:
        path=config['intermediate_results_models']+'/emission_flux.china.{frequency}.{sdate}-{edate}.nc'
    output:
        outpath=config['flexdust_results']+'/emission_flux.time_series.{region}.{frequency}.{sdate}-{edate}.nc'
    wildcard_constraints:
        region='|'.join(config['source_regions'].keys())+'|total'
    run:
        
        if wildcards.region=='total':
            ts=create_timeseries(input.path)
        else:
            region=config['source_regions'][wildcards.region]
            ts=create_timeseries(input.path,region['lon0'], region['lon1'], region['lat0'], region['lat1'])
            ts=xr.decode_cf(ts)
#             from IPython import embed; embed()
            t0=str(ts.time[0].dt.strftime('%Y').values)
            t1=str(ts.time[-1].dt.strftime('%Y').values)
            ts.attrs['title']='FLEXDUST simulated dust emissions {}-{}, '.format(t0,t1) + wildcards.region
        ts.to_netcdf(output.outpath)