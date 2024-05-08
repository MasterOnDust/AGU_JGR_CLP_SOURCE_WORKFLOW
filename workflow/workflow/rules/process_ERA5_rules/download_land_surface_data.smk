
rule download_land_surface_data:
    input:
        expand(rules.seasonal_average.output.outpath, level="single_level", 
                    varName=['sea_ice_cover', 'snow_albedo', 'snow_depth',
            'soil_temperature_level_1', 'soil_temperature_level_2', 'volumetric_soil_water_layer_1',
            'volumetric_soil_water_layer_2'], season=['DJF','MAM'],
            sdate=1979, edate=2019)