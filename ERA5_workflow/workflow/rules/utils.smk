

rule setup_folder_structure:
    output:
        data_folder=config['intermediate_files'],
        figs=config['figs_path'],
        monthly_data=config['download_path'],
        nao=config['nao_output_path']
    shell:
        """
        [[ -d {output.data_folder} ]] || mkdir {output.data_folder}
        [[ -d {output.monthly_data} ]] || mkdir {output.monthly_data}
        [[ -d {output.figs} ]] || mkdir {output.figs}
        [[ -d {output.nao} ]] || mkdir {output.nao}
        """
rule geopotential_to_geopotential_height:
    input: DOWNLOADS + "era5.{plevel}hPa.Geopotential.monthly.{sdate}-{edate}.nc"
    output: DOWNLOADS + "era5.{plevel}hPa.GeopotHeight.monthly.{sdate}-{edate}.nc"
    params: g=g0
    shell:
        """
        ncap2 -v -s'z=z/{params.g}f*0.1f;' {input} {output} 
        ncatted -a units,z,o,c,'dam' {output}
        ncatted -a long_name,z,o,c,'Geopotential height' {output}  
        ncatted -a standard_name,z,o,c,'geopotential height' {output}
        ncrename -v z,Z {output}
        """


rule clean:
    params: 
        data_folder=config['intermediate_files'],
        figs=config['figs_path']
    shell:
        """
        find {params.data_folder} -type f \! -name "*monthly*.nc" -exec rm {{}} \;
        rm {params.data_folder}*GeopotHeight*.nc
        rm -r {params.figs}
        """
