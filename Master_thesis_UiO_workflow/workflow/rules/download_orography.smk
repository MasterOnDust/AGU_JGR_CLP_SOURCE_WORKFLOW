rule download_orography:
    output:
        outpath = 'downloads/ERA5_orography.nc'
    run:
        import cdsapi

        c = cdsapi.Client()

        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': 'orography',
                'year': '2019',
                'month': '03',
                'day': '10',
                'time': '09:00',
            },
            output.outpath)
