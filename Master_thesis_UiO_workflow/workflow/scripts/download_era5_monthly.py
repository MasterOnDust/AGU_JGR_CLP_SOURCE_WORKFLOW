import cdsapi
import argparse


def download_era5_monthly_pressure_level(pressure_level, variable,sdate,edate,filename):
    pressure_level=str(pressure_level)
    c = cdsapi.Client()
    years = [str(year) for year in range(sdate, edate+1)]
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variable,
            'pressure_level': pressure_level,
            'year': years
            ,
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
        },
        filename)

def download_era5_monthly_single_level(variable,sdate,edate,filename):
    c = cdsapi.Client()
    years = [str(year) for year in range(sdate, edate+1)]
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variable,
            'year': years
            ,
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
        },
        filename)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Download monthly ERA5 data from cds")
    parser.add_argument('pressure_level', 
            help="Which pressure level to download, if = -1, download single level data")
    parser.add_argument("variable", help="variable to download")
    parser.add_argument('sdate', help="Starting year to download datafrom", type=int)
    parser.add_argument("edate", help='last year to download data from', type=int)
    parser.add_argument("--filename", "--fn", default="download.nc", 
                        help="which filename to save the data to")
    
    args=parser.parse_args()
    pressure_level=args.pressure_level
    sdate=args.sdate
    edate=args.edate
    variable=args.variable
    filename=args.filename
    if pressure_level=='-1':
        download_era5_monthly_single_level(variable,sdate,edate,filename)
    else:
        download_era5_monthly_pressure_level(pressure_level,variable,sdate,edate,filename)
    