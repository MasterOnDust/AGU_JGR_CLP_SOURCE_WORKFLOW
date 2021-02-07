if __name__=='__main__':
    parser=ap.ArgumentParser(description='Command line interface for plotting ERA5 data')
    parser.add_argument('path', help='path to ERA5 dataset')
    parser.add_argument('--method', '--m', default='pcolormesh', help='which ploting function to use')
    parser.add_argument('--title', '--t', default=None, help='title of plot')
    parser.add_argument('--filename', '--fn', default=None, help='Filename to which plot should be stored')
    parser.add_argument('--variable', '--v', default=None, 
            help='If there are more than one variable in the dataset variable name has to be specified')
    args = parser.parse_args()
    path=args.path
    title=args.title
    file_name=args.filename
    var=args.variable
    method=args.method
    dset=read_data(path,var)
    da=dset[dset.varName]
    fig,ax =plt.subplots(subplot_kw={'projection':ccrs.Robinson(central_longitude=90.0)})
    print(method)
    ax = draw_map(ax)
    if method=='pcolormesh':
        ax=plot_pcolormesh(da,ax=ax)
    elif method=='contour':
        ax=plot_contour(da, ax=ax)
    elif method=='contourf':
        ax=plot_contourf(da,ax=ax)
    else:
        raise(ValueError('{} is not a valid plotting method choose either pcolormesh, contour, contourf'.format(method)))
    
    plt.show()
    if file_name.endswith('.pdf'):
        plt.savefig(file_name, bbox_inches='tight')
    elif file_name.endswith('png'):
        plt.savefig(file_name, dpi=300,bbox_inches='tight')
    else: 
        raise(ValueError("invalid file format {}".format(file_name)))