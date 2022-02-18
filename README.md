# CDSAPItools

These functions are designed to download ERA5 data per month. In other words, the function will split the required duration into 1-month intervals. The idea behind this module is to prevent having to keep a Python instance running during the entire duration of the queue process. This is implemented in two steps:

1. Request the data to the server, which will queue all requests. All requests are stored in `out_path/ERA5props.csv`, but no download is made yet. 
2. Manually check what requests have been completed and download those that have. Again, this updates `ERA5props.csv` so only new files are downloaded


## Requirements

`CDSAPItools` require the `pandas` and `cdsapi` modules:

```bash
pip install cdsapi pandas
````
 
## References:

### CDS background

- Homepage: https://cds.climate.copernicus.eu/#!/home
- Datasets: https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset


### Method implementation:

- https://github.com/ecmwf/cdsapi/issues/2
- https://github.com/ecmwf/cdsapi/blob/master/examples/example-era5-update.py


## To do:

- [ ] 2022-02-17: Still need to adapt the queries in `submitIt` to match different datasets - for now only `reanalysis-era5-pressure-levels`
- [ ] 2022-02-18: Find why counter is not updatin