# ------------------------------------------------------------------------- #
# IMPORT LIBRARIES                                                          #
# ------------------------------------------------------------------------- #

import geopandas as gpd
import xarray as xr
from shapely.geometry import Polygon
from shapely.geometry import Point
import numpy as np
import os
import glob
import h5py
from collections import OrderedDict
#import matplotlib.pyplot as plt

# ------------------------------------------------------------------------- #
# BEGIN FUNCTION DEFITION                                                   #
# ------------------------------------------------------------------------- #


def bin_IS2(grid_file, data_file, reduce='mean', itd=False, thres=50):
    '''
    "bin_IS2" aggregates point data (satellite measurements, buoys, anything
    with a (lat,lon), into bins that correspond to locations on a grid; the
    function also reduces the data into one measurement per grid cell.

    Inputs:
    (1) grid_file: a netcdf file containing at minimum the four vertices of
        every point on the grid.
    (2) data_file: a netcdf file containing the (lat,lon) time, and data 
        values for each point measurement. There are no limits on the number
        of type of data included in the file (though strings may get tricky
        when using the reduction function (means, medians, etc). Pre-process
        as you dare.
    (3) reduce: the type of aggregate function used to reduce the data. Options
        include any possible "aggfunc" argument to geopandas.dissolve.
    (4) itd: flag for whether the "ITD" should be computed. ITD is here used
        to represent any histogram type distribution, it does not have to be
        related to ice thickness. Default is False. THIS IS NOT YET IMPLEMENTED.
    (5) thres: threshold for the number of points within a grid cell considered
        that is considered feasible for aggregation into a distribution.
        Default is 50.

    Outputs:
    (1) reduced_data: a GeoPandas DataFrame of containing the Polygon
        information for each gridcell, the reduced measurements included in
        data_file, and the locations and time of each measurement.

    Author: Molly M. Wieringa (mmw906@uw.edu, Univ. Washington ICG)
    '''
    # open the grid file
    grid = xr.open_dataset(grid_file)

    # open the data file
    ds = xr.open_dataset(data_file)

    # cycle through all grid indices ...
    grid_cells = []
    for x in range(0, np.shape(grid.TLAT)[0]):
        for y in range(0, np.shape(grid.TLON)[1]):

            # collect the latitude (y) & longitude (x) of each of the vertices
            LONT = grid.lont_bounds[x, y, :].values
            LATT = grid.latt_bounds[x, y, :].values

            # if any of the grid corner data is NaN, pass on it
            if np.any(np.isnan(LONT, LATT)) > 0:
                pass

            # otherwise, generate a Polygon from the corners
            else:
                grid_cells.append(Polygon([(LONT[0], LATT[0]),
                                           (LONT[1], LATT[1]),
                                           (LONT[2], LATT[2]),
                                           (LONT[3], LATT[3])]))

    # turn the Polygons into a GeoDataFrame with the standard WGS84 CRS
    cell_data = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                                 crs='epsg:3411')

    # convert the point data data into Points
    df = ds.to_dataframe().reset_index()

    # generate Points with an iterator
    points = [Point(x, y) for x, y in zip(df['longitude'], df['latitude'])]

    # convert to a GeoDataFrame
    point_data = gpd.GeoDataFrame(df, geometry=points, crs=cell.crs)

    # combine the grid data and the point data into on GeoDF and drop NaNs
    binned_data = gpd.sjoin(cell_data, point_data, how='left',
                            predicate='contains').dropna()

    # pull and assign each bin a "cell_number". This is identical to the index,
    # but accessible to operate over for grouping functions
    index_list = list(binned_data.index.values)
    binned_data["cell_number"] = index_list

    # reduce the data with the a given function
    reduced_data = binned_data.dissolve(by='cell_number', aggfunc=reduce)

    if itd is True:
        raise NotImplementedError

        # now do the distribution: group by cell_number and aggregate in lists
        # temp = binned_data.groupby('cell_number').agg(lambda x: list(x))

        # determine the number of freeboard values per index
        # counts = [len(x) for x in temp["freeboard"]]

        # save counts as part of the GeoDataFrame
        # reduced_data["count"] = counts

        # filter out cells that do not have sufficent satellite returns
        # reduced_data = reduced_data.where(reduced_data["count"] > thres).dropna()

    return reduced_data


def read_ATL10(data_in, data_out,
               tracks=['gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l']):
    '''
    "read_IS2" preprocesses IceSat-2 ATL10, converting it from hdf5 files to
    netcdf, and extracting the freeboard, uncertainty and surface type flags.

    Inputs:
    (1) data_in: the file path to the top-level directory where the data is
        stored
    (2) data_out: the file path the the directory where the processed data
        should be saved
    (3) tracks: the IS2 beams that to be included in the dataset

    Outputs:
    (1) null: saved files are written to the designated directory. A friendly
        message indicates when the function finishes processing each directory.

    Author: Molly M. Wieringa (mmw906@uw.edu, Univ. Washington ICG)
    '''

    dir_list = [name for name in sorted(os.listdir(data_in))]
    for dir_day in dir_list:
        files = sorted(glob.glob(data_in + dir_day + '/ATL10-01*[0-9]*.h5'))
        if len(files) == 0:
            print('No files in '+dir_day+'.')
        else:
            data_list = []
            for file in files:
                with h5py.File(file, mode='r') as f:
                    for track in tracks:
                        if len(f['%s/freeboard_beam_segment/' % track].keys()) < 2:
                            pass
                        else:
                             # access the latitude
                            latvar = f['/%s/freeboard_beam_segment/beam_freeboard/latitude' % track]
                            latitude = latvar[:]

                            # access the longitude
                            lonvar = f['/%s/freeboard_beam_segment/beam_freeboard/longitude' % track]
                            longitude = lonvar[:]

                            # access the freeboard values 
                            datavar = f['/%s/freeboard_beam_segment/beam_freeboard/beam_fb_height' % track]
                            data = datavar[:]

                            # access the freeboard uncertainty
                            unc_var = f['/%s/freeboard_beam_segment/beam_freeboard/beam_fb_sigma' % track]
                            unc = unc_var[:]

                            # access the open ocean mask
                            mask_var = f['/%s/freeboard_beam_segment/height_segments/height_segment_ssh_flag' % track]
                            mask = mask_var[:]
                            
                            # access the surface type classification
                            surf_var = f['/%s/freeboard_beam_segment/height_segments/height_segment_type' % track]
                            surf = surf_var[:]

                            # handle FillValue
                            _FillValue = datavar.attrs['_FillValue']
                            data[data == _FillValue] = np.nan
                            unc[unc == _FillValue] = np.nan

                            # collect time information
                            timevar = f['/%s/freeboard_beam_segment/delta_time' % track]
                            time = timevar[:]

                            # make dataset, decode time to cftime object, and save in a list 
                            data_list.append(xr.decode_cf(xr.Dataset({"freeboard": (["time"], data),
                                                                      "uncertainty": (["time"], unc),
                                                                      "mask": (["time"], mask),
                                                                      "surface_type": (["time"], surf),
                                                                      "longitude": (["time"], longitude),
                                                                      "latitude": (["time"], latitude)
                                                                      },
                                                                      coords={"time": (['time'], time, {"units": "seconds since 2018-01-01"}),
                                                                              "track":track,
                                                                              "granule":file[-18:-10]
                                                                              }
                                                                     )
                                                          )
                                             )

            day_data = xr.concat(data_list, dim = 'time', compat = 'no_conflicts')
            filename = 'IS2_'+dir_day+'.nc'
            day_data.to_netcdf(data_out+filename)
            print(dir_day + ' done!')