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
import scipy.interpolate.lagrange as lagrange
from numpy.polynomial.polynomial import Polynomial

# ------------------------------------------------------------------------- #
# BEGIN FUNCTION DEFINITION                                                   #
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

    Adding in floe-chord detection and pairing
              level vs deformed ice (gradient filter, if gradient in FB exceeds
              a threshold)
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
                            ssh_flag_var = f['/%s/freeboard_beam_segment/height_segments/height_segment_ssh_flag' % track]
                            ssh_flag = ssh_flag_var[:]
                            

                            # access the surface type classification
                            surf_var = f['/%s/freeboard_beam_segment/height_segments/height_segment_type' % track]
                            surf = surf_var[:]
                            
                            # access the segment length
                            lenseg_var = f['/%s/freeboard_beam_segment/height_segments/height_segment_length_seg' % track]
                            lenseg = lenseg_var[:]

                            # handle FillValue
                            _FillValue = datavar.attrs['_FillValue']
                            data[data == _FillValue] = np.nan
                            unc[unc == _FillValue] = np.nan

                            # collect time information
                            timevar = f['/%s/freeboard_beam_segment/beam_freeboard/delta_time' % track]
                            time = timevar[:]


                            # detect_level_ice
                            level_ice_mask = detect_level_ice(datavar)

                            # make dataset, decode time to cftime object, and save in a list 
                            data_list.append(xr.decode_cf(xr.Dataset({"freeboard": (["time"], data),
                                                                      "uncertainty": (["time"], unc),
                                                                      "ssh_flag": (["time"], ssh_flag),
                                                                      "level_ice":(["time"], level_ice_mask),
                                                                      "surface_type": (["time"], surf),
                                                                      "segment_length": (["time"], lenseg),
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
                            
                            # Run additional processing functions
                            data = compute_floe_chords(data)
                            
                            # Save in a list 
                            data_list.append(data)

            day_data = xr.concat(data_list, dim = 'time', compat = 'no_conflicts')
            filename = 'IS2_'+dir_day+'.nc'
            day_data.to_netcdf(data_out+filename)
            print(dir_day + ' done!')

            

            
def compute_floe_chords(data):
    """
    "compute_floe_chords" generates an ice-ocean mask for the ICESat-2 data in
    xarray format. Floe chords are measured from the ice-ocean mask and each 
    ice segment in ICESat-2 is assigned with the length of the floe chord it
    is part of.
    
    Inputs:
    (1) xarray.DataSet: ICESat-2 data
    
    Outputs:
    (1) xarray.DataSet: ICESat-2 data including a floe_chord variable

    Author: Nils Hutter (nhutter@uw.edu, Univ. Washington ICG)
    """
    
    # Initialise field to output floe chords
    data['floe_chord'] = data.freeboard.copy()
                                             
    
    # Create ice mask
    non_nan = np.isfinite(data.freeboard)
    ice = (~np.logical_or(np.logical_and(data.freeboard==0,data.surface_type>=6),
                                            data.ssh_flag>0)).where(np.isfinite(data.freeboard))

    # Check for series of valid segments
    valid_segs = contiguous_regions(non_nan.data)

    # Loop through all series of valid segments
    for ivseg in valid_segs:
        # Detect floes in the ice-ocean mask
        segs = contiguous_regions(ice.data[ivseg[0]:ivseg[1]])
        for iseg in segs:
            if iseg[0]!=0 and iseg[1]!=ivseg[1]-ivseg[0]:
                data['floe_chord'][iseg[0]+ivseg[0]] = np.sum(data.segment_length[iseg[0]+ivseg[0]])
                
    return data

    
    
            
def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero() 

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size-1] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx



def detect_level_ice(freeboard_data):
    '''
    Implements a gradient test to detect level or presumably deformed ice. For
    now, this test is applied to freeboard data and follows von Albedyll et al.
    (2022). The approach was to anchor the search at each subsequent index and
    scan forward through the satellite points until a distance between points
    of more than 100m is reached. This segment is then treated according to
    Rabenstein et al. (2010) to determine whether any sufficient part of the
    segment (including the beginning portion that may continue from the
    previous segment) is considered level. The function then steps forward to
    the next index and interates.

    Inputs:
    (1) freeboard data

    Outputs:
    (1) mask array of "level" or "deformed" ice

    From Rabenstein et al. (2010): Level ice was identified using two criteria.
    First the numerical differentiation of sea-ice thickness along the profile
    using a three-point Langrangian interpolator must be <0.04, and second,
    level-ice sections must extend at least 100m in length...
    '''

    # initialize an empty flags object, set index step to 1, and set found to False
    flags = np.zeros(len(freeboard_data))
    count = 1
    found = False

    # cycle through all the data...
    for i in range(0, len(freeboard_data)):
        # use the found flag to keep searching while dist < 100 m
        x_interp = []
        while found is not True:
            # using the anchor point of the cycle... 
            x1 = freeboard_data.latitude[i]
            y1 = freeboard_data.longitude[i]
            # grab the next point along the satellite track to consider
            # based on the count...
            x2 = freeboard_data.latitude[i+count]
            y2 = freeboard_data.longitude[i+count]
            # and calculate the haversine distance between the two points.
            dist = 1000 * haversine_distance(x1, x2, y1, y2)

            # if the distance is < 100 meters...
            if dist < 100:
                # up the count to look one point further...
                count += 1
                # and store the distance for the interpolation later.
                x_interp.append(dist)
            # if the distance is > 100 meters
            else:
                # set the found flag True to break the loop
                found = True
                # set the index to the count reached
                index = count
                # set the count back to 1 for next while loop
                count = 1

        # for the 100 meters-ish segment that's been found...
        data_to_interpolate = freeboard_data.freeboard[i:index]

        # interpolate between the points...
        poly = lagrange(x_interp, data_to_interpolate)
        curve = Polynomial(poly.coef[::-1])*(x_interp)

        # ... and calculate the gradient of the interpolated curve
        gradient = np.gradient(curve)

        # assign a temporary flag object
        temp_flags = flags[i:index]

        # find indicies in gradient for which points do not satisfy reqs
        bad_points = []
        for x in range(0, len(gradient)):
            if gradient[x] > 0.04:
                bad_points.append(x)

        # find the bad point indices pairs that have level ice > 100 m between
        level_lengths = []
        prev_bp = bad_points[0]
        for bp in bad_points[1:]:
            # if the index before the bad point is not equal to the previous bad point... 
            if bp-1 != prev_bp:
                # calculate the length between the point before the bad one and the previous bad one
                sub_length = x_interp[bp-1] - x_interp[prev_bp]
                # if the length is >= 100m...
                if sub_length >= 100:
                    # store the indices within the gradient curve for those segments
                    level_lengths.append((prev_bp+1, bp-1))
            prev_bp = bp

        # first, if this ISN'T the first segment and the last flag of the previous one was 1...
        if i != 0 and flags[i-1] == 1:
            # ... set the flags up to the first bad point to 1
            temp_flags[0:bad_points[0]] = 1
        # if this IS the first segment or the last flag ISN'T 1...
        else:
            # go through the indices pairs for level_lengths > 100 m...
            for pair in level_lengths:
                # ... and replace flags with 1 between the two indices
                temp_flags[pair[0]:pair[1]] = 1

        # check that new flags and old flags don't conflict
        if i == 0:
            prev_flags = temp_flags

        # determine conflicts between overlapping flag segments
        temp_overlap = temp_flags[:-1]
        prev_overlap = prev_flags[1:]
        conflicts = temp_overlap - prev_overlap

        # collect locations of conflicts
        conflict_locations = []
        for x in range(0, len(conflicts)):
            if conflicts[x] != 0:
                conflict_locations.append(x)

        # for each conflict location...
        for c_loc in conflict_locations:
            # ... if the previous set of flags saved the conflict location as level...
            if prev_overlap[c_loc] == 1:
                # ... save the new flag at that overlap as level as well.
                temp_overlap[c_loc] = 1
            # if the previous set didn't save as level...
            else:
                # ... resolve the conflict by saving as not level
                temp_overlap[c_loc] = 0

        # save the corrected conflicts back into the temporary flag object
        temp_flags[:-1] = temp_overlap

        # once the corrections have been applied save a for the next set of prev_flags
        prev_flags = temp_flags

        # save current temp_flags back into list of all flags
        flags[i:index] = temp_flags

    # save all calculated flags to an array and return mask
    mask = np.array(flags)

    return mask


def haversine_distance(x1, x2, y1, y2):

    dist = 2*np.arcsin(np.sqrt(np.sin((x1 - y1)/2)**2 +
                               np.cos(x1)*np.cos(y1)*np.sin((x2-y2)/2)**2))
    return dist

