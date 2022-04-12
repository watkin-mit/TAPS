#/usr/bin/env python3

import numpy as np

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in radians)
    """
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371*c
    return km

def gen_gc_xmat(hrz_grid,lat_idx=None):
    lon_vec = (np.pi/180.0)*hrz_grid['lon'].values
    lat_vec = (np.pi/180.0)*hrz_grid['lat'].values
    n_lon = len(lon_vec)
    n_lat = len(lat_vec)

    # Starting location
    src_lon = lon_vec[0]

    # Generate mesh of lons and lats
    lat_mesh, lon_mesh = np.meshgrid(lat_vec,lon_vec,indexing='ij')

    # gc_dist[j,k,l] will be the distance from location [lat = j, lon = i] to location [lat = k, lon = i +/- l]
    if lat_idx is None: 
        # Calculate for all latitudes
        gc_dist = np.zeros((n_lat,n_lat,n_lon))
        # Don't need to explicitly calculate values past the equator
        mid_lat = np.ceil(n_lat/2)
        for i_lat, src_lat in enumerate(lat_vec):
            # i_lat and src_lat refer to the "starting point" of the calculation
            if i_lat > mid_lat:
                # Use the equivalent point in the other hemisphere
                x_lat = (n_lat - 1) - i_lat
                gc_dist[i_lat,...] = np.flip(gc_dist[x_lat,...],axis=0)
            else:
                gc_dist[i_lat,...] = haversine(src_lat,src_lon,lat_mesh,lon_mesh)
    else:
        # Only calculate for the current latitude
        gc_dist = haversine(lat_vec[lat_idx],src_lon,lat_mesh,lon_mesh)
    return gc_dist
