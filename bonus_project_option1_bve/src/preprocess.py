"""
Preprocessing module for reading and processing local ERA5 NetCDF data.

Supports NetCDF files from data/raw/.  Adapts to common variable name
conventions and handles conversion of geopotential to geopotential height.
"""

import os
import glob
import numpy as np
import xarray as xr


# Common variable name patterns
VAR_NAMES = {
    'geopotential': ['z', 'geopotential', 'gh', 'hgt', 'phi'],
    'u_wind': ['u', 'u_component_of_wind', 'UGRD', 'eastward_wind'],
    'v_wind': ['v', 'v_component_of_wind', 'VGRD', 'northward_wind'],
    'geopotential_height': ['z', 'geopotential_height', 'gh', 'hgt'],
}

COORD_NAMES = {
    'latitude': ['latitude', 'lat'],
    'longitude': ['longitude', 'lon'],
    'pressure': ['pressure_level', 'level', 'isobaricInhPa', 'plev', 'isobaric'],
    'time': ['time', 'valid_time', 't'],
}


def _find_name(ds, candidates, category):
    """Find a variable or coordinate name from a list of candidates."""
    all_names = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for name in candidates:
        if name in all_names:
            return name
    available = ', '.join(sorted(all_names))
    raise KeyError(
        f"Cannot find {category}. Tried: {candidates}. Available: {available}"
    )


def find_data_files(data_dir='data/raw', search_parent=True):
    """
    Search for NetCDF files in data_dir.

    Also searches the parent directory if search_parent is True and
    data_dir yields no results.
    """
    patterns = [os.path.join(data_dir, '*.nc'),
                os.path.join(data_dir, '*.nc4'),
                os.path.join(data_dir, '*.h5')]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    if not files and search_parent:
        parent = os.path.dirname(data_dir)
        for pat in [os.path.join(parent, '*.nc'), os.path.join(parent, '*.nc4')]:
            files.extend(glob.glob(pat))

    # Deduplicate
    files = sorted(set(f for f in files if os.path.isfile(f)))
    return files


def open_dataset(filepath):
    """Open a NetCDF dataset with xarray and print diagnostics."""
    print(f"\nOpening: {filepath}")
    ds = xr.open_dataset(filepath)

    print(f"  Dimensions: {dict(ds.sizes)}")
    print(f"  Coordinates: {list(ds.coords)}")
    print(f"  Data variables: {list(ds.data_vars)}")

    # Print time info
    time_name = _find_name(ds, COORD_NAMES['time'], 'time coordinate')
    time_vals = ds[time_name].values
    print(f"  Time range: {time_vals[0]} to {time_vals[-1]}")
    print(f"  Time steps: {len(time_vals)}")

    # Print pressure level info
    try:
        p_name = _find_name(ds, COORD_NAMES['pressure'], 'pressure coordinate')
        p_vals = ds[p_name].values
        print(f"  Pressure levels: {p_vals}")
    except KeyError:
        print("  No pressure-level coordinate found (single-level data?)")

    # Print lat/lon info
    lat_name = _find_name(ds, COORD_NAMES['latitude'], 'latitude coordinate')
    lon_name = _find_name(ds, COORD_NAMES['longitude'], 'longitude coordinate')
    print(f"  Lat range: {ds[lat_name].values[0]:.2f} to {ds[lat_name].values[-1]:.2f}"
          f"  ({len(ds[lat_name])} pts)")
    print(f"  Lon range: {ds[lon_name].values[0]:.2f} to {ds[lon_name].values[-1]:.2f}"
          f"  ({len(ds[lon_name])} pts)")

    return ds


def detect_variable_names(ds):
    """Detect geopotential, u, v variable names in the dataset."""
    z_name = _find_name(ds, VAR_NAMES['geopotential'], 'geopotential')
    u_name = _find_name(ds, VAR_NAMES['u_wind'], 'u wind')
    v_name = _find_name(ds, VAR_NAMES['v_wind'], 'v wind')
    print(f"  Detected: Z={z_name}, U={u_name}, V={v_name}")
    return z_name, u_name, v_name


def select_pressure_level(ds, level_hpa=500.0):
    """Select a specific pressure level if the coordinate exists."""
    try:
        p_name = _find_name(ds, COORD_NAMES['pressure'], 'pressure coordinate')
    except KeyError:
        print("  Single-level data — skipping pressure selection.")
        return ds

    p_vals = ds[p_name].values
    if len(p_vals) == 1:
        # Already single-level
        p_actual = float(p_vals[0])
        tolerance = max(1.0, 0.1 * p_actual)
        if abs(p_actual - level_hpa) < tolerance:
            print(f"  Data is at {p_actual} hPa — using as-is.")
            return ds.sel({p_name: p_vals[0]})
        else:
            print(f"  WARNING: data at {p_actual} hPa, expected {level_hpa} hPa")

    if level_hpa in p_vals:
        ds_sel = ds.sel({p_name: level_hpa})
        print(f"  Selected {level_hpa} hPa.")
        return ds_sel

    # Nearest neighbour
    idx = np.argmin(np.abs(p_vals - level_hpa))
    nearest = float(p_vals[idx])
    print(f"  Nearest level to {level_hpa} hPa: {nearest} hPa (selected)")
    return ds.isel({p_name: idx})


def select_times(ds, times):
    """
    Select specific times from the dataset.

    Parameters
    ----------
    ds : xarray.Dataset
    times : list of str
        ISO 8601 datetime strings.

    Returns
    -------
    datasets : list of xarray.Dataset (one per time)
    """
    time_name = _find_name(ds, COORD_NAMES['time'], 'time coordinate')
    time_vals = ds[time_name].values

    results = []
    for t in times:
        t_np = np.datetime64(t)
        if t_np in time_vals:
            ds_t = ds.sel({time_name: t_np}, method=None)
            print(f"  Selected time: {t}")
            results.append(ds_t)
        else:
            # Nearest
            idx = np.argmin(np.abs(time_vals - t_np))
            nearest = time_vals[idx]
            print(f"  Nearest time to {t}: {nearest}")
            results.append(ds.isel({time_name: idx}))

    return results


def subset_region(ds, lat_range=(15.0, 65.0), lon_range=(60.0, 170.0)):
    """
    Subset the dataset to a latitude–longitude box.

    Parameters
    ----------
    ds : xarray.Dataset
    lat_range : tuple (lat_min, lat_max)
    lon_range : tuple (lon_min, lon_max)

    Returns
    -------
    ds_subset : xarray.Dataset
    """
    lat_name = _find_name(ds, COORD_NAMES['latitude'], 'latitude coordinate')
    lon_name = _find_name(ds, COORD_NAMES['longitude'], 'longitude coordinate')

    lat = ds[lat_name].values

    # Determine ordering
    if lat[0] > lat[-1]:
        # Latitudes decreasing (N to S) — typical ERA5
        lat_slice = slice(lat_range[1], lat_range[0])
    else:
        lat_slice = slice(lat_range[0], lat_range[1])

    ds_sub = ds.sel(**{lat_name: lat_slice})

    # Longitude — handle 0–360 wrapping
    lon = ds[lon_name].values
    if lon[0] >= 0 and lon[-1] > 180:
        # 0–360 convention
        lon_min = lon_range[0] % 360
        lon_max = lon_range[1] % 360
        if lon_min <= lon_max:
            ds_sub = ds_sub.sel(**{lon_name: slice(lon_min, lon_max)})
        else:
            # Wrapping across 0
            raise NotImplementedError("Longitude wrapping not implemented; adjust domain.")
    else:
        # -180–180 convention
        ds_sub = ds_sub.sel(**{lon_name: slice(lon_range[0], lon_range[1])})

    lat_vals = ds_sub[lat_name].values
    lon_vals = ds_sub[lon_name].values
    print(f"  Subset region: {lat_vals[0]:.1f}–{lat_vals[-1]:.1f}°N, "
          f"{lon_vals[0]:.1f}–{lon_vals[-1]:.1f}°E")
    print(f"  Subset shape: {ds_sub.sizes}")

    return ds_sub


def coarsen_to_resolution(ds, target_dx_deg=1.0, target_dy_deg=1.0):
    """
    Coarsen a dataset to approximately target resolution via block averaging.

    Parameters
    ----------
    ds : xarray.Dataset
    target_dx_deg, target_dy_deg : float
        Target grid spacing in degrees.

    Returns
    -------
    ds_coarse : xarray.Dataset
    """
    lat_name = _find_name(ds, COORD_NAMES['latitude'], 'latitude coordinate')
    lon_name = _find_name(ds, COORD_NAMES['longitude'], 'longitude coordinate')

    lat = ds[lat_name].values
    lon = ds[lon_name].values

    dlat = abs(lat[1] - lat[0]) if len(lat) > 1 else target_dy_deg
    dlon = abs(lon[1] - lon[0]) if len(lon) > 1 else target_dx_deg

    if dlat >= target_dy_deg and dlon >= target_dx_deg:
        print(f"  Resolution ({dlat:.3f}° x {dlon:.3f}°) already ≥ target "
              f"({target_dy_deg}° x {target_dx_deg}°) — skipping coarsen.")
        return ds

    factor_lat = max(1, int(round(target_dy_deg / dlat)))
    factor_lon = max(1, int(round(target_dx_deg / dlon)))

    print(f"  Coarsening by factor {factor_lat} (lat) × {factor_lon} (lon)")

    # Use xarray coarsen with mean
    ds_coarse = ds.coarsen(
        **{lat_name: factor_lat, lon_name: factor_lon},
        boundary='trim'
    ).mean()

    print(f"  Coarsened shape: {ds_coarse.sizes}")
    return ds_coarse


def convert_geopotential(ds, z_name):
    """
    Convert geopotential Φ (m²/s²) to geopotential height Z (m).

    Returns the variable as Z in metres.  If the variable is already in
    geopotential metres (height), return as-is.
    """
    G = 9.80665
    var = ds[z_name]
    units = var.attrs.get('units', '').lower()

    # Check if data is geopotential (m²/s²) vs height (m)
    if 'm**2' in units or 'm2' in units or 'm^2' in units:
        print(f"  Converting geopotential ({units}) to height (Z = Φ/g)")
        ds = ds.copy()
        ds[z_name] = var / G
        ds[z_name].attrs['units'] = 'm'
        ds[z_name].attrs['long_name'] = 'Geopotential height'
        return ds, z_name

    # Check magnitude to guess
    vals = var.values
    mean_val = np.nanmean(np.abs(vals))
    if mean_val > 10000:
        # Likely geopotential in m²/s² (~50000 at 500hPa)
        print(f"  Variable '{z_name}' has mean |value| = {mean_val:.0f}. "
              f"Assuming geopotential (m²/s²) → converting to height (m).")
        ds = ds.copy()
        ds[z_name] = var / G
        ds[z_name].attrs['units'] = 'm'
        ds[z_name].attrs['long_name'] = 'Geopotential height'
    else:
        print(f"  Variable '{z_name}' appears to already be geopotential height (m).")

    return ds, z_name


def extract_field(ds, var_name):
    """Extract a 2-D numpy array from a dataset, squeezing singleton dims."""
    data = ds[var_name].values
    # Squeeze any singleton dimensions
    while data.ndim > 2:
        # Remove the first singleton dim
        for ax, s in enumerate(data.shape):
            if s == 1:
                data = data.squeeze(axis=ax)
                break
        else:
            # No singleton dim found
            if data.ndim == 3:
                data = data[0]  # Take first if 3-D
            elif data.ndim > 3:
                data = data[0, 0]
            break
    return np.asarray(data, dtype=np.float64)
