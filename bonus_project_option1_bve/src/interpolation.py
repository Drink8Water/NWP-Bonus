"""Small interpolation helpers for regular latitude-longitude source grids."""

import numpy as np


def interp_regular_latlon(field, lon2d, lat2d, lon_target, lat_target):
    """Bilinearly interpolate a regular lat-lon field to target lon/lat points."""
    lon = np.asarray(lon2d[0, :], dtype=float)
    lat = np.asarray(lat2d[:, 0], dtype=float)
    src = np.asarray(field, dtype=float)

    if lon[0] > lon[-1]:
        lon = lon[::-1]
        src = src[:, ::-1]
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        src = src[::-1, :]

    if not np.any(np.isfinite(src)):
        return np.full_like(lon_target, np.nan, dtype=float)

    lon_t = np.clip(lon_target, lon[0], lon[-1])
    lat_t = np.clip(lat_target, lat[0], lat[-1])

    ix1 = np.searchsorted(lon, lon_t, side='right')
    iy1 = np.searchsorted(lat, lat_t, side='right')
    ix1 = np.clip(ix1, 1, len(lon) - 1)
    iy1 = np.clip(iy1, 1, len(lat) - 1)
    ix0 = ix1 - 1
    iy0 = iy1 - 1

    x0 = lon[ix0]
    x1 = lon[ix1]
    y0 = lat[iy0]
    y1 = lat[iy1]
    wx = np.divide(lon_t - x0, x1 - x0, out=np.zeros_like(lon_t), where=(x1 != x0))
    wy = np.divide(lat_t - y0, y1 - y0, out=np.zeros_like(lat_t), where=(y1 != y0))

    f00 = src[iy0, ix0]
    f10 = src[iy0, ix1]
    f01 = src[iy1, ix0]
    f11 = src[iy1, ix1]
    interp = ((1 - wx) * (1 - wy) * f00
              + wx * (1 - wy) * f10
              + (1 - wx) * wy * f01
              + wx * wy * f11)

    nearest = src[np.where(wy < 0.5, iy0, iy1), np.where(wx < 0.5, ix0, ix1)]
    return np.where(np.isfinite(interp), interp, nearest)
