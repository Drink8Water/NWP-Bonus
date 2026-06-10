"""Lambert conformal model grid for the finite-area BVE experiments."""

import numpy as np
import cartopy.crs as ccrs

from .grid import A, OMEGA, G


def _haversine(lon1, lat1, lon2, lat2):
    """Great-circle distance in metres."""
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    h = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * A * np.arcsin(np.sqrt(np.maximum(h, 0.0)))


class LambertGrid:
    """Regular Lambert conformal grid in projected metres."""

    A = A
    OMEGA = OMEGA
    G = G

    def __init__(self, center_lon=115.0, center_lat=35.0,
                 central_lon=115.0, central_lat=35.0,
                 standard_parallels=(25.0, 45.0),
                 width_m=8.0e6, height_m=5.0e6, dx=150_000.0):
        self.center_lon = center_lon
        self.center_lat = center_lat
        self.central_lon = central_lon
        self.central_lat = central_lat
        self.standard_parallels = tuple(standard_parallels)
        self.dx = float(dx)
        self.dy = float(dx)
        self.width_m = float(width_m)
        self.height_m = float(height_m)

        self.proj = ccrs.LambertConformal(
            central_longitude=central_lon,
            central_latitude=central_lat,
            standard_parallels=standard_parallels,
        )
        pc = ccrs.PlateCarree()
        cx, cy = self.proj.transform_point(center_lon, center_lat, pc)

        self.nx = int(np.round(width_m / dx)) + 1
        self.ny = int(np.round(height_m / dx)) + 1

        x0 = cx - 0.5 * (self.nx - 1) * dx
        y0 = cy - 0.5 * (self.ny - 1) * self.dy
        self.x_1d = x0 + np.arange(self.nx) * dx
        self.y_1d = y0 + np.arange(self.ny) * self.dy
        self.X, self.Y = np.meshgrid(self.x_1d, self.y_1d)

        lonlat = pc.transform_points(self.proj, self.X, self.Y)
        self.lon2d = lonlat[..., 0]
        self.lat2d = lonlat[..., 1]
        self.f = 2.0 * OMEGA * np.sin(np.radians(self.lat2d))
        self.f0 = 2.0 * OMEGA * np.sin(np.radians(center_lat))
        self.m = self._compute_map_factor()

        print(f"Lambert grid: {self.nx} x {self.ny} points")
        print(f"  dx = dy = {self.dx/1000:.1f} km")
        print(f"  lon range: {np.nanmin(self.lon2d):.1f}–{np.nanmax(self.lon2d):.1f}°E")
        print(f"  lat range: {np.nanmin(self.lat2d):.1f}–{np.nanmax(self.lat2d):.1f}°N")
        print(f"  map factor m: {np.nanmin(self.m):.3f}–{np.nanmax(self.m):.3f}")

    @classmethod
    def from_lonlat_bounds(cls, lon_min, lon_max, lat_min, lat_max,
                           central_lon=115.0, central_lat=35.0,
                           standard_parallels=(25.0, 45.0),
                           dx=150_000.0, buffer_m=450_000.0,
                           samples=300):
        """Create a Lambert grid that fully covers a lat-lon rectangle."""
        proj = ccrs.LambertConformal(
            central_longitude=central_lon,
            central_latitude=central_lat,
            standard_parallels=standard_parallels,
        )
        pc = ccrs.PlateCarree()

        bottom_lons = np.linspace(lon_min, lon_max, samples)
        bottom_lats = np.full(samples, lat_min)
        right_lons = np.full(samples, lon_max)
        right_lats = np.linspace(lat_min, lat_max, samples)
        top_lons = np.linspace(lon_max, lon_min, samples)
        top_lats = np.full(samples, lat_max)
        left_lons = np.full(samples, lon_min)
        left_lats = np.linspace(lat_max, lat_min, samples)

        lons = np.concatenate([bottom_lons, right_lons, top_lons, left_lons])
        lats = np.concatenate([bottom_lats, right_lats, top_lats, left_lats])
        xy = proj.transform_points(pc, lons, lats)
        x = xy[:, 0]
        y = xy[:, 1]
        valid = np.isfinite(x) & np.isfinite(y)
        x = x[valid]
        y = y[valid]

        x_min = x.min() - buffer_m
        x_max = x.max() + buffer_m
        y_min = y.min() - buffer_m
        y_max = y.max() + buffer_m
        center_x = 0.5 * (x_min + x_max)
        center_y = 0.5 * (y_min + y_max)
        center_lon, center_lat = pc.transform_point(center_x, center_y, proj)

        width_m = np.ceil((x_max - x_min) / dx) * dx
        height_m = np.ceil((y_max - y_min) / dx) * dx
        return cls(center_lon=float(center_lon), center_lat=float(center_lat),
                   central_lon=central_lon, central_lat=central_lat,
                   standard_parallels=standard_parallels,
                   width_m=float(width_m), height_m=float(height_m), dx=dx)

    def _compute_map_factor(self):
        dist_x = np.empty((self.ny, self.nx - 1))
        dist_y = np.empty((self.ny - 1, self.nx))
        dist_x[:, :] = _haversine(
            self.lon2d[:, :-1], self.lat2d[:, :-1],
            self.lon2d[:, 1:], self.lat2d[:, 1:]
        )
        dist_y[:, :] = _haversine(
            self.lon2d[:-1, :], self.lat2d[:-1, :],
            self.lon2d[1:, :], self.lat2d[1:, :]
        )
        mx_edge = self.dx / np.maximum(dist_x, 1.0)
        my_edge = self.dy / np.maximum(dist_y, 1.0)

        mx = np.empty((self.ny, self.nx))
        my = np.empty((self.ny, self.nx))
        mx[:, 1:-1] = 0.5 * (mx_edge[:, :-1] + mx_edge[:, 1:])
        mx[:, 0] = mx_edge[:, 0]
        mx[:, -1] = mx_edge[:, -1]
        my[1:-1, :] = 0.5 * (my_edge[:-1, :] + my_edge[1:, :])
        my[0, :] = my_edge[0, :]
        my[-1, :] = my_edge[-1, :]
        return 0.5 * (mx + my)

    def to_npz_kwargs(self):
        return {
            'X': self.X,
            'Y': self.Y,
            'lon2d': self.lon2d,
            'lat2d': self.lat2d,
            'm': self.m,
            'f': self.f,
            'dx_m': self.dx,
            'dy_m': self.dy,
            'nx': self.nx,
            'ny': self.ny,
            'center_lon': self.center_lon,
            'center_lat': self.center_lat,
            'central_lon': self.central_lon,
            'central_lat': self.central_lat,
            'standard_parallels': np.array(self.standard_parallels),
            'width_m': self.width_m,
            'height_m': self.height_m,
            'f0': self.f0,
        }


def load_lambert_grid(npz_path):
    """Load a lightweight grid object from grid_info_lambert.npz."""
    data = np.load(npz_path)
    grid = object.__new__(LambertGrid)
    grid.X = data['X']
    grid.Y = data['Y']
    grid.lon2d = data['lon2d']
    grid.lat2d = data['lat2d']
    grid.m = data['m']
    grid.f = data['f']
    grid.dx = float(data['dx_m'])
    grid.dy = float(data['dy_m'])
    grid.nx = int(data['nx'])
    grid.ny = int(data['ny'])
    grid.f0 = float(data['f0'])
    grid.G = G
    grid.A = A
    grid.OMEGA = OMEGA
    return grid
