"""
Beta-plane Cartesian grid for the barotropic vorticity equation model.

Creates a local Cartesian grid on the beta-plane centred at a reference
latitude phi0.  The grid is regular in (x, y) but the physical lat-lon
locations are used for interpolation of reanalysis data.
"""

import numpy as np

# Earth and physical constants
A = 6.371e6          # Earth radius (m)
OMEGA = 7.292e-5     # Earth rotation rate (s^-1)
G = 9.80665          # Gravity (m s^-2)


class BetaPlaneGrid:
    """Regular Cartesian beta-plane grid centred at (lat0, lon0)."""

    # Expose constants as class attributes for convenient access
    A = A
    OMEGA = OMEGA
    G = G

    def __init__(self, lat0=40.0, lon0=115.0,
                 lat_range=(15.0, 65.0), lon_range=(60.0, 170.0),
                 dx_deg=1.0, dy_deg=1.0):
        """
        Parameters
        ----------
        lat0, lon0 : float
            Reference point for the beta-plane (degrees).
        lat_range : tuple
            (lat_min, lat_max) in degrees.
        lon_range : tuple
            (lon_min, lon_max) in degrees.
        dx_deg, dy_deg : float
            Grid spacing in degrees longitude / latitude.
        """
        self.lat0 = np.radians(lat0)
        self.lon0 = np.radians(lon0)
        self.lat0_deg = lat0
        self.lon0_deg = lon0

        # 1-D coordinate arrays in degrees
        self.lat_1d = np.arange(lat_range[0], lat_range[1] + dy_deg/2, dy_deg)
        self.lon_1d = np.arange(lon_range[0], lon_range[1] + dx_deg/2, dx_deg)

        self.ny = len(self.lat_1d)
        self.nx = len(self.lon_1d)
        self.dy_deg = dy_deg
        self.dx_deg = dx_deg

        # 2-D meshgrid (degrees)
        self.lon2d, self.lat2d = np.meshgrid(self.lon_1d, self.lat_1d)

        # Convert to radians
        self.lat_rad = np.radians(self.lat2d)
        self.lon_rad = np.radians(self.lon2d)

        # Beta-plane Cartesian coordinates (m)
        #   x = a cos(phi0) * (lon - lon0)
        #   y = a * (lat - phi0)
        self.x = A * np.cos(self.lat0) * (self.lon_rad - self.lon0)
        self.y = A * (self.lat_rad - self.lat0)

        # Grid spacings in metres
        self.dx = A * np.cos(self.lat0) * np.radians(dx_deg)
        self.dy = A * np.radians(dy_deg)

        # Coriolis parameters
        self.f0 = 2.0 * OMEGA * np.sin(self.lat0)        # constant f0
        self.beta = 2.0 * OMEGA * np.cos(self.lat0) / A  # df/dy

        # Full f field on the grid
        self.f = 2.0 * OMEGA * np.sin(self.lat_rad)

        print(f"Beta-plane grid: {self.nx} x {self.ny} points")
        print(f"  dx = {self.dx/1000:.1f} km, dy = {self.dy/1000:.1f} km")
        print(f"  f0 = {self.f0:.2e} s^-1, beta = {self.beta:.2e} m^-1 s^-1")
        print(f"  lat: {lat_range[0]:.0f}–{lat_range[1]:.0f}°N")
        print(f"  lon: {lon_range[0]:.0f}–{lon_range[1]:.0f}°E")

    def latlon_to_xy(self, lat, lon):
        """Convert lat/lon (degrees) to beta-plane (x, y) in metres."""
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        x = A * np.cos(self.lat0) * (lon_rad - self.lon0)
        y = A * (lat_rad - self.lat0)
        return x, y
