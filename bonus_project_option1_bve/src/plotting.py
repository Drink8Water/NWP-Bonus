"""
Plotting module for the BVE forecast project.

Generates maps of geopotential height, vorticity, forecast verification,
and summary score charts.

Plotting philosophy:
    Regional Lambert conformal maps should use the projected outline of the
    source latitude-longitude domain as the map frame.  The visible boundary is
    therefore a curved quadrilateral: constant-longitude side edges and
    constant-latitude top/bottom edges after Lambert projection.  Filled
    fields, coastlines, borders, and gridlines are clipped to this boundary.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
from matplotlib.patches import PathPatch

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    ccrs = None
    cfeature = None
    print("cartopy not installed — using plain lat-lon contour plots.")


# ═══════════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════════

MAP_EXTENT = [60, 170, 15, 65]
VERIF_EXTENT = [90, 145, 25, 55]

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 11,
    'axes.labelsize': 9, 'figure.titlesize': 13,
})


# ═══════════════════════════════════════════════════════════════════════════
#  Lambert projection
# ═══════════════════════════════════════════════════════════════════════════

def get_lambert_projection():
    if not HAS_CARTOPY:
        return None
    return ccrs.LambertConformal(
        central_longitude=105,
        central_latitude=35,
        standard_parallels=(25, 47),
    )


LAMBERT_PROJ = get_lambert_projection()


# ═══════════════════════════════════════════════════════════════════════════
#  Rectangular Lambert plotting domain
# ═══════════════════════════════════════════════════════════════════════════

def build_lambert_rectangular_domain(center_lon=115.0, center_lat=38.0,
                                      width_m=8.0e6, height_m=5.0e6,
                                      nx=240, ny=150):
    """
    Build a regular rectangular Lambert x-y plotting domain.

    The domain is centred at (center_lon, center_lat) in PlateCarree,
    with given width and height in metres.  Returns the Lambert X, Y
    grids, their inverse-transformed lon/lat, and the x/y axis limits.

    Returns
    -------
    X, Y : ndarray (ny, nx) — Lambert coordinates (m)
    lon_target, lat_target : ndarray (ny, nx)
    xlim : (x_min, x_max)
    ylim : (y_min, y_max)
    """
    lambert = get_lambert_projection()
    pc = ccrs.PlateCarree()

    centre_xy = lambert.transform_point(center_lon, center_lat, pc)
    x0 = centre_xy[0] - width_m / 2.0
    x1 = centre_xy[0] + width_m / 2.0
    y0 = centre_xy[1] - height_m / 2.0
    y1 = centre_xy[1] + height_m / 2.0

    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    X, Y = np.meshgrid(x, y)

    lonlat = pc.transform_points(lambert, X, Y)
    lon_target = lonlat[..., 0]
    lat_target = lonlat[..., 1]

    return X, Y, lon_target, lat_target, (x0, x1), (y0, y1)


def interpolate_to_lambert_rectangle(field, lon2d, lat2d, lon_target, lat_target):
    """
    Interpolate a lat-lon source field to target lon/lat points
    (which define a rectangular Lambert domain).

    Uses bilinear interpolation on the regular source lat-lon grid.
    No footprint mask — the full rectangle is returned.
    """
    return _interp_regular_latlon(field, lon2d, lat2d, lon_target, lat_target)


def _interp_regular_latlon(field, lon2d, lat2d, lon_target, lat_target):
    """Bilinear interpolation from a regular lat-lon grid to target points."""
    lon = lon2d[0, :]
    lat = lat2d[:, 0]
    src = np.asarray(field)

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


def interpolate_fields_to_lambert_rectangle(fields_dict, lon2d, lat2d,
                                              lon_target, lat_target):
    """
    Interpolate multiple fields to the same target lon/lat points.
    Returns a dict of {name: field_lcc}.
    """
    result = {}
    for name, field in fields_dict.items():
        result[name] = _interp_regular_latlon(
            field, lon2d, lat2d, lon_target, lat_target
        )

    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Axes creation
# ═══════════════════════════════════════════════════════════════════════════

def set_lambert_curved_boundary(ax, proj, lon_min, lon_max, lat_min, lat_max, n=300):
    """
    Replace a Lambert GeoAxes frame with the projected lat-lon domain boundary.

    The boundary is constructed in geographic coordinates, transformed to the
    Lambert projection, then applied in projected data coordinates. This yields
    a curved quadrilateral with constant-longitude side edges and
    constant-latitude top/bottom edges.
    """
    bottom_lons = np.linspace(lon_min, lon_max, n)
    bottom_lats = np.full(n, lat_min)

    right_lons = np.full(n, lon_max)
    right_lats = np.linspace(lat_min, lat_max, n)

    top_lons = np.linspace(lon_max, lon_min, n)
    top_lats = np.full(n, lat_max)

    left_lons = np.full(n, lon_min)
    left_lats = np.linspace(lat_max, lat_min, n)

    lons = np.concatenate([bottom_lons, right_lons, top_lons, left_lons])
    lats = np.concatenate([bottom_lats, right_lats, top_lats, left_lats])

    xy = proj.transform_points(ccrs.PlateCarree(), lons, lats)
    x = xy[:, 0]
    y = xy[:, 1]

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if len(x) < 4:
        return None

    verts = np.column_stack([x, y])
    codes = np.full(len(verts), Path.LINETO)
    codes[0] = Path.MOVETO

    verts = np.vstack([verts, verts[0]])
    codes = np.append(codes, Path.CLOSEPOLY)
    boundary_path = Path(verts, codes)

    ax.set_boundary(boundary_path, transform=ax.transData)
    if 'geo' in ax.spines:
        ax.spines['geo'].set_visible(False)

    patch = PathPatch(
        boundary_path,
        transform=ax.transData,
        facecolor='none',
        edgecolor='black',
        linewidth=1.2,
        zorder=20,
    )
    ax.add_patch(patch)

    pad_x = 0.0 * (x.max() - x.min())
    pad_y = 0.0 * (y.max() - y.min())
    ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
    ax.set_ylim(y.min() - pad_y, y.max() + pad_y)

    ax._lambert_boundary_path = boundary_path
    ax._lambert_boundary_patch = patch
    return boundary_path


def _clip_to_lambert_boundary(ax, artists=None):
    """Clip plotted artists to the curved Lambert boundary if one is set."""
    boundary_path = getattr(ax, '_lambert_boundary_path', None)
    if boundary_path is None:
        return
    clip_patch = PathPatch(boundary_path, transform=ax.transData)
    if artists is None:
        artists = list(ax.collections) + list(ax.patches)
    elif hasattr(artists, 'collections'):
        artists = list(artists.collections)
    elif not isinstance(artists, (list, tuple)):
        artists = [artists]
    for artist in artists:
        if artist is getattr(ax, '_lambert_boundary_patch', None):
            continue
        try:
            artist.set_clip_path(clip_patch)
        except Exception:
            pass


def make_lambert_rect_axes(fig, subplot_spec, xlim=None, ylim=None,
                           lon_min=None, lon_max=None, lat_min=None, lat_max=None,
                           proj=None):
    """
    Create a Cartopy GeoAxes with a curved Lambert lat-lon boundary.
    """
    proj = proj or LAMBERT_PROJ
    lon_min = MAP_EXTENT[0] if lon_min is None else lon_min
    lon_max = MAP_EXTENT[1] if lon_max is None else lon_max
    lat_min = MAP_EXTENT[2] if lat_min is None else lat_min
    lat_max = MAP_EXTENT[3] if lat_max is None else lat_max

    ax = fig.add_subplot(subplot_spec, projection=proj)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_facecolor('white')

    boundary_path = set_lambert_curved_boundary(
        ax, proj, lon_min, lon_max, lat_min, lat_max, n=300)

    coastline = ax.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor='#333333')
    borders = ax.add_feature(cfeature.BORDERS, linewidth=0.35, edgecolor='#999999')

    # Subtle unlabelled gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=0.4, color='gray', alpha=0.45,
                      linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    _draw_verification_box(ax)
    if boundary_path is not None:
        _clip_to_lambert_boundary(ax, [coastline, borders])
        _clip_to_lambert_boundary(ax)
    return ax


def _draw_verification_box(ax):
    import matplotlib.patches as mpatches
    lon0, lon1, lat0, lat1 = VERIF_EXTENT
    rect = mpatches.Rectangle(
        (lon0, lat0), lon1 - lon0, lat1 - lat0,
        linewidth=0.8, edgecolor='#333333', facecolor='none',
        linestyle='--', alpha=0.7, transform=ccrs.PlateCarree())
    ax.add_patch(rect)


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _save_and_close(fig, filepath, dpi=300):
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {filepath}")


def _symmetric_levels(data, interval, for_error=False):
    if for_error:
        vmax = np.ceil(np.nanmax(np.abs(data)) / interval) * interval
        if vmax < interval:
            vmax = interval
        return np.arange(-vmax, vmax + interval, interval)
    vmin = np.floor(np.nanmin(data) / interval) * interval
    vmax = np.ceil(np.nanmax(data) / interval) * interval
    return np.arange(vmin, vmax + interval, interval)


def _add_colorbar(fig, mappable, ax, label, orientation='horizontal',
                  pad=0.06, shrink=0.82):
    cbar = fig.colorbar(mappable, ax=ax, orientation=orientation,
                        pad=pad, shrink=shrink)
    cbar.set_label(label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    return cbar


# ═══════════════════════════════════════════════════════════════════════════
#  Plot a gridded field on Lambert axes
# ═══════════════════════════════════════════════════════════════════════════

def plot_lambert_gridded_field(ax, X, Y, field_lcc, levels, cmap='RdBu_r',
                                title='', contour=True, clabel=False):
    """Plot a field on a regular Lambert x-y grid."""
    cf = ax.contourf(X, Y, field_lcc, levels=levels, cmap=cmap,
                     transform=LAMBERT_PROJ, extend='both')
    _clip_to_lambert_boundary(ax, cf)
    if contour:
        cs = ax.contour(X, Y, field_lcc, levels=levels,
                        colors='k', linewidths=0.25,
                        transform=LAMBERT_PROJ)
        _clip_to_lambert_boundary(ax, cs)
        if clabel:
            ax.clabel(cs, fmt='%d', fontsize=6)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=6)
    return cf


def plot_lambert_latlon_field(ax, lon2d, lat2d, field, levels, cmap='RdBu_r',
                              title='', contour=True, clabel=False):
    """Plot a lon-lat field on Lambert axes and clip it to the curved boundary."""
    cf = ax.contourf(
        lon2d, lat2d, field, levels=levels, cmap=cmap,
        transform=ccrs.PlateCarree(), extend='both')
    _clip_to_lambert_boundary(ax, cf)
    if contour:
        cs = ax.contour(
            lon2d, lat2d, field, levels=levels,
            colors='k', linewidths=0.25,
            transform=ccrs.PlateCarree())
        _clip_to_lambert_boundary(ax, cs)
        if clabel:
            ax.clabel(cs, fmt='%d', fontsize=6)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=6)
    return cf


def _lonlat_extent(lon2d, lat2d):
    return (float(np.nanmin(lon2d)), float(np.nanmax(lon2d)),
            float(np.nanmin(lat2d)), float(np.nanmax(lat2d)))


# ═══════════════════════════════════════════════════════════════════════════
#  Figure-building functions
# ═══════════════════════════════════════════════════════════════════════════

# ── Shared Lambert domain (built once, reused) ──

_DOMAIN_CACHE = None  # (X, Y, lon_t, lat_t, xlim, ylim)


def _get_shared_domain(nx=240, ny=150):
    """Return the cached Lambert rectangular domain."""
    global _DOMAIN_CACHE
    if _DOMAIN_CACHE is None or _DOMAIN_CACHE[0].shape != (ny, nx):
        _DOMAIN_CACHE = build_lambert_rectangular_domain(
            center_lon=115.0, center_lat=42.0,
            width_m=7.5e6, height_m=4.4e6,
            nx=nx, ny=ny,
        )
    return _DOMAIN_CACHE


def plot_geopotential_height(Z, lon2d, lat2d, title, filepath,
                              contour_interval=60, cmap='RdYlBu_r'):
    """Plot 500 hPa geopotential height (full field)."""
    vmin = np.floor(np.nanmin(Z) / contour_interval) * contour_interval
    vmax = np.ceil(np.nanmax(Z) / contour_interval) * contour_interval
    levels = np.arange(vmin, vmax + contour_interval, contour_interval)

    if HAS_CARTOPY:
        lon_min, lon_max, lat_min, lat_max = _lonlat_extent(lon2d, lat2d)
        fig = plt.figure(figsize=(13, 7.5))
        gs = GridSpec(1, 1, figure=fig)
        ax = make_lambert_rect_axes(
            fig, gs[0, 0],
            lon_min=lon_min, lon_max=lon_max,
            lat_min=lat_min, lat_max=lat_max)

        cf = plot_lambert_latlon_field(
            ax, lon2d, lat2d, Z, levels, cmap=cmap,
            title=title, contour=True, clabel=True)
        _add_colorbar(fig, cf, ax, 'Geopotential height (m)')
    else:
        fig, ax = plt.subplots(figsize=(12, 7))
        cf = ax.contourf(lon2d, lat2d, Z, levels=levels, cmap=cmap,
                         extend='both')
        cs = ax.contour(lon2d, lat2d, Z, levels=levels, colors='k',
                        linewidths=0.5)
        ax.clabel(cs, fmt='%d', fontsize=7)
        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
        ax.set_xlim(MAP_EXTENT[0], MAP_EXTENT[1])
        ax.set_ylim(MAP_EXTENT[2], MAP_EXTENT[3])
        _add_colorbar(fig, cf, ax, 'Geopotential height (m)')
        ax.set_title(title, fontsize=12, fontweight='bold')

    _save_and_close(fig, filepath)


def make_verification_figure(zf_anom, za_anom, lon2d, lat2d, title_prefix,
                              filepath, anomaly_interval=30,
                              error_interval=15):
    """Build a 3-panel verification figure with Lambert maps."""
    error_field = zf_anom - za_anom

    if not HAS_CARTOPY:
        _make_verification_figure_nocartopy(
            zf_anom, za_anom, error_field, lon2d, lat2d,
            title_prefix, filepath, anomaly_interval, error_interval)
        return

    # Symmetric anomaly levels
    max_abs = max(np.nanmax(np.abs(zf_anom)),
                  np.nanmax(np.abs(za_anom)))
    max_anom = np.ceil(max_abs / anomaly_interval) * anomaly_interval
    if max_anom < anomaly_interval:
        max_anom = anomaly_interval
    anomaly_levels = np.arange(-max_anom, max_anom + anomaly_interval,
                                anomaly_interval)

    max_abs_err = np.nanmax(np.abs(error_field))
    max_err = np.ceil(max_abs_err / error_interval) * error_interval
    if max_err < error_interval:
        max_err = error_interval
    error_levels = np.arange(-max_err, max_err + error_interval, error_interval)

    fig = plt.figure(figsize=(18, 5.8))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[10, 1],
                  hspace=0.15, wspace=0.12,
                  left=0.03, right=0.97, bottom=0.12, top=0.88)

    titles = [f'{title_prefix} Forecast anomaly',
              f'{title_prefix} Analysis anomaly',
              f'{title_prefix} Error (fcst − anal)']
    fields = [zf_anom, za_anom, error_field]
    level_lists = [anomaly_levels, anomaly_levels, error_levels]
    cbar_labels = ['Height anomaly (m)', 'Height anomaly (m)', 'Error (m)']
    lon_min, lon_max, lat_min, lat_max = _lonlat_extent(lon2d, lat2d)

    cfs = []
    for idx in range(3):
        ax = make_lambert_rect_axes(
            fig, gs[0, idx],
            lon_min=lon_min, lon_max=lon_max,
            lat_min=lat_min, lat_max=lat_max)
        cf = plot_lambert_latlon_field(
            ax, lon2d, lat2d, fields[idx],
            levels=level_lists[idx], cmap='RdBu_r',
            title=titles[idx], contour=True, clabel=False)
        cfs.append(cf)

    for idx in range(3):
        cax = fig.add_subplot(gs[1, idx])
        cbar = fig.colorbar(cfs[idx], cax=cax, orientation='horizontal')
        cbar.set_label(cbar_labels[idx], fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        cax.set_facecolor('none')
        for s in cax.spines.values():
            s.set_visible(False)
        cax.set_xticks([]); cax.set_yticks([])

    fig.suptitle(f'{title_prefix} BVE Forecast Verification',
                 fontsize=14, fontweight='bold', y=0.96)
    _save_and_close(fig, filepath)


def _make_verification_figure_nocartopy(zf_anom, za_anom, error_field,
                                         lon2d, lat2d, title_prefix,
                                         filepath, anomaly_interval,
                                         error_interval):
    anomaly_levels = _symmetric_levels(za_anom, anomaly_interval, for_error=False)
    error_levels = _symmetric_levels(error_field, error_interval, for_error=True)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    titles = [f'{title_prefix} Forecast anomaly',
              f'{title_prefix} Analysis anomaly',
              f'{title_prefix} Error (fcst − anal)']
    fields = [zf_anom, za_anom, error_field]
    level_lists = [anomaly_levels, anomaly_levels, error_levels]
    cbar_labels = ['Height anomaly (m)', 'Height anomaly (m)', 'Error (m)']
    for idx, ax in enumerate(axes):
        cf = ax.contourf(lon2d, lat2d, fields[idx], levels=level_lists[idx],
                         cmap='RdBu_r', extend='both')
        ax.contour(lon2d, lat2d, fields[idx], levels=level_lists[idx],
                   colors='k', linewidths=0.3)
        ax.set_title(titles[idx], fontsize=11, fontweight='bold')
        ax.set_xlabel('Longitude (°E)', fontsize=8)
        if idx == 0:
            ax.set_ylabel('Latitude (°N)', fontsize=8)
        ax.set_xlim(MAP_EXTENT[0], MAP_EXTENT[1])
        ax.set_ylim(MAP_EXTENT[2], MAP_EXTENT[3])
        _add_colorbar(fig, cf, ax, cbar_labels[idx])
    fig.suptitle(f'{title_prefix} BVE Forecast Verification',
                 fontsize=14, fontweight='bold')
    _save_and_close(fig, filepath)


def plot_forecast_verification(zf_anom, za_anom, lon2d, lat2d, title_prefix,
                                filepath, contour_interval=30,
                                error_interval=15):
    make_verification_figure(zf_anom, za_anom, lon2d, lat2d, title_prefix,
                              filepath, anomaly_interval=contour_interval,
                              error_interval=error_interval)


def plot_vorticity(zeta, lon2d, lat2d, title, filepath, cmap='RdBu_r'):
    """Plot relative vorticity on a Lambert map."""
    zeta_plot = zeta * 1e5
    vmax = np.nanmax(np.abs(zeta_plot))
    vmax = np.ceil(vmax * 2) / 2
    if vmax < 1:
        vmax = 2
    levels = np.linspace(-vmax, vmax, 21)

    if HAS_CARTOPY:
        lon_min, lon_max, lat_min, lat_max = _lonlat_extent(lon2d, lat2d)
        fig = plt.figure(figsize=(13, 7.5))
        gs = GridSpec(1, 1, figure=fig)
        ax = make_lambert_rect_axes(
            fig, gs[0, 0],
            lon_min=lon_min, lon_max=lon_max,
            lat_min=lat_min, lat_max=lat_max)

        cf = plot_lambert_latlon_field(
            ax, lon2d, lat2d, zeta_plot, levels, cmap=cmap,
            title=title, contour=True, clabel=False)
        _add_colorbar(fig, cf, ax, 'Relative vorticity (10⁻⁵ s⁻¹)')
    else:
        fig, ax = plt.subplots(figsize=(12, 7))
        cf = ax.contourf(lon2d, lat2d, zeta_plot, levels=levels, cmap=cmap,
                         extend='both')
        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
        ax.set_xlim(MAP_EXTENT[0], MAP_EXTENT[1])
        ax.set_ylim(MAP_EXTENT[2], MAP_EXTENT[3])
        _add_colorbar(fig, cf, ax, 'Relative vorticity (10⁻⁵ s⁻¹)')
        ax.set_title(title, fontsize=12, fontweight='bold')

    _save_and_close(fig, filepath)


def plot_scores(scores_12h, scores_24h, filepath):
    """Bar chart and summary table of verification scores."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    metrics_rmse = ['height_rmse_m', 'height_debiased_rmse_m', 'vorticity_rmse']
    labels_rmse = ['Height Anom.\nRMSE (m)', 'Debiased Height\nRMSE (m)',
                   'Vorticity\nRMSE (10⁻⁵ s⁻¹)']
    vals_12 = [scores_12h.get(m, scores_12h.get('vorticity_rmse_s-1', 0))
               if m == 'vorticity_rmse' else scores_12h.get(m, 0)
               for m in metrics_rmse]
    vals_24 = [scores_24h.get(m, scores_24h.get('vorticity_rmse_s-1', 0))
               if m == 'vorticity_rmse' else scores_24h.get(m, 0)
               for m in metrics_rmse]
    vals_12[2] *= 1e5; vals_24[2] *= 1e5
    x = np.arange(len(labels_rmse)); w = 0.35
    axes[0].bar(x - w/2, vals_12, w, label='+12 h', color='#4472C4',
                edgecolor='k', linewidth=0.5)
    axes[0].bar(x + w/2, vals_24, w, label='+24 h', color='#ED7D31',
                edgecolor='k', linewidth=0.5)
    axes[0].set_xticks(x); axes[0].set_xticklabels(labels_rmse, fontsize=9)
    axes[0].set_ylabel('RMSE'); axes[0].set_title('RMSE Scores', fontweight='bold', fontsize=12)
    axes[0].legend(fontsize=8); axes[0].grid(axis='y', alpha=0.3)

    metrics_acc = ['height_acc', 'vorticity_acc']
    labels_acc = ['Height Anom.\nCorrelation', 'Vorticity\nCorrelation']
    vals_12_acc = [scores_12h.get(m, 0) for m in metrics_acc]
    vals_24_acc = [scores_24h.get(m, 0) for m in metrics_acc]
    x_acc = np.arange(len(labels_acc))
    axes[1].bar(x_acc - w/2, vals_12_acc, w, label='+12 h', color='#4472C4',
                edgecolor='k', linewidth=0.5)
    axes[1].bar(x_acc + w/2, vals_24_acc, w, label='+24 h', color='#ED7D31',
                edgecolor='k', linewidth=0.5)
    axes[1].set_xticks(x_acc); axes[1].set_xticklabels(labels_acc, fontsize=9)
    axes[1].set_ylabel('Correlation'); axes[1].set_ylim(0, 1.05)
    axes[1].set_title('Pattern Correlation', fontweight='bold', fontsize=12)
    axes[1].legend(fontsize=8); axes[1].grid(axis='y', alpha=0.3)

    axes[2].axis('off')
    table_data = [
        ['Metric', '+12 h', '+24 h'],
        ['Height anom. RMSE (m)', f"{scores_12h.get('height_rmse_m',0):.1f}", f"{scores_24h.get('height_rmse_m',0):.1f}"],
        ['Debiased RMSE (m)', f"{scores_12h.get('height_debiased_rmse_m',0):.1f}", f"{scores_24h.get('height_debiased_rmse_m',0):.1f}"],
        ['Height anom. corr.', f"{scores_12h.get('height_acc',0):.4f}", f"{scores_24h.get('height_acc',0):.4f}"],
        ['Vorticity corr.', f"{scores_12h.get('vorticity_acc',0):.4f}", f"{scores_24h.get('vorticity_acc',0):.4f}"],
        ['Height bias (m)', f"{scores_12h.get('height_bias_m',0):.1f}", f"{scores_24h.get('height_bias_m',0):.1f}"],
    ]
    tab = axes[2].table(cellText=table_data, cellLoc='center',
                        loc='center', colWidths=[0.32, 0.18, 0.18])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.2, 1.6)
    for j in range(3):
        tab[0, j].set_facecolor('#4472C4')
        tab[0, j].set_text_props(color='white', fontweight='bold')
    axes[2].set_title('Score Summary', fontweight='bold', fontsize=12, y=0.98)
    fig.suptitle('BVE Forecast Verification Scores', fontsize=14, fontweight='bold')
    _save_and_close(fig, filepath)


def plot_dummy(lon2d, lat2d, filepath):
    ny, nx = lon2d.shape
    Z = 5500 + 200 * np.sin(2*np.pi*np.arange(ny)[:,None]/ny) * \
                np.cos(3*np.pi*np.arange(nx)[None,:]/nx)
    plot_geopotential_height(Z, lon2d, lat2d,
                             'Test: Dummy 500 hPa Geopotential Height', filepath)
    return True
