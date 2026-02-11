import matplotlib
import numpy as np
import cartopy.crs as ccrs
import cartopy.util as cutil
from .grid_utils import get_lon_name, get_lat_name


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    if isinstance(cmap, str):
        cmap = matplotlib.cm.get_cmap(cmap)
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_contour_map(ax, data, **kwargs):
    ax.coastlines(lw=.8)
    cdata, clon, clat = cutil.add_cyclic(data, data[get_lon_name(data)], data[get_lat_name(data)])
    cont = ax.contour(clon, clat, cdata, **kwargs)
    return cont


def plot_contourf_map(ax, data, **kwargs):
    ax.coastlines(lw=.8)
    cdata, clon, clat = cutil.add_cyclic(data, data[get_lon_name(data)], data[get_lat_name(data)])
    cont = ax.contourf(clon, clat, cdata, **kwargs)
    cont.set_edgecolor('face')
    return cont


def plot_pcolormesh_map(ax, data, **kwargs):
    ax.coastlines(lw=.8)
    cdata, clon, clat = cutil.add_cyclic(data, data[get_lon_name(data)], data[get_lat_name(data)])
    mesh = ax.pcolormesh(clon, clat, cdata, **kwargs)
    # cont.set_edgecolor('face')
    return mesh
