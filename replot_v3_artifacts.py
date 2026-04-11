import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import flopy
import flopy.utils.binaryfile as bf
import json

base_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(base_dir, "output")
mf6_dir = os.path.join(base_dir, "mf6_workspace")
mp7_dir = os.path.join(base_dir, "mp7_workspace")

# Load geom
with open(os.path.join(output_dir, "modflow_domain.json"), 'r') as f:
    domain = json.load(f)
nrow, ncol = domain['nrow'], domain['ncol']
xmin, xmax, ymin, ymax = domain['xmin'], domain['xmax'], domain['ymin'], domain['ymax']

top = np.loadtxt(os.path.join(output_dir, "model_top.csv"), delimiter=',')
ibound = np.loadtxt(os.path.join(output_dir, "model_ibound.csv"), delimiter=',')
botm = np.zeros((4, nrow, ncol))
botm[0] = top - 50.0
botm[1] = botm[0] - 50.0
botm[2] = botm[1] - 150.0
botm[3] = botm[2] - 450.0

# Load Heads
head_file = os.path.join(mf6_dir, "MF_FIJI_v3.hds")
hds = bf.HeadFile(head_file)
heads = hds.get_data(totim=hds.get_times()[-1])

# Load Obs Data
river_obs = pd.read_csv(os.path.join(output_dir, "river_obs_v3.csv"))
gw_obs = pd.read_csv(os.path.join(output_dir, "gw_obs_v3.csv"))

river_gdf = gpd.GeoDataFrame(river_obs, geometry=gpd.points_from_xy(river_obs.lon, river_obs.lat), crs='EPSG:4326').to_crs('EPSG:32760')
gw_gdf = gpd.GeoDataFrame(gw_obs, geometry=gpd.points_from_xy(gw_obs.lon, gw_obs.lat), crs='EPSG:4326').to_crs('EPSG:32760')

fig, ax = plt.subplots(figsize=(14, 10))

# Plan View Plotting
extent = [xmin, xmax, ymin, ymax]
heads_ma = np.ma.masked_where((heads[0] > 1e10) | (heads[0] < -100) | (ibound == 0), heads[0])

# To ensure Flopy works properly, create a dummy model for spatial reference
gwf = flopy.mf6.ModflowGwf(flopy.mf6.MFSimulation(), modelname='dummy')
flopy.mf6.ModflowGwfdis(gwf, nlay=4, nrow=nrow, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm, idomain=np.broadcast_to(ibound, (4, nrow, ncol)), xorigin=xmin, yorigin=ymin)

pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
im = pmv.plot_array(heads_ma, cmap='viridis', alpha=0.8, vmin=0, vmax=150)

# Plot Head Contours!
contours = pmv.contour_array(heads[0], masked_values=[1e30], levels=np.arange(0, 160, 10), colors='white', linewidths=0.5, alpha=0.8)
plt.clabel(contours, fmt='%.0f', colors='white', fontsize=8)

p_file = os.path.join(mp7_dir, "MF_FIJI_v3_mp7.mppth")
if os.path.exists(p_file):
    plines = flopy.utils.PathlineFile(p_file).get_alldata()
    pmv.plot_pathline(plines, layer='all', colors='white', lw=1.2, alpha=0.8, label='Pathlines')

# Plot Points MANUALLY to guarantee they are HUGE
if not river_gdf.empty:
    ax.scatter(river_gdf.geometry.x, river_gdf.geometry.y, c='red', s=350, marker='s', edgecolors='black', linewidths=2.0, zorder=10, label='River Obs (FR)')
if not gw_gdf.empty:
    ax.scatter(gw_gdf.geometry.x, gw_gdf.geometry.y, c='aqua', s=350, marker='^', edgecolors='black', linewidths=2.0, zorder=10, label='GW Obs (FG)')

ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])

ax.set_title("MODFLOW 6 v3 - Layer 1 Hydraulic Head & Pathlines", fontsize=18)
ax.set_xlabel("UTM Easting (m)", fontsize=14)
ax.set_ylabel("UTM Northing (m)", fontsize=14)
plt.colorbar(im, ax=ax, label="Head (m)")
ax.legend(loc='lower right', fontsize=14)

plt.savefig(os.path.join(output_dir, 'modflow_v3_plan_view.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Plan view generated.")

# Cross Section
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
col_idx = ncol // 2
xs = np.arange(nrow) * 100.0 + ymin

# Dem is index 0 for top, but xs corresponds to y traversing rows!
# Let's plot y vs z instead of x vs z, since column specifies X, we are changing Y (Northing).
# row 0 is ymax, row -1 is ymin. So xs needs to map to actual Y coordinates:
ys = ymax - np.arange(nrow) * 100.0 - 50.0 # Cell centers

ax1.plot(ys, top[:, col_idx], 'k-', linewidth=3, label='DEM')
for l in range(4):
    ax1.plot(ys, botm[l, :, col_idx], 'k--', linewidth=1, alpha=0.5)

colors = ['blue', 'cyan', 'green', 'goldenrod']
for l in range(4):
    heads_ma_xs = np.ma.masked_where((heads[l, :, col_idx] > 1e10) | (heads[l, :, col_idx] < -100) | (ibound[:, col_idx] == 0), heads[l,:,col_idx])
    ax1.plot(ys, heads_ma_xs, label=f'L{l+1} Head', linewidth=2.5, color=colors[l])

ax1.set_ylim([-600, np.nanmax(top[:, col_idx]) + 50])
ax1.set_ylabel("Elevation (m)", fontsize=14)
ax1.set_title(f"Cross-section (Northing profile) at Column {col_idx} (Easting ~ {xmin + col_idx*100}m)", fontsize=16)
ax1.legend(loc='upper right', fontsize=12)
ax1.grid(True)

# Just copy the array K plot so it's not empty
# Rough K mock for brevity (or load actual K from NPF, but it's hard without NPF created)
ax2.plot(ys, np.where(ys < 0, 0, 2), label='K (Dummy)', linewidth=0)
ax2.set_xlabel("Northing (m)", fontsize=14)
ax2.set_visible(False) # Hide second plot to just show the clean cross-section!

plt.savefig(os.path.join(output_dir, 'modflow_v3_cross_section.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Cross-section generated.")
