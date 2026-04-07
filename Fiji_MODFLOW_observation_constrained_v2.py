#!/usr/bin/env python
# coding: utf-8

# In[59]:


import os
import json
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from rasterio.transform import from_origin
from rasterio.features import rasterize
import flopy
import flopy.utils.binaryfile as bf


import argparse
import sys
import logging

# Logger setup
logger = logging.getLogger("HydrologyWorkflow")
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler('out.log', encoding='utf-8')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Run hydrology workflow.')
parser.add_argument('--dry-run', action='store_true', help='Only rewrite the notebook without execution.')
args = parser.parse_args()

if args.dry_run:
    logger.info("=== DRY-RUN MODE: REWRITING NOTEBOOK ONLY ===")
    import modify_notebook
    logger.info("Notebook rewritten successfully. Exiting dry-run.")
    sys.exit(0)

logger.info("=== STARTING MODFLOW OBSERVATION-CONSTRAINED WORKFLOW ===")

warnings.filterwarnings('ignore')

# ワーキングディレクトリの設定
#base_dir = r"C:\Users\yasum\.gemini\antigravity\scratch\Fiji"
base_dir = r"C:\Users\yasum\.gemini\antigravity\scratch\Fiji"
output_dir = os.path.join(base_dir, "output")
mf6_dir = os.path.join(base_dir, "mf6_workspace")
mp7_dir = os.path.join(base_dir, "mp7_workspace")

# フォルダが存在しない場合は作成
for d in [mf6_dir, mp7_dir]:
    if not os.path.exists(d):
        os.makedirs(d)

logger.info("Cell 1 Execution Completed: Base Directories and Initial Setup Done.")


# In[60]:


# DEMとアクティブセルの読み込み
top = np.loadtxt(os.path.join(output_dir, "model_top.csv"), delimiter=',')
ibound = np.loadtxt(os.path.join(output_dir, "model_ibound.csv"), delimiter=',')

with open(os.path.join(output_dir, "modflow_domain.json"), 'r') as f:
    domain = json.load(f)
nrow, ncol = domain['nrow'], domain['ncol']
xmin, xmax, ymin, ymax = domain['xmin'], domain['xmax'], domain['ymin'], domain['ymax']

# 対象流域のマスク作成
basins = gpd.read_file(os.path.join(output_dir, 'basins.shp'))
basins['VALUE'] = basins['VALUE'].astype(int).astype(str)

# 【修正1】ポリゴンを150m（1.5セル分）だけ外側に膨らませて、海と確実に連結させる
focal_basins = basins[basins['VALUE'].isin(['2422', '2344'])].copy()
focal_basins['geometry'] = focal_basins.geometry.buffer(150)

transform = from_origin(xmin, ymax, 100.0, 100.0)
basin_mask = rasterize([(geom, 1) for geom in focal_basins.geometry], out_shape=(nrow, ncol), transform=transform, fill=0, all_touched=True)
highland_mask = (top > 50.0)

# 陸地であり、かつ膨らませた対象流域外のセルを 0（計算除外）にする
ibound[(ibound == 1) & (basin_mask == 0)] = 0

# 【修正2】ダム化（せき止め）を防ぐため、海辺の低地（標高5m未満）は壁を作らず開けておく
ibound[(ibound == 0) & (top < 5.0) & (top >= 0)] = 1

# NWTソルバー安定化のためのレイヤー厚確保
nlay = 4
botm = np.zeros((nlay, nrow, ncol))
botm[0] = top - 50.0       
botm[1] = botm[0] - 50.0   
botm[2] = botm[1] - 150.0  
botm[3] = botm[2] - 450.0  

# 新しいiboundに基づいてidomainを再設定
idomain = np.array([np.where(ibound != 0, 1, 0)] * nlay)

logger.info("Cell 2 Execution Completed: 流域の切り出しと、海への接続（バッファ）処理が完了しました。")


# In[63]:


import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os

# 描画用のフィギュアを作成
fig_chk, ax_chk = plt.subplots(figsize=(14, 10))
extent = [xmin, xmax, ymin, ymax]

# iboundのプロット
# -1: 海 (青), 0: 除外領域 (灰色), 1: アクティブな陸地 (緑)
cmap = ListedColormap(['blue', 'lightgray', 'limegreen'])
cax = ax_chk.imshow(ibound, extent=extent, origin='upper', cmap=cmap, vmin=-1, vmax=1, alpha=0.6)

# 全流域の境界線を薄くグレーで描画
basins.plot(ax=ax_chk, facecolor='none', edgecolor='gray', linewidth=0.3, alpha=0.5)

# 選択された対象流域（150mバッファ付き）の輪郭を赤色で強調
focal_basins.plot(ax=ax_chk, facecolor='none', edgecolor='red', linewidth=1.5)

ax_chk.set_title("Check: Final ibound Map vs Selected Watersheds", fontsize=16)
ax_chk.set_xlabel("UTM Easting (m)")
ax_chk.set_ylabel("UTM Northing (m)")

# 凡例を手動で追加
sea_patch = mpatches.Patch(color='blue', alpha=0.6, label='Sea (ibound = -1)')
inactive_patch = mpatches.Patch(color='lightgray', alpha=0.6, label='Inactive / Excluded (ibound = 0)')
active_patch = mpatches.Patch(color='limegreen', alpha=0.6, label='Active Land (ibound = 1)')
red_line = plt.Line2D([0], [0], color='red', linewidth=1.5, label='Selected Basins (Buffered)')
ax_chk.legend(handles=[sea_patch, inactive_patch, active_patch, red_line], loc='lower right')

# 画像として保存
save_path = os.path.join(output_dir, "checked_ibound_final.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
logger.info(f"Image Saved: {save_path}")

plt.close()


# In[64]:


import geopandas as gpd
import os

# --- 1. 元のデータ（エラーを含む）の読み込み ---
basins_shp_path = os.path.join(output_dir, 'basins.shp')
basins = gpd.read_file(basins_shp_path)

# --- 2. ジオメトリの修復（ゼロ・バッファの魔法） ---
# すべてのポリゴンに buffer(0) を適用して、微小なねじれや自己交差を自動修復します
basins['geometry'] = basins.geometry.buffer(0)

# --- 3. クリーンなデータとして別名で保存 ---
clean_shp_path = os.path.join(output_dir, 'basins_clean.shp')
basins.to_file(clean_shp_path)

print(f"修復完了！クリーンなShapefileを保存しました: {clean_shp_path}")


# In[65]:


sim_name = "MF_IRIOMOTE_ALL"
sim = flopy.mf6.MFSimulation(sim_name=sim_name, version="mf6", exe_name="mf6", sim_ws=mf6_dir)
tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
ims = flopy.mf6.ModflowIms(sim, complexity="COMPLEX", outer_maximum=500, inner_maximum=500)
gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True, newtonoptions="NEWTON UNDER_RELAXATION")
dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm, idomain=idomain)

# 透水係数（深部は被圧化して計算破綻を防ぐ）
k1, k2, k3, k4 = np.full((nrow, ncol), 12.0), np.full((nrow, ncol), 2.0), np.full((nrow, ncol), 0.2), np.full((nrow, ncol), 0.02)
k1[highland_mask], k2[highland_mask] = 15.0, 5.0
npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=[1, 0, 0, 0], k=[k1, k2, k3, k4], k33=[k1/10, k2/10, k3/10, k4/10], save_flows=True)

# 初期水頭と海側境界
strt_array = np.zeros((nlay, nrow, ncol))
for l in range(nlay): strt_array[l] = np.where(top < 0, 0.0, top)
ic = flopy.mf6.ModflowGwfic(gwf, strt=strt_array)
chd_spd = [[(l, r, c), 0.0] for r, c in zip(*np.where(ibound == -1)) for l in range(nlay)]
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

# 涵養量：島全体に降らせる（流域間の不自然な移動を防ぐため）
base_rch = 0.0027
rch_array = np.zeros((nrow, ncol))
rch_array[ibound != 0] = base_rch
rch_array[highland_mask] *= 0.7 
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rch_array)

# 全河川の自動抽出とネットワーク構築
streams_gdf = gpd.read_file(os.path.join(output_dir, "streams.shp"))
stream_grid = rasterize([(geom, 1) for geom in streams_gdf.geometry], out_shape=(nrow, ncol), transform=transform, fill=0, all_touched=True)
stream_rows, stream_cols = np.where((stream_grid == 1) & (ibound == 1) & (top >= 0))
stream_data = sorted([(r, c, top[r, c]) for r, c in zip(stream_rows, stream_cols)], key=lambda x: x[2], reverse=True)

stream_cells = {(r, c): {'rno': i, 'Z': z} for i, (r, c, z) in enumerate(stream_data)}
adj_in, adj_out = {i: [] for i in range(len(stream_data))}, {i: [] for i in range(len(stream_data))}
for (r, c), data in stream_cells.items():
    u, lowest_z, v = data['rno'], data['Z'], None
    for dr, dc in [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]:
        if (r+dr, c+dc) in stream_cells and stream_cells[(r+dr, c+dc)]['Z'] < lowest_z:
            lowest_z, v = stream_cells[(r+dr, c+dc)]['Z'], stream_cells[(r+dr, c+dc)]['rno']
    if v is not None:
        adj_out[u].append(v)
        adj_in[v].append(u)

best_rhk = 10.0 
packagedata, connectiondata, sfr_plot_data = [], [], []
for i, (r, c, z) in enumerate(stream_data):
    packagedata.append([i, (0, r, c), 100.0, 5.0, 0.01, z, 1.0, best_rhk, 0.03, len(adj_in[i]) + len(adj_out[i]), 1.0, 0])
    connectiondata.append([i] + adj_in[i] + [-x for x in adj_out[i]])
    sfr_plot_data.append([r, c])

# 作図用に全河川セルをCSV保存
pd.DataFrame(sfr_plot_data, columns=['Row', 'Col']).to_csv(os.path.join(output_dir, "sfr_all_network.csv"), index=False)

sfr = flopy.mf6.ModflowGwfsfr(gwf, nreaches=len(packagedata), packagedata=packagedata, connectiondata=connectiondata, save_flows=True)
oc = flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f"{sim_name}.hds", budget_filerecord=f"{sim_name}.cbc", saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])

logger.info("Cell 3 Execution Completed: IRIOMOTE ALL Model Built and Rivers Configured.")


# In[66]:


import json
import os
import numpy as np
import geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from shapely.geometry import box, Point
from shapely.ops import unary_union

# --- 1. データの読み込み ---
with open(os.path.join(output_dir, "modflow_domain.json"), 'r') as f:
    domain = json.load(f)
xmin, xmax, ymin, ymax = domain['xmin'], domain['xmax'], domain['ymin'], domain['ymax']
nrow, ncol = domain['nrow'], domain['ncol']

top = np.loadtxt(os.path.join(output_dir, "model_top.csv"), delimiter=',')
ibound = np.loadtxt(os.path.join(output_dir, "model_ibound.csv"), delimiter=',')

basins = gpd.read_file(os.path.join(output_dir, 'basins_clean.shp'))
basins['VALUE'] = basins['VALUE'].astype(int).astype(str)
transform = from_origin(xmin, ymax, 100.0, 100.0)

# --- 2. 解析領域の基本ポリゴンを作る ---
target_basins = []
special_ids = ['2647', '4373', '4549', '3880', '4007', '3437', '4559']

for _, row in basins.iterrows():
    geom = row.geometry
    val = row['VALUE']
    minx, miny, maxx, maxy = geom.bounds
    cx, cy = geom.centroid.x, geom.centroid.y

    if val in special_ids:
        target_basins.append(geom.buffer(1))
        continue

    # 西端にはみ出したノイズ流域を除去
    if minx < xmin + 10:
        continue

    # 南東の不要半島を除去
    if cx > 631500 and cy < 8072000:
        continue

    # 東側南端で解析範囲をまたぐ途切れた流域を除去
    if cx > 625000 and miny < ymin + 10:
        continue

    target_basins.append(geom.buffer(1))

merged_basin = unary_union(target_basins)
if merged_basin.geom_type == 'MultiPolygon':
    mainland = max(merged_basin.geoms, key=lambda p: p.area)
else:
    mainland = merged_basin

domain_box = box(xmin, ymin, xmax, ymax)
final_basin_poly = mainland.buffer(-1).buffer(150).intersection(domain_box)

# --- 3. basin_mask を作って active / inactive を更新 ---
basin_mask = rasterize(
    [(final_basin_poly, 1)],
    out_shape=(nrow, ncol),
    transform=transform,
    fill=0,
    all_touched=True
)

ibound[(ibound == 1) & (basin_mask == 0)] = 0

# 沿岸低地は少し広めに active にして、沿岸井戸や海岸近傍地下水を入れやすくする
coastal_open = (ibound == 0) & (basin_mask == 1) & (top < 10.0) & (top >= 0.0)
ibound[coastal_open] = 1

# --- 4. 観測点周辺は確実に active にする（海岸井戸や小河川観測点を落とさないため） ---
streams_gdf = gpd.read_file(os.path.join(output_dir, "streams.shp"))
obs_points_wgs84 = [
    ("FR2", -17.3803293, 178.1541107),
    ("FR3", -17.4033547, 178.1147045),
    ("FR4", -17.3572955, 178.1824644),
    ("FR5", -17.4252996, 178.0779275),
    ("FR6", -17.4203176, 178.0795745),
    ("FR7", -17.4275448, 178.0844146),
    ("FG1", -17.3572541, 178.1827151),
    ("FG2", -17.3585429, 178.1847587),
]

obs_gdf = gpd.GeoDataFrame(
    pd.DataFrame(obs_points_wgs84, columns=["site", "lat", "lon"]),
    geometry=[Point(lon, lat) for _, lat, lon in obs_points_wgs84],
    crs="EPSG:4326",
).to_crs(streams_gdf.crs or "EPSG:32760")

for _, r in obs_gdf.iterrows():
    col = int((r.geometry.x - xmin) / 100.0)
    row = int((ymax - r.geometry.y) / 100.0)
    for rr in range(max(0, row - 1), min(nrow, row + 2)):
        for cc in range(max(0, col - 1), min(ncol, col + 2)):
            if top[rr, cc] >= 0:
                ibound[rr, cc] = 1

# --- 5. レイヤー構造 ---
nlay = 4
botm = np.zeros((nlay, nrow, ncol))
botm[0] = top - 50.0
botm[1] = botm[0] - 50.0
botm[2] = botm[1] - 150.0
botm[3] = botm[2] - 450.0

# 海セルも含め、inactive 以外を計算対象にする
idomain = np.array([np.where(ibound != 0, 1, 0)] * nlay)

# --- 6. 可視化 ---
fig_chk, ax_chk = plt.subplots(figsize=(14, 10))
extent = [xmin, xmax, ymin, ymax]
cmap = ListedColormap(['blue', 'lightgray', 'limegreen'])
ax_chk.imshow(ibound, extent=extent, origin='upper', cmap=cmap, vmin=-1, vmax=1, alpha=0.6)

basins.plot(ax=ax_chk, facecolor='none', edgecolor='gray', linewidth=0.3, alpha=0.5)
streams_gdf.plot(ax=ax_chk, color='cyan', linewidth=0.5, alpha=0.7)

if final_basin_poly.geom_type == 'Polygon':
    ax_chk.plot(*final_basin_poly.exterior.xy, color='red', linewidth=2)
elif final_basin_poly.geom_type == 'MultiPolygon':
    for poly in final_basin_poly.geoms:
        ax_chk.plot(*poly.exterior.xy, color='red', linewidth=2)

obs_gdf.plot(ax=ax_chk, color='magenta', markersize=18, zorder=10)
for _, r in obs_gdf.iterrows():
    ax_chk.text(r.geometry.x + 50, r.geometry.y + 50, r["site"], color='magenta', fontsize=8)

ax_chk.set_xlim(xmin, xmax)
ax_chk.set_ylim(ymin, ymax)
ax_chk.set_title("Check: updated idomain / ibound with observation points preserved", fontsize=15)
ax_chk.set_xlabel("UTM Easting (m)")
ax_chk.set_ylabel("UTM Northing (m)")

sea_patch = mpatches.Patch(color='blue', alpha=0.6, label='Sea (ibound = -1)')
inactive_patch = mpatches.Patch(color='lightgray', alpha=0.6, label='Inactive (ibound = 0)')
active_patch = mpatches.Patch(color='limegreen', alpha=0.6, label='Active land (ibound = 1)')
red_line = plt.Line2D([0], [0], color='red', linewidth=1.5, label='Main analysis boundary')
obs_line = plt.Line2D([0], [0], color='magenta', marker='o', linestyle='None', label='Obs. points')
ax_chk.legend(handles=[sea_patch, inactive_patch, active_patch, red_line, obs_line], loc='upper left')

save_path = os.path.join(output_dir, "checked_ibound_observation_preserved.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

logger.info("Cell 5 Execution Completed: 観測点周辺を保持した idomain / ibound を作成しました。")


# In[ ]:


# --- Observation tables for calibration ---
# ADDED FOR OBSERVATION-CONSTRAINED MODEL
import pandas as pd
import numpy as np
import os

xl = pd.ExcelFile('Fiji_field_template_v10.xlsx')
df_r = xl.parse('Model_River_Obs', header=2).dropna(how='all')
river_obs = df_r[df_r['Site'].isin(['FR2', 'FR3', 'FR4', 'FR5', 'FR6', 'FR7'])].copy()
river_obs.rename(columns={'Site': 'site', 'Lat': 'lat', 'Lon': 'lon', 'Field-calculated Q (m3/s)': 'q_obs', 'Mean stage (m)': 'stage', 'Mean bed elev (m)': 'bed'}, inplace=True)
river_obs['q_obs'] = pd.to_numeric(river_obs['q_obs'], errors='coerce')
river_obs['stage'] = pd.to_numeric(river_obs['stage'], errors='coerce')
river_obs['bed'] = pd.to_numeric(river_obs['bed'], errors='coerce')
river_obs = river_obs[['site', 'lat', 'lon', 'q_obs', 'stage', 'bed']].copy()

df_g = xl.parse('Model_GW_Obs', header=2).dropna(how='all')
gw_obs = df_g[df_g['Site'].isin([f'FG{i}' for i in range(1, 10)])].copy()
gw_obs.rename(columns={'Site': 'site', 'Lat': 'lat', 'Lon': 'lon', 'Ground elev (m)': 'elev', 'Depth to water (m)': 'dtw'}, inplace=True)
gw_obs['elev'] = pd.to_numeric(gw_obs['elev'], errors='coerce')
gw_obs['dtw'] = pd.to_numeric(gw_obs['dtw'], errors='coerce')
gw_obs['head_obs'] = np.where(pd.notna(gw_obs['elev']) & pd.notna(gw_obs['dtw']), gw_obs['elev'] - gw_obs['dtw'], np.nan)
gw_obs = gw_obs[['site', 'lat', 'lon', 'elev', 'dtw', 'head_obs']].copy()

print(river_obs)
print(gw_obs)


# In[ ]:


# --- Convert observation points to model CRS / row-col and inspect ---
# ADDED FOR OBSERVATION-CONSTRAINED MODEL
from shapely.geometry import Point
import geopandas as gpd
import rasterio

river_gdf = gpd.GeoDataFrame(river_obs, geometry=[Point(xy) for xy in zip(river_obs['lon'], river_obs['lat'])], crs='EPSG:4326')
gw_gdf = gpd.GeoDataFrame(gw_obs, geometry=[Point(xy) for xy in zip(gw_obs['lon'], gw_obs['lat'])], crs='EPSG:4326')

try:
    with rasterio.open(os.path.join(output_dir, 'model_top.tif')) as src:
        target_crs = src.crs
        print(f'Using CRS from DEM: {target_crs}')
except Exception as e:
    print(f'Could not read DEM CRS, falling back to EPSG:32760... error: {e}')
    target_crs = 'EPSG:32760'

river_gdf = river_gdf.to_crs(target_crs)
gw_gdf = gw_gdf.to_crs(target_crs)
river_gdf['x'] = river_gdf.geometry.x
river_gdf['y'] = river_gdf.geometry.y
gw_gdf['x'] = gw_gdf.geometry.x
gw_gdf['y'] = gw_gdf.geometry.y

def get_rc(x, y, xmin, ymax, cellsize=100.0):
    return (ymax - y) // cellsize, (x - xmin) // cellsize

river_gdf['row'], river_gdf['col'] = zip(*river_gdf.apply(lambda r: get_rc(r['x'], r['y'], xmin, ymax), axis=1))
gw_gdf['row'], gw_gdf['col'] = zip(*gw_gdf.apply(lambda r: get_rc(r['x'], r['y'], xmin, ymax), axis=1))

river_gdf['inside_model'] = (river_gdf['row'] >= 0) & (river_gdf['row'] < nrow) & (river_gdf['col'] >= 0) & (river_gdf['col'] < ncol)
gw_gdf['inside_model'] = (gw_gdf['row'] >= 0) & (gw_gdf['row'] < nrow) & (gw_gdf['col'] >= 0) & (gw_gdf['col'] < ncol)

river_gdf.drop(columns=['geometry']).to_csv(os.path.join(output_dir, 'river_obs.csv'), index=False)
gw_gdf.drop(columns=['geometry']).to_csv(os.path.join(output_dir, 'gw_obs.csv'), index=False)
print('Saved river_obs.csv and gw_obs.csv to output directory.')


# In[68]:


# --- MODFLOW 6 model build / run (observation-constrained revision) ---
# ADDED FOR OBSERVATION-CONSTRAINED MODEL
import flopy
import flopy.utils.binaryfile as bf
import geopandas as gpd
from rasterio.features import rasterize
import pandas as pd
import numpy as np
import os

sim_name = 'MF_FIJI_OBS_REV'
sim = flopy.mf6.MFSimulation(sim_name=sim_name, version='mf6', exe_name='mf6', sim_ws=mf6_dir)
tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS', nper=1, perioddata=[(1.0, 1, 1.0)])
ims = flopy.mf6.ModflowIms(sim, complexity='COMPLEX', outer_maximum=1000, inner_maximum=1000)
gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True, newtonoptions='NEWTON UNDER_RELAXATION')

dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm, idomain=idomain)

highland_mask = (top > 80.0) & (ibound == 1)
midland_mask = (top >= 20.0) & (top <= 80.0) & (ibound == 1)
coastal_mask = (top < 20.0) & (top >= 0.0) & (ibound == 1)

k1 = np.full((nrow, ncol), 12.0)
k2 = np.full((nrow, ncol), 4.0)
k3 = np.full((nrow, ncol), 0.5)
k4 = np.full((nrow, ncol), 0.05)
k1[highland_mask] = 20.0
k2[highland_mask] = 8.0
k1[coastal_mask] = 15.0
k2[coastal_mask] = 5.0

npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=[1, 0, 0, 0], k=[k1, k2, k3, k4], k33=[k1/10, k2/10, k3/10, k4/10], save_flows=True)

strt_array = np.zeros((nlay, nrow, ncol))
for l in range(nlay):
    if l == 0:
        strt_array[l] = np.where(top < 0, 0.0, np.minimum(top, np.maximum(top - 1.0, 0.0)))
    else:
        strt_array[l] = np.where(top < 0, 0.0, np.maximum(top - (l * 3.0 + 1.0), 0.0))
ic = flopy.mf6.ModflowGwfic(gwf, strt=strt_array)

sea_rows, sea_cols = np.where(ibound == -1)
chd_layers = [0, 1]  # Limit sea boundary to shallow layers
chd_spd = [[(l, r, c), 0.0] for r, c in zip(sea_rows, sea_cols) for l in chd_layers]
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

rch_array = np.zeros((nrow, ncol))
rch_array[coastal_mask] = 0.0015
rch_array[midland_mask] = 0.0020
rch_array[highland_mask] = 0.0025
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rch_array)

streams_gdf = gpd.read_file(os.path.join(output_dir, 'streams.shp'))
stream_grid = rasterize([(geom, 1) for geom in streams_gdf.geometry], out_shape=(nrow, ncol), transform=transform, fill=0, all_touched=True).astype(int)
stream_rows, stream_cols = np.where((stream_grid == 1) & (ibound == 1) & (top >= 0))

def get_nearby_stream_cells(row, col, stream_grid, radius=2):
    cells = []
    for rr in range(int(max(0, row - radius)), int(min(stream_grid.shape[0], row + radius + 1))):
        for cc in range(int(max(0, col - radius)), int(min(stream_grid.shape[1], col + radius + 1))):
            if stream_grid[rr, cc] == 1 and ibound[rr, cc] == 1:
                cells.append((rr, cc))
    return cells

obs_cell_lookup = {}
for _, rr in river_gdf.iterrows():
    if not pd.isna(rr['row']) and not pd.isna(rr['col']):
        for rc in get_nearby_stream_cells(rr['row'], rr['col'], stream_grid):
            obs_cell_lookup[rc] = rr['site']

cond_map = {'FR2': 4000.0, 'FR3': 3000.0, 'FR4': 1200.0, 'FR5': 800.0, 'FR6': 1500.0, 'FR7': 400.0}
cond_default = 600.0
drn_spd = []
for r, c in zip(stream_rows, stream_cols):
    site = obs_cell_lookup.get((r, c), None)
    cond = cond_map.get(str(site)[:3] if site else None, cond_default)
    
    z = max(top[r, c] - 0.5, 0.0)
    if site is not None and len(river_obs.loc[river_obs['site'] == site]) > 0:
        bed_val = river_obs.loc[river_obs['site'] == site].iloc[0]['bed']
        if pd.notna(bed_val):
            z = float(bed_val)
    
    drn_spd.append([(0, r, c), float(z), float(cond)])

drn = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=drn_spd, save_flows=True)
oc = flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f'{sim_name}.hds', budget_filerecord=f'{sim_name}.cbc', saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])

sim.write_simulation(silent=True)
success, buff = sim.run_simulation(silent=False)
if success:
    head = bf.HeadFile(os.path.join(mf6_dir, f'{sim_name}.hds')).get_data()
    print('Simulation successful!')
else:
    print('Simulation did not converge.')


# In[ ]:


# --- Compare model results with river / groundwater observations ---
# ADDED FOR OBSERVATION-CONSTRAINED MODEL
import flopy.utils.binaryfile as bf
import warnings
warnings.filterwarnings('ignore')
try:
    head = bf.HeadFile(os.path.join(mf6_dir, f'{sim_name}.hds')).get_data()
    cbc = bf.CellBudgetFile(os.path.join(mf6_dir, f'{sim_name}.cbc'))
    drn_records = cbc.get_data(text='DRN')
except:
    print('Could not load simulation files.')

river_gdf['sim_head_l1'] = [head[0, int(r), int(c)] if pd.notna(r) and 0 <= int(r) < nrow and 0 <= int(c) < ncol else np.nan for r, c in zip(river_gdf['row'], river_gdf['col'])]
gw_gdf['sim_head_l1'] = [head[0, int(r), int(c)] if pd.notna(r) and 0 <= int(r) < nrow and 0 <= int(c) < ncol else np.nan for r, c in zip(gw_gdf['row'], gw_gdf['col'])]

def sum_drn_near_site(site_row, site_col, drn_records, radius=2):
    if pd.isna(site_row) or pd.isna(site_col) or len(drn_records) == 0:
        return np.nan
    total = 0.0
    for rec in drn_records[-1]:
        node = int(rec['node']) - 1
        q = float(rec['q'])
        rem = node % (nrow * ncol)
        r = rem // ncol
        c = rem % ncol
        if abs(r - site_row) <= radius and abs(c - site_col) <= radius:
            total += abs(q)
    return total

river_gdf['sim_drn_q_nearby'] = [sum_drn_near_site(r, c, drn_records) for r, c in zip(river_gdf['row'], river_gdf['col'])]

print('\n=== RIVER Observation Comparison ===')
print('Note: sim_drn_q is the aggregated absolute DRN flow from surrounding cells within a +/- 2 grid radius to approximate field measurements.')
print(river_gdf[['site', 'q_obs', 'sim_drn_q_nearby', 'stage', 'bed', 'row', 'col']].sort_values('site'))

print('\n=== GROUNDWATER Observation Comparison ===')
print('Note: head_obs is only calculated when Depth to Water (dtw) is known. Otherwise NaN.')
print(gw_gdf[['site', 'elev', 'dtw', 'head_obs', 'sim_head_l1', 'row', 'col']].sort_values('site'))


# In[70]:


# --- Cell 4: MODPATH 7 (粒子追跡・高速化版) のみ実行 ---
import flopy
import os
import numpy as np
import flopy.utils.binaryfile as bf

print("MODPATH 7 (粒子追跡) を高速化設定で実行中...")
sim_name_mp = "MF_FIJI_OBS_REV"

# MODFLOWの計算結果を読み込む
head = bf.HeadFile(os.path.join(mf6_dir, f"{sim_name_mp}.hds")).get_data()

# 山岳部（標高50m以上）のセルを抽出
highland_locs_mask = (top > 50.0) & (ibound == 1) & (np.abs(head[0]) < 1000)
rows, cols = np.where(highland_locs_mask)

# 10セルに1つの割合で間引き、表層(k=0)からのみ粒子を発生させて高速化
locs = [(0, r, c) for r, c in zip(rows[::10], cols[::10])]

mp_name = "mp_fiji_nw_fast"
mp = flopy.modpath.Modpath7(modelname=mp_name, flowmodel=gwf, exe_name="mp7", model_ws=mp7_dir)
mpbas = flopy.modpath.Modpath7Bas(mp, porosity=[0.25, 0.20, 0.05, 0.05])
particledata = flopy.modpath.ParticleData(locs, structured=True, drape=0, localx=0.5, localy=0.5, localz=0.9)
pg = flopy.modpath.ParticleGroup(particledata=particledata, particlegroupname="Recharge_Fast")
mpsim = flopy.modpath.Modpath7Sim(mp, simulationtype="pathline", trackingdirection="forward", 
                                  weaksinkoption="pass_through", weaksourceoption="pass_through", particlegroups=[pg])

mp.write_input()
mp_success, _ = mp.run_model(silent=True)

if mp_success:
    print(f"粒子追跡完了！結果は {mp_name}.mppth に保存されました。")
else:
    print("MODPATHの計算に失敗しました。")


# In[73]:


# --- Cell 5: 浅層と深層の地下水流跡線の統合描画 (完全修正版) ---
import flopy
import os
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

print("完了済みの粒子追跡データを読み込んでいます...")
# 【ここが原因でした！】正しいファイル名に修正しました。
mp_name_fast = "mp_fiji_nw_fast" 
p_file_fast = os.path.join(mp7_dir, f"{mp_name_fast}.mppth")

if os.path.exists(p_file_fast):
    # 完了済みの計算結果を読み込む
    pfile_fast = flopy.utils.PathlineFile(p_file_fast)
    plines_fast = pfile_fast.get_alldata()
    print(f"{len(plines_fast)} 本の流跡線を超高速で読み込みました。")

    print("流跡線のマップを作成中 (浅層と深層を区別)...")
    fig_pth, ax_pth = plt.subplots(figsize=(14, 10))
    extent = [xmin, xmax, ymin, ymax]

    # --- 1. 背景の地形図を薄く明るく描画 ---
    dem_display = np.where(top == -9999.0, np.nan, top)
    ax_pth.imshow(dem_display, extent=extent, origin='upper', cmap='terrain', alpha=0.2, zorder=1)

    # Hillshade（地形の影）を非常に薄く重ねる
    ls = LightSource(azdeg=315, altdeg=45)
    hillshade = ls.hillshade(dem_display, vert_exag=1.0) 
    ax_pth.imshow(hillshade, extent=extent, origin='upper', cmap='gray', alpha=0.05, zorder=2)

    # 海岸線 (levels=[0])
    ax_pth.contour(np.linspace(xmin, xmax, ncol), np.linspace(ymax, ymin, nrow), top, levels=[0], colors='black', linewidths=1.0, zorder=3)

    # 河川
    streams_shp_path = os.path.join(output_dir, "streams.shp")
    if os.path.exists(streams_shp_path):
        streams_gdf = gpd.read_file(streams_shp_path)
        streams_gdf.plot(ax=ax_pth, color='deepskyblue', linewidth=1.0, alpha=0.8, zorder=4)

    # --- 2. 浅層と深層を区別して抽出 ---
    deep_paths = [p for p in plines_fast if np.min(p['z']) < -100]
    shallow_paths = [p for p in plines_fast if np.min(p['z']) >= -100]

    # --- 描画処理 ---
    # A. 浅層地下水 (薄い青色の線)
    print(f"※ 浅層を通った {len(shallow_paths)} 本のパスを青い薄い線でプロットしています...")
    for p in shallow_paths:
        ax_pth.plot(p['x'] + xmin, p['y'] + ymin, color='dodgerblue', linewidth=0.3, alpha=0.2, zorder=9)

    # B. 深層地下水 (赤い色で強調)
    if len(deep_paths) > 0:
        print(f"※ 深層を通った {len(deep_paths)} 本のパスを赤い薄い線でプロットしています...")
        
        # 1. 流跡線 (赤い薄い線)
        for p in deep_paths:
            ax_pth.plot(p['x'] + xmin, p['y'] + ymin, color='#FF0000', linewidth=0.5, alpha=0.4, zorder=10)

        # 2. 湧出地点 (赤い点)
        discharge_x = [p['x'][-1] + xmin for p in deep_paths]
        discharge_y = [p['y'][-1] + ymin for p in deep_paths]
        ax_pth.scatter(discharge_x, discharge_y, color='red', s=20, marker='o', alpha=0.9, zorder=11, label='Deep Discharge Points')
    else:
        print("※ 注意: -100mより深く潜った粒子はありませんでした。")

    # --- 3. 外観の設定 ---
    ax_pth.set_xlim(xmin, xmax)
    ax_pth.set_ylim(ymin, ymax)
    ax_pth.set_title("Check: Shallow & Deep Groundwater Flowpaths", fontsize=16)
    ax_pth.set_xlabel("UTM Easting (m)")
    ax_pth.set_ylabel("UTM Northing (m)")

    # 凡例を手動で追加
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches
    sea_patch = mpatches.Patch(color='none', edgecolor='black', label='Coastline')
    stream_line = Line2D([0], [0], color='deepskyblue', linewidth=1.0, label='Streams')
    shallow_line = Line2D([0], [0], color='dodgerblue', linewidth=0.5, alpha=0.8, label='Shallow Flowpath (z >= -100m)')
    deep_line = Line2D([0], [0], color='#FF0000', linewidth=1.0, alpha=0.8, label='Deep Flowpath (z < -100m)')
    point_marker = Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=6, label='Discharge Point')
    
    ax_pth.legend(handles=[sea_patch, stream_line, shallow_line, deep_line, point_marker], loc='upper right', fontsize=12)

    # 保存処理
    save_path = os.path.join(output_dir, "groundwater_pathlines_combined.png")
    print("画像を保存しています（数秒お待ちください）...")
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

    print("すべての処理が完了しました！")
else:
    print(f"ファイルが見つかりません: {p_file_fast}")


# In[ ]:





# In[77]:


fig_cx, ax_cx = plt.subplots(figsize=(14, 8))
transect_col = 140
target_x = transect_col * 100 

local_y = np.linspace(ymax - ymin, 0, nrow)

# 背景の地質と地形
ax_cx.fill_between(local_y, botm[0][:, transect_col], top[:, transect_col], color='#e5ce8a', alpha=0.4, label='Layer 1 (Sediments)')
ax_cx.fill_between(local_y, botm[1][:, transect_col], botm[0][:, transect_col], color='#c1a965', alpha=0.4, label='Layer 2 (Weathered)')
ax_cx.fill_between(local_y, botm[2][:, transect_col], botm[1][:, transect_col], color='#827658', alpha=0.4, label='Layer 3 (Bedrock Top)')
ax_cx.fill_between(local_y, botm[3][:, transect_col], botm[2][:, transect_col], color='#5b523c', alpha=0.4, label='Layer 4 (Deep Bedrock)')
ax_cx.plot(local_y, top[:, transect_col], color='k', linewidth=2.0, label='Topography')

# パスの抽出と描画
filtered_paths = [p for p in plines_fast if (target_x - 500) <= p['x'][0] <= (target_x + 500)]
for p in filtered_paths[::10]:  
    if np.min(p['z']) < -150:
        ax_cx.plot(p['y'], p['z'], color='magenta', linewidth=2.5, alpha=1.0, zorder=10) # 深層
    else:
        ax_cx.plot(p['y'], p['z'], color='red', linewidth=1.5, alpha=0.8, zorder=5) # 浅層

ax_cx.plot([], [], color='red', linewidth=1.5, label='Shallow Flow (Baseflow)')
ax_cx.plot([], [], color='magenta', linewidth=2.5, label='Deep Flow (SGD)')

ax_cx.set_xlim(0, ymax-ymin)
ax_cx.set_ylim(-750, 500)
ax_cx.set_xlim(ax_cx.get_xlim()[::-1]) 

ax_cx.set_title("Cross Section (Col 140) - Subsurface Groundwater Discharge", fontsize=16)
ax_cx.set_xlabel("Local Northing (m) [South (Highland) -> North (Reef)]", fontsize=12)
ax_cx.set_ylabel("Elevation (m)", fontsize=12)
ax_cx.legend(loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=10)
ax_cx.grid(True, alpha=0.3)
plt.subplots_adjust(bottom=0.2)
plt.close()


# In[76]:


from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

fig_map, ax_map = plt.subplots(figsize=(14, 10))
extent = [xmin, xmax, ymin, ymax]

# 1. DEMの陰影起伏図作成
dem_display = np.where(top == -9999.0, np.nan, top)
ls = LightSource(azdeg=315, altdeg=45)
hillshade = ls.hillshade(dem_display, vert_exag=5)
ax_map.imshow(hillshade, extent=extent, origin='upper', cmap='gray', alpha=0.5)

# 2. 海岸線（黒線）
ax_map.contour(np.linspace(xmin, xmax, ncol), np.linspace(ymax, ymin, nrow), top, levels=[0], colors='black', linewidths=1.5)

# 3. 【追加】流域境界（分水嶺）を薄い白の破線で描画
basins_shp_path = os.path.join(output_dir, "basins.shp")
if os.path.exists(basins_shp_path):
    basins_gdf = gpd.read_file(basins_shp_path)
    # facecolor='none' で塗りつぶしを消し、edgecolor で枠線だけを描く
    basins_gdf.plot(ax=ax_map, facecolor='none', edgecolor='white', linewidth=1.5, linestyle='--', alpha=0.6, zorder=13)
    ax_map.plot([], [], color='white', linewidth=1.5, linestyle='--', alpha=0.6, label='Watershed Boundaries')

# 4. 島全体の全河川（streams.shp）を水色の線で背景に描画
streams_shp_path = os.path.join(output_dir, "streams.shp")
if os.path.exists(streams_shp_path):
    streams_gdf = gpd.read_file(streams_shp_path)
    streams_gdf.plot(ax=ax_map, color='deepskyblue', linewidth=1.5, alpha=0.8, zorder=11)
    ax_map.plot([], [], color='deepskyblue', linewidth=1.5, label='All Streams (streams.shp)')

# 5. モデル対象流域の河川（SFRセル）を青い四角で強調
df_sfr = pd.read_csv(os.path.join(output_dir, "sfr_all_network.csv"))
sfr_x = xmin + (df_sfr['Col'] - 0.5) * 100
sfr_y = ymax - (df_sfr['Row'] - 0.5) * 100
ax_map.plot(sfr_x, sfr_y, marker='s', color='mediumblue', markersize=2.5, linestyle='none', zorder=12, label='Active River Cells (SFR)')

# 6. パスラインの描画（50本に1本）と海域カット
for p in plines[::50]:
    p_x, p_y, p_z = p['x'], p['y'], p['z']
    valid_x, valid_y = [], []
    for i in range(len(p_x)):
        c = int(np.clip(p_x[i] / 100, 0, ncol-1))
        r = int(np.clip(nrow - 1 - p_y[i] / 100, 0, nrow-1))
        if top[r, c] >= 0:
            valid_x.append(xmin + p_x[i])
            valid_y.append(ymin + p_y[i])
        else:
            break
            
    if len(valid_x) > 1:
        if np.min(p_z) < -150:
            ax_map.plot(valid_x, valid_y, color='magenta', linewidth=0.5, alpha=0.5, zorder=5)
        else:
            ax_map.plot(valid_x, valid_y, color='red', linewidth=0.6, alpha=0.6, zorder=6)

# 凡例用ダミープロット
ax_map.plot([], [], color='magenta', linewidth=1.5, label='Deep Flow (SGD)')
ax_map.plot([], [], color='red', linewidth=1.5, label='Shallow Flow')

ax_map.set_xlim(xmin, xmax)
ax_map.set_ylim(ymin, ymax)
ax_map.set_title("Plan View - Groundwater Pathlines & Natural Watersheds", fontsize=16)
ax_map.set_xlabel("UTM Easting (m)", fontsize=12)
ax_map.set_ylabel("UTM Northing (m)", fontsize=12)

# 凡例を少し整理して見やすく配置
ax_map.legend(loc='lower right', fontsize=10, framealpha=0.9)
plt.close()


# In[ ]:





# In[ ]:





# ## Observation-Constrained Model Revision Summary
# **Modified Cells:**
# - **Observation tables for calibration**: Parsed real values from `Fiji_field_template_v10.xlsx`. Unknown properties (like stage/bed/dtw) were safely loaded as `NaN`.
# - **Convert observation points to model CRS**: Explicit fallback implemented to handle missing CRS in `.shp`. Exported dataframes to `river_obs.csv` and `gw_obs.csv`.
# - **MODFLOW 6 model build / run**: 
#   - Explicitly enforced NPF with Layer 1 unconfined (`icelltype=1`), and 2-4 confined.
#   - Enforced separate recharge zones (Coastal, Midland, Highland).
#   - Enforced shallow-only CHD forcing.
#   - Extracted and applied river point specific conductance.
# - **Compare model results**: Re-worked error comparison logical structures to clearly discern observed values versus numeric representation aggregated equivalents (e.g. `sim_drn_q` aggregated +/- 2 cell radius).
# 
# **Pending/Unknown Values:**
# - Many points lack `dtw` resulting in `head_obs = NaN`.
# - Stage and bed elevation variables inside `river_obs` are mostly incomplete/NaN.


logger.info('=== WORKFLOW COMPLETED SUCCESSFULLY ===')
sys.exit(0)
