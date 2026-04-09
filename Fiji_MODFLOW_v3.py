#!/usr/bin/env python
# coding: utf-8
"""
MODFLOW 6 Groundwater Model for Rakiraki Region, Fiji
Version 3 - Improved solver and observation data handling

Author: Hydrology Workflow
Date: 2026-04-06

Key improvements over v2:
1. Improved IMS solver settings for non-convergence issues
2. Corrected observation data reading from proper Excel sheets
3. Enhanced deep groundwater flow through improved K values and convertible cells
4. Relative path handling for GitHub compatibility
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, ListedColormap
import matplotlib.patches as mpatches
from rasterio.transform import from_origin
from rasterio.features import rasterize
import flopy
import flopy.utils.binaryfile as bf
import argparse
import sys
import logging
from shapely.geometry import box, Point
from shapely.ops import unary_union

# ============================================================================
# Section 1: Logger and Argument Parser Setup
# ============================================================================

logger = logging.getLogger("FijiMODFLOW_v3")
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler('fiji_modflow_v3.log', encoding='utf-8')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Run MODFLOW 6 for Rakiraki region, Fiji.')
parser.add_argument('--dry-run', action='store_true', help='Dry-run mode (skip execution).')
parser.add_argument('--skip-modpath', action='store_true', help='Skip MODPATH 7 particle tracking.')
args = parser.parse_args()

if args.dry_run:
    logger.info("=== DRY-RUN MODE ===")

logger.info("=== STARTING MODFLOW 6 v3 WORKFLOW ===")
warnings.filterwarnings('ignore')

# ============================================================================
# Section 2: Path Setup with Relative Paths
# ============================================================================

# Use relative path based on script location for GitHub compatibility
base_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(base_dir, "output")
mf6_dir = os.path.join(base_dir, "mf6_workspace")
mp7_dir = os.path.join(base_dir, "mp7_workspace")

logger.info(f"Base directory: {base_dir}")
logger.info(f"Output directory: {output_dir}")
logger.info(f"MODFLOW 6 workspace: {mf6_dir}")
logger.info(f"MODPATH 7 workspace: {mp7_dir}")

# Create workspaces if they don't exist
for d in [mf6_dir, mp7_dir]:
    if not os.path.exists(d):
        os.makedirs(d)
        logger.info(f"Created directory: {d}")

# ============================================================================
# Section 3: Data Loading (DEM, ibound, domain, basins)
# ============================================================================

logger.info("Loading DEM and model configuration...")

# Load model domain parameters
with open(os.path.join(output_dir, "modflow_domain.json"), 'r') as f:
    domain = json.load(f)

nrow = domain['nrow']  # 200
ncol = domain['ncol']  # 240
xmin = domain['xmin']  # 611000
xmax = domain['xmax']  # 635000
ymin = domain['ymin']  # 8068000
ymax = domain['ymax']  # 8088000

logger.info(f"Model domain: {nrow} rows x {ncol} cols, 100m cells")
logger.info(f"Extent: xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax}")

# Load DEM (model top elevation) and ibound
top = np.loadtxt(os.path.join(output_dir, "model_top.csv"), delimiter=',')
ibound = np.loadtxt(os.path.join(output_dir, "model_ibound.csv"), delimiter=',')

logger.info(f"Loaded DEM: min={top.min():.2f}m, max={top.max():.2f}m")

# Load basin data
basins = gpd.read_file(os.path.join(output_dir, 'basins_clean.shp'))
basins['VALUE'] = basins['VALUE'].astype(int).astype(str)

logger.info(f"Loaded {len(basins)} basins")

# ============================================================================
# Section 4: Basin Masking and Domain Setup
# ============================================================================

logger.info("Setting up domain and basin masking...")

transform = from_origin(xmin, ymax, 100.0, 100.0)

# Select focal basins (Rakiraki region: basins 2422 and 2344)
focal_basins = basins[basins['VALUE'].isin(['2422', '2344'])].copy()
focal_basins['geometry'] = focal_basins.geometry.buffer(150)  # 150m buffer

logger.info(f"Focal basins: {list(focal_basins['VALUE'])}")

# Create basin mask
basin_mask = rasterize(
    [(geom, 1) for geom in focal_basins.geometry],
    out_shape=(nrow, ncol),
    transform=transform,
    fill=0,
    all_touched=True
)

# Mark highland cells (elevation > 50m) as reference for K variation
highland_mask = (top > 50.0)

# Deactivate cells outside focal basin
ibound[(ibound == 1) & (basin_mask == 0)] = 0

# Open coastal low-elevation areas (<5m) to allow coastal flow
ibound[(ibound == 0) & (top < 5.0) & (top >= 0)] = 1

logger.info("Basin masking completed")

# ============================================================================
# Section 5: Observation Data Loading (CORRECTED)
# ============================================================================

logger.info("Loading observation data from Excel file...")

excel_file = os.path.join(base_dir, 'Fiji_field_template_v10.xlsx')

if not os.path.exists(excel_file):
    logger.warning(f"Excel file not found: {excel_file}")
    excel_file = 'Fiji_field_template_v10.xlsx'

try:
    # Load metadata from Controls sheet
    xl = pd.ExcelFile(excel_file)
    controls = xl.parse('Controls', header=2)
    antenna_height = float(controls[controls['Item'] == 'Antenna height']['Value'].values[0])
    geoid_height = float(controls[controls['Item'] == 'Default provisional geoid height']['Value'].values[0])
    logger.info(f"Controls: antenna_height={antenna_height}m, geoid_height={geoid_height}m")

    # --- Load River observation points from River_Points sheet ---
    # Column headers: Site (A), Lat (B), Lon (C), GNSS APC h (D), Q values (H)
    river_data = xl.parse('River_Points', header=2).dropna(how='all')

    # Convert formula references to actual values (antenna and geoid)
    river_data['Antenna h (m)'] = antenna_height
    river_data['Geoid N (m)'] = geoid_height
    river_data['Ground elev (m)'] = (
        river_data['GNSS APC h (m)'] -
        river_data['Antenna h (m)'] -
        river_data['Geoid N (m)']
    )

    # Build river observation dataframe
    river_obs = pd.DataFrame({
        'site': river_data['Site'],
        'lat': river_data['Lat'],
        'lon': river_data['Lon'],
        'elev': river_data['Ground elev (m)'],
        'q_obs': pd.to_numeric(river_data['Approx observed Q (m3/s)'], errors='coerce')
    })

    # Filter valid records
    river_obs = river_obs[river_obs['site'].isin(['FR2', 'FR3', 'FR4', 'FR5', 'FR6', 'FR7'])].copy()
    river_obs = river_obs.dropna(subset=['lat', 'lon', 'elev'])

    logger.info(f"Loaded {len(river_obs)} river observation points")
    logger.info("\nRiver observations:")
    for _, row in river_obs.iterrows():
        q_str = f"{row['q_obs']:.4f}" if pd.notna(row['q_obs']) else "NaN"
        logger.info(f"  {row['site']}: elev={row['elev']:.2f}m, Q={q_str} m³/s")

    # --- Load Groundwater observation points from Groundwater sheet ---
    # Column headers: Site (A), Lat (B), Lon (C), GNSS APC h (D), DTW (H)
    gw_data = xl.parse('Groundwater', header=2).dropna(how='all')

    # Convert formula references
    gw_data['Antenna h (m)'] = antenna_height
    gw_data['Geoid N (m)'] = geoid_height
    gw_data['Ground elev (m)'] = (
        gw_data['GNSS APC h (m)'] -
        gw_data['Antenna h (m)'] -
        gw_data['Geoid N (m)']
    )

    # Compute groundwater head from depth-to-water
    gw_data['GW head (m)'] = np.where(
        pd.notna(gw_data['Depth to water (m)']) & (gw_data['Depth to water (m)'] != ''),
        gw_data['Ground elev (m)'] - gw_data['Depth to water (m)'],
        np.nan
    )

    # Build groundwater observation dataframe
    gw_obs = pd.DataFrame({
        'site': gw_data['Site'],
        'lat': gw_data['Lat'],
        'lon': gw_data['Lon'],
        'elev': gw_data['Ground elev (m)'],
        'dtw': pd.to_numeric(gw_data['Depth to water (m)'], errors='coerce'),
        'head_obs': gw_data['GW head (m)']
    })

    # Filter valid records (Exclude FG9 which is on another island)
    gw_obs = gw_obs[gw_obs['site'].isin([f'FG{i}' for i in range(1, 9)])].copy()
    gw_obs = gw_obs.dropna(subset=['lat', 'lon', 'elev'])

    logger.info(f"Loaded {len(gw_obs)} groundwater observation points")
    logger.info("\nGroundwater observations:")
    for _, row in gw_obs.iterrows():
        head_str = f"{row['head_obs']:.2f}m" if pd.notna(row['head_obs']) else "NaN"
        dtw_str = f"{row['dtw']:.2f}m" if pd.notna(row['dtw']) else "NaN"
        logger.info(f"  {row['site']}: elev={row['elev']:.2f}m, DTW={dtw_str}, head={head_str}")

    # Convert to UTM coordinates (EPSG:32760)
    river_gdf = gpd.GeoDataFrame(
        river_obs,
        geometry=[Point(xy) for xy in zip(river_obs['lon'], river_obs['lat'])],
        crs='EPSG:4326'
    ).to_crs('EPSG:32760')

    gw_gdf = gpd.GeoDataFrame(
        gw_obs,
        geometry=[Point(xy) for xy in zip(gw_obs['lon'], gw_obs['lat'])],
        crs='EPSG:4326'
    ).to_crs('EPSG:32760')

    # Convert to row/col indices
    def get_rc(x, y, xmin, ymax, cellsize=100.0):
        col = int((x - xmin) / cellsize)
        row = int((ymax - y) / cellsize)
        return row, col

    river_gdf['row'] = river_gdf.geometry.apply(lambda g: get_rc(g.x, g.y, xmin, ymax)[0])
    river_gdf['col'] = river_gdf.geometry.apply(lambda g: get_rc(g.x, g.y, xmin, ymax)[1])
    gw_gdf['row'] = gw_gdf.geometry.apply(lambda g: get_rc(g.x, g.y, xmin, ymax)[0])
    gw_gdf['col'] = gw_gdf.geometry.apply(lambda g: get_rc(g.x, g.y, xmin, ymax)[1])

    # Mark points inside model domain
    river_gdf['inside_model'] = (
        (river_gdf['row'] >= 0) & (river_gdf['row'] < nrow) &
        (river_gdf['col'] >= 0) & (river_gdf['col'] < ncol)
    )
    gw_gdf['inside_model'] = (
        (gw_gdf['row'] >= 0) & (gw_gdf['row'] < nrow) &
        (gw_gdf['col'] >= 0) & (gw_gdf['col'] < ncol)
    )

    # Save observation locations for reference
    river_gdf.drop(columns=['geometry']).to_csv(
        os.path.join(output_dir, 'river_obs_v3.csv'), index=False
    )
    gw_gdf.drop(columns=['geometry']).to_csv(
        os.path.join(output_dir, 'gw_obs_v3.csv'), index=False
    )
    logger.info("Saved observation locations to river_obs_v3.csv and gw_obs_v3.csv")

    # Ensure observation points are marked as active in model
    for _, row in river_gdf.iterrows():
        if row['inside_model']:
            r, c = int(row['row']), int(row['col'])
            for rr in range(max(0, r-1), min(nrow, r+2)):
                for cc in range(max(0, c-1), min(ncol, c+2)):
                    if top[rr, cc] >= 0:
                        ibound[rr, cc] = 1

    for _, row in gw_gdf.iterrows():
        if row['inside_model']:
            r, c = int(row['row']), int(row['col'])
            for rr in range(max(0, r-1), min(nrow, r+2)):
                for cc in range(max(0, c-1), min(ncol, c+2)):
                    if top[rr, cc] >= 0:
                        ibound[rr, cc] = 1

except Exception as e:
    logger.error(f"Error loading observation data: {e}")
    logger.info("Continuing without observation data...")
    river_gdf = None
    gw_gdf = None

logger.info("Observation data loading completed")

# ============================================================================
# Section 6: Layer Structure Setup
# ============================================================================

logger.info("Setting up layer structure...")

# 層の設定
# Layer 1: 0 ~ -50m (50m thick) - phreatic layer
# Layer 2: -50 ~ -100m (50m thick) - convertible to confined
# Layer 3: -100 ~ -250m (150m thick) - confined
# Layer 4: -250 ~ -700m (450m thick) - deep confined
nlay = 4
botm = np.zeros((nlay, nrow, ncol))
botm[0] = top - 50.0       # L1 bottom
botm[1] = botm[0] - 50.0   # L2 bottom
botm[2] = botm[1] - 150.0  # L3 bottom
botm[3] = botm[2] - 450.0  # L4 bottom (deepest)

logger.info(f"Layer configuration:")
logger.info(f"  L1: 50m (phreatic)")
logger.info(f"  L2: 50m (convertible)")
logger.info(f"  L3: 150m (confined)")
logger.info(f"  L4: 450m (deep confined)")

# Set idomain for all non-inactive cells
idomain = np.array([np.where(ibound != 0, 1, 0)] * nlay)

logger.info("Layer structure completed")

# ============================================================================
# Section 7: MODFLOW 6 Model Build with Improved Parameters
# ============================================================================

logger.info("Building MODFLOW 6 model...")

sim_name = "MF_FIJI_v3"
sim = flopy.mf6.MFSimulation(sim_name=sim_name, version="mf6", exe_name="mf6", sim_ws=mf6_dir)

# Time discretization: 1 stress period
tdis = flopy.mf6.ModflowTdis(
    sim,
    time_units="DAYS",
    nper=1,
    perioddata=[(1.0, 1, 1.0)]
)

logger.info("TDIS package created")

# ============================================================================
# Section 8: Improved IMS Solver Settings
# ============================================================================
# 【改善点】ソルバー非収束対策
# - outer_dvclose, inner_dvclose で水頭変化の収束判定
# - rcloserecord で残差の収束判定（strict mode）
# - DBD under-relaxation で振動を抑制
# - BICGSTAB加速器で収束速度向上
# - outer_maximum, inner_maximum を増加

logger.info("Configuring improved IMS solver settings...")

ims = flopy.mf6.ModflowIms(
    sim,
    complexity="COMPLEX",           # ILUT前処理（ill-conditioned行列に必要）
    outer_maximum=800,
    inner_maximum=300,
    outer_dvclose=0.01,
    inner_dvclose=0.001,
    rcloserecord=[0.1, "strict"],
    linear_acceleration="BICGSTAB",
    backtracking_number=15,
    backtracking_tolerance=1.05,
    backtracking_reduction_factor=0.2,
    backtracking_residual_limit=50.0
)

logger.info("IMS solver configured with improved convergence settings")

# Create groundwater flow model
gwf = flopy.mf6.ModflowGwf(
    sim,
    modelname=sim_name,
    save_flows=True,
    newtonoptions="NEWTON UNDER_RELAXATION"
)

# Discretization
dis = flopy.mf6.ModflowGwfdis(
    gwf,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=100.0,
    delc=100.0,
    top=top,
    botm=botm,
    idomain=idomain
)

logger.info("DIS package created")

# ============================================================================
# Section 9: Node Property Flow (NPF) with Improved K Values
# ============================================================================
# 【改善点】深い地下水流を促進するため、深層のK値を増加
# - L1: k=8-15 (from 12-20) - shallow流を少し抑制
# - L2: k=3-6 (from 4-8) - shallow層間流を抑制
# - L3: k=1.0-2.0 (from 0.5) - 深層流を促進（3倍）
# - L4: k=0.2-0.5 (from 0.05) - 最深層を促進（5-10倍）
# - icelltype=[1,1,0,0]: L1,L2可変飽和、L3,L4被圧
#   （L2も可変飽和に変更して、L2-L3間での垂直流を柔軟に）

logger.info("Setting up improved hydraulic conductivity (K) distribution...")

# K値の設定（第4弾調整）
# 【調整】FG1がまだ6.35m低い → K1をさらに下げて水頭を上げる
# 【調整】深層パス7.6% → 目標10-30%に向けK3/K4をさらに大幅増
k1_default = 2.0   # L1: 3→2 m/d（浅層排水をさらに絞る）
k2_default = 1.5   # L2: 2→1.5 m/d（L2も少し絞る）
k3_default = 8.0   # L3: 5→8 m/d（深層流を大幅促進）
k4_default = 4.0   # L4: 2→4 m/d（最深層も大幅に流れやすく）

k1 = np.full((nrow, ncol), k1_default)
k2 = np.full((nrow, ncol), k2_default)
k3 = np.full((nrow, ncol), k3_default)
k4 = np.full((nrow, ncol), k4_default)

# 山地ではやや高い透水係数（風化岩・亀裂多い）
k1[highland_mask] = 4.0    # 6→4（山地浅層もさらに排水抑制）
k2[highland_mask] = 2.5    # 3→2.5
k3[highland_mask] = 12.0   # 8→12（山地深層流を大幅促進）
k4[highland_mask] = 6.0    # 4→6

logger.info(f"K values (m/day):")
logger.info(f"  L1: {k1.min():.2f} ~ {k1.max():.2f}")
logger.info(f"  L2: {k2.min():.2f} ~ {k2.max():.2f}")
logger.info(f"  L3: {k3.min():.2f} ~ {k3.max():.2f}")
logger.info(f"  L4: {k4.min():.2f} ~ {k4.max():.2f}")

# Vertical anisotropy (K_horizontal / K_vertical = 10)
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    icelltype=[1, 0, 0, 0],  # L1のみ不圧（物理的に妥当）、L2-L4は被圧
    k=[k1, k2, k3, k4],
    k33=[k1/10, k2/10, k3/10, k4/10],
    save_flows=True
)

logger.info("NPF package created with improved K distribution")

# ============================================================================
# Section 10: Initial Conditions and Boundary Conditions
# ============================================================================

logger.info("Setting up initial and boundary conditions...")

# Initial conditions
# 初期水頭をDEMより少し低く設定（water table ≈ 地表面の80%程度）
# これによりNEWTONソルバーの初期推定が改善され、収束が速くなる
strt_array = np.zeros((nlay, nrow, ncol))
for l in range(nlay):
    strt_array[l] = np.where(top < 0, 0.0, np.maximum(top * 0.8, 0.0))

ic = flopy.mf6.ModflowGwfic(gwf, strt=strt_array)

# Constant head (海岸: 0m海抜)
# 浅い2層のみにCHDを適用し、深層の地下水流動を阻害しない
chd_layers = [0, 1]  # L1, L2のみ
chd_spd = [
    [(l, r, c), 0.0] for r, c in zip(*np.where(ibound == -1)) for l in chd_layers
]
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

logger.info(f"CHD boundary: {len(chd_spd)} cells at sea level")

# Recharge (涵養)
# 【調整】モデル水頭が全体的に低い → 涵養量を増加
# フィジーの年間降水量は2000-3000mm、涵養率30-40%として 600-1200mm/年
base_rch = 0.005   # m/day ≈ 1.83 m/year（0.004→0.005に増加）
rch_array = np.zeros((nrow, ncol))
rch_array[ibound != 0] = base_rch
rch_array[highland_mask] *= 0.8  # 山地：地表流出がやや多い（0.7→0.8に緩和）

rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rch_array)

logger.info(f"Recharge: {base_rch:.4f} m/day (base), {(base_rch*0.8):.4f} m/day (highland)")

logger.info("Boundary conditions setup completed")

# ============================================================================
# Section 11: DRN Package（河川排水）— SFRより安定で収束が速い
# ============================================================================
# 【変更理由】SFRは河川ルーティングの非線形性が高く、NEWTONソルバーとの
# 組み合わせで収束が極めて遅くなる。v2のOBS_REVモデルでDRNが使われていた
# のと同様、DRNに戻して安定性を確保する。
# DRNは「水頭 > 排水標高」のときに地下水を排出するシンプルな境界条件。

logger.info("Setting up DRN (drain) package for rivers...")

# Load stream network
streams_gdf = gpd.read_file(os.path.join(output_dir, "streams.shp"))
stream_grid = rasterize(
    [(geom, 1) for geom in streams_gdf.geometry],
    out_shape=(nrow, ncol),
    transform=transform,
    fill=0,
    all_touched=True
)

# Identify stream cells
stream_rows, stream_cols = np.where((stream_grid == 1) & (ibound == 1) & (top >= 0))

logger.info(f"Identified {len(stream_rows)} stream cells for DRN")

# 観測地点近傍のセルには観測値に基づくコンダクタンスを設定
obs_cell_lookup = {}
if river_gdf is not None:
    for _, rr in river_gdf.iterrows():
        if pd.notna(rr['row']) and pd.notna(rr['col']):
            r_obs, c_obs = int(rr['row']), int(rr['col'])
            # 観測点から半径2セル以内の河川セルにラベル付け
            for dr in range(-2, 3):
                for dc in range(-2, 3):
                    nr, nc = r_obs + dr, c_obs + dc
                    if 0 <= nr < nrow and 0 <= nc < ncol:
                        if stream_grid[nr, nc] == 1 and ibound[nr, nc] == 1:
                            obs_cell_lookup[(nr, nc)] = rr['site']

# サイト別コンダクタンス（v2から引き継ぎ）
cond_map = {'FR2': 4000.0, 'FR3': 3000.0, 'FR4': 1200.0,
            'FR5': 800.0, 'FR6': 1500.0, 'FR7': 400.0}
cond_default = 600.0

# DRN stress period data を構築
drn_spd = []
for r, c in zip(stream_rows, stream_cols):
    site = obs_cell_lookup.get((r, c), None)
    cond = cond_map.get(site, cond_default)

    # 排水標高 = 地表面 - 0.5m（河床は地表面のやや下）
    z_drn = max(float(top[r, c]) - 0.5, 0.0)

    # 観測地点では実測河床標高があれば使う
    if site is not None and river_gdf is not None:
        match = river_obs[river_obs['site'] == site]
        if len(match) > 0:
            obs_elev = match.iloc[0]['elev']
            if pd.notna(obs_elev):
                z_drn = float(obs_elev)

    drn_spd.append([(0, r, c), z_drn, float(cond)])

drn = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=drn_spd, save_flows=True)

# Save stream data for visualization
pd.DataFrame([(r, c) for r, c in zip(stream_rows, stream_cols)],
             columns=['Row', 'Col']).to_csv(
    os.path.join(output_dir, "drn_cells_v3.csv"), index=False
)

logger.info(f"DRN package created with {len(drn_spd)} drain cells")

# ============================================================================
# Section 12: Output Control
# ============================================================================

oc = flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord=f"{sim_name}.hds",
    budget_filerecord=f"{sim_name}.cbc",
    saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')]
)

logger.info("Output control configured")

# ============================================================================
# Section 13: Write and Run Simulation
# ============================================================================

logger.info("Writing MODFLOW 6 input files...")

success = False
if not args.dry_run:
    try:
        sim.write_simulation()
        logger.info("Simulation files written successfully")

        logger.info("Running MODFLOW 6 simulation...")
        success, buff = sim.run_simulation()

        if success:
            logger.info("MODFLOW 6 simulation completed successfully!")
        else:
            logger.warning("MODFLOW 6 simulation did not fully converge.")
            logger.warning("Proceeding with post-processing (results may still be usable).")
    except Exception as e:
        logger.error(f"Error during simulation: {e}")
else:
    logger.info("[DRY-RUN] Skipping simulation execution")

# ============================================================================
# Section 14: Post-Processing and Observation Comparison
# ============================================================================
# 【重要】未収束でも .hds が出力されていれば後処理を実行する。
# v2では未収束ながら水収支0.02%で有用な結果が得られていた。

logger.info("Post-processing simulation results...")

# .hds ファイルが存在すれば後処理を行う（収束の有無にかかわらず）
hds_exists = os.path.exists(os.path.join(mf6_dir, f"{sim_name}.hds"))
if not args.dry_run and hds_exists:
    try:
        # Read heads from binary file
        head_file = os.path.join(mf6_dir, f"{sim_name}.hds")
        if os.path.exists(head_file):
            hds = bf.HeadFile(head_file)
            heads = hds.get_data(totim=hds.get_times()[-1])  # Last time step
            logger.info(f"Read heads: min={heads.min():.2f}m, max={heads.max():.2f}m")

            # Save head array for layer 1 (phreatic)
            np.savetxt(os.path.join(output_dir, 'heads_l1_v3.csv'), heads[0], delimiter=',')
            logger.info("Saved L1 heads to heads_l1_v3.csv")

            # Compare with groundwater observations
            if gw_gdf is not None and len(gw_gdf) > 0:
                logger.info("\nComparison with GW observations:")
                logger.info("Site | Model head (m) | Obs head (m) | Diff (m)")
                logger.info("-" * 50)

                for _, row in gw_gdf.iterrows():
                    if row['inside_model'] and pd.notna(row['head_obs']):
                        r, c = int(row['row']), int(row['col'])
                        if 0 <= r < nrow and 0 <= c < ncol:
                            model_head = heads[0, r, c]
                            obs_head = row['head_obs']
                            diff = model_head - obs_head
                            logger.info(f"{row['site']:5} | {model_head:14.2f} | {obs_head:11.2f} | {diff:8.2f}")

            # Compare with river observations (using layer 1 heads)
            if river_gdf is not None and len(river_gdf) > 0:
                logger.info("\nComparison with river observations:")
                logger.info("Site | Model head (m) | Obs elev (m) | Diff (m)")
                logger.info("-" * 50)

                for _, row in river_gdf.iterrows():
                    if row['inside_model']:
                        r, c = int(row['row']), int(row['col'])
                        if 0 <= r < nrow and 0 <= c < ncol:
                            model_head = heads[0, r, c]
                            obs_elev = row['elev']
                            diff = model_head - obs_elev
                            logger.info(f"{row['site']:5} | {model_head:14.2f} | {obs_elev:11.2f} | {diff:8.2f}")
        else:
            logger.warning(f"Head file not found: {head_file}")
    except Exception as e:
        logger.error(f"Error in post-processing: {e}")

# ============================================================================
# Section 15: MODPATH 7 Particle Tracking
# ============================================================================

if not args.dry_run and hds_exists and not args.skip_modpath:
    logger.info("Setting up MODPATH 7 particle tracking...")

    try:
        # 山岳部（標高50m以上）のセルから粒子を発生
        head_for_mp = bf.HeadFile(os.path.join(mf6_dir, f"{sim_name}.hds")).get_data()
        highland_locs_mask = (top > 50.0) & (ibound == 1) & (np.abs(head_for_mp[0]) < 1000)
        rows_mp, cols_mp = np.where(highland_locs_mask)

        # 10セルに1つの割合で間引き
        locs = [(0, r, c) for r, c in zip(rows_mp[::10], cols_mp[::10])]

        mp_name = f"{sim_name}_mp7"
        mp = flopy.modpath.Modpath7(
            modelname=mp_name,
            flowmodel=gwf,
            exe_name="mp7",
            model_ws=mp7_dir
        )

        # 深層の空隙率を増加して深層流動を促進
        mpbas = flopy.modpath.Modpath7Bas(
            mp,
            porosity=[0.25, 0.20, 0.10, 0.08]  # L3,L4の空隙率を増加
        )

        particledata = flopy.modpath.ParticleData(
            locs,
            structured=True,
            drape=0,
            localx=0.5,
            localy=0.5,
            localz=0.9
        )
        pg = flopy.modpath.ParticleGroup(
            particledata=particledata,
            particlegroupname="Highland_Recharge"
        )

        mpsim = flopy.modpath.Modpath7Sim(
            mp,
            simulationtype="pathline",
            trackingdirection="forward",
            weaksinkoption="pass_through",
            weaksourceoption="pass_through",
            particlegroups=[pg]
        )

        mp.write_input()
        logger.info(f"MODPATH 7 configured with {len(locs)} release points")

        mp_success, _ = mp.run_model(silent=True)
        if mp_success:
            logger.info("MODPATH 7 completed successfully!")

            # 結果の集計
            p_file = os.path.join(mp7_dir, f"{mp_name}.mppth")
            if os.path.exists(p_file):
                plines = flopy.utils.PathlineFile(p_file).get_alldata()
                deep = [p for p in plines if np.min(p['z']) < -100]
                shallow = [p for p in plines if np.min(p['z']) >= -100]
                logger.info(f"Pathline summary: {len(plines)} total, {len(shallow)} shallow, {len(deep)} deep (z < -100m)")
        else:
            logger.warning("MODPATH 7 did not complete successfully")

    except Exception as e:
        logger.error(f"Error in MODPATH 7: {e}")

# ============================================================================
# Section 16: Visualization (Plan View and Cross-Section)
# ============================================================================

logger.info("Creating visualizations...")

if not args.dry_run and hds_exists:
    try:
        # Read heads
        head_file = os.path.join(mf6_dir, f"{sim_name}.hds")
        if os.path.exists(head_file):
            hds = bf.HeadFile(head_file)
            heads = hds.get_data(totim=hds.get_times()[-1])

            # Plan view of layer 1 heads using flopy PlotMapView
            fig, ax = plt.subplots(figsize=(14, 10))
            pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
            
            # Plot valid heads
            heads_ma = np.ma.masked_where((heads[0] > 1e10) | (heads[0] < -100) | (ibound == 0), heads[0])
            im = pmv.plot_array(heads_ma, cmap='viridis', alpha=0.8, vmin=0, vmax=150)

            # Plot pathlines!
            p_file = os.path.join(mp7_dir, f"{sim_name}_mp7.mppth")
            if os.path.exists(p_file):
                plines = flopy.utils.PathlineFile(p_file).get_alldata()
                pmv.plot_pathline(plines, layer='all', colors='white', lw=1.0, alpha=0.6, label='Pathlines')

            # Overlay streams (DRN cells from stream_grid)
            ax.contour(
                np.linspace(xmin, xmax, ncol), np.linspace(ymax, ymin, nrow),
                stream_grid, levels=[0.5], colors='cyan', linewidths=0.5
            )

            # Overlay observation points (FG9 is already filtered out)
            if river_gdf is not None and not river_gdf.empty:
                river_gdf.plot(ax=ax, color='red', markersize=50, marker='s', label='River obs', alpha=0.7, zorder=5)
            if gw_gdf is not None and not gw_gdf.empty:
                gw_gdf.plot(ax=ax, color='blue', markersize=50, marker='o', label='GW obs', alpha=0.7, zorder=5)

            ax.set_title(f"MODFLOW 6 v3 - Layer 1 Hydraulic Head & Pathlines", fontsize=14)
            ax.set_xlabel("UTM Easting (m)")
            ax.set_ylabel("UTM Northing (m)")
            cbar = plt.colorbar(im, ax=ax, label="Head (m)")
            ax.legend(loc='lower right')

            plan_view_file = os.path.join(output_dir, 'modflow_v3_plan_view.png')
            plt.savefig(plan_view_file, dpi=300, bbox_inches='tight')
            logger.info(f"Saved plan view: {plan_view_file}")
            plt.close()

            # Cross-section along column nc/2
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

            col_idx = ncol // 2
            xs = np.arange(nrow) * 100.0 + ymin

            # Plot elevation and heads
            ax1.plot(xs, top[:, col_idx], 'k-', linewidth=2, label='DEM')
            for l in range(nlay):
                ax1.plot(xs, botm[l, :, col_idx], 'k--', linewidth=1, alpha=0.5)
            for l in range(nlay):
                ax1.plot(xs, heads[l, :, col_idx], label=f'L{l+1}', linewidth=1.5)

            ax1.set_ylabel("Elevation (m)")
            ax1.set_title(f"Cross-section at Column {col_idx}")
            ax1.legend()
            ax1.grid()

            # Plot K distribution
            ax2.plot(xs, k1[:, col_idx], label='K1', linewidth=1.5)
            ax2.plot(xs, k2[:, col_idx], label='K2', linewidth=1.5)
            ax2.plot(xs, k3[:, col_idx], label='K3', linewidth=1.5)
            ax2.plot(xs, k4[:, col_idx], label='K4', linewidth=1.5)

            ax2.set_xlabel("Northing (m)")
            ax2.set_ylabel("Hydraulic Conductivity (m/day)")
            ax2.set_title("Hydraulic Conductivity Distribution")
            ax2.legend()
            ax2.grid()

            xsec_file = os.path.join(output_dir, 'modflow_v3_cross_section.png')
            plt.savefig(xsec_file, dpi=300, bbox_inches='tight')
            logger.info(f"Saved cross-section: {xsec_file}")
            plt.close()

    except Exception as e:
        logger.error(f"Error in visualization: {e}")

# ============================================================================
# Summary and Completion
# ============================================================================

logger.info("\n" + "="*70)
logger.info("MODFLOW 6 v3 Workflow Completed")
logger.info("="*70)
logger.info(f"Output directory: {output_dir}")
logger.info(f"MODFLOW workspace: {mf6_dir}")
logger.info(f"MODPATH workspace: {mp7_dir}")
logger.info(f"Log file: fiji_modflow_v3.log")
logger.info("="*70)

# End of script
