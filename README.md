# Rakiraki Groundwater MODFLOW 6 Model

Steady-state groundwater flow model for the Rakiraki region in northern Viti Levu, Fiji, built with MODFLOW 6 and FloPy.

## Overview

This project simulates groundwater flow and submarine groundwater discharge (SGD) in the Rakiraki coastal watershed using MODFLOW 6 with MODPATH 7 particle tracking. The model is constrained by field observations including river discharge measurements and groundwater level data collected during the 2025-2026 field campaigns.

### Model Specifications

| Parameter | Value |
|-----------|-------|
| Domain | 24 km x 20 km (611,000-635,000 E / 8,068,000-8,088,000 N, UTM Zone 60S) |
| Grid | 200 rows x 240 columns, 100 m cell size |
| Layers | 4 (L1: 50m unconfined, L2: 50m confined, L3: 150m confined, L4: 450m deep confined) |
| Solver | NEWTON with UNDER_RELAXATION, COMPLEX IMS (ILUT + BICGSTAB) |
| Packages | DIS, NPF, IC, CHD, RCH, DRN, OC + MODPATH 7 |
| CRS | EPSG:32760 (WGS 84 / UTM Zone 60S) |

## Project Structure

```
.
├── Fiji_MODFLOW_v3.py              # Main model script (improved v3)
├── Fiji_MODFLOW_observation_constrained_v2.py  # Previous version (v2)
├── Fiji_field_template_v10.xlsx    # Field observation data (Excel)
├── RakiRaki_FieldSurve_20250929.xlsx  # Field survey data
├── output/                         # Model input/output data
│   ├── model_top.csv              # DEM (model top elevation)
│   ├── model_ibound.csv           # Active/inactive cell mask
│   ├── modflow_domain.json        # Domain configuration
│   ├── basins.shp                 # Watershed boundaries
│   ├── basins_clean.shp           # Cleaned watershed boundaries
│   ├── streams.shp                # Stream network
│   └── *.png                      # Visualization outputs
├── mf6_workspace/                  # MODFLOW 6 input/output files
├── mp7_workspace/                  # MODPATH 7 input/output files
└── .gitignore
```

## Requirements

- Python 3.8+
- FloPy (`pip install flopy`)
- MODFLOW 6 and MODPATH 7 executables (install via `python -m flopy.utils.get_modflow`)
- NumPy, Pandas, GeoPandas, Matplotlib, Rasterio, Shapely

### Installation

```bash
pip install flopy numpy pandas geopandas matplotlib rasterio shapely openpyxl
python -m flopy.utils.get_modflow :flopy
```

## Usage

### Run the full workflow

```bash
python Fiji_MODFLOW_v3.py
```

### Dry run (no simulation, just file writing)

```bash
python Fiji_MODFLOW_v3.py --dry-run
```

### Skip MODPATH particle tracking

```bash
python Fiji_MODFLOW_v3.py --skip-modpath
```

## Data Sources

- **DEM**: ALOS 30m DSM (JAXA), resampled to 100m
- **Bathymetry**: GEBCO 2025
- **Hydrology**: PySheds-derived stream network and watersheds
- **Field data**: GNSS survey, river cross-section measurements, groundwater level observations (2025-2026 fieldwork)

## Calibration Results (v3)

| Metric | Result | Target |
|--------|--------|--------|
| FG1 head error | -2.10 m | < 5 m |
| Deep flowpaths (z < -100m) | 28.2% (51/181) | 10-30% |
| Water budget discrepancy | 6.11% | < 1% (not yet achieved) |
| Convergence | Not converged (800 iterations) | Normal termination |

The model successfully reproduces observed groundwater heads and deep groundwater flow paths. Convergence remains an ongoing optimization target due to the ill-conditioned coefficient matrix from large K contrasts and layer thickness variations.

## Key Features in v3

1. **DRN package**: Replaced SFR (Streamflow Routing) with DRN (Drain) for stable river-aquifer interaction with the NEWTON solver
2. **Observation data**: Reads directly from River_Points, Groundwater, and Controls sheets with GNSS elevation calculation (APC_h - antenna_height - geoid_height)
3. **Deep groundwater flow**: Calibrated K values to promote deep flow paths (K3=8.0, K4=4.0 m/d) with reduced shallow K (K1=2.0 m/d)
4. **MODPATH 7 particle tracking**: Particles released from highland areas (elevation > 50m) to trace shallow and deep groundwater flowpaths
5. **Portable paths**: Relative paths for cross-platform compatibility

## License

This project is part of the Fiji groundwater research initiative.

## Contact

Jun Yasumoto (junysmt@gmail.com)
