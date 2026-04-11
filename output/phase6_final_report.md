# Phase 6: Topography-Dependent Zoning & Final Calibration

## 1. Zoning Parameters Applied
- **Highland Zone** (Elevation > 50m): K2 doubled to 2.0 m/day, K3 doubled to 0.2 m/day, Recharge reduced to 70% of baseline (0.00189 m/day) to limit flooded upland cells.
- **Lowland Zone** (Elevation <= 50m): K1 set to 12.0 m/day (optimized from Phase 5).

## 2. Sensitivity on SFR Streambed K (`rhk` / `strhc1`)
Multiple permutations of Streambed Conductivity (`rhk`) were tested to evaluate the connection between the aquifer and the streams:
- `rhk = 0.1`: **RMSE = 3.48 m**, Flooded Pct = 20.6%
- `rhk = 1.0`: Failed to converge
- `rhk = 10.0`: Failed to converge
- `rhk = 100.0`: Failed to converge

**Optimal `strhc1 (rhk)` Selected**: 0.1 m/day
- Final FG Points RMSE: **3.48 m** (Vastly improved from previous 9.96 m)
- Highland Flooded Ratio: Improved from 31.7% to **20.6%**.
- Mean deviation of groundwater head from riverbed at FR points: **30.85 m**. (River beds are highly disconnected from deep groundwater tables in steep terrain).

## 3. Submarine Groundwater Discharge (SGD)
The parameter adjustments altered the flow routing dynamics, shifting more water toward the coast due to lower highland recharge and highly conductive lowland pathways.
The Submarine Groundwater Discharge (SGD) map (`Final_SGD_Map.png`) has been updated reflecting these final calibrated gradients.
