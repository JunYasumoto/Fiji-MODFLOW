# Phase 4: MODFLOW 6 Results & Calibration Report

## 1. Validation vs Observed Points (FR/FG)
- **Points Assessed**: 14 (Excluding offshore No1-15)
- **RMSE**: 45.68 m
- **MAE**: 29.23 m
*(Note: Due to lack of measured piezometric heads, ground Z-value was used as proxy for shallow groundwater/river heads.)*

## 2. Submarine Groundwater Discharge (SGD) Analysis
 SGD and coastal groundwater discharge flows were extracted from the `.cbc` file as a combination of CHD outflows at the northern boundary and upward vertical flux (`FLOW LOWER FACE` Layer 1) in marine topography (`TOP < 0`).
**Malake Island Channel Focus**: The upward flux in the shallow bathymetric channel between the main island and Malake Island totals approx **0.00 m³/day**.

## 3. Global Water Budget Information
```text
  VOLUME BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   1
  ---------------------------------------------------------------------------------------------------

     CUMULATIVE VOLUME      L**3       RATES FOR THIS TIME STEP      L**3/T          PACKAGE NAME    
     ------------------                 ------------------------                     ----------------

           IN:                                      IN:
           ---                                      ---
                RCHA =      300348.0000                  RCHA =      300348.0000     RCHA_0                           
                 CHD =           0.0000                   CHD =           0.0000     CHD_0                            
                 SFR =     5144586.7641                   SFR =     5144586.7641     SFR_0                            

            TOTAL IN =     5444934.7641              TOTAL IN =     5444934.7641

          OUT:                                     OUT:
          ----                                     ----
                RCHA =           0.0000                  RCHA =           0.0000     RCHA_0                           
                 CHD =       57559.3883                   CHD =       57559.3883     CHD_0                            
                 SFR =     5384909.9740                   SFR =     5384909.9740     SFR_0                            

           TOTAL OUT =     5442469.3623             TOTAL OUT =     5442469.3623

            IN - OUT =        2465.4018              IN - OUT =        2465.4018

 PERCENT DISCREPANCY =           0.05     PERCENT DISCREPANCY =           0.05




         TIME SUMMARY AT END OF TIME STEP    1 IN STRESS PERIOD    1
                    SECONDS     MINUTES      HOURS       DAYS        YEARS
                    -----------------------------------------------------------
   TIME STEP LENGTH  86400.      1440.0      24.000      1.0000     2.73785E-03
 STRESS PERIOD TIME  86400.      1440.0      24.000      1.0000     2.73785E-03
         TOTAL TIME  86400.      1440.0      24.000      1.0000     2.73785E-03


end timestep


```
