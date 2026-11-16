# NOW-Project

This repository contains the R processing and analysis scripts associated with the manuscript **“The North Water Polynya in Transition: Multi-Decadal Changes of Hydrographic and Nutrient Dynamics”** (Schiffrine & Tremblay; submitted to *Journal of Geophysical Research: Oceans*).

The North Water Polynya (NOW), at the confluence of Nares Strait and Baffin Bay, is a productive Arctic region where Pacific- and Atlantic-origin waters interact. The analysis concerns hydrographic and biogeochemical observations from ArcticNet expeditions, with the study transect selected by latitude and longitude bounds. The scripts describe coverage over 2005–2023. The workflow standardizes multi-year nutrient data, extracts and clusters transect stations, calculates TEOS-10 seawater properties, classifies water masses in potential-density/spicity space, and derives nutrient and water-mass tracers.

Five candidate endmembers are defined in the classification code: Polar Mixed Water (PMW), Upper Halocline Water (UHW), Canadian Basin Atlantic Water (CBAW), Baffin Bay Polar Water (BBPW), and Subpolar Mode Water (SPMW). **The current master, `20261006_NOW.R`, explicitly removes CBAW before classification**, so its active classifier uses four endmembers. The final factor-level list also omits CBAW. See the code notes below for further implementation details and limitations.

## README — `20261006_NOW.R` (master pipeline)

### Purpose

`20261006_NOW.R` is the current master/orchestration script among the supplied versions. It creates the in-memory analysis object `NOW`. It is not a self-contained command-line pipeline as uploaded: it sets a personal absolute working directory, expects source and raw-data files at particular relative paths, and the final `select()` is commented out. No save/export of `NOW` is performed in this script.

### File inventory

The inventory below covers the 19 distinct analysis files/data files in the supplied set. Repeated upload copies with suffixes such as `(1)` and `(2)` are duplicates or alternate revisions, not additional pipeline modules. Names in the table are the canonical filenames referenced by the code.

| File | Role / status |
|---|---|
| `20261006_NOW.R` | **Current master.** Loads packages and helper scripts, loads and filters data, clusters stations, calculates TEOS-10 properties, classifies water masses and derives tracers. |
| `20260209_NOW.R` | **Legacy/superseded master.** Earlier pipeline; uses the older salinity/temperature spicity workflow and contains an active final `select()` that the current master has commented out. |
| `20261006_NOW_Figure.R` | **Current supplied figure script.** Builds the overview, water-mass diagnostic, vertical structure and temporal-trend figures (labelled manuscript Figures 1–4 in the script), plus additional/supplementary plots and summaries. Requires `NOW` and other objects in the R session. |
| `20260209_NOW_Figure.R` | **Legacy/superseded figure script** relative to `20261006_NOW_Figure.R`. Earlier plotting version. |
| `20260209_Package.R` | Defines `required_packages` and `load_now_packages()`; sourcing the file also calls the loader immediately. |
| `20260209_DataLoading.R` | Defines `load_data()`, which reads and standardizes the three external ArcticNet tab-delimited datasets listed below. |
| `20260209_Extract_Transect.R` | Defines `extract_transect()`; applies global `NOW_BOUNDS`. |
| `20260209_ClusterStation.R` | Defines `cluster_stations()` and `rename_clusters()` for hierarchical station clustering and cluster naming. |
| `20260209_sw_pspi.r` | Defines the older/custom `sw_pspi()` spicity implementation. Sourced by both masters, but the current master calculates `spicity` with `gsw_pspi()` instead. |
| `20260209_WaterMass_Classif.R` | Defines `create_endmembers()` and `classify_watermass()`. |
| `20260209_PacificWater.R` | Defines `ANP()`, `fpw()` and `POs_star()`. |
| `20260209_ratio_calc.R` | Defines `LOD` and detection-limit-aware `ratio_calc()`. |
| `20260209_WaterMass_Sensitivity.R` | Optional sensitivity functions: `test_ratios()`, `test_thresholds()`, `test_parameter_grid()` and `calculate_performance_metrics()`. Not sourced by the current master; its example call is commented out there. |
| `20260209_SensPlot.R` | Optional `plot_trend()` trend-plot function. Sourced by the current figure script, not the master. |
| `20260209_query_wm.R` | Optional summary/query helpers `query_wm_stat()` and `query_wm_range_global()`. |
| `20260209_Depth-layer_test.R` | Optional depth-strata and station pairwise Dunn-test/heat-map analysis code; not sourced by the master. |
| `20260209_pw_regime_calc.R` | Optional annual Pacific Water regime analysis function `pw_regime_calc()`; sourced by the figure script, not the master. |
| `20260209_PlotWaterMassSection.R` | Optional water-mass section plot function `plot_watermass_section()`; sourced by the figure script, not the master. |
| `NOW_bathymetry.csv` | Bathymetric coordinates (`decimalLongitude`, `decimalLatitude`, `Depth`) used by figure-generation code. Not an input to the master processing pipeline. |

The figure script additionally sources `Script_R/sen_slope.R` and `Script_R/plot_section.R`. Those files are **not** among the 19 analysis files listed above, so the supplied figure script will require those dependencies to be present in the repository/session.

### Required R packages

The exact vector in `20260209_Package.R` is:

```r
c(
  "broom", "colorspace", "factoextra", "fields", "forcats",
  "GGally", "ggcorrplot", "ggh4x", "ggnewscale", "ggoce",
  "ggOceanMaps", "ggplot2", "ggpmisc", "ggpubr", "ggrepel",
  "ggspatial", "grid", "gridExtra", "gsw", "lubridate",
  "MBA", "metR", "pals", "RColorBrewer", "rstatix",
  "see", "tibble", "tidyverse", "trend", "legendry"
)
```

Sourcing `20260209_Package.R` calls `load_now_packages(required_packages)`. By default, `load_now_packages()` attempts to install packages not already available (`install_missing = TRUE`). The figure script and optional analyses may also need packages or sourced helpers beyond what the master pipeline itself uses.

### Data sources expected by `load_data()`

`load_data(parameters = parameter)` reads these files using paths relative to the current working directory:

| Expected path | Coverage indicated in source comments |
|---|---|
| `DataBase_Nut/Data/R_Database_dec2021.txt` | Historical data through 2021 (commented as 2005–2021) |
| `DataBase_Nut/Data/R_Data_2022.txt` | 2022 expedition data |
| `DataBase_Nut/Data/R_Data_2023.txt` | 2023 expedition data |

For each file, `load_data()` reads a tab-delimited file with empty strings as `NA`, forward-fills `Cruise:Bot.Depth`, parses `eventDate`, extracts `year`, `month`, `day`, and `DayOfYear`, subtracts 360 from `decimalLongitude`, selects metadata plus `Strate` (if present) and requested parameters, sorts and removes exact duplicates, then row-binds the files. It does **not** explicitly coerce each nutrient field to numeric in the function shown.

The master requests these parameters:

```r
parameter <- c(
  "Temperature", "Salinity", "Nitrate", "Nitrite", "Ammonium", "Oxconc",
  "Phosphate", "Silicate"
)
nut_vars <- c("Nitrate", "Nitrite", "Ammonium", "Oxconc", "Phosphate", "Silicate")

NOW_BOUNDS <- list(
  longitude = c(-80, -65),
  latitude  = c(76.2, 76.5)
)
```

`extract_transect(data = raw_data)` uses inclusive `dplyr::between()` filters against those longitude and latitude limits.

### Running the pipeline

The master currently contains a machine-specific `setwd()` pointing to `/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/`. Its `source()` statements expect helper scripts under `DataBase_Nut/20260209/`, and `load_data()` expects the raw text files under `DataBase_Nut/Data/`. These paths are not portable as written. The upload set also does not contain the three raw data files, so it is insufficient on its own to run the data-loading step.

For a repository-root setup, place the source files and external datasets in the paths expected above and replace the personal `setwd()` with the local repository root (or launch R with that root as the working directory). Alternatively, edit the `source()` and data-file paths to match the actual directory layout. Do not assume the uploaded flat file bundle already has the expected folder hierarchy. The figure script has a separate personal absolute path for `NOW_bathymetry.csv`; change it to a repository-relative path (for example, `DataBase_Nut/Data/NOW_bathymetry.csv`) before running that script.

After path and data setup, run `20261006_NOW.R` in an R session. It clears the current workspace with `rm(list = ls())`; execute it in a dedicated session if you need to preserve other objects. The script leaves `NOW` in memory but does not write an output file. The figure script is a separate downstream script and assumes `NOW`, `wm_ref` and supporting objects have already been created in the session; it does not source the master or call `load_data()` itself.

### Workflow — actual steps in `20261006_NOW.R`

#### STEP 0 — Environment and paths

```r
rm(list = ls())
setwd(dir = "/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/")
```

The first line removes all objects from the workspace. The second is a personal absolute path and must be adapted outside the author's environment.

#### STEP 1 — Inline configuration

The master defines `parameter`, `nut_vars`, and `NOW_BOUNDS` inline before sourcing helpers. The bounding box is −80 to −65 longitude and 76.2 to 76.5 latitude. `nut_vars` identifies the six nutrient/oxygen fields later preserved as `_L` columns and converted by division by `rho`.

#### STEP 2 — Package loader and function sources

The master sources, in order:

```r
source("DataBase_Nut/20260209/20260209_Package.R")
source("DataBase_Nut/20260209/20260209_DataLoading.R")
source("DataBase_Nut/20260209/20260209_Extract_Transect.R")
source("DataBase_Nut/20260209/20260209_ClusterStation.R")
source("DataBase_Nut/20260209/20260209_sw_pspi.r")
source("DataBase_Nut/20260209/20260209_WaterMass_Classif.R")
source("DataBase_Nut/20260209/20260209_PacificWater.R")
source("DataBase_Nut/20260209/20260209_ratio_calc.R")
```

The package file both defines and invokes its package loader. `20260209_sw_pspi.r` is sourced, but its `sw_pspi()` function is not used in the current master.

#### STEP 3 — Load and extract observations

```r
raw_data <- load_data()
raw_now <- extract_transect(data = raw_data)
```

`raw_data` is the combined, standardized data returned by `load_data()`. `raw_now` is the geographic subset. `extract_transect()` prints a count and retention-rate summary; its code refers to `raw_data` in that message even though the formal argument is `data`, so this reporting line depends on an object named `raw_data` being present in the global environment (as it is when called in the master).

#### STEP 4 — Cluster stations and attach labels

```r
now_clusters <- cluster_stations(
  data = raw_now,
  n_clusters = 15,
  min_stations_per_cluster = 7
)
now_renamed <- rename_clusters(now_clusters)

now_geo <- raw_now %>%
  left_join(now_renamed,
    by = c("Station", "decimalLongitude", "decimalLatitude", "year", "month", "day")
  ) %>%
  select(Cruise, Station, Station_Cluster, everything(), -Cluster) %>%
  filter(!is.na(Station_Cluster))
```

`cluster_stations()` extracts distinct station/coordinate/date combinations, computes `dist()` on latitude and longitude values, applies `hclust(..., method = "ward.D2")`, cuts the tree at 15, and sets `Cluster` to `NA` when its count is below 7. The coordinates are treated as ordinary Euclidean coordinates by `dist()` (not geodesic distances). `rename_clusters()` chooses each cluster's most frequent station label and resolves duplicate dominant labels using longitude-based `W`/`E` or numeric suffixes. It does not itself assign the fixed numeric west-to-east identifiers claimed in the old README. The join then drops observations without a `Station_Cluster` assignment.

#### STEP 5 — Endmembers and TEOS-10 hydrography

The master calculates median longitude and latitude of `now_geo` as `ref_coords`, creates endmembers, and removes CBAW:

```r
ref_coords <- list(
  longitude = median(now_geo$decimalLongitude, na.rm = TRUE),
  latitude = median(now_geo$decimalLatitude, na.rm = TRUE)
)

wm_ref <- create_endmembers(
  reference_lon = ref_coords$longitude,
  reference_lat = ref_coords$latitude
) %>%
  filter(!watermass %in% "CBAW")
```

`create_endmembers()` defines five sets of `SA`, `CT` and their standard deviations in `20260209_WaterMass_Classif.R`. The master then removes the `CBAW` row. `now_hydro` is constructed with the following TEOS-10 calculations:

```r
Pressure = gsw_p_from_z(z = -Depth, latitude = decimalLatitude)
SA       = gsw_SA_from_SP(SP = Salinity, p = Pressure,
                          longitude = decimalLongitude, latitude = decimalLatitude)
CT       = gsw_CT_from_t(SA = SA, t = Temperature, p = Pressure)
pt       = gsw_pt_from_CT(SA = SA, CT = CT)
sigma0   = gsw_sigma0(SA = SA, CT = CT)
rho      = gsw_rho(SA = SA, CT = CT, p = Pressure) / 1000
spicity  = gsw_pspi(SA = SA, CT = CT, pr = 0)
```

Here `z` is supplied as negative depth in metres. `pt` is calculated but oxygen solubility/AOU lines using it are commented out. After the calculations, rows are filtered to retain nonmissing `Salinity`, `Temperature`, `sigma0`, and `spicity`; a `stopifnot()` checks for missing values in `SA`, `CT`, `sigma0`, and `spicity` in the retained data. `rho` is in kg L⁻¹ according to the code's division by 1000.

#### STEP 6 — Water-mass classification

The active call is:

```r
now_wm <- classify_watermass(
  now_hydro,
  wm_ref %>% filter(!watermass %in% "CBAW"),
  threshold = 0.6,
  comparable_ratio = 1.3,
  return_long = F
)
```

The CBAW exclusion is present both when creating `wm_ref` and in the classifier argument. The classification function uses Euclidean distances in `(sigma0, spicity)` coordinates to each supplied endmember; the code does not scale by uncertainty ellipses. With `return_long = FALSE`, it appends `classification` and the nearest-endmember `distance` to the valid input rows. Its code labels a sample `Mixed` if either (a) at least two endmembers have comparable distances according to `comparable_ratio` or (b) all endmember distances exceed `threshold`. Otherwise a nearest endmember is assigned only when its distance is strictly less than `threshold`; in the remaining case the result is `Unclassified`. Thus the actual logic is not simply “nearest wins unless the ratio is small.”

The master also prints a frequency/percentage summary of `classification_summary`. Since the active endmember table excludes CBAW, CBAW cannot be returned by this call. The available endmember names are PMW, UHW, BBPW, and SPMW; `Mixed` and potentially `Unclassified` are classification outcomes.

#### STEP 7 — Derived variables and final `NOW`

The master mutates `now_wm` into `NOW` as follows (formulae abbreviated only for line wrapping):

```r
across(any_of(nut_vars), ~ .x, .names = "{.col}_L")
across(any_of(nut_vars), ~ .x / rho)
Nstar = 0.87 * (Nitrate - 16 * Phosphate + 2.9)
TIN = Nitrate + replace_na(Nitrite, 0) + replace_na(Ammonium, 0)
NP  = ratio_calc(Nitrate, Phosphate, LOD["Nitrate"], LOD["Phosphate"])
SiN = ratio_calc(Silicate, Nitrate, LOD["Silicate"], LOD["Nitrate"])
SiP = ratio_calc(Silicate, Phosphate, LOD["Silicate"], LOD["Phosphate"])
ANP = ANP(Phosphate = Phosphate, Nitrogen = TIN)
fPW = fpw(Phosphate = Phosphate, Nitrogen = TIN,
          atlantic_ref = "Jones2008", pacific_ref = "Yamamoto.Kawai2008")
NO = (9 * Nitrate) + Oxconc
PO = (135 * Phosphate) + +Oxconc
NO_PO = NO / PO
POs_star = POs_star(phosphate = Phosphate, oxygen = Oxconc, salinity = Salinity)
```

The ratio helper uses `LOD = c(Nitrate = 0.05, Phosphate = 0.05, Silicate = 0.10)` and returns `NA` when either input is below its corresponding LOD. `fpw()` uses the Jones **2008** Atlantic reference and Yamamoto-Kawai 2008 Pacific reference in this master, with its default `clip = TRUE`; it is not using Jones 1998. `ANP()` uses its own default Atlantic/Pacific line coefficients. `POs_star()` applies the salinity-dependent blend of its 175:1 and 125:1 oxygen-to-phosphate forms as implemented in `20260209_PacificWater.R`.

The six fields in `nut_vars` are first copied to `Nitrate_L`, `Nitrite_L`, `Ammonium_L`, `Oxconc_L`, `Phosphate_L`, and `Silicate_L`, then the original six columns are divided by `rho`. The `_L` suffix therefore denotes preserved pre-division values in the current implementation; the source labels are not documented as units. The code then factors `classification` with levels `PMW`, `UHW`, `BBPW`, `SPMW`, `Mixed` (CBAW is commented out) and `Station_Cluster` with levels `100`, `101`, `103`, `105`, `107`, `108`, `110`, `111`, `113`, `115`, `116`. Any classification or cluster labels outside those levels become `NA` upon factor conversion.

Finally, `filter(!year %in% c(1997, 1998, 1999, 2007, 2008))` removes those years. It does not restrict the data to a stated start/end range such as 2005–2023; it only excludes this explicit set. **There is no active final `select()` in `20261006_NOW.R`**: the `select()` block on lines 260–266 is commented out. Consequently, `NOW` retains all upstream columns, as well as the added columns and changed nutrient fields described above.

### Script dependency map

```text
20261006_NOW.R (master; expects raw data and path hierarchy)
  ├── 20260209_Package.R          → required_packages, load_now_packages()
  ├── 20260209_DataLoading.R      → load_data()
  ├── 20260209_Extract_Transect.R → extract_transect()
  ├── 20260209_ClusterStation.R   → cluster_stations(), rename_clusters()
  ├── 20260209_sw_pspi.r          → sw_pspi() [sourced but unused by current master]
  ├── 20260209_WaterMass_Classif.R → create_endmembers(), classify_watermass()
  ├── 20260209_PacificWater.R     → ANP(), fpw(), POs_star()
  └── 20260209_ratio_calc.R       → LOD, ratio_calc()

20261006_NOW_Figure.R (separate downstream script; assumes NOW and wm_ref)
  ├── 20260209_PlotWaterMassSection.R → plot_watermass_section()
  ├── 20260209_SensPlot.R             → plot_trend()
  ├── 20260209_pw_regime_calc.R       → pw_regime_calc()
  ├── Script_R/sen_slope.R            → external, not in supplied analysis-file set
  └── Script_R/plot_section.R         → external, not in supplied analysis-file set

Optional / not sourced by current master:
  ├── 20260209_WaterMass_Sensitivity.R
  ├── 20260209_query_wm.R
  └── 20260209_Depth-layer_test.R
```

The custom `sw_pspi()` implementation belongs to the legacy workflow; the current master calls `gsw_pspi()` from package `gsw`. `20260209_PlotWaterMassSection.R`, `20260209_SensPlot.R`, and `20260209_pw_regime_calc.R` are not sourced by the master, but are sourced by the current figure script.

### Output object structure

The current master constructs `NOW` without an active final column selection, so the definitive schema is inherited from `load_data()` and the intervening mutations rather than a fixed select list. The core expected groups are:

| Group | Columns added or expected in `NOW` |
|---|---|
| Expedition / station metadata | `Cruise`, `Station`, `Station_Cluster`, `Type`, `eventDate`, `year`, `month`, `day`, `DayOfYear`, `decimalLongitude`, `decimalLatitude`, `Depth`; `Strate` where present in input. `Cluster` is dropped when `now_geo` is assembled. |
| Original requested observations | `Temperature`, `Salinity`, and the six `nut_vars`: `Nitrate`, `Nitrite`, `Ammonium`, `Oxconc`, `Phosphate`, `Silicate`. In `NOW`, each of the six nutrient/oxygen originals is divided by `rho`; corresponding `_L` copies are also retained. |
| TEOS-10 / physical properties | `Pressure`, `SA`, `CT`, `pt`, `sigma0`, `rho`, `spicity`. |
| Classification | `classification`, `distance` (nearest supplied endmember distance). Classification factor levels exclude CBAW and include `Mixed`. |
| Derived nutrient indices | `Nstar`, `TIN`, `NP`, `SiN`, `SiP`, `ANP`, `fPW`, `NO`, `PO`, `NO_PO`, `POs_star`. |
| Preserved input-unit fields | `Nitrate_L`, `Nitrite_L`, `Ammonium_L`, `Oxconc_L`, `Phosphate_L`, `Silicate_L` (pre-division copies as coded). |

`NOW` is an R object only; the current master does not save it to disk. The older master has an active `select()` and a different final schema, so that older selection should not be used to infer current output columns.

### What's new (Oct/Nov 2026 update)

Compared with `20260209_NOW.R`, the current master:

- Adds `nut_vars` and retains pre-division versions of six nutrient/oxygen variables under `_L` names, while dividing the original columns by the calculated `rho`.
- Replaces the earlier simplified `sigma0`/custom `sw_pspi()` pathway with pressure, Absolute Salinity, Conservative Temperature, potential temperature, in-situ density, and `gsw_pspi()` calculations using TEOS-10 functions.
- Filters samples missing salinity, temperature, density anomaly or spicity and runs an explicit `stopifnot()` check on `SA`, `CT`, `sigma0`, and `spicity`.
- Excludes CBAW from the reference endmembers used to classify samples; it also activates and prints the classification count/percentage summary.
- Keeps the same threshold values (`0.6` and `1.3`) and same `Jones2008` / `Yamamoto.Kawai2008` references for `fPW` as actually used in the newer master.
- Comments out the former final `select()` rather than replacing it with another selection; more columns therefore remain in `NOW`.

Compared with `20260209_NOW_Figure.R`, `20261006_NOW_Figure.R` switches the principal hydrographic diagram to `SA`/`CT` and corresponding TEOS-10 labels/units, suppresses CBAW from the displayed water-mass palette/order and endmember overlays, adds water-mass distance/statistics export/table code and a Pacific-water reference-method comparison, and has some newer supplementary figure code. Both supplied figure scripts identify the principal manuscript products as Figure 1 (transect overview), Figure 2 (water-mass diagnostics), Figure 3 (vertical structure) and Figure 4 (temporal trends). Most figure-script `ggsave()` calls are commented out, but two are active in the supplied script: they write the Pacific-water boxplot and the Si/P supplementary figure (to the configured figure directory). Other plot objects are available for interactive review unless their export calls are uncommented.

### Code and documentation caveats

- The current master’s `setwd()` and the figure script’s bathymetry read path are hard-coded to the author's OneDrive directory. Source/data paths also assume a specific folder tree.
- The three raw text datasets expected by `load_data()` are not included in the supplied analysis-file set. The supplied `NOW_bathymetry.csv` is for plotting and does not replace those inputs.
- `extract_transect()`'s diagnostic uses the global name `raw_data`, rather than its argument `data`.
- `create_endmembers()` defines five candidates but the current master filters out CBAW for classification. The code uses direct Euclidean distances in `(sigma0, spicity)` and not endmember uncertainty-ellipse scaling.
- The current master sets `Station_Cluster` factor levels to numeric-looking labels, while `rename_clusters()` constructs labels from dominant station names and conflict suffixes. Verify the actual input station labels before relying on the fixed factor conversion; unmatched labels turn into `NA`.
- The code comments describe `fPW` as a fraction and the function returns clipped values between 0 and 1 by default. The `mutate()` block does not multiply it by 100.
- The current master contains `PO = (135 * Phosphate) + +Oxconc`; in R the repeated unary `+` is equivalent to adding `Oxconc`, but is likely a typographical artifact.
- `gsw_pt_from_CT()` produces a property that is calculated but unused in the active pipeline; oxygen solubility and AOU calculations remain commented out.
