# NOW-Project
This repository contains the data processing and analysis pipeline supporting the manuscript "The North Water Polynya in Transition: Multi-Decadal Changes of Hydrographic and Nutrient Dynamics" (Schiffrine & Tremblay, submitted to Journal of Geophysical Research: Oceans).
The North Water Polynya (NOW), located at the confluence of Nares Strait and Baffin Bay, is a biologically productive Arctic region where Pacific-origin and Atlantic-origin waters interact. Understanding how upstream hydrographic changes propagate through this system is essential for predicting shifts in nutrient supply, primary production, and carbon cycling under ongoing Arctic transformation.

This project implements a complete analytical workflow for processing multi-decadal (2005–2023) hydrographic and biogeochemical observations collected during ArcticNet expeditions along a zonal transect across the NOW (~76.2–76.5°N). The pipeline includes: (1) standardization of nutrient databases across multiple expedition years, (2) geographic extraction and hierarchical station clustering for spatial consistency, (3) water mass classification using a σ₀–spicity methodology adapted from Huang et al. (2024), which provides more robust discrimination than traditional θ–S diagrams in high-latitude systems, and (4) derivation of geochemical tracers including N* (Gruber & Sarmiento), Pacific Water fraction (fPW), and Arctic N–P indices.

Five water masses are identified and characterized: Polar Mixed Water (PMW), Upper Halocline Water (UHW), Baffin Bay Polar Water (BBPW), Subpolar Mode Water (SPMW), and Canadian Basin Atlantic Water (CBAW). The classification framework includes sensitivity analysis to optimize discrimination thresholds and ensure robust detection of transitional waters.

Je vois maintenant le fichier réel. Voici la documentation corrigée :

---

# README — `20260209_NOW.R`

## Purpose

`20260209_NOW.R` is the master orchestration script that constructs the final analysis-ready dataset `NOW`. Crucially, **there is no separate `parameters.R` file** — all global constants are defined inline at the top of this script before any function is sourced.

---

## Workflow

---

### STEP 0 — Environment Setup

```r
rm(list = ls())
setwd(dir = "/Users/nicolas/.../")
```

Clears the workspace and sets the working directory to the OneDrive root containing the `DataBase_Nut/` folder.

---

### STEP 1 — Inline Configuration

Defined directly in `20260209_NOW.R`, before any `source()` call:

```r
parameter <- c("Temperature", "Salinity", "Nitrate", "Nitrite", 
               "Ammonium", "Oxconc", "Phosphate", "Silicate")

NOW_BOUNDS <- list(
  longitude = c(-80, -65),
  latitude  = c(76.2, 76.5)
)
```

`parameter` defines the biogeochemical variables extracted from raw files. `NOW_BOUNDS` defines the geographic bounding box delimiting the NOW transect (76.2–76.5°N, 80–65°W). Both objects are global constants available to all downstream functions.

**Note:** Clustering parameters (`n_clusters = 15`, `min_stations_per_cluster = 7`) and classification parameters (`threshold = 0.6`, `comparable_ratio = 1.3`) are passed directly as arguments in the relevant function calls rather than stored as named constants.

---

### STEP 2 — Package Loading

```r
source("DataBase_Nut/20260209/20260209_Package.R")
```

Loads `load_now_packages()` and immediately calls it. This function iterates over `required_packages` — a vector of ~30 packages including `tidyverse`, `gsw`, `ggplot2`, `ggOceanMaps`, `metR`, `lubridate`, `trend` — and installs any missing packages automatically when `install_missing = TRUE`.

---

### STEP 3 — Function Sourcing

```r
source("DataBase_Nut/20260209/20260209_DataLoading.R")
source("DataBase_Nut/20260209/20260209_Extract_Transect.R")
source("DataBase_Nut/20260209/20260209_ClusterStation.R")
source("DataBase_Nut/20260209/20260209_sw_pspi.r")
source("DataBase_Nut/20260209/20260209_WaterMass_Classif.R")
source("DataBase_Nut/20260209/20260209_PacificWater.R")
```

All project-specific functions are loaded into the global environment. No computation occurs at this stage.

---

### STEP 4 — Data Loading

```r
raw_data <- load_data()
```

`load_data()` is sourced from `20260209_DataLoading.R`. It reads and standardizes three tab-delimited ArcticNet source files:

| File | Coverage |
|------|----------|
| `R_Database_dec2021.txt` | 2005–2021 |
| `R_Data_2022.txt` | 2022 |
| `R_Data_2023.txt` | 2023 |

For each file the function:
1. Reads raw data with empty strings coerced to `NA`
2. Forward-fills station metadata to resolve orphan rows
3. Converts Western-positive longitudes to negative-East convention
4. Extracts year, month, day, and day-of-year from date strings
5. Coerces biogeochemical fields to numeric (non-parseable values → `NA`)
6. Selects only the variables listed in `parameter`

The three processed data frames are row-bound and deduplicated into a single standardized database.

**Output:** `raw_data`

---

### STEP 5 — Transect Extraction

```r
raw_now <- extract_transect(data = raw_data)
```

`extract_transect()` is sourced from `20260209_Extract_Transect.R`. It applies a `between()` filter on `decimalLongitude` and `decimalLatitude` using the bounds defined in `NOW_BOUNDS`, retaining only observations within the NOW transect spatial domain. A console report logs original count, retained count, and retention rate.

**Output:** `raw_now`

---

### STEP 6 — Station Clustering

```r
now_clusters <- cluster_stations(
  data = raw_now,
  n_clusters = 15,
  min_stations_per_cluster = 7
)

now_renamed <- rename_clusters(now_clusters)

now_geo <- raw_now %>%
  left_join(now_renamed, by = c("Station","decimalLongitude",
                                 "decimalLatitude","year","month","day")) %>%
  select(Cruise, Station, Station_Cluster, everything(), -Cluster) %>%
  filter(!is.na(Station_Cluster))
```

Two functions are sourced from `20260209_ClusterStation.R`:

**`cluster_stations()`**
Addresses positional offsets between nominally identical stations across expedition years. The function:
1. Extracts unique (Station, lon, lat) combinations
2. Computes a pairwise geographic Euclidean distance matrix
3. Applies Ward's minimum variance agglomerative clustering (`ward.D2`)
4. Cuts the dendrogram at `n_clusters = 15`
5. Discards clusters with fewer than `min_stations_per_cluster = 7` temporal instances

**`rename_clusters()`**
Assigns canonical west-to-east station identifiers:
1. Computes median longitude per cluster
2. Ranks clusters geographically and assigns IDs (100, 101, 103, …, 116)
3. Detects and resolves naming conflicts when multiple historical station names map to a single cluster

Cluster assignments are joined back to `raw_now`; observations outside any valid cluster are discarded via `filter(!is.na(Station_Cluster))`.

**Output:** `now_geo`

---

### STEP 7 — Oceanographic Property Calculation

```r
now_hydro <- now_geo %>%
  mutate(
    sigma0  = gsw_sigma0(SA = Salinity, CT = Temperature),
    spicity = sw_pspi(S = Salinity, temp = Temperature,
                      temp_unit = "conservative", sal_unit = "SA",
                      longitude = ref_coords$longitude,
                      latitude  = ref_coords$latitude, pr = 0)
  )
```

Reference coordinates are first computed as the median longitude and latitude of `now_geo`, then passed to both functions.

`gsw_sigma0()` from the TEOS-10 Gibbs SeaWater toolbox computes potential density anomaly (σ₀, kg m⁻³) referenced to surface pressure.

`sw_pspi()` is sourced from `20260209_sw_pspi.R`. It computes potential spicity (π₀, kg m⁻³) using the 41-coefficient polynomial of Huang et al. (2011). Spicity varies along isopycnals and provides enhanced water mass discrimination in Arctic stratified systems where contrasting T–S combinations share similar densities. Together, σ₀ and π₀ define the two-dimensional classification space used in Step 8.

**Output:** `now_hydro`

---

### STEP 8 — Water Mass Classification

```r
wm_ref <- create_endmembers(
  reference_lon = ref_coords$longitude,
  reference_lat = ref_coords$latitude
)

now_wm <- classify_watermass(
  now_hydro, wm_ref,
  threshold        = 0.6,
  comparable_ratio = 1.3,
  return_long      = FALSE
)
```

Two functions are sourced from `20260209_WaterMass_Classif.R`:

**`create_endmembers()`**
Defines five Arctic water mass endmembers by their canonical T–S properties and converts each to (σ₀, π₀) coordinates using the reference position:

| Water Mass | Code |
|-----------|------|
| Polar Mixed Water | PMW |
| Upper Halocline Water | UHW |
| Baffin Bay Polar Water | BBPW |
| Subpolar Mode Water | SPMW |
| Canadian Basin Atlantic Water | CBAW |

**`classify_watermass()`**
For each observation the algorithm:
1. Computes scaled Euclidean distances to all five endmembers in (σ₀, π₀) space, normalized by each endmember's uncertainty ellipse
2. Identifies the nearest (d₁) and second-nearest (d₂) endmember
3. Applies a two-rule decision:
   - d₁ < `threshold` (0.6) → assign to nearest endmember
   - d₁/d₂ < `comparable_ratio` (1.3) → label **Mixed** (observation lies in a transitional zone equidistant between two water masses)
   - Otherwise → assign to nearest endmember

**Output:** `now_wm` with `classification` column (PMW | UHW | BBPW | SPMW | CBAW | Mixed)

---

### STEP 9 — Geochemical Tracer Derivation and Final Assembly

```r
NOW <- now_wm %>%
  mutate(
    Nstar    = 0.87 * (Nitrate - 16 * Phosphate + 2.9),
    TIN      = Nitrate + replace_na(Nitrite, 0) + replace_na(Ammonium, 0),
    ANP      = ANP(Phosphate, TIN),
    fPW      = fpw(Phosphate, TIN, "Jones1998", "Yamamoto.Kawai2008"),
    NO       = 9   * Nitrate   + Oxconc,
    PO       = 135 * Phosphate + Oxconc,
    NO_PO    = NO / PO,
    POs_star = POs_star(Phosphate, Oxconc, Salinity),
    classification  = factor(classification, 
                             c("PMW","UHW","BBPW","SPMW","CBAW","Mixed")),
    Station_Cluster = factor(Station_Cluster, 
                             c("100","101","103","105","107","108",
                               "110","111","113","115","116"))
  ) %>%
  select(Cruise:Salinity, sigma0:spicity, Nitrate:Silicate,
         TIN, Nstar, NO:POs_star, classification, distance, ANP, fPW) %>%
  filter(!year %in% c(1997, 1998, 1999, 2007, 2008))
```

Three functions are sourced from `20260209_PacificWater.R`:

**`fpw()`**
Estimates the Pacific Water fraction (fPW, 0–1) using a two-endmember PO₄–N mixing model. The Atlantic endmember follows Jones et al. (1998) and the Pacific endmember follows Yamamoto-Kawai et al. (2008). Values approaching 1 indicate predominantly Pacific-origin waters; values near 0 indicate Atlantic-origin dominance.

**`ANP()`**
Computes the Arctic N–P tracer index following Newton et al. (2013). The function calculates orthogonal distances from each observation in (PO₄, TIN) space to both the Atlantic and Pacific N–P regression lines, returning the normalized ratio d_AW / (d_AW + d_PW). ANP = 0 indicates pure Atlantic affinity; ANP = 1 indicates pure Pacific affinity. This approach discriminates water mass origin independently of absolute concentration levels.

**`POs_star()`**
Computes a salinity-normalized phosphate–oxygen combined tracer adapted for Arctic waters. It combines phosphate and dissolved oxygen signals with a salinity normalization to account for freshwater dilution effects characteristic of Arctic shelf and halocline environments.

The inline `mutate()` block additionally computes:

| Tracer | Formula | Physical Meaning |
|--------|---------|-----------------|
| `TIN` | NO₃ + NO₂ + NH₄ | Total inorganic nitrogen; reduced forms set to 0 when missing via `replace_na()` |
| `Nstar` | 0.87 × (NO₃ − 16×PO₄ + 2.9) | Deviation from Redfield N:P stoichiometry after Gruber & Sarmiento (1997); negative values flag Pacific-origin denitrified waters |
| `NO` | 9 × NO₃ + O₂ | Broecker (1974) semiconservative tracer; 9:1 weighting reflects stoichiometric O₂ consumption per mole of NO₃ produced during remineralization |
| `PO` | 135 × PO₄ + O₂ | Analogous semiconservative tracer using Redfield O₂:P ratio of 135:1 |
| `NO_PO` | NO ÷ PO | Deviations from unity signal non-Redfield processes or mixing of water masses with contrasting preformed nutrient ratios |

**Temporal filtering** excludes years 1997–1999 (pre-ArcticNet data with heterogeneous protocols) and 2007–2008 (anomalous sea-ice conditions precluding standard transect occupation).

**Factor ordering** enforces:
- `classification`: PMW → UHW → BBPW → SPMW → CBAW → Mixed (surface-to-deep vertical structure)
- `Station_Cluster`: 100 → 101 → … → 116 (west-to-east geographic ordering)

**Output:** `NOW` — the final analysis-ready data frame consumed by `20260209_NOW_Figure.R`

---

## Script Dependency Map

```
20260209_NOW.R
  │
  ├── 20260209_Package.R          → load_now_packages()
  ├── 20260209_DataLoading.R      → load_data()
  ├── 20260209_Extract_Transect.R → extract_transect()
  ├── 20260209_ClusterStation.R   → cluster_stations(), rename_clusters()
  ├── 20260209_sw_pspi.R          → sw_pspi()
  ├── 20260209_WaterMass_Classif.R → create_endmembers(), classify_watermass()
  └── 20260209_PacificWater.R     → fpw(), ANP(), POs_star()
```

---

## Output Object Structure

| Group | Variables |
|-------|-----------|
| Metadata | Cruise, Station, Station_Cluster, year, month, day, doy, decimalLatitude, decimalLongitude, Depth |
| Physical | Temperature, Salinity, sigma0, spicity |
| Nutrients | Nitrate, Phosphate, Silicate, Oxconc, TIN |
| Tracers | Nstar, NO, PO, NO_PO, POs_star, ANP, fPW |
| Classification | classification, distance |
