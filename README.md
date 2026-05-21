# Data and code for *Henosepilachna diekei* reproductive isolation analyses

This repository contains the data and R code associated with the manuscript:

**The evolution of reproductive isolation beyond a strong first barrier in speciation between micro-allopatric host races of a phytophagous ladybird beetle, *Henosepilachna diekei***

The repository includes datasets, analysis scripts, visualisation scripts, separated RI estimates + bootstrap analyses script, and figure icons used for the manuscript.

---

## Repository structure

```text
.
├── data/
│   ├── 00_2026_EE_data_host_preference.csv
│   ├── 01_2026_EE_data_nonchoice_mating.csv
│   ├── 01_2026_EE_vis_nonchoice_mating.csv
│   ├── 02_2026_EE_data_egg.csv
│   ├── 02_2026_EE_vis_egg.csv
│   ├── 03_2026_EE_data_larval_performance.csv
│   ├── 03_2026_EE_vis_larval_performance.csv
│   ├── 04_2026_EE_data_fidelity.csv
│   ├── 04_2026_EE_vis_fidelity.csv
│   ├── 05_2026_EE_vis_dispersal.csv
│   ├── 05a_2026_EE_data_dispersal.csv
│   ├── 05b_2026_EE_data_migration.csv
│   ├── 06_2026_EE_data_choice_mating.csv
│   ├── 07a_2026_EE_vis_ri_calculation1.csv
│   └── 07b_2026_EE_vis_ri_calculation2.csv
├── icon/
│   ├── dic_d.png
│   ├── dic_m.png
│   ├── dic_dd.png
│   ├── dic_mm.png
│   ├── dic_dm.png
│   └── dic_md.png
├── 2026_EE_script_analysis.R
├── 2026_EE_script_visualisation.R
├── 2026_EE_script_RI_estimate.R
├── LICENSE
└── README.md
```

---

## Data files

The `data/` folder contains raw data and processed data used for analyses and figure generation.

| File | Description |
|---|---|
| `00_2026_EE_data_host_preference.csv` | Raw data for host preference and habitat isolation analysis |
| `01_2026_EE_data_nonchoice_mating.csv` | Raw data for non-choice mating trials |
| `01_2026_EE_vis_nonchoice_mating.csv` | Processed data for visualising non-choice mating results |
| `02_2026_EE_data_egg.csv` | Raw data for egg hatchability and postmating prehatching isolation analysis |
| `02_2026_EE_vis_egg.csv` | Processed data for visualising egg hatchability |
| `03_2026_EE_data_larval_performance.csv` | Raw data for larval performance and hybrid inviability analysis |
| `03_2026_EE_vis_larval_performance.csv` | Processed data for visualising larval performance |
| `04_2026_EE_data_fidelity.csv` | Raw data for host fidelity analysis |
| `04_2026_EE_vis_fidelity.csv` | Processed data for visualising host fidelity |
| `05_2026_EE_vis_dispersal.csv` | Processed data for visualising dispersal distance |
| `05a_2026_EE_data_dispersal.csv` | Raw data for dispersal analysis |
| `05b_2026_EE_data_migration.csv` | Raw data for migration pattern analysis |
| `06_2026_EE_data_choice_mating.csv` | Raw data for choice mating trials and realised sexual isolation |
| `07a_2026_EE_vis_ri_calculation1.csv` | Processed data for reproductive isolation contribution modelling |
| `07b_2026_EE_vis_ri_calculation2.csv` | Processed data for visualising reproductive isolation asymmetry |

---

## Icon files

The `icon/` folder contains image files used in the figures.

| File | Description |
|---|---|
| `dic_d.png` | Icon for *H. diekei* D-race |
| `dic_m.png` | Icon for *H. diekei* M-race |
| `dic_dd.png` | Icon for DD-pair |
| `dic_mm.png` | Icon for MM-pair |
| `dic_dm.png` | Icon for DM-pair |
| `dic_md.png` | Icon for MD-pair |

---

## R scripts

| Script | Purpose |
|---|---|
| `2026_EE_script_analysis.R` | Main statistical analyses |
| `2026_EE_script_visualisation.R` | Figure generation and visualisation |
| `2026_EE_script_RI_estimate.R` | Reproductive isolation strength and contribution estimates |

---

## How to reproduce the analyses

Clone or download this repository, then open the repository folder in R or RStudio.

```bash
git clone https://github.com/diekei/2026_EE_diekei_dicliptera_speciation.git
cd 2026_EE_diekei_dicliptera_speciation
```

Run the scripts from the repository root directory.

Recommended order:

```r
source("2026_EE_script_analysis.R")
source("2026_EE_script_visualisation.R")
source("2026_EE_script_RI_estimate.R")
```

The scripts assume that the `data/` and `icon/` folders remain in their original locations.

---

## Citation

If you use the data or code from this repository, please cite the associated manuscript:

Maulana, A. et al. The evolution of reproductive isolation beyond a strong first barrier in speciation between micro-allopatric host races of a phytophagous ladybird beetle, *Henosepilachna diekei*. Manuscript submitted to *Ecology and Evolution*.

A full citation will be added after publication.
