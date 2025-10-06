# glp1r-reward-circuit
Analysis code and example datasets from the Güler lab for behavioral neuroscience research.

## Overview
This repository contains analysis code and toy datasets used to reproduce key steps of our pipeline for studying how GLP-1R agonists affect home cage behavior (Figure 2).  
The core pipeline is written in Python and Jupyter notebooks, with toy datasets provided for reproducibility.  

Additional statistical analyses (e.g., PCA, MANOVA, permutation testing) are implemented in R for downstream interpretation of behavioral data. 
Separate MATLAB scripts are included for fiber photometry data processing and calcium signal analysis from a different experiment within the same study.

## Installation
This project requires Python 3.11 and the packages listed in `environment.yml`.  
To set up a conda environment with all dependencies:

```bash
# Create the environment
conda env create -f environment.yml

# Activate it
conda activate glp1r-reward-circuit
```

All required packages, including numpy, pandas, matplotlib, networkx, h5py, and opencv-python, will be installed automatically.
After activating the environment, you can run the notebooks as described below.

## Usage
Run the notebooks in order with the provided toy datasets:

1. **locations.ipynb**  
   - Identifies whether a mouse is in a region of interest (ROI) using keypoint coordinates from SLEAP.  
   - Loads the following toy datasets:  
     ```
     saline_orfo_arena2_2785_re_15s.h5    # SLEAP keypoint coordinates
     saline_orfo_arena2_2785_re_15s.mp4   # Cropped video corresponding to the .h5 file
     ```

2. **location_behavior.ipynb**  
   - Groups/names Keypoint-MoSeq syllables using ROI information.  
   - Loads the following toy datasets:  
     ```
     dan_arena1_1001_re.csv            # Keypoint-MoSeq syllables
     dan_arena2_1002_re.csv
     veh_arena1_1001_re.csv
     veh_arena2_1002_re.csv
     location_dan_arena1_1001_new.csv  # ROI assignments
     location_dan_arena2_1002_new.csv
     location_veh_arena1_1001_new.csv
     location_veh_arena2_1002_new.csv
     ```
     *Comments indicate the type of data in each file to help users understand how they are used in the notebook.*

3. **build_dataset.ipynb**  
   - Performs transition and bout analyses on the processed data.

4. **visualization.ipynb**  
   - Generates select plots and figures from the processed dataset.
  
## Additional Analyses (R, MATLAB)

The file `PCA.R` contains analyses and visualizations performed on the processed behavioral dataset.  
This script complements the Python pipeline.

Key steps include:
- Principal Component Analysis (PCA) on standardized behavioral features
- Visualizations: Scree plot, PCA scatter plots, loadings heatmap

To run this script, open `PCA.R` in R (tested with R ≥ 4.1.2) and ensure required packages are installed:
```r
install.packages(c("tidyverse", "plotly", "scales"))
```

Fiber Photometry Signal Processing
MATLAB scripts for analyzing fiber photometry data from experiments in this study. The provided scripts perform analysis of calcium and 
dopamine signals.

Available Scripts
- `GCaMP Calcium Events.m` – Analysis of calcium dynamics (related to Figure 4q–t).
- `Dopamine dLight Z Score.m` – Z-scored analysis of dopamine signals (related to Figure 5i–v).

To run these scripts open the desired script in MATLAB (tested with version R2024b) and ensure that any required MATLAB 
toolboxes/packages are installed.
