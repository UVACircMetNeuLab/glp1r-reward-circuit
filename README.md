# glp1r-reward-circuit
Analysis code and example datasets from the Güler lab for behavioral neuroscience research.

## Overview
This repository contains analysis code and toy datasets used to reproduce key steps of our pipeline for studying how GLP-1R agonists affect home cage behavior.  
The code is written in Python and Jupyter notebooks, and toy datasets are provided so users can run the workflow without large raw data files.

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

