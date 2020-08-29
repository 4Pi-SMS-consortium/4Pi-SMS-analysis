# 4Pi-SMS-analysis
This package includes the Matlab code for analyzing the dataset recorded on the 4Pi-SMS microscope.

# Requirements
  - Microsoft Windows 7 or newer, 64-bit
  - CUDA capable graphics card, minimum Compute Capability 3.0
  - CUDA 8 compatible graphics driver (for GeForce products 378.66 or later)
  - Matlab R2017b or newer  
    - Curve Fitting Toolbox
    - Optimization Toolbox
  - DIPimage Toolbox (http://www.diplib.org/download)
  - Some files from SMAP (included in this package) (https://github.com/jries/SMAP) 
    
# Installation
  - Download and install Matlab of the right version 
  - Download and install DIPimage
  - Download the code to your computer
  - Download example dataset
  - Start Matlab and run the code

# Data required
An example dataset of immunolabeled microtubules of a COS-7 cell is provided. 

Download the dataset from: https://drive.google.com/open?id=1I3vBu9UyKaMJ6h1eGCE8DSasIcNGrG9X

Additonal information of the dataset
  - Cell images: 20 cycles included (3000 frames per cycle)
    - Labeling: anti-a-tubulin primary antibody and CF660C conjugated secondary antibody
    - Imaging conditions: 100 fps with a 642 nm laser at 7.5 kW/cm2 
  - Images of fluorescence beads for channel alignment (10 frames)
  - Images of fluorescence beads with 20 nm step sizes for phase shift estimation (61 steps)

# How to run
  - Start Matlab program
    - In Matlab, run "dipstart" to initialize DIPimage toolbox
    - In Matlab, run "SMS-4Pi" to open the GUI interface 
    - Click "Main Folder" and select the folder of the example dataset
  - Generate calibration file 1
    - On the GUI Menu, click "Channel Alignment" -> "Load", select the file in the folder "Beads_align"
    - A calibration file named "align_642_FMTtransform_datestring.mat" will be generated
  - Generate calibration file 2
    - On the GUI Menu, click "Find Phaseshift" -> "Load", select all the file in the folder "Beads_stack"
    - Two calibration files named "bead_642_Astfit_datestring.mat" and "bead_642_dphi_cali_datestring.mat" will be generated
  - Step 1: get positions of the single molecules
    - On the GUI file dispaly window, select the example dataset folder "\4Pi-SMS-Example-Dataset\Cell04", click "Get Positions"
    - A result file named "Cell04_642_tmpresult_datestring.mat" will be genrated
  - Step 2: phase unwrapping, drift correction and stitching
    - Click radio button "Reconstruct" to switch the GUI display
    - On the GUI file display window, select "Cell04_642_tmpresult_datestring.mat", click "Reconstruction"
  - Step 3: visualization
    - A result file named "Cell04_642v20_60.mat" will be generated, which includes 3D positions of localized molecules and can be loaded in PYME for visualization
    - A folder named "Cell04_ll" will be genreated with a file named "particles.csv", which can be loaded in Vutara SRX Software (Bruker) for 3D visualization and rendering

# Contact
For any questions / comments about this software, please contact [Bewersdorf Lab](http://www.bewersdorflab.org/).

# Copyright and Software License
Copyright (c) 2020 Bewersdorf Lab, Yale Univeristy School of Medcine, USA.

The package is licenced under the [GNU GPL](https://www.gnu.org/licenses/). 
