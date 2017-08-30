# Quick Start Manual:

OBJECT: This software calculates: (a) Effective solid angle under specimen tilt and shift condition, (b) Counts ratio, (c) Absolute counts of two elements. (d) Composition analysis based on counts/ratio method. It is based on numeric approach taking effects from multiple-detector geometry, Be holder absorption, holder frame shadowing, and specimen absorption in flat or complex geometry into consideration.

NEW UPDATE: INCORPORATE COMPLEX SPECIMEN GEOMETRY (2017)
Effective solid angle, counts and counts ratio prediction under complex specimen geometry are now available in SuperAngle. This helps make better EDS prediction in systems with strong absorption issues.

HOW TO RUN: The program is coded in MATLAB. To start, find file superangle.m in the software folder and run. Then a graphic user interface (GUI) will show up, as illustrated below. 

** To conduct calculation under complex specimen geometry, please Run SuperAngle_script_specimen_geometry_*.m. GUI version is currently under construction. 

You will see parameters of specimen information be loaded from the file. The software also pre-loaded default Super-X detector parameters and FEI low background holder setting for calculation. These parameters can be further changed in the GUI interface via direct screen input or file input.
(1)	Detector setting: detector.xlsx (or other names)
(2)	Specimen information: specimen_startup.xlsx (or other names)
(3)	Holder information: holder_FEI_LB.xlsx
For example, click Browse bottom in the specimen information panel, choose file “specimen_Ni3Al_demo.xlsx“ as specimen information input.


To calculate, click functional bottoms in the RUN panel

(1) Single spot column
Solid Angle: calculate effective solid angle of A, B elements at single point
Counts: calculate absolute counts of A, B elements at single point
Composition: For unknown composition sample, composition will be calculated based on iteration method.

(2) Tilt series column
Solid Angle: calculate effective solid angle of A, B elements along X-tilt (alpha-tilt) or Y-tilt (Beta-tilt) series. X/Y-tilt can be setup in “Experimental Comparison” panel. 
Counts: calculate absolute counts of A, B elements along X-tilt (alpha-tilt) or Y-tilt (Beta-tilt) series. To compare with experimental results, dataset must be input in the “Experimental Comparison” panel.
Composition: If composition is unknown, composition will be calculated based on iteration method. For known composition condition, comparison between experimental counts/ratio and model prediction will be analyzed.

(3) 2D map (X-Y-tilt)
Solid Angle: calculate 2D map of effective solid angle of A, B elements along X-Y-tilt. Tilt range and step size must be set in the “RUN” panel. The calculation may take long time to finish.
Counts: calculate 2D map of absolute counts.

(4) 2D map (X-Y-shift)
Solid Angle: calculate 2D map of effective solid angle of A, B elements along X-Y-shift. Shift range and size must be set in the “RUN” panel. The calculation may take long time to finish.
Counts: calculate 2D map of absolute counts.

Tips: To save time, choose Display figures from calculated data to re-draw the calculated 1D/2D figures. Click Clear All Figures to close all open figures. You can save figures into files when choose Figure output. Click Update Info to get updates from MATLAB command window.



## Something beyond: How to change settings for your sample?

###1. Detector configuration

dAngle control the accuracy of the calculation. For typical calculation, 0.2 degree is recommended for convergent data. 0.5 degree is also chosen for fast calculation with sacrificing much accuracy. 

Super-X four-quadrant detector configuration has been set as default configuration in the code. To further modify, open detector.xlsx and change settings.  Note: the code support one to multiple detectors setting. 

Detector efficiencies for windowless SDD detector in Titan are incorporated in SDD_windowless_efficiency.xlsx. 

Below is the orientation setting of the detector for your reference

=======Coordinate========
   Titan G2 ChemiSTEM
           Y+
            |
       D1   |     D2  
     135deg |   45deg
           (Z)-----X+
            
       D4         D3
     225deg     315deg
            
=========================


### 2. Sample information

Absorption coefficient is fully incorporated in the software. The value can be found in Absorption coefficient_K.xlsx, Absorption coefficient_L.xlsx and Absorption coefficient_M.xlsx. However, ionization cross-section and fluorescence yield data for elements are not incorporated. To calculate counts/counts ratio, these info must be added into Atomic_info.xlsx first. Note: the effective solid can still be calculated without knowing these values.

It is recommended to use a well-known shape sample (FIB grid) to calibrated center of specimen shift from the readings in TIA software. 

The calculation does not require wedge angle as input, but does require the knowledge of inclination angle between specimen (top surface) and the primary axis of the holder in the absorption correction calculation. For FIB sample with simple geometry, the inclination angle can be considered close to zero degree. For wedge shaped sample, the inclination angle is different from wedge angle and cannot be directly measured as it changes with grid mounting and specimen local distortion specifically at thin areas. Strictly speaking, the inclination angle is an immeasurable (or unknown) parameter for wedge shaped sample from external measurement. Alternatively, an internal way is recommended here to estimate the value of inclination angle by searching the smallest deviation between calculation and experimental results from both counts and counts ratios at different tilt angles series. The inclination angle can be considered as an internally calibrated parameter via tilt series prior to any further EDS quantification at local region of interest.


### 3. Holder information

The current version supports FEI low background holder. For other types of holder, code modification might be needed for accurate prediction. The tilt limit of the holder should be less than 35 degree. Considering the possibility that sample may mount on the grid. The grid beam-blocking scenario is also considered in the calculation as long as the grid setting is properly set.

### 4. Experiment information

Live acquisition time will vary with dwell time setting. It must be read from raw data or EDS software. It is not the total time of data acquisition in experimental.


### 5. Complex Sample Geometry
Run through script based matlab file SuperAngle_script_specimen_geometry_*.m to calculate. The setting is similar to the standard script version with additional info about the 3D specimen object input and location of the probe with this model. Demo files are provided for example. GUI version is under construction. Matlab 2016a and later version is needed.

To make a model, please run main_convert_CAD2MAT_*.m in /3D_model_build folder to load model from stl format to Matlab. Then, run main_model_build_*.m to build final dual-models for SuperAngle calculation. Detailed comments of usage can be found in these demo files. 

## Citation Reference

The authors request that any published work or images created using SuperAngle include the following reference:

W. Xu, J.H. Dycus, X. Sang, J.M. Lebeau, A Numerical Model for Multiple Detector Energy Dispersive X-ray Spectroscopy in the Transmission Electron Microscope, Ultramicroscopy, 2016, 164 (2016) 51-61.

W. Xu, J.H. Dycus, J.M. Lebeau, Numerical modeling of specimen geometry for quantitative energy dispersive X-ray spectroscopy, Ultramicroscopy, 2017, accepted. (arXiv:1708.04565)