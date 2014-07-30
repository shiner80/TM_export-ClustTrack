TM_export-ClustTrack
====================

This is the Matlab code that has been used in the pubblication (Mossuto M. et al., 2014) for importing TrackMate xlm files and calcuating mean square displacements and particle sizes from the movies. 

The routines require the Image processing toolbox and the Optimization toolbox from matlab and require Matlab 2011a or higher. The folder containing the Matlab files of this repository need to be added to the Matlab path. Also the Matlab scripts from Fiji need to be included in the Matlab path. 

To run the code open the file MSD_scriptForTrackMate.m in matlab and run by typing 

run MSD_scriptForTrackMate

in the command window
The script will open a dialog to load the trackmate export xlm file, and a second dialog to import the .tif file containing the corresponding timeseries. After running the script 3 outputs will be placed in the workspace:

AA_tlist: vector with the timepoints for which the MSD has been calculated

AA_OUT_COEF: coefficients of the analysis: 

AA_OUT_COEF (1,:) = Diffusion constants for each track;

AA_OUT_COEF (2,:) = Anomalous coefficient for each track;

AA_OUT_COEF(3,:) = sigma of gaussian fit to the particle in pixels in the first channel of the movie;

AA_OUT_COEF(4,:) = total intensity of the particle in the first channel of the movie;

AA_OUT_COEF(5,:) = sigma of gaussian fit to the particle in pixels in the second channel of the movie;

AA_OUT_COEF(6,:) = total intensity of the particle in the second channel of the movie;

AA_OUT_MSD: Mean squared displacements curves for each of the analyzed tracks.

=======================================================================
2014. Davide Mazza. San Raffaele Scientific Institute. Milan (Italy).

