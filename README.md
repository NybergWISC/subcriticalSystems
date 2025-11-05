# OpenMC Files for Subcritical Systems Inc. Meeting

Initial Models for Subcritical Systems Homogenized ADS Pb-Cooled System

__Files__: 
- generalRuns.py: Base model with most of the model implemented including both eigenvalue and fixed source settings.
- searchForKeff.py: Example of a `search_for_keff` model with a target of 0.95 with the base implementation of the model.
- *.xml: Model files based on the OpenMC xml structure
- tallies.out: Example of an OpenMC tallies output file
- *.h5 OpenMC output files that I should probably clean/clear if this is used as anything besides a storage repo (DELETE)
- trackNeutronsNoNorm/*: Outputs from a generalRuns.py fixed source run without any normalization
- sweepOfEigenvalue/*: Inputs and outputs for a sweep over an eigenvalue calculation. Run through defining the `paramArray` in the bash script loop_coreRadius.sh and running.
- printImportant_11_4/*: Directory with a copy of generalRuns.py which includes extra print statements for questions from the 11/4 interview and does not run a full calculation (commented)(DELETE)
- meshTally/*: Some scripts to produce and postprocess a cylindrical mesh tally for visualization.
- kappaFiss/*: Extra outputs checking the difference between heating and kappFiss (DELETE)
