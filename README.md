# DINOS Model Code

## Overview of DINOs

This codebase includes the code used to generate results for DINOs, a configurable theoretical framework that can assemble arbitrary messages of information that are to be written into DNA. The model in this code assumes a overhang based assembly model, an assembly method commonly used in Golden Gate Assembly. The DINOs model implemented in this code models 3 different methods for removing redundant reactions for a given data set: an ideal method whose only constraint is that reactions can be removed if they assemble the same data, a realistic optimization method that is constrained by both data and the bookending overhangs of an assembled strand, and finally a optimization method that tries to improve the opportunity of redundant reactions by using the rotation properties of DINOs assembly.

## Requirements 

This module is relatively self contained, but at least needs Matplotlib to plot results after analyzing a data set. Requirements should be installable through the use of `make init` in the top level directory of this project. As of now, the project is still only compatible with Python 2.7.

## Installing 

After requirements are installed simply use `make install` to install the module so that it can be loaded.

## Running DINOs Experiments over a Data Set

With the requirements installed and the module installed, the analysis performed in the DINOs paper can be run on a data set by first `cd PATH_TO_DINOS_PROJECT/tools`, then launching the analysis on a data set with `python tree_analysis.py --out_dir OUT_DIR_PATH --w_dir DATA_SET_PATH --1_bit`. In this command `OUT_DIR_PATH` is top level output directory that you would like the results to be dumped out to, and `DATA_SET_PATH` is the top level path to the data set that is going to be analyzed. When data is dumped out, both raw data and figures used in the DINOs paper will be generated. It should be noted that this analysis uses a lot of memory, for the data set analyzed in DINOS we used 64 GB of memory.
