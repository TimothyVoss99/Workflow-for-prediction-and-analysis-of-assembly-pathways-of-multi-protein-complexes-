# Assembly Pathway Prediction Workflow

This repository provides a Jupyter notebook-based workflow for predicting the assembly pathway of protein complexes, based on the spatial proximity of their subunits.

## Features

- Predicts assembly pathways using 3D spatial relationships between subunits.
- Implemented as an interactive Jupyter Notebook.
- Uses PTGLTools to generate structural graphs from annotated protein structures.
- Supports five clustering techniques to model and visualize the assembly process.
## Installation

To clone this repository **including the required submodule**, use the following command:

git clone --recurse-submodules https://github.com/TimothyVoss99/Workflow-for-prediction-and-analysis-of-assembly-pathways-of-multi-protein-complexes-.git

If you already cloned the repository without --recurse-submodules, you can initialize and update the submodule manually:

git submodule update --init --recursive

## User Guide

### 1. Prepare the Structure File

Generate a DSSP-annotated structure file of your protein complex using the online tool:  
https://pdb-redo.eu/dssp

### 2. Place the File

Copy the DSSP-annotated structure file into the following directory: \PTGLtools\PTGLgraphComputation\dist
### 3. Run the Workflow

Open the provided Jupyter Notebook and follow the instructions step by step:

- Load the structure and generate the graph using PTGLTools.
- Choose from five clustering methods to model the assembly.
- Visualize and interpret the resulting pathway.

## Contact

For questions or suggestions, feel free to open an issue or contact the maintainer.
