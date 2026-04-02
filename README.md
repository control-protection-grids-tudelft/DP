![alt text](IEPG_logo.jpg?raw=true) $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$  ![alt text](cresym.png?raw=true)    

[![DOI](https://zenodo.org/badge/951306773.svg)](https://doi.org/10.5281/zenodo.18544531)

# DQsym modelling library

DQsym: A Dynamic Phasor-Based library for Analysis of Modern Power Systems
This repository presents the Dynamic Phasor library for Analysis of Modern Power Systems (DQsym), a Matlab-based library for dynamic phasor modeling and SSA of hybrid AC/DC power systems.

DQsym was developed within the Intelligent Electrical Power Grids research group at TU Delft.

![alt text](DQsym_logo_v2.png?raw=true)

# Installation:

Follow these steps to install the **DQsym** library in MATLAB/Simulink:

1. Download or clone the repository
2. Open MATLAB
3. Navigate to the toolbox folder
4. Locate the DQsym toolbox directory
5. Open it in MATLAB (using the Current Folder panel or cd command)
6. Run the setup script addlib.m
7. When prompted, confirm adding the folder to the MATLAB path
8. The Simulink Library Browser will refresh automatically

It is recommended to have MATLAB R2024a or a later version installed. Simulink, SimPowerSystems, the Signal Processing Toolbox, and the DSP System Toolbox 2024a or a later version are required for the toolbox's features.


# Usage:

The framework assumes that users can obtain a discrete-time state-space representation of the system.



# Features: 
Version 1.0 includes:

Support for Modular Multilevel Converters (MMC) using an averaged model . 

The simulations models a Modular Multilevel Converter (MMC), in a single-station and point-to-point HVDC transmission setup, using both EMT, Matlab/Simulink Specialized power systems, and DQ-sym library.

Control blocks govern both models, with a sequencer managing the timing of control events.

Scopes and measurement blocks are provided for detailed tracing of outputs and internal states.

Current limitations:
More detailed converter models require additional derivation and implementation
Generator models with mechanical dynamics are not supported for  Version 1.0, these are represented as ideal sources.

# Contributing:

This project is under active development. Future updates will include a complete expansion into a full Multi-Terminal DC (MTDC) network


# License:

DQsym is licensed under the BSD 3-Clause License.

# Original Authors:

Saif Alsarayreh, Robert Dimitrovski, Aleksandra Lekić.

Copyright (c) 2025 TU DELFT 

# Citation:

If you use DQsym in your research, please cite the following preprint:

```
Robert Dimitrovski, Saif Alsarayreh, Aleksandra Lekić. A Novel Dynamic Phasor-based Mathematical Framework
for Hybrid AC/DC Power System Simulation. TechRxiv. September 15, 2025.
```

# Acknowledgement:

This work was funded by the CRESYM project Harmony (https://cresym.eu/harmony/).


# Contact:

For any type of inquiry, please contact: S.T.S.Alsarayreh@tudelft.nl, R.Dimitrovski@tudelft.nl, A.Lekic@tudelft.nl.  
