![alt text](IEPG_logo.jpg?raw=true) $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$  ![alt text](cresym.png?raw=true)    

# DQsym modelling library

DQsym: A Dynamic Phasor-Based library for Analysis of Modern Power Systems
This repository presents the Dynamic Phasor library for Analysis of Modern Power Systems (DQsym), a Matlab-based library for dynamic phasor modeling and SSA of hybrid AC/DC power systems.

DQsym was developed within the Intelligent Electrical Power Grids research group at TU Delft.

# Installation:

Interested users should download a local version of the repository. Aside from adding the entire DQsym library to the MATLAB path, no other installation action is required.

It is recommended to have MATLAB R2024a or a later version installed. Simulink, SimPowerSystems, the Signal Processing Toolbox, and the DSP System Toolbox are required for the toolbox's features.


# Usage:

To run the a test case:

1- Add the #lib folder to Matlab path.

2- Navigate to test case folder and execute , e.g. testcase.m in MATLAB. This will load .mat and other necessary parameters.

3- Selected test case .slx file will open and run.


# Features: 

The simulations models a Modular Multilevel Converter (MMC), in a single-station and point-to-point HVDC transmission setup, using both EMT, Matlab/Simulink Specialized power systems, and DQ-sym library.

Control blocks govern both models, with a sequencer managing the timing of control events.

Scopes and measurement blocks are provided for detailed tracing of outputs and internal states.

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
