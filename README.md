![alt text](IEPG_logo.jpg?raw=true) $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$  ![alt text](cresym.png?raw=true)    

[![DOI](https://zenodo.org/badge/951306773.svg)](https://doi.org/10.5281/zenodo.18544531)

# DQsym modelling library

- DQsym: A Dynamic Phasor-Based library for Analysis of Modern Power Systems
This repository presents the Dynamic Phasor library for Analysis of Modern Power Systems (DQsym), a Matlab-based library for dynamic phasor modeling and SSA of hybrid AC/DC power systems.

- DQsym was developed within the Intelligent Electrical Power Grids research group at TU Delft.

![alt text](DQsym_logo_v2.png?raw=true)

# Installation:

Follow these steps to install the **DQsym** library in MATLAB/Simulink:

1. Download or clone the repository
2. Open MATLAB
3. Navigate to the toolbox folder (DQsym Toolbox)
4. Locate the DQsym toolbox directory
5. Open it in MATLAB (using the Current Folder panel or cd command)
6. Run the setup script addlib.m
7. When prompted, confirm adding the folder to the MATLAB path
8. The Simulink Library Browser will refresh automatically

It is recommended to have MATLAB R2024a or a later version installed. Simulink, SimPowerSystems, the Signal Processing Toolbox, and the DSP System Toolbox 2024a or a later version are required for the toolbox's features and validations.


# Usage:
- The library supports the representation of higher-order harmonics using dynamic phasors.To capture these harmonics, the corresponding system states must be included in the model formulation.

- Example models are provided and implemented using both **MATLAB/Simulink Specialized Power Systems (SPS)** and the **DQsym** library. Users can enable or disable the SPS-based models as needed for comparison or validation.

- Scope and measurement blocks are included to support detailed observation of output signals and internal states during simulation.

# Features: 

Version 1.0 includes:

- A set of masked blocks, including:
  - State-space block  
  - Multiplication block  
  - Integration block  
  - Mux/Demux blocks  
  - DQsym-to-ABC transformation block  

- Full access to internal block implementations through masked subsystems

- Support for **Modular Multilevel Converters (MMC)** using an averaged model in a state-space formulation.

- A unified **multi-harmonic state-space framework** for modelling and analysis.


# Creating Custom Test Cases

- Custom test cases can be developed by defining the system in a discrete-time state-space form. The DQsym framework operates directly on system matrices, allowing users to represent networks, converters, and control systems within a unified multi-harmonic formulation.

- For standard network components such as AC sources, transformers, and RLC elements, models can be implemented directly by deriving the corresponding state-space matrices and connecting them using the provided blocks.

- The framework is also compatible with averaged models of power converters. In Version 1.0, the Modular Multilevel Converter (MMC) averaged model has been implemented and validated. Other converter topologies can be incorporated using the same approach; however, their implementation and validation are left to the user.

- This design provides a high degree of modelling flexibility, but it assumes familiarity with state-space modelling and system formulation.

# limitations:
- The framework assumes that users can obtain a **discrete-time state-space representation** of the system to be built.

- More detailed converter models require additional derivation and implementation.

- Generator models with mechanical dynamics are not supported in Version 1.0. These are currently represented as ideal sources.


# Contributing:

- This project is under active development. Future updates will include a complete expansion into a full Multi-Terminal DC (MTDC) network


# License:

- DQsym is licensed under the BSD 3-Clause License.

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
