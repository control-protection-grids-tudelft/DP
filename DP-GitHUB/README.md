Dynamic Phasor-Based Harmonic Modelling for Modular Multilevel Converters

- Installation
  
MATLAB/Simulink 2024a is required to run the simulations.

- Usage
  
To run the single MMC simulation :

Execute dqnMMCsingle.m in MATLAB. This will load Spq.mat and other necessary parameters.

Open and run Dqn_C_PQ_HVsMc.slx in Simulink. The simulation demonstrates power sharing between the MMC station and the grid.

To run the point-to-point HVDC MMC simulation:

Execute P2PHVDCMMC.m in MATLAB. This will load SFu.mat and other necessary parameters.

Open and run P2PHVDCMMC.slx in Simulink. The simulation demonstrates power sharing between the  two MMC stations.

- Features
  
The simulation models a Modular Multilevel Converter (MMC) using both EMT and DQ-sym approaches.

Control blocks govern both models, with a sequencer managing the timing of control events.

Scopes and measurement blocks are provided for detailed tracing of outputs and internal states.

- Contributing
  
This project is under active development. Future updates will include a complete expansion into a full Multi-Terminal DC (MTDC) network.
