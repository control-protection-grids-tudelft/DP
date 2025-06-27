Dynamic Phasor-Based Harmonic Modelling for Modular Multilevel Converters
1. Installation
MATLAB/Simulink 2024a is required to run the simulations.

2. Usage
To run the single MMC simulation :

Execute dqnMMCsingle.m in MATLAB. This will load Spq.mat and other necessary parameters.

Open and run Dqn_C_PQ_HVsMc.slx in Simulink.

To run the point-to-point HVDC MMC simulation:

Execute P2PHVDCMMC.m in MATLAB. This will load SFu.mat and other necessary parameters.

Open and run P2PHVDCMMC.slx in Simulink.

3. Features
The simulation models a Modular Multilevel Converter (MMC) using both EMT and DQ-sym approaches.

Control blocks govern both models, with a sequencer managing the timing of control events.

The simulation demonstrates power sharing between the MMC station and the grid.

Scopes and measurement blocks are provided for detailed tracing of outputs and internal states.

4. Contributing
This project is under active development. Future updates will include a complete expansion into a full Multi-Terminal DC (MTDC) network.
