Dynamic Phasor-Based Harmonic Modelling for Modular Multilevel Converters

1-Installation 
  - Matlab/Simulink 2024a version is needed to run the simulation. 
2- Usage 
  - Run the [dqnMMCsingle.m] file in Matlab to upload Spq.mat and other needed parameters.
3- Features  
  - The simulation showcases the MMC model in both EMT and DQsym approaches; both models are controlled from each control block,
    with a sequencer to time the control events. Power will be shared between the station and the grid. The scopes and measurements blocks 
    are available to trace the simulation and the outputs.
4- Contributing
  This model will be expanded to a complete point-to-point model with two MMCS and later into an MTDC network.
  
