# DQsym vs HSS vs DP: open-loop MMC

MATLAB/Simulink code that reproduces the simulation results of the paper

> S. Alsarayreh, R. Dimitrovski, A. Lekić, *Characterisation of the DQsym Framework Against Harmonic State-Space and Dynamic Phasor Models for Modular Multilevel Converters.*

The code builds four models of the same open-loop MMC and compares them:

- **abc**: the nonlinear arm-averaged time-domain model (abc ODE), used as the reference
- **HSS**: the harmonic state-space model
- **DP**: the real-valued dynamic phasor model
- **DQsym**: the symmetrical dq dynamic phasor model, simulated in Simulink

All harmonic-domain models use the truncation order K = 4.

## Files

| File | Purpose |
|---|---|
| `OL_TD_DP_HSSv2_cos.m` | Builds the abc, HSS, DP and DQsym models. Produces Figs. 2–5, Fig. 9 and Table 3. Does not need Simulink. |
| `mmc_dp_hss_vs_simulink_cos.m` | Simulates the DQsym Simulink model and compares it with the abc, DP and HSS models. Produces Figs. 6–8 and Table 4. |
| `OL_MMC_dqn_cos_V2.slx` | DQsym Simulink model of the open-loop MMC. |
| `SOL.mat` | Solver matrices (`SOL`) of the `DQsym_SS` block. |

## Requirements

- MATLAB R2024a (tested). The second script also runs the first one.
- Symbolic Math Toolbox
- Optimization Toolbox (`fsolve`)
- Simulink and Simscape Electrical (Specialized Power Systems), for the second script only
- The DQsym Simulink library (`DQsym_Lib`) on the MATLAB path, for the second script only

## Usage

Run the scripts from the folder that contains them.

Model-based results only (Figs. 2–5, 9 and Table 3):

```matlab
OL_TD_DP_HSSv2_cos
```

All results (Figs. 2–9 and Tables 3–4):

```matlab
mmc_dp_hss_vs_simulink_cos
```

The full run takes about 25 minutes on a standard PC. Most of that time is the symbolic DP derivation and the two Simulink simulations.

## Outputs

The figures are saved as vector PDFs in `./figs`, drawn at their printed size:

| File | Paper figure |
|---|---|
| `fig2_sparsity_abc.pdf` | Fig. 2: sparsity of the linearised abc model |
| `fig3_sparsity_hss.pdf` | Fig. 3: sparsity of the HSS model |
| `fig4_sparsity_dp.pdf` | Fig. 4: sparsity of the DP model |
| `fig5_sparsity_dqsym.pdf` | Fig. 5: sparsity of the DQsym solver (pnz coordinates) |
| `fig6_steady_state.pdf` | Fig. 6: periodic steady state over one fundamental period |
| `fig7_step_response.pdf` | Fig. 7: response over three periods and error w.r.t. the abc ODE |
| `fig8_harmonics.pdf` | Fig. 8: harmonic magnitudes of the phase-a upper arm current and capacitor voltage |
| `fig9_eigenvalues.pdf` | Fig. 9: eigenvalues of DP, HSS and DQsym |

Tables 3 and 4 are printed in the Command Window as `table3` and `table4`.

## Benchmark parameters (Table 2)

| Parameter | Value |
|---|---|
| Pole-to-pole dc voltage V<sub>DC</sub> | 640 kV |
| Arm resistance, inductance R<sub>arm</sub>, L<sub>arm</sub> | 1.07 Ω, 50 mH |
| Submodule capacitance, submodules per arm C<sub>SM</sub>, N<sub>sm</sub> | 10 mF, 400 |
| Star-point resistance, inductance R<sub>star</sub>, L<sub>star</sub> | 10 Ω, 5 H |
| AC phase voltage, peak V<sub>g</sub> | 271.9 kV |
| Modulation amplitude, angle | 0.85, −2° |
| Fundamental frequency | 50 Hz |

The parameters are set in Section 1 of `OL_TD_DP_HSSv2_cos.m`, and the truncation order `K` in Section 0.

## Citation

If you use this code, please cite the paper above.
