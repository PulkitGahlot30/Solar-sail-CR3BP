# Solar Sail CR3BP — PINN-Based Stability Analysis of Triangular Equilibrium Points

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-orange)](https://www.mathworks.com/)
[![Status](https://img.shields.io/badge/Status-Active%20Research-green)]()
[![License](https://img.shields.io/badge/License-MIT-blue)]()

## Overview

This repository contains the MATLAB simulation and Physics-Informed Neural Network (PINN) pipeline for analyzing the **stability of triangular equilibrium points (L4/L5)** in the **Solar Sail Circular Restricted Three-Body Problem (CR3BP)**.

The work extends the classical CR3BP framework by incorporating:
- **Solar sail radiation pressure** (McInnes flat-sail SRP model)
- **Sail lightness number (β)** as a key control parameter
- **Pitch and clock angles** for sail orientation
- **Artificial Equilibrium Points (AEPs)** near L4/L5
- **PINN-based stability classification** using a three-output `dlnetwork` architecture

---

## Repository Structure

```
solar-sail-CR3BP/
│
├── README.md
│
├── MATLAB/
│   ├── main_simulation.m          # Main entry point
│   ├── equations_of_motion.m      # Solar sail CR3BP EOM
│   ├── equilibrium_points.m       # L4/L5 AEP computation
│   ├── stability_analysis.m       # Eigenvalue / Lyapunov analysis
│   └── plots/                     # Output figures
│
├── PINN/
│   ├── pinn_training.m            # PINN training loop
│   ├── network_architecture.m     # dlnetwork (3-output) definition
│   └── loss_functions.m           # Physics + data loss functions
│
└── results/
    └── figures/                   # Saved plots and results
```

---

## Physical Model

### Solar Sail CR3BP Equations of Motion

The non-dimensional equations incorporating solar sail SRP acceleration:

$$\ddot{x} - 2\dot{y} = \frac{\partial \Omega}{\partial x} + a_x^{sail}$$

$$\ddot{y} + 2\dot{x} = \frac{\partial \Omega}{\partial y} + a_y^{sail}$$

$$\ddot{z} = \frac{\partial \Omega}{\partial z} + a_z^{sail}$$

where $\Omega$ is the effective pseudo-potential and $\mathbf{a}^{sail} = \beta \frac{\mu_1}{r_1^2}(\hat{n} \cdot \hat{r}_1)^2 \hat{n}$ is the sail acceleration.

### Key Parameters

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| Sail lightness number | β | Ratio of SRP to gravitational force |
| Mass ratio | μ | Primary mass ratio |
| Critical mass ratio | μ_c ≈ 0.0385 | Stability boundary |
| Pitch angle | α | Sail orientation w.r.t. Sun-line |
| Clock angle | δ | Out-of-plane sail orientation |

---

## PINN Architecture

- **Input:** State vector `[x, y, ẋ, ẏ, β, α]`
- **Hidden layers:** Fully connected, tanh activation
- **Output:** 3 nodes → `[x(t), y(t), stability_index]`
- **Loss function:** Physics residual loss + Initial condition loss + Data loss
- **Framework:** MATLAB Deep Learning Toolbox (`dlnetwork`, `dlarray`)

---

## Results

> *(Figures and results will be updated as research progresses)*

- Stability maps in (β, μ) parameter space
- PINN-predicted trajectories vs numerical integration
- Eigenvalue analysis at AEPs near L4/L5

---

## Target Publication

- **Advances in Space Research** (Elsevier, Q2)
- **Chaos, Solitons & Fractals** (Elsevier, Q1)

---

## Dependencies

- MATLAB R2021b or later
- Deep Learning Toolbox
- Optimization Toolbox

---

## Author

**Pulkit**  
Researcher, Astrodynamics & Celestial Mechanics  
*Assistant Professor Applicant — India*

---

## Citation

> *(To be updated upon publication)*

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
