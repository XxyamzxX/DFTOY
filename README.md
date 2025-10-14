# 🧩 DFToy: A Pedagogical Self-Consistent Quantum Simulator

**DFToy** is an open-source educational software designed to illustrate the **self-consistent quantum problem** at the heart of *Density Functional Theory (DFT)* and related mean-field methods.  
It allows users to explore how a quantum system reaches equilibrium by iteratively updating the **density**, **effective potential**, and **wavefunction** until self-consistency is achieved.

> DFToy is not intended as a high-precision computational package, but as a **didactic tool** to help students, educators, and researchers visualize the logic and physics behind quantum self-consistency.

---

## 🧠 Concept

The program models a one-dimensional quantum system of interacting particles confined by an external potential \( V_{\text{ext}}(x) \).  
It solves the **effective Schrödinger equation**:

\[
\left[-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V_{\text{eff}}(x; \rho)\right] \psi(x) = \mu \, \psi(x)
\]

where:
- $( \rho(x) = 2|\psi(x)|^2)$ is the particle density (with spin degeneracy 2),
- $( V_{\text{eff}}(x; \rho) )$ includes both external and density-dependent (self-consistent) terms,
- $( \mu )$ is the chemical potential of the ground state.

Through successive iterations, DFToy converges to a stationary density where \( V_{\text{eff}} \) and \( \rho(x) \) are mutually consistent — just like in full DFT calculations.

---

## 🧰 Features

- **Interactive Graphical Interface (GUI):**  
  Intuitive controls for running and visualizing self-consistent simulations.

- **Multiple Potential Types:**  
  Choose among harmonic, double-well, and linear potentials.

- **Dynamic Visualization:**  
  Real-time plots of density \( \rho(x) \), effective potential \( V_{\text{eff}}(x) \), and convergence progress.

- **Numerical Stability:**  
  Damped density mixing and adaptive convergence control.

- **Educational Focus:**  
  Clearly structured equations and comments throughout the code to promote understanding of self-consistent quantum models.

---
