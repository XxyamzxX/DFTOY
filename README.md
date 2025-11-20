# DFTOY
This repository implements a simplified 1D Density Functional Theory (DFT) solver for teaching, concept verification, and numerical tests. It performs a Kohn–Sham self-consistent cycle: builds the effective potential, assembles the Hamiltonian, solves for the ground state, updates the density, and iterates until convergence.



#Instructions


If you want to run this program in the console you have to write this:

python3 -m gui_app.gui


If you wanna run pytest you have to write this:

pytest tests/
