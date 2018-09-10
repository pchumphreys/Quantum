## Summary 
This directory contains simulations made to understand the dynamics of a single NV under optical excitation. 

### Contents
* [Pulsed optical excitation simulation.ipynb](Pulsed optical excitation simulation.ipynb) - This notebook contains a simple quantum jump type code used to calculate tthe chance that an NV will be excited twice during one excitation pulse.
* [NVpop_monteCarlo](NVpop_monteCarlo) - This module is a fairly general code for simulating probabilistic transitions between a manifold of states. It has been mildly optimised using cython.
* [Repumping Monte Carlo](Repumping Monte Carlo) - This directory contains a classical simulation built on NVpop_monteCarlo, and used to simulate the dynamics of the NV centre during repumping.
