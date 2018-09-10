## Electron_nuclear_sim.py - A module to calculate 13C nuclear and NV spin dynamics.
### By Peter Humphreys, 2016-2018 (inspiration and some early code from a previous approach by Tim Taminiau)

This module runs simulations of the nitrogen-vacancy centre electronic state and its interactions with carbon 13 nuclear spins, along with, if desired, the nitrogen nuclear spin of the NV. 

It is built to simulate a range of imperfections, including finite microwave pulses, amplitude errors, and detuning of the MW field from the NV frequency.

A set of examples of the use of the code are given in [NV_C13_simulation_examples.ipynb](NV_C13_simulation_examples.ipynb)

Here is the basic code structure:
1) Define an NV_system, this is the base class that holds all the physics

`nvs = noisy_NV_system(mw_duration=180e-9,carbon_params = [],inc_nitrogen=False,pulse_shape='Hermite')`

2) Make an experiment object, this holds the state of the sytem, the desired gate sequence, and allows for quick and easy measurements

`nv_expm = NV_experiment(nvs)`

3) Make a gate_sequence, here with one simple gate

`desr_seq = nv_expm.gate_sequence()
desr_seq.re(theta = 0,phi=1.0)`

4) Now one can, for example, sweep the NV mw detuning and measure the e state 

`for i,freq in enumerate(freq_range):\
    noisy_NV_system.set_NV_detuning(freq)\
    noisy_NV_system.recalculate()\
    nv_expm.apply_gates(desr_seq)\
    results[i] = nv_expm.measure_e()\
    nv_expm.reset_output_state()`\

