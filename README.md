# ExactResonance
This package of codes provides the kinematic model simulating the energy transfer by exact resonances for surface gravity waves in a finite periodic domain.

Reference: Zhou Zhang, Yulin Pan, Do surface gravity waves have a frozen turbulence state? (2022)

## Dependencies
MATLAB (versions newer than 2018b) is needed to run the codes.

## Usage
The codes consist of 3 parts:
* __preComputaion__: Computation of the scale-resonant quartets in the finite wavenumber domain (`exactResonance.m`).
* __cascade__: Initialization of active region (`initialModes.m`); Generation of new active modes from exact resonances (`cascade.m`); Summary of results with different initial conditions (`kmaxOut.m`,`kminIn.m`).
* __quartetStructure__: Calculation of the quantities manifesting the structure of resonant quartets (`resQuart.m`).

## Contact
Zhou Zhang \
Email: joezhang@umich.edu
