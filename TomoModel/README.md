# TomoModel

## Overview

ReadS20RTS takes in models in the S20RTS/S40RTS format and prints out to an ascii file the spherical harmonic coefficients for the model at a discrete set of radii. 

Inputs: [Name of spherical reference model in deck format] [Name of 3D model file in S20RTS format] [name for the output file]

Within the data directory examples of the two input files can be found (specifically, PREM and S20RTS). Each line of the output file contains the following:

radius (m)  | reference density (kg / m^{3}) | reference shear velocity (m / s) | spherical harmonic coefficients for dlnvs 

Radii can be repeated to indicate discontinuities within the parameters.

For each radius the spherical harmonic coefficients are listed in the following order, denoting the function by $f$ for simplicity:


$$
\mathrm{Re} f_{00} \quad \mathrm{Re} f_{1-1} \quad \mathrm{Re} f_{10} \quad  \mathrm{Re} f_{11} \quad  \mathrm{Re} f_{2-2} \quad \dots \mathrm{Re} f_{LL} \quad \mathrm{Im} f_{00} \quad \mathrm{Im} f_{1-1} \quad \dots
$$

where $L$ is the maximum degree of the model (e.g. 20 for S20RTS). Note that the coefficients are for complex fully normalised spherical harmonics as defined in Appendix B of Dahlen \& Tromp (1998). Such spherical harmonics
do include to Condon-Shortley Phase within their definition. 

## Installation

Clone the repository with:

git clone https://github.com/da380/TomoModel.git

Configure the build using:

cmake -S [Souce directory] -B [build directory]

Compile using:

cmake --build [build directory]

The executable is ReadS20RTS which is located in "bin" within the build directory. Example input files are found in "data" within the build directory. 
