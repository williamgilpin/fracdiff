# fracdiff

A set of tools for simulation of fractional Brownian diffusion subject to a harmonic potential and a finite "reaction" zone.

## Dependencies

+ Python 3.2 or greater (this code may work with 2.X, but it has not been tested)
+ Standard scientific Python packages: numpy, scipy, matplotlib
+ The notebook files ending in  .ipynb require iPython notebook and its dependencies.

## Contents

The main working directory has a main iPython notebook, **frac_integrate.ipynb**, that calls other Python .py files (which must be in the working directory). This code walks through the process of setting up an integrator with parameter values of interest, solving it for hte full distribution $p(r, t)$, and then manipulating this distribution in various ways to generate survival probability curves.


The folder **old_code** has iPython notebooks and pure python libaries associated with the following operations:

+ Analytic solutions of the fractional diffusion equation using the methods of time-dependent perturbation theory. This code computes the perturbation matrix elements for various types of reaction well and guiding potential

+ Numerical integration of the fractional diffusion equation using Runge-Kutta with a Crank-Nicolson method, as well as a "split operator" method.