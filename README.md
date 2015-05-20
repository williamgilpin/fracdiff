# fracdiff

A set of tools for simulation of fractional Brownian diffusion subject to a harmonic potential and a finite "reaction" zone.

## Dependencies

+ Python 3.2 or greater (this code may work with 2.X, but it has not been tested)
+ Standard scientific Python packages: numpy, scipy, matplotlib
+ The notebook files ending in  .ipynb require iPython notebook and its dependencies.
+ PyPDF2 for the use of the optional fig_annotate function

## Contents

The main working directory has a main iPython notebook, **frac_integrate.ipynb**, that calls other Python .py files (which must be in the working directory). This code walks through the process of setting up an integrator with parameter values of interest, solving it for hte full distribution $p(r, t)$, and then manipulating this distribution in various ways to generate survival probability curves.

The files **fig_annotate.py** and **plt_fmt.py** are non-essential and the code can by stripped of their use if they are causing problems. plt_fmt can be disabled by commenting out the lines at the top of hte iPython ipynb files in shich it gets imported. fig_annotate.py allows simulation parameter values to be saved as metadata in exported PDF figures, and it can be disabled by commenting out any line on which it appears.


The folder **old_code** has iPython notebooks and pure python libaries associated with the following operations:

+ Analytic solutions of the fractional diffusion equation using the methods of time-dependent perturbation theory. This code computes the perturbation matrix elements for various types of reaction well and guiding potential

+ Numerical integration of the fractional diffusion equation using Runge-Kutta with a Crank-Nicolson method, as well as a "split operator" method.