# fracdiff

A set of tools for simulation of fractional Brownian diffusion subject to a harmonic potential and a finite "reaction" zone.

This code is developed and maintained by the Spakowitz Group at Stanford University

## Dependencies

+ Python 3.2 or greater (this code may work with 2.X, but it has not been tested)
+ Standard scientific Python packages: numpy, scipy, matplotlib
+ The notebook files ending in  .ipynb additionally require iPython notebook and its dependencies.
+ *(Optional)* PyPDF2 for the use of the fig_annotate function

*All of these packages are available on PyPI via 'pip install'*

## Contents

The main working directory has a main iPython notebook, **main_rxndiff_solver.ipynb**, that calls other Python .py files (which must be in the working directory). This code walks through the process of setting up an integrator with parameter values of interest, solving it for the full distribution p(r, t), and then manipulating this distribution in various ways to generate survival probability curves. the integration scheme currently being used calls the standard LSODA Fortran library (which is installed automatically in numpy), which switches between stiff and non-stiff solving depending on how the solution is behaving. There are additional options to try an exlicitly variable-step integrator, as well as my own implementation of Crank-Nicolson integrator with a "split operator." LSODA performs the best (and fastest) of all of these so it is used by default.

The file **rxn_diffusion_nondim.ipynb** contains a similar structure to the main notebook, but it represents a reparametricization of the diffusion equation so that the fractional time term appears on the reaction part, not the diffusion part.

The files **brownian_integrator.py**,**diffusion_integrator_funcs.py**, and **frac_brown.py** generally do not need to be edited unless serious under-the-hood changes to the code are necessary; they represent libraries of functions that get called in the main notebooks.

The file **rk4_demo.py** is just a proof-of-concept Runge-Kutta integration of the fractional diffusion equation. This goes unstable for most interesting parameter values, and so it mainly exists as a sanity check to make sure that it agrees with the fancier integrator for the range of aprameters in which it is stable.

The file **1D_reaction_diffusion.ipynb** is just a proof-of-concept showing that this integration seems to work just fine for one-dimensional diffusion in a harmonic potential.

The files **fig_annotate.py** and **plt_fmt.py** are non-essential and the code can by stripped of their use if they are causing problems. plt_fmt can be disabled by commenting out the lines at the top of hte iPython ipynb files in shich it gets imported. fig_annotate.py allows simulation parameter values to be saved as metadata in exported PDF figures, and it will just do nothing if it runs into any problems with dependencies.


The folder **old_code** has iPython notebooks and pure python libaries associated with the following operations:

+ Analytic solutions of the fractional diffusion equation using the methods of time-dependent perturbation theory. This code computes the perturbation matrix elements for various types of reaction well and guiding potential

+ Numerical integration of the fractional diffusion equation using Runge-Kutta with a Crank-Nicolson method, as well as a "split operator" method.