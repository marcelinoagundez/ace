# ACE (Atmospheric Chemical Equilibrium) #

  Fortran code that computes the chemical equilibrium composition for a given
   - elemental composition
   - pressure-temperature profile

  The code can be typically used to compute the chemical equilibrium composition in atmospheres of
  AGB stars, brown dwarfs, and exoplanets.

  The code is based on the algorithm implemented in the NASA/CEA program (Gordon & McBride 1994, NASA Reference Publication 1311, 1).

# If you use ACE for your work, I would appreciate if you cite the papers
   - Agundez et al. 2020, A&A, 637, A59
   - Agundez 2025, A&A, submitted

# Inquiries can be sent to marcelino.agundez@csic.es

# Compilation:

  Download and uncompress the zip file. This will create three directories:

   - `src`
   - `examples`

  You need to have gfortran installed in your system. Compilation has been tested with gfortran.
  Other fortran compilers such as ifort may work but they have not been tested.
  To compile simply type:

`
$ cd src
$ make clean       (only needed if the code has been compiled previously)
$ make
`

  This should have created a binary file named `ace`.
  It is convenient to move the binary file to the bin folder of your system where binaries are located.
  This allows you to run the code simply typing `ace` in any folder without the need to duplicate binaries.

# Usage:

  In the directory `examples` you have two sub-directories named `agb` and `exoplanet`.
  These sub-directories contain specific input files to run ACE for AGB and exoplanet atmospheres.
  The file `ace.inp`contains the input file names of the models to be run, e.g., `model_xx.inp`.
  Each input file `model_xx.inp` contains in turn three file names with extensions:
   - `.apt`        (altitude-pressure-temperature profile, the altitude is a dummy variable)
   - `.spec`       (species file)
   - `.therm`      (thermochemical data file)
 
  To run ACE, move to corresponding sub-directory and simply type:
`
$ ace              (assuming you have the binary in the bin folder of your system)
`
