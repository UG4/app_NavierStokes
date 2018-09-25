# navier_stokes #

**UG4-App** providing example scripts for Navier-Stokes based problems.

Copyright 2010-2017 Goethe Center for Scientific Computing, University Frankfurt

Please install/clone this repository through UG4's package manager
[ughub](https://github.com/UG4/ughub):

    ughub install navier_stokes


## How to use the navier_stokes scripts ##
Once you compiled UG4 (with the cmake option -DNavierStokes=ON)
and after sourcing 'ugbash'
(cf. https://github.com/UG4/ughub#compilation-of-ug4),
you may execute a script like this (from any folder you like):

    ugshell -ex navier_stokes/drivencavity.lua
