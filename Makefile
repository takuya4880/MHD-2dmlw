main: main.f90
	gfortran -O3 defstruct.f90 initial_jet.f90 boundary_parker.f90 detdt.f90 flux.f90 source.f90 rev_lw1.f90 rev_lw2.f90 artvis.f90 pressure.f90 step.f90 output.f90 main.f90 
