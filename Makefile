main: main.f90
	gfortran -O3 -fopenmp defstruct.f90 initial_jet.f90 boundary_parker.f90 detdt.f90 lw.f90 pressure.f90 step.f90 cansio.f output.f90 main.f90 
