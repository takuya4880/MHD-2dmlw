main: main.f90
	gfortran -O3 -fopenmp defstruct.f90 initial_jet.f90 boundary_parker.f90 detdt.f90 flux.f90 source.f90 lw1.f90 lw2.f90 artvis.f90 pressure.f90 step.f90 output.f90 dacdef0s.f dacdef2s.f dacputparamd.f dacdef1d.f dacdefparam.f dacputparami.f main.f90 
