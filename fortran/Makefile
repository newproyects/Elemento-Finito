ejecutable: main.f90
	gfortran -c draw_structure.f90 elem_distance.f90 load_data.f90 perifea2d_linear.f90
	gfortran main.f90 draw_structure.o elem_distance.o load_data.o perifea2d_linear.o -o Peridynamics
	rm *.o
