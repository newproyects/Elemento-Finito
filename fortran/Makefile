##
## Makefile
##  
##
## Started on  Noviembre 2022 Grupo Simulacion Peridinanimca
## Last update Time-stamp: <18-nov-2022 15:33 Voronoi>
##
#
# This makefile is intended for the program Peridynamics Simulation Software
#
#===============================================
# Fortran Parameters
#===============================================
FC = gfortran #Compiler
FFLAGS = -O3  #Optimization Flags (maybe tinker with these?)
FSRC = main.f90 load_data.f90 elem_distance.f90 draw_structure.f90 perifea2d_linear.f90 #SGESV.f
FOBJECTS = $(FSRC:.f90=.o) #Objets
#===============================================
# C Parameters
#===============================================
CC = gcc #Compiler
CFLAGS = -O   #-ffixed-line-length-132 #Optimization Flags (maybe tinker with these?)
CSRC = 
COBJECTS = $(CSRC:.c=.u) # C Objets 
CHEAD = 

#===============================================
# General Parameters.  
#===============================================
OBJECTS =	$(FOBJECTS) $(COBJECTS)
EXEC =  Peridynamics

#===============================================
# Compilation rules
#===============================================
default: $(EXEC)

clean:	
	rm -f *~ $(OBJECTS) *.pro *.ascii *.nat *.mod *.Gh

spotless:	
	rm -f *~ $(OBJECTS) $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) $(FLFLAGS) -o $(EXEC) \
		 $(OBJECTS) -lm

%.o:  %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.u:	%.c $(CHEAD)	
	$(CC) $(CFLAGS) -c $< -o $@ 
