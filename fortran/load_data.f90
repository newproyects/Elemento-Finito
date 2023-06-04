SUBROUTINE load_data(file_name,coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,dims,stat)

!!.................................................................................
!! Load_data file to a peridynamics simulation software. Here we define arrays
!! that correspond with the data of our system.
!!
!! Made in MATLAB by:		Nicolas Sau
!! Translated into Fortran by:  Byron Encinas
!!				Luis Gutierrez
!! 				Rebeca Arredondo
!!				Carlos Munoz
!!				Kevin Moreno
!!The data given (in "prueba.csv" file) to plot the structure, is read and stored
!!in arrays previously defined in the main program:
!!	coords		node coordinates
!!	bc_vals		force, torque & its direction at every node
!!	bcs		    free degrees of freedom (translational & rotational)
!!	elem_con	indicates how nodes are connected to form every element
!!	elem_props	material's physical properties
!!The variables and parameters used in this load_data file are:
!!	file_name	input data file
!!	i,j,k		loop counters
!!	stat,msg	state variables
!!	num_nodes	number of nodes						
!!	num_elems	number of elements					
!!	dims		dimensions						
!!	Formato		printing format
!!
!!.................................................................................

!!.................................................................................
!!
!! El modulo iso_fortran_env se utiliza para rescatar el kind de los números reales
!!
USE iso_fortran_env, ONLY: real64
!!
!!.................................................................................

IMPLICIT NONE

SAVE

INTEGER,INTENT(INOUT):: num_nodes,num_elems,dims
REAL(real64), INTENT(INOUT):: coords(num_nodes,dims), bcs_vals(num_nodes,2*dims)
REAL(real64), INTENT(INOUT):: bcs(num_nodes,2*dims)
REAL(real64), INTENT(INOUT):: elem_con(num_nodes,num_elems), elem_props(12,num_elems)
INTEGER:: i,j,k
INTEGER, INTENT(INOUT) :: stat
CHARACTER(LEN=50) :: msg
CHARACTER(LEN=20), INTENT(IN) :: file_name


!!Nulificar arrays

coords(:,:) 	= 0.0
bcs_vals(:,:) 	= 0.0
bcs(:,:) 	= 0.0
elem_con(:,:) 	= 0.0
elem_props(:,:) = 0.0


!! Initial conditions are pulled into the code from a .csv file

OPEN(UNIT=13, FILE = File_name, ACTION = 'READ' )
	
	IF (stat == 0) THEN 
	
	!! We gotta pull data into arrays Coords, bcs_vals(forces), bcs(fixity points)
	
		WRITE(*,*) ' Data from disk has been loaded (load_data funciona sin errores 1)' 
          
		DO i=1,num_nodes,1
	
			READ(13,*) (coords(i,k), k=1, dims),(bcs_vals(i,k), k=1, 2*dims),(bcs(i,k), k=1, 2*dims),(elem_con(i,k), k=1, num_elems)
	
		END DO
	
	ELSE
	
		WRITE(*,*) 'There has been an error loading data from disk' 
	
	ENDIF
	
CLOSE(13)


!! Finally we pull the elem_props, which are macro properties of the material
!! in this case concrete
elem_props(1,1) = 3			!! Tipo de elemento "membranas"
elem_props(2,1) = 4 			!! Numero de nodos que hacen un elemento 
elem_props(3,1) = 100 			!! 
elem_props(4,1) = 0			!!
elem_props(5,1) = 1			!!
elem_props(6,1) = 3			!!
elem_props(7,1) = 5			!!
elem_props(8,1) = 0.8			!!
elem_props(9,1) = 0.001			!!
elem_props(10,1) = 0.001		!!
elem_props(11,1) = 1			!!
elem_props(12,1) = 1			!!
elem_props(:,2) = elem_props(:,1)	!!

RETURN 

END SUBROUTINE load_data
