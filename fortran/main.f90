PROGRAM main
!!.................................................................................
!! Main file to a peridynamics simulation software.
!!
!! Made in MATLAB by:		Nicolas Sau
!! Translated into Fortran by:  Byron Encinas
!!				Luis Gutierrez
!! 				Rebeca Arredondo
!!				Carlos Muñoz
!!				Kevin Moreno
!!
!!The data given (in "prueba.csv" file) to plot the structure are:
!!	coords		node coordinates
!!	bc_vals		force, torque & its direction at every node
!!	bcs		free degrees of freedom (translational & rotational)
!!	elem_con	indicates how nodes are connected to form every element
!!	elem_props	material's physical properties
!!The variables and parameters used in this main file are:
!!	file_name	input data file
!!	i,j,k		loop counters
!!	stat		state variable
!!	num_nodes	number of nodes						
!!	num_elems	number of elements					*
!!	dims		dimensions						*
!!	time		computing time
!! It is important to note that this program is able to graph and calculate
!! displacements for structures consisting of 4-node elements, in 3 dimensions only*.
!!.................................................................................

!!.................................................................................
!!
!! El modulo iso_fortran_env se utiliza para rescatar el kind de los números reales
!!
USE iso_fortran_env, ONLY: real64, stdin=>input_unit
!!
!!.................................................................................

	IMPLICIT NONE
	
	REAL(real64),ALLOCATABLE:: coords(:,:), bcs_vals(:,:), bcs(:,:)
	REAL(real64),ALLOCATABLE:: elem_con(:,:), elem_props(:,:)
	REAL(real64),ALLOCATABLE:: nodal_disp(:,:), nodal_int_forces(:,:), new_elem_props(:,:)
	CHARACTER(LEN=20):: file_name
	INTEGER:: i,j,k, stat=0, num_nodes,num_elems,dims,deform
	REAL::time
	

	file_name = 'prueba.csv'	!Nombre del archivo de entrada
	num_nodes = 0				!Número de nodos
	num_elems = 2				!Número de elementos
	dims	  = 3				!Dimensiones

!!Las líneas a continuación permiten obtener el número de nodos a partir de lectura del archivo
!!de datos, generalización que se hará más adelante
!!------ ------	------	------	------	------	------	------	------	------	------
!!Contador de nodos

OPEN(UNIT = 12, FILE = File_name, action = 'read' , IOSTAT = stat)
	DO
	    READ(12,*, END = 10) 
	    num_nodes = num_nodes + 1
	END DO
10 CLOSE(12)

!! Number of nodes/dof

PRINT*, ""
PRINT*, " La red posee ", num_nodes, "nodos (main funciona sin errores 1)"
PRINT*, ""

!!------ ------	------	------	------	------	------	------	------	------	------

ALLOCATE(coords(num_nodes,dims))
ALLOCATE(bcs_vals(num_nodes,2*dims))
ALLOCATE(bcs(num_nodes,2*dims))
ALLOCATE(elem_con(num_nodes,num_elems))
ALLOCATE(elem_props(12,num_elems))

!!------ ------	------	------ ------ ------ ------	------ ------

	WRITE(*,*) '(main funciona sin errores 2)'

CALL load_data(file_name,coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,dims,stat)

	WRITE(*,*) '(main funciona sin errores 3)'

deform=0	
CALL draw_structure(coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,dims,nodal_disp,deform)

CALL CPU_TIME(time)			!Indicar tiempo
	PRINT*, ""
	WRITE(*,*)'Elapsed time', time
	PRINT*, ""

	WRITE(*,*)'Presione enter para continuar'

READ(stdin,*)

CALL perifea2d_linear(coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,&
dims,nodal_disp,nodal_int_forces,new_elem_props)


ALLOCATE(nodal_disp(num_nodes,dims))

nodal_disp(1,1) = 0
nodal_disp(1,2) = 0
nodal_disp(2,1) = 0.044782
nodal_disp(2,2) = -0.058702
nodal_disp(3,1) = 0
nodal_disp(3,2) = 0
nodal_disp(4,1) = -0.041105
nodal_disp(4,2) = -0.0598043
nodal_disp(5,1) = 0.054565
nodal_disp(5,2) = -0.16239
nodal_disp(6,1) = -0.032196
nodal_disp(6,2) = -0.18043
	
deform=1
CALL draw_structure(coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,dims,nodal_disp,deform)

CALL CPU_TIME(time)			!Indicar tiempo
	PRINT*, ""
!	WRITE(*,*)'Elapsed time', time
	PRINT*, ""

END PROGRAM main
