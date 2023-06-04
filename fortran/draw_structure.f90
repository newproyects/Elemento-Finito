SUBROUTINE draw_structure(coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,dims,nodal_disp,deform)

!!.................................................................................
!! draw_structure file to a peridynamics simulation software. Here we make a dat file used for 
!! gnuplot component
!!
!! Made in MATLAB by:		Nicolas Sau
!! Translated into Fortran by:  Byron Encinas
!!				Luis Gutierrez
!! 				Rebeca Arredondo
!!				Carlos Munoz
!!				Kevin Moreno
!!
!! The data loaded by load_data is the input and the output is a dat file with only coordinates
!! for the lines we want to plot:
!!	coords		node coordinates
!!	bc_vals		force, torque & its direction at every node
!!	bcs		free degrees of freedom (translational & rotational)
!!	elem_con	indicates how nodes are connected to form every element
!!	elem_props	material's physical properties
!!	num_nodes	total number of nodes
!!	num_elems	total number of elements
!!	dims		total number of dimensions of interest
!!The variables and parameters used locally are:
!!	end_coorda	an extra set of coordinates used for plotting line segments
!!	i,j,k		loop counters
!!	node		current node identifier						
!!	elem_nodes	number of nodes in the current element
!!
!!.................................................................................
	
!! El modulo iso_fortran_env se utiliza para rescatar el kind de los numeros reales

USE iso_fortran_env, ONLY: real64

  IMPLICIT NONE
 
  ! Dummy variables brought from the main program
  INTEGER, INTENT(IN) :: num_nodes,num_elems,dims,deform
  REAL(real64), INTENT(IN):: coords(num_nodes,dims), bcs_vals(num_nodes,2*dims), bcs(num_nodes,2*dims)
  REAL(real64), INTENT(IN):: elem_con(num_nodes,num_elems), elem_props(12,num_elems)
  REAL(real64), INTENT(IN):: nodal_disp(num_nodes,dims)
  
  ! Local variables
  REAL (real64) :: end_coords(num_nodes,dims)
  INTEGER:: i,j,k, node, elem_nodes
  
  
!We generate the dat file for the coordinates, which will later be used by gnuplot
OPEN(UNIT=25, FILE='draw_structure.dat', STATUS='REPLACE', ACTION = 'WRITE')
	
	! We identify the purpose of the file, it will be separated by double blank spaces to use the index function in gnuplot
	WRITE(25,*) '# Este documento es usado por gnuplot para graficar los nodos y los elementos'
	WRITE(25,*) '# La informacion contenida son solo puntos y no contienen informacion fisica sobre el sistema'
	WRITE(25,*) ''
	
	! write the coordinates of every node
	WRITE(25,*) '##### Nodos'
		DO i=1,num_nodes,1
			WRITE(25,*) (coords(i,k), k=1, dims)
		END DO
	WRITE(25,*) ''
	WRITE(25,*) ''
	
	! write coordinates of the nodes in every element to join them in the appropiate order, adding a single blank space between each element 
	WRITE(25,*) '##### Elementos'
		DO j=1,num_elems,1
			WRITE(25,*) '### Elemento', j
			elem_nodes=INT(elem_props(2,j))
			DO i=1,elem_nodes,1
				node = elem_con(i,j)
				WRITE(25,*) (coords(node,k), k=1, dims)
			END DO
			node = elem_con(1,j)
			WRITE(25,*) (coords( node ,k), k=1, dims)
			WRITE(25,*) ''
		END DO
		WRITE(25,*) ''
		
	! write bcs where they are appropiate. scale can be changed by the numbre multipliying bcs array
	WRITE(25,*) '##### Fijeza de los nodos (bcs)'
	end_coords(:,:)=coords(:,:)
		WRITE(25,*) '### Traslacionales'
		DO i=1,num_nodes,1
			WRITE(25,*) '# Traslacional', i
				DO j=1,dims,1
					WRITE(25,*) (coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)+0.3*bcs(i,j)
					WRITE(25,*) (end_coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)
					WRITE(25,*) ''
				END DO
		END DO
		WRITE(25,*) ''
		WRITE(25,*) '### Rotacionales'
		DO i=1,num_nodes,1
			WRITE(25,*) '# Rotacional', i
				DO j=1,dims,1
					WRITE(25,*) (coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)+0.3*bcs(i,j+3)
					WRITE(25,*) (end_coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)
					WRITE(25,*) ''
				END DO
		END DO
	WRITE(25,*) ''
	
	! write bcs vals
	WRITE(25,*) '##### Esfuerzos y momentos (bc values)'
		WRITE(25,*) '### Fuerzas traslacionales'
		DO i=1,num_nodes,1
			WRITE(25,*) '# Fuerzas', i
				DO j=1,dims,1
					WRITE(25,*) (coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)+1.0*bcs_vals(i,j)
					WRITE(25,*) (end_coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)
					WRITE(25,*) ''
				END DO
		END DO
		WRITE(25,*) ''
		WRITE(25,*) '### Torcas'
		DO i=1,num_nodes,1
			WRITE(25,*) '# Momentos', i
				DO j=1,dims,1
					WRITE(25,*) (coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)+1.0*bcs_vals(i,j+3)
					WRITE(25,*) (end_coords(i,k), k=1, dims)
					end_coords(i,j)=coords(i,j)
					WRITE(25,*) ''
				END DO
		END DO
		WRITE(25,*) ''
		
	! write coordinates of the nodes in every deformed element, only on the second plotting
	IF (deform==1) THEN 
	WRITE(25,*) '##### Elementos deformados'
		DO j=1,num_elems,1
			WRITE(25,*) '### Elemento', j
			elem_nodes=INT(elem_props(2,j))
			DO i=1,elem_nodes,1
				node = elem_con(i,j)
				WRITE(25,*) (coords(node,k)+nodal_disp(node,k), k=1, dims)
			END DO
			node = elem_con(1,j)
			WRITE(25,*) (coords( node ,k)+nodal_disp(node,k), k=1, dims)
			WRITE(25,*) ''
		END DO
		WRITE(25,*) ''
	END IF
CLOSE(25)

! here we call the gp file to plot all the data above
	
CALL SYSTEM("gnuplot -p draw_structure.gp")

!WRITE(*,*) '(draw_structure funciona sin errores 1)'


RETURN 

END SUBROUTINE draw_structure

