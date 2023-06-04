SUBROUTINE perifea2d_linear(coords,bcs_vals,bcs,elem_con,elem_props,num_nodes,num_elems,&
  dims,nodal_disp,nodal_int_forces,new_elem_props)
!!.................................................................................
!! perifea2d_lineal file to a peridynamics simulation software.
!!
!! Made in MATLAB by:		        		Nicolas Sau

!! Translated into Fortran by:  			Byron Encinas
!!				                       Luis Gutierrez
!!                       				Rebeca Arredondo
!!                      				Carlos Muñoz
!!                      				Kevin Moreno
!!
!!  The data given (in "prueba.csv" file) to plot the structure are:
!!	coords		node coordinates
!!	bc_vals		force, torque & its direction at every node
!!	bcs		free degrees of freedom (translational & rotational)
!!	elem_con	indicates how nodes are connected to form every element
!!	elem_props	material's physical properties
!!The variables and parameters used in this file are:
!!
!!.................................................................................
!! THIS CODE ASSUMES 2D SYSTEM WITH (x,y) COORDINATES
!!.................................................................................
!!
!! El modulo iso_fortran_env se utiliza para rescatar el kind de los números reales
!!
USE iso_fortran_env, ONLY: real64
!!
!!.................................................................................
    
    IMPLICIT NONE

    SAVE
!!------------     DUMMY VARIABLES BROUGHT FROM MAIN 	-----------------  
    INTEGER,INTENT(IN):: num_nodes,dims,num_elems
    REAL(real64), INTENT(IN):: coords(num_nodes,dims), bcs_vals(num_nodes,2*dims)
    REAL(real64), INTENT(IN):: bcs(num_nodes,2*dims)
    REAL(real64), INTENT(IN):: elem_con(num_nodes,num_elems), elem_props(12,num_elems)

    REAL(real64),INTENT(OUT):: nodal_disp(dims), nodal_int_forces(dims), new_elem_props(dims)

!!-------VARIABLES USADAS EN LA ASIGNACION DE F_free y D_fixed----------------
   
    INTEGER:: elems_load_factor
    INTEGER:: worst_elem, num_free_dof, num_fixed_dof 
    REAL(real64):: fixed_dof(num_nodes,dims), free_dof(num_nodes,dims) !! This is supposed to be allocated
    INTEGER:: node = 0

    REAL(real64), ALLOCATABLE:: kff(:,:), kfs(:,:), F_free(:), D_fixed(:)
    
!!-----------VARIABLES USADAS EN LA CONSTRUCCION DE kff y kfs-------------  

    INTEGER:: ielem_type, ielem, ielem_num_nodes, jelem_type, jelem, jelem_num_nodes

    REAL(real64), ALLOCATABLE::ielem_coords( : ), jelem_coords( : )
    INTEGER:: elem_row, el_row_nod, el_row_node, el_row, glob_row
    INTEGER:: elem_col, el_col_nod, el_col_node, el_col, glob_col

!!------------------------------------------------------------------------------  

    CHARACTER(LEN=50) :: Frmt1,Frmt2,Frmt3,Frmt4 	!! No se necesita
    INTEGER:: i=0,j=0,k=0,w=0,elem,elem_nodes, aux,w1,w2			!! No se necesita "stat"
    REAL::time	 					!! Opcional
    
!!------------------K_INDEXES RELATEX VARIABLES------------------------  
    
    INTEGER:: k_elem, INFO
    
    INTEGER(real64), ALLOCATABLE:: IPIV(:)

!!------------------------------------------------------------------------------  
    

    Frmt1 = '(f9.3,f9.3,f9.3,f9.3,f9.3,f9.3,f9.3,f9.3)'
    Frmt2 = ''
    Frmt3 = ''
    Frmt4 = ''

!!------ ------	------	------	------	------	------	------	------	------
  
  elems_load_factor = 0          !! No se usa en el primer segmento del codigo
  worst_elem = 0                 !! No se usa en el primer segmento del codigo

  free_dof(:,:) = 0

  num_free_dof = count(bcs(:,1:dims) == 0)
  num_fixed_dof = count(bcs(:,1:dims) == 1)

  ALLOCATE(kff(num_free_dof,num_free_dof))
  ALLOCATE(kfs( num_free_dof , num_fixed_dof ) )

kfs(:,:) = 0
kff(:,:) = 0  


  !! -------------------------TEST SECTION-----------------------------------------
  !WRITE(*,*) 'Should be 6x6 = ',size(bcs)
  !WRITE(*,*) 'Should be 8, 10',num_free_dof, num_fixed_dof
  !WRITE(*,*) 'Should be 2' , '',num_elems
  !WRITE(*,*) 'Should be 8^2, 8*10' , size(kff),size(kfs)
  !WRITE(*,*)''

  


!! -------------------------TEST SECTION-----------------------------------------


  ALLOCATE(F_free(num_free_dof))
  ALLOCATE(D_fixed(num_fixed_dof))
  
  D_fixed(:) = 0
  F_free(:) = 0
  
!! Vamos a contar la cantidad de zeros que tiene bcs(N,6)
  
WRITE(*,*) '(perifea2d_linear funciona sin errores 1) '

!! -------- Number the free degrees of freedom and fill the forcing vector ------  

  
  DO i = 1, num_nodes, 1
    DO j = 1, dims, 1
      IF (bcs(i,j) == 0) THEN
  
        !num_free_dof = num_free_dof + 1       !! Cuenta las entradas nulas en Ftranslacional de bcs
        free_dof(i,j) = num_free_dof          !! Almacena ese numero para una entrada especifica de bcs
        F_free(num_free_dof) = bcs_vals(j,i)  !! Almacena el vector de fuerza sobre los nodos libres
        
      ENDIF
    END DO
  END DO


!! -------- Number the fixed degrees of freedom and fill the forcing vector ------

  DO i = 1, num_nodes, 1 !! Prueba de 6 

    DO j = 1, dims, 1

      IF (bcs(i,j) /= 0) THEN 
	
        !num_fixed_dof = num_fixed_dof + 1       
        fixed_dof(i,j) = num_fixed_dof          
        D_fixed(num_fixed_dof) = bcs_vals(j,i)  
       
      ENDIF
    END DO
  END DO

WRITE(*,*) '(perifea2d_linear funciona sin errores 2)'

!PRINT*, num_fixed_dof


!! -----------------Assemble the global stiffness matrix -----------------------

WRITE(*,*) '(perifea2d_linear funciona sin errores 3)'


DO ielem = 1, 1,num_elems

  ielem_type = elem_props(1,ielem)
  ielem_num_nodes = elem_props(2,ielem)

  ALLOCATE(ielem_coords( INT( elem_props(2,ielem) ) * (dims-1))) ! 2D

!! -----------------ielem_coords OBTAIN -----------------------

  WRITE(*,*) '(perifea2d_linear funciona sin errores 4)'
  
  ielem_coords(:) = 0

  WRITE(*,*) '(perifea2d_linear funciona sin errores 5)'
  k = 1

            elem_nodes = INT(elem_props(2,elem))               ! one element corresponds with "elem_props(2,ielem)" nodes

            DO i=1,elem_nodes,1                             ! so we save up each node 
            
                node = elem_con(i,elem)                        ! each node is saved up in elem_con(ielem,:)
             
                DO j=1,dims-1,1

                  IF (coords(node,j) == 0) THEN
                  
                    CONTINUE
                  
                  ELSE
                                          
                    ielem_coords(k) = coords(node,j)
                    !ielem_coords(k) = 0
                     
                    !WRITE(*,*) node, ielem_coords(k)

                    k = k + 1

                  ENDIF
                
                END DO
                                              
            END DO


WRITE(*,*) '(perifea2d_linear funciona sin errores 6)'



!! ----------------------------------------------------------- -----------------------


  IF (ielem_num_nodes == 1) THEN 
    
    STOP

  END IF

  DEALLOCATE(ielem_coords)

  WRITE(*,*) '(perifea2d_linear funciona sin errores 7)'

  DO jelem = 1, 1,num_elems
    
    jelem_type = elem_props(1,jelem)
    jelem_num_nodes = elem_props(2,jelem)

    ALLOCATE(jelem_coords( INT( elem_props(2,jelem) ) * (dims-1) ) ) ! 2D

    
    !! -----------------ielem_coords OBTAIN -----------------------
  jelem_coords(:) = 0

  k = 1

WRITE(*,*) '(perifea2d_linear funciona sin errores 8)'

            elem_nodes = INT(elem_props(2,elem))               ! one element corresponds with "elem_props(2,ielem)" nodes

            !WRITE(25,*) '### Elemento', elem 
            
            !WRITE(*,*) ''

            !WRITE(*,*) '### Elemento', elem 

            DO i=1,elem_nodes,1                             ! so we save up each node 
            
                node = elem_con(i,elem)                        ! each node is saved up in elem_con(ielem,:)
                
                DO j=1,dims-1,1

                  IF (coords(node,j) == 0) THEN
                  
                    CONTINUE
                  
                  ELSE
                                          
                    jelem_coords(k) = coords(node,j)
                    !ielem_coords(k) = 0
                    
                    !WRITE(*,*) node, jelem_coords(k)

                    k = k + 1

                  ENDIF
                
                END DO
                                              
            END DO



!! ----------------------------------------------------------- -----------------------

WRITE(*,*) '(perifea2d_linear funciona sin errores 9)'

    IF (jelem_num_nodes == 1) THEN 
    
      STOP
  
    END IF

    DEALLOCATE(jelem_coords)


!! -------------------------------ROW CALCULATION MATLAB - 114 and on-----------------------


    IF (ielem /= jelem) THEN


      !CALL k_ij

      i = 0 !! in  MATLAB line 109 i = el_row_node
      j = 0 !! in  MATLAB line 117 j = el_row_node_dof
      k = 0 !! in  MATLAB line 123 k = el_col_node
      w = 0 !! in  MATLAB line 135 w = el_col_node_dof



      aux = INT(ielem_num_nodes + jelem_num_nodes)

      DO i = 1, aux, 1 !!el_row_node

        IF (i <= ielem_num_nodes) THEN

          elem_row = ielem
          el_row_nod = i !! el_row_node

        ELSE

          elem_row = jelem
          el_row_nod = i - ielem_num_nodes
          
        END IF

        DO j = 1, 3, 1 !! el_row_node_dof

          el_row =  3*(i  - 1) + j !! el_row_node - el_row_node_dof

          glob_row = free_dof(INT(elem_con(elem_row,el_row_nod) ), j) !! el_row_node_dof
        
          IF (glob_row /= 0)THEN

            DO k = 1 , aux, 1 !! el_col_node

              IF (k <= ielem_num_nodes)THEN

                elem_col = ielem
                el_col_nod = k !!el_col_node

              ELSE

                elem_col = jelem
                el_col_node = k - ielem_num_nodes

              END IF

              DO w = 1, 3, 1 !! el_col_node_dof

                el_col = 3*(k - 1) + w
                glob_col = free_dof(INT(elem_con(elem_col,el_col_nod) ), w)

                IF(glob_col /= 0)THEN

                  kff(glob_row,glob_col) = kff(glob_row,glob_col) !! + k_ijelem(el_row,el_col)

                ELSE

                  glob_col = fixed_dof( INT(elem_con(elem_col,el_col_nod)) ,w)

                  kfs(glob_row,glob_col - num_free_dof) = kfs(glob_row,glob_col - num_free_dof) !! + k_ijelem(el_row,el_col)



                END IF

              END DO

            END DO

          ENDIF

        END DO
        
      END DO

      i = 0 !! in  MATLAB line 109 i = el_row_node
      j = 0 !! in  MATLAB line 117 j = el_row_node_dof
      k = 0 !! in  MATLAB line 123 k = el_col_node
      w = 0 !! in  MATLAB line 135 w = el_col_node_dof
      
    ELSE

      !k_elem = k_ii_2D(ielem_coords,elem_props(ielem,:))

      DO i = 1, ielem_num_nodes, 1 !! el_row_node

        DO j = 1,3,1 !!el_row_node_dof

          el_row = 3*(i - 1) + j

          glob_row = free_dof(  INT(elem_con(ielem , i) ) ,j )      
          
          IF(glob_row /= 0) THEN

            DO k = 1, ielem_num_nodes, 1 !!el_col_node

              DO w = 1, 3, 1 !!el_col_node_dof

                el_col = 3*(k - 1) + w
                glob_col = free_dof(INT(elem_con(ielem,w)),k)

                IF(glob_col /= 0) THEN 
                
                  kff(glob_row,glob_col) = kff(glob_row,glob_col) ! + k_elem(el_row,el_col)

                ELSE

                  glob_col = fixed_dof( INT(elem_con(ielem,k)) ,w)
                  kfs(glob_row,glob_col) = kff(glob_row,glob_col) ! + k_elem(el_row,el_col)

                ENDIF

              END DO

            END DO

          END IF

        END DO        

      END DO

    END IF
!! ------------------------------------------------------------------------------------------

  END DO

END DO

!!--------------------------SPRANK FUNCTION-----------------------------

i = 0 
j = 0 
k = 0 
w1 = 0
w2 = 0 

!! --  Count null columns in kff and kfs


DO j = 1, num_free_dof,1

  !PRINT*, k,w
  
  DO i = 1, num_free_dof, 1
    
    IF( kff(i,j) == 0)THEN

    k = k + 1

    !PRINT*, kff(i,j)
    END IF

    IF (k == num_free_dof )THEN

      
      !PRINT*, kff(i,j)

      w1 = w1 + 1
    
    ENDIF
  
  ENDDO

  k = 0

ENDDO

i = 0 
j = 0 
k = 0 

DO j = 1, num_fixed_dof,1

  !PRINT*, k,w2
  
  DO i = 1, num_free_dof, 1
    
    IF( kfs(i,j) == 0)THEN

    k = k + 1

    !PRINT*, kfs(i,j)
    END IF

    IF (k == num_free_dof )THEN

      
      !PRINT*, kff(i,j)

      w2 = w2 + 1
    
    ENDIF
  
  ENDDO

  k = 0

ENDDO

!PRINT*, w1,w2 !! SPRANK(kff),SPRANK(kfs)

!! --  Count null columns in kff and kfs

!WRITE(*,*) k

!!----------------------------------------------------------------------------------------
!!--------------------------NODAL DISPLACEMENTS CALCULATION-----------------------------

!WRITE(*,*) size(kfs),size(D_fixed) !! 8x10 - 5x1
!WRITE(*,*) size(kff),size(F_free) 
!WRITE(*,*) size(Matmul(kfs,D_fixed))
!WRITE(*,*) size(Matmul(kff,F_free))

ALLOCATE(IPIV(num_free_dof))

IF(w1 == 0)THEN
  
  ! disp_free = 0
  WRITE(*,*) 'there are no null columns in kff (line 515 perifea2d_linear.f90)'

ELSEIF (w1 == num_free_dof)THEN

  !CALL SGESV( num_free_dof, num_free_dof, kff , num_free_dof, &
  !IPIV, F_free - Matmul(kfs,D_fixed), num_free_dof, INFO	)

  !disp_free = x : kff * x = (F_free - Matmul(kfs,D_fixed))
  WRITE(*,*) 'there are no null columns in kff (line 515 perifea2d_linear.f90)'
  
ELSE 

  WRITE(*,*) 'the structure is unstable'
  !nodal_disp =

ENDIF


!!----------------------------------------------------------------------------------------




RETURN
  
END SUBROUTINE perifea2d_linear
