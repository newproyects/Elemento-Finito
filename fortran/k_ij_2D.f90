SUBROUTINE k_ij_2D(icoords,jcoords,ielem_props,jelem_props)

  IMPLICIT NONE

  !ielem_coords = icoords
  !jelem_coords = jcoords
  !elem_props(ielem) = ielem_props
  !elem_props(jelem) = jelem_props

  !   Variables de entrada
  !icoords, jcoords

  !   Variables de salida
  !el_stiff_c1 

  !   Cosas
  INTEGER :: , num_nodes_i, c_i, d_i, t_i, n_i, delta_i, as_i, num_nodes_i, c_j, d_j, t_j, n_j, delta_j, as_j, c, d , del_i, del_j

  INTEGER :: i1 , i2, j1, j2

  REAL(real64) ::ielem_props(12,2) , jelem_props(12,2)

  REAL(real64), INTENT(OUT) :: xi_i1 , xi_i2,xi_j1 , xi_j2

  num_nodes_i = ielem_props(2,1)

  num_nodes_j = jelem_props(2,2)
  
  c_i = ielem_props(3,1)
  d_i = ielem_props(4,1)
  t_i = ielem_props(5,1)
  n_i = ielem_props(6,1)
  delta_i = ielem_props(7,1)
  
  as_i = ielem_props(5,1)
  as_j = jelem_props(5,2)

  c_j = jelem_props(3,2)
  d_j = jelem_props(4,2)
  t_j = jelem_props(5,2)
  n_j = jelem_props(6,2)
  delta_j = jelem_props(7,2)

  c = (c_i + c_j)/2
  d = (d_i + d_j)/2
  mat_horizon = min(delta_i, delta_j)

  el_dist = elem_distance(icoords, j coords)

  !IF (el_dist > mat_horizon) THEN

!K

! Calcula las xi_i1 y *2 para mandarla al periquants 
DO i1 = 1,n_i,1
   DO i2 = 1,n_i,1 
      del_i = 2/n_i
      xi_i1 = -1 + del_i*i1 - del_i/2
      xi_i2 = -1 + del_i*i2 - del_i/2

      CALL peri_quants( xi_i1 , xi_i2  , icoords)
      !Usando un dato salida de periquants y ielemprop 5
      dvol_i = det_J_i*t_1*(del_i**2)

 !Calcula las xi_j1 y *2 para mandarla al periquants 
      DO j1 = 1,n_j,1
         DO j2 = 1,n_j,1 
            del_j = 2/n_j
            xi_j1 = -1 + del_j*j1 - del_j/2
            xi_j2 = -1 + del_j*j2 - del_j/2

            CALL peri_quants( xi_j1 , xi_j2  , jcoords)
            !Usando dato salida de periquants y t_j  
            dvol_j = det_J_j*t_j*(del_j**2)

            !Usando coor de las dos salidas de periquants
            del_coor = coor_j - coor_i 
            L = norm(del_coor)


   END DO
END DO

    !falta
!dvol_j
!det_J
!del_coor
!coor_j
!coor_j
!L
  

  


  

end SUBROUTINE k_ij_2D


