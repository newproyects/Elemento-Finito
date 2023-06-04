SUBROUTINE peri_quantsxi_i1, xi_i2 , elcoords )

  !saca xypos, N, det_J

  IMPLICIT NONE

  INTEGER :: row,col,node, num_nodes, i, j

  REAL(real64), INTENT(IN) :: elcoords(8,1), xi_i1, xi,i2

  REAL(real64) :: J(2,2) ,N(3,8) , N_derivs(4,2) ,xi(1,2), xypos(1,2), J_adj(2,2),  inv_J(2,2)

  REAL(real64) :: ms, ps, mt, pt, det_J


  xi(:,:) = (0.0)
  
  num_nodes = size(elcoords)/2
  
 ! elcoords(1) = 4.2543
 ! elcoords(2) = 7.1502
 ! elcoords(3) = 4.3416
 ! elcoords(4) = 4.9074
 ! elcoords(5) = 7.2016
 ! elcoords(6) = 5.0514
 ! elcoords(7) = 7.0370
 ! elcoords(8) = 7.5617

    ms = 1-xi(1,1)
    ps = 1+xi(1,1)
    mt = 1-xi(1,2)
    pt = 1+xi(1,2)
    N1 = ms*mt/4

    N(:,:) = 0.0

    N(1,1) = (ms*mt)/4

    N(1,3) = (ps*mt)/4

    N(1,5) = (ps*pt)/4

    N(1,7) = (ms*pt)/4
    

    N(2,2) = N(1,1)

    N(2,4) = N(1,3)

    N(2,6) = N(1,5)

    N(2,8) = N(1,7)

    xypos(:,:) = 0.0

    DO i= 1,7,2
       N_1(1,i) = N(1,i)
       elcoords_1(i,1) = elcoords(i,1)

   END DO

   DO i= 2,8,2
      N_2(1,i) = N(1,i)
      elcoords_2(i,1) = elcoords(i,1)
   END DO

   DO i=1,3,1
      xypos(1,1) = N_1(1,i) * elcoords_1(i,1)
      xypos(1,2) = N_2(1,i) * elcoords_2(i,1)
   END DO

    
    N_derivs(:,:) = 0.0

    N_derivs(1, 1) = -0.25*(1 - xi(1,2))
    N_derivs(2, 1) =  0.25*(1 - xi(1,2))
    N_derivs(3, 1) =  0.25*(1 + xi(1,2))
    N_derivs(4, 1) = -0.25*(1 + xi(1,2))
    N_derivs(1, 2) = -0.25*(1 - xi(1,1))
    N_derivs(2, 2) = -0.25*(1 + xi(1,1))
    N_derivs(3, 2) =  0.25*(1 + xi(1,1))
    N_derivs(4, 2) =  0.25*(1 - xi(1,1))

     
  J(:,:) = 0.0
  
  IF (num_nodes > 2) THEN
     DO row = 1,2,1
        DO col = 1,2,1
           DO node = 1,num_nodes,1
              J(row, col) = J(row,col) + N_derivs(node,row)*elcoords(2*(node-1) +col)
           END DO

        END DO

     END DO     
           
  END IF


   adj_J(:,:) = 0.0
  
  DO row = 1,2,1
     DO col = 1,2,1
      adj_J(col,row) = J(row,col)
     END DO
  END DO

  det_J = J(1,1)*J(2,2) - J(1,2)*J(2,1)


  DO row = 1,2,1
     DO col = 1,2,1
      inv_J(row,col) = (1/det_J) * adj_J(row,col)
     END DO
  END DO
  

  

END SUBROUTINE peri_quants



  
