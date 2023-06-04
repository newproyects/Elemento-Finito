SUBROUTINE elem_distance(vec1, vec2, dis_ij) 

	USE iso_fortran_env, ONLY: real64	
		
	IMPLICIT NONE

  REAL(real64),INTENT(IN)::vec1(3),vec2(3)
  REAL(real64),INTENT(OUT)::dis_ij
  REAL(real64)::a
  
	 a = (  vec1(1)  -  vec2(1)  )**2 & 
	   + (  vec1(2)  -  vec2(2)  )**2 & 
	   + (  vec1(3)  -  vec2(3)  )**2
	 
	dis_ij = Sqrt( a )

!! COPY IN main.f90 to test

!	REAL(real64):: distance
!	WRITE(*,'(F9.3,F9.3,F9.3)') coords(1,:)
!	WRITE(*,'(F9.3,F9.3,F9.3)') coords(2,:)
!	Call elem_distance(coords(1,:),coords(2,:), distance)
!	WRITE(*,*) distance**2

END SUBROUTINE elem_distance


