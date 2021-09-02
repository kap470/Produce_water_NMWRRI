!Author: Kevin Perez
!Email: kap470@nmsu.edu
!New Mexico Water Resources Research Institute 

MODULE FUNCTIONS
    
	contains
	
	!This function creates a matrix from a points database
	!through the assignation of 1 values to those pixels containing the points
	!can be used to create a mask raster layer related to a points database
	
	!x is an array with X points coordinates
	!Y is an array with Y points coordinates
	!n is dimension of the arrays
	SUBROUTINE CREATE_MATRIX_INT(Ncols,Nrows,Xcor,Ycor,delta,x,y,n,A)
	
		IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(IN)                              	:: Ncols,Nrows,n
        REAL,INTENT(IN)                                 	:: Xcor,Ycor,delta
		REAL,DIMENSION(:),ALLOCATABLE,INTENT(IN)   		    :: x,y
        INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: A
		!Local variables section
		REAL												::xi,xf,yi,yf
		INTEGER												::i,j,k
	
		DO i=1,Nrows
			yi=Ycor+(Nrows-1)*delta-delta/2-(i-1)*delta
			yf=Ycor+(Nrows-1)*delta+delta/2+-(i-1)*delta
			DO j=1,Ncols
				xi=Xcor-delta/2+(j-1)*delta
				xf=Xcor+delta/2+(j-1)*delta
				A(i,j)=0
				DO k=1,n
					IF (xi<x(k).and.x(k)<xf) THEN 
						IF (yi<y(k).and.y(k)<yf) THEN 
							A(i,j)=1
						END IF
					END IF
				END DO
			END DO
		END DO
	
	END SUBROUTINE CREATE_MATRIX_INT
	
	SUBROUTINE CREATE_RANDOM_NUMBERS(owf,mean,sd,array)
		 
		IMPLICIT NONE
		 
		INTEGER, INTENT(IN) :: owf
		REAL, INTENT (IN)  ::mean, sd 
		REAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: array
		REAL :: pi, temp 
		INTEGER :: i,n
		
		n=owf
		 
		pi = 4.0*ATAN(1.0)
		CALL RANDOM_NUMBER(array) ! Uniform distribution
		 
		! Now convert to normal distribution
		DO i = 1, n-1, 2
			temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
			array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
			array(i) = temp
			IF (array(i)<0) THEN
				array(i) = array(i)*(-1)
			END IF
			IF (array(i+1)<0) THEN
				array(i+1) = array(i+1)*(-1)
			END IF
		END DO
	END SUBROUTINE CREATE_RANDOM_NUMBERS
	
	SUBROUTINE CREATE_RANDOM_LOCATION(Ncols,Nrows,Xcor,Ycor,n,delta,x,y)
		
		IMPLICIT NONE
		
		REAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: x,y
		REAL, INTENT(IN) :: delta,Xcor,Ycor
		INTEGER, INTENT(IN) :: Ncols,Nrows,n
		INTEGER :: i
		
		DO i=1,n
            !Random locations in (0,1) interval
            CALL RANDOM_NUMBER(X(i))
            CALL RANDOM_NUMBER(Y(i))
            !Random locations in real spatial domain
            X(i)=X(i)*delta*Ncols+Xcor
            Y(i)=Y(i)*delta*Nrows+Ycor
        END DO
        
	END SUBROUTINE CREATE_RANDOM_LOCATION
	
END MODULE

