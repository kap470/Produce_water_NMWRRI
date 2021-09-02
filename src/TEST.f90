PROGRAM Random
 
  INTEGER, PARAMETER :: n = 300
  INTEGER :: i
  REAL :: array(n), pi, temp, mean = 600.0, sd = 500.0
 
  pi = 4.0*ATAN(1.0)
  CALL RANDOM_NUMBER(array) ! Uniform distribution
 
! Now convert to normal distribution
  DO i = 1, n-1, 2
    temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
    array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
    array(i) = temp
	do while (array(i)<0) 
		CALL RANDOM_NUMBER(array(i))
		temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
		array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
		array(i) = temp
	end do
	WRITE(*,*)array(i)
  END DO
 
! Check mean and standard deviation
  mean = SUM(array)/n
  sd = SQRT(SUM((array - mean)**2)/n)
 
  WRITE(*, "(A,F8.6)") "Mean = ", mean
  WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd
 
END PROGRAM Random