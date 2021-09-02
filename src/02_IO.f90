!Author: Kevin Perez
!Email: kap470@nmsu.edu
!New Mexico Water Resources Research Institute 

MODULE IO

    contains

    !OWR=Oil Water Ratio
    !This routine reads the inputs contained in a txt file
    SUBROUTINE read_input_file(ow,gw,iw,tini,tend,dt,Ncols,Nrows,Xcor,Ycor,delta,nan,RT,cost,OWS0&
	&,INJS0,NTP,TPCAP,RTPT,RTPO,TPS0)
    
        IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(OUT)                             :: Ncols,Nrows,nan,ow,gw,iw,NTP
        REAL,INTENT(OUT)                                :: Xcor,Ycor,delta,cost,tini,tend,dt
        REAL,INTENT(OUT)                                :: OWS0,INJS0,TPCAP,RTPT,RTPO,TPS0
        CHARACTER(200),INTENT(IN)                       :: RT
        !Local variables section
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        
        OPEN(24, FILE=RT, STATUS='OLD',ACTION='READ')
        !Reading values
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,ow
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,gw
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,iw
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.3))')aux,tini
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.3))')aux,tend 
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.3))')aux,dt
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,Ncols
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,Nrows
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.3))')aux,Xcor
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.3))')aux,Ycor
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.2))')aux,delta
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,nan
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.2))')aux,cost
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.2))')aux,OWS0
		READ(24, '(A6)')aux
        READ(24, '((A6),(F10.2))')aux,INJS0
		READ(24, '(A6)')aux
        READ(24, '((A6),(I10))')aux,NTP
		READ(24, '(A6)')aux
		READ(24, '((A6),(F10.2))')aux,TPCAP
		READ(24, '(A6)')aux
		READ(24, '((A6),(F10.2))')aux,RTPT
		READ(24, '(A6)')aux
		READ(24, '((A6),(F10.2))')aux,RTPO
		READ(24, '(A6)')aux
		READ(24, '((A6),(F10.2))')aux,TPS0

        CLOSE(24)
        
        WRITE(*,*)"******* The input file has been read succesfully *******"

    END SUBROUTINE read_input_file
        
	!*******************************************************
	!****Function to read Oil Wells outflow coefficients****
	!*******************************************************
	!ny number of years
    SUBROUTINE READ_OW_COEFF(ny,RT,years,a,b,c)

        USE TYPES

        IMPLICIT NONE
        !Global variables
        INTEGER,INTENT(IN)                              :: ny
        CHARACTER(200),INTENT(IN)                       :: RT
        REAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT)     :: a,b,c,years
        !Local variables
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(12)                                   ::intfile
        
        OPEN(24, FILE=RT, STATUS='OLD',ACTION='READ')
        !Lectura de metadatos
        READ(24, '((A14),(I10))')aux
        !Se crea internal file para descriptor de formato
        WRITE(intfile,'(I12)')4
        !Lectura de matriz
		
        DO i=1,ny
            READ(24,'('//intfile//'((ES12.2),(1X),(ES12.2),(1X),(ES12.2),(1X),(ES12.2))))')years(i),a(i),b(i),c(i)
        END DO
        CLOSE(24)
		
		WRITE(*,*)"******* The wells coeff file "//TRIM(ADJUSTL(RT))//" has been read succesfully *******"
		
    END SUBROUTINE READ_OW_COEFF
	
    !Function to read injection wells
    SUBROUTINE READ_IW_FILE(iw,RT,SDW)

        USE TYPES

        IMPLICIT NONE
        !Global variables
        INTEGER,INTENT(IN)                              :: iw
        CHARACTER(200),INTENT(IN)                       :: RT
        TYPE(SALT_WATER_DISPOSAL),INTENT(INOUT)             :: SDW
        !Local variables
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(12)                                   ::intfile
        
        OPEN(24, FILE=RT, STATUS='OLD',ACTION='READ')
        !Lectura de metadatos
        READ(24, '((A14),(I10))')aux
        !Se crea internal file para descriptor de formato
        WRITE(intfile,'(I12)')4
        !Lectura de matriz
        DO i=1,iw
            READ(24,'('//intfile//'((A12),(1X),(ES12.2),(1X),(ES12.2),(1X),(ES12.2))))')SDW%ID(i),SDW%TD(i),SDW%X(i),SDW%Y(i)
        END DO
        CLOSE(24)
		
		WRITE(*,*)"******* The Injection Wells file has been read succesfully *******"
		
    END SUBROUTINE READ_IW_FILE

    !Function to read oil wells
    SUBROUTINE READ_OW_FILE(RT,Oilw,ow)

        USE TYPES

        IMPLICIT NONE
        !Global variables
        INTEGER,INTENT(IN)                              :: ow
        CHARACTER(200),INTENT(IN)                       :: RT
        TYPE(OG_WELLS),INTENT(INOUT)                   :: Oilw
        !Local variables
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(12)                                   ::intfile,a
		REAL											::b,c
		INTEGER::d
		CHARACTER(1)::e
        
        OPEN(24, FILE=RT, STATUS='OLD',ACTION='READ')
        !Lectura de header
        READ(24, '((A14),(I10))')aux
        !Se crea internal file para descriptor de formato
        WRITE(intfile,'(I12)')3
        !Lectura de matriz
		
        DO i=1,ow
            READ(24,'('//intfile//'((A12),(1X),(ES12.2),(1X),(ES12.2),(1X),(I4),(1X),(1A))))')&
			&Oilw%ID(i),Oilw%X(i),Oilw%Y(i),Oilw%In_Year(i),Oilw%Dir(i)
        END DO
        CLOSE(24)
		
		WRITE(*,*)"******* The wells file "//TRIM(ADJUSTL(RT))//" has been read succesfully *******"
		
    END SUBROUTINE READ_OW_FILE
	
	!This routine writes an INTEGER raster layer in ".asc" format
    SUBROUTINE WRITE_INT_ASC(Ncols,Nrows,Xcor,Ycor,delta,nan,A,RT)
    
        IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(IN)                              :: Ncols,Nrows,nan
        REAL,INTENT(IN)                                 :: Xcor,Ycor,delta
        INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)   :: A
        CHARACTER(200),INTENT(IN)                       :: RT
        !Local variables section
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(10)                                   ::intfile
        
        OPEN(25, FILE=RT, STATUS='REPLACE')
        !Writing values
        WRITE(25, '((A14),(I10))')'ncols         ',Ncols
        WRITE(25, '((A14),(I10))')'nrows         ',Nrows
        WRITE(25, '((A14),(F10.3))')'xllcorner     ',Xcor
        WRITE(25, '((A14),(F11.3))')'yllcorner     ',Ycor
        WRITE(25, '((A14),(F10.2))')'cellsize      ',delta
        WRITE(25, '((A14),(I10))')'NODATA_Value  ',nan
        !An internal file its created
        WRITE(intfile,'(I1)')Ncols
    
        DO i=1,Nrows
            WRITE(25,'('//intfile//'((I1),(1X))))')(A(i,j),j=1,Ncols)
        END DO
        CLOSE(25)

        WRITE(*,*)"******* The file "//TRIM(ADJUSTL(RT))//" has been written successfully *******"
        
    END SUBROUTINE WRITE_INT_ASC

        !This routine writes a REAL raster layer in ".asc" format
    SUBROUTINE WRITE_REAL_ASC(Ncols,Nrows,Xcor,Ycor,delta,nan,A,RT)
    
        IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(IN)                              :: Ncols,Nrows,nan
        REAL,INTENT(IN)                                 :: Xcor,Ycor,delta
        REAL,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)      :: A
        CHARACTER(200),INTENT(IN)                       :: RT
        !Local variables section
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(10)                                   ::intfile
        
        OPEN(25, FILE=RT, STATUS='REPLACE')
        !Writing values
        WRITE(25, '((A14),(I10))')'ncols         ',Ncols
        WRITE(25, '((A14),(I10))')'nrows         ',Nrows
        WRITE(25, '((A14),(F10.3))')'xllcorner     ',Xcor
        WRITE(25, '((A14),(F10.3))')'yllcorner     ',Ycor
        WRITE(25, '((A14),(F10.2))')'cellsize      ',delta
        WRITE(25, '((A14),(I10))')'NODATA_Value  ',nan
        !An internal file its created
        WRITE(intfile,'(I10)')Ncols
    
        DO i=1,Nrows
            WRITE(25,'('//intfile//'((F10.3),(1X))))')(A(i,j),j=1,Ncols)
        END DO
        CLOSE(25)

        WRITE(*,*)"******* The file "//TRIM(ADJUSTL(RT))//" has been written successfully *******"
        
    END SUBROUTINE WRITE_REAL_ASC

            !This routine writes a REAL raster layer in ".asc" format
    SUBROUTINE WRITE_CHAR_MATRIX(t,ow,A,X,RT)
    
        IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(IN)                              :: t,ow
        CHARACTER (LEN=12),ALLOCATABLE,INTENT(IN)       :: X(:)
        CHARACTER (LEN=12),ALLOCATABLE,INTENT(IN)       :: A(:,:)
        CHARACTER(200),INTENT(IN)                       :: RT
        !Local variables section
        INTEGER                                         ::i,j
        CHARACTER(15)                                   ::aux
        CHARACTER(12)                                   ::intfile
        
        OPEN(25, FILE=RT, STATUS='REPLACE')
        !An internal file its created
        
        WRITE(intfile,'(I12)')ow
        
        WRITE(25,'('//intfile//'((A12),(A1))))')(X(j),",",j=1,ow)

        DO i=1,t
            WRITE(25,'('//intfile//'((A12),(A1))))')(A(i,j),",",j=1,ow)
        END DO
        CLOSE(25)

        WRITE(*,*)"******* The file "//TRIM(ADJUSTL(RT))//" has been written successfully *******"
        
    END SUBROUTINE WRITE_CHAR_MATRIX

    SUBROUTINE WRITE_TXT_RESULTS(Nwells,t,A,X,RT)

        USE TYPES
    
        IMPLICIT NONE
        !Global variables section
        INTEGER,INTENT(IN)                              :: Nwells,t
        CHARACTER (LEN=12),ALLOCATABLE,INTENT(IN)       :: X(:)
        REAL,DIMENSION(:,:),ALLOCATABLE,INTENT(IN)      :: A
        CHARACTER(200),INTENT(IN)                       :: RT
        !Local variables section
        INTEGER                                         ::i,j
        CHARACTER(12)                                   ::intfile
        
        OPEN(25, FILE=RT, STATUS='REPLACE')

        !An internal file its created
        WRITE(intfile,'(I12)')Nwells

        WRITE(25,'('//intfile//'((A12),(A1))))')(X(j),",",j=1,Nwells)
        
        DO i=0,t
            WRITE(25,'('//intfile//'((ES12.5),(A1))))')(A(i,j),",",j=1,Nwells)
        END DO
        CLOSE(25)

        WRITE(*,*)"******* The file "//TRIM(ADJUSTL(RT))//" has been written successfully *******"
        
    END SUBROUTINE WRITE_TXT_RESULTS

END MODULE