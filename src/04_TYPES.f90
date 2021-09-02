!Author: Kevin Perez
!Email: kap470@nmsu.edu
!New Mexico Water Resources Research Institute 

MODULE TYPES
	!OG=Oil and Gas
    !OPPD=Oil production per day
    TYPE OG_WELLS
		INTEGER, DIMENSION(:),ALLOCATABLE :: In_Year
		CHARACTER (LEN=1),ALLOCATABLE :: Dir(:)
        REAL, DIMENSION(:),ALLOCATABLE :: X,Y
        CHARACTER (LEN=12),ALLOCATABLE :: ID(:)
        REAL, DIMENSION(:,:),ALLOCATABLE :: OUTF,STOCK,PW
    END TYPE

    !RC=remaining capacity
    !TD=Total depth
	!TV=Total volume
    !STATE is a matrix with oil wells in rows and IJ in columns
    TYPE SALT_WATER_DISPOSAL
        REAL, DIMENSION(:),ALLOCATABLE :: X,Y,TD,TV
        CHARACTER (LEN=12),ALLOCATABLE :: ID(:)
        REAL, DIMENSION(:,:),ALLOCATABLE :: INF,OUTF,STOCK,RC
        LOGICAL, DIMENSION(:),ALLOCATABLE :: STATE
    END TYPE    
    
	!TC=total capacity
    TYPE TREATMENT_PLANT
		REAL::TC
        REAL, DIMENSION(:),ALLOCATABLE :: X,Y
        CHARACTER (LEN=12),ALLOCATABLE :: ID(:)
        REAL, DIMENSION(:,:),ALLOCATABLE :: INF,OUTF,STOCK,RC
    END TYPE    
	
END MODULE

