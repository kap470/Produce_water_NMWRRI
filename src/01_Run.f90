!Author: Kevin Perez
!Email: kap470@nmsu.edu
!New Mexico Water Resources Research Institute 

PROGRAM main

    USE TYPES
	USE FUNCTIONS
    USE IO

    IMPLICIT NONE

    INTEGER                                        :: Ncols,Nrows,nan,ow,gw,t,i,j,k,iw,NTP,aux,m,ny,acm,n,tf,owf,gwf,CountE
    REAL                                           :: Xcor,Ycor,delta,xi,xf,yi,yf,tini,tend,dt,RTPT,RTPO,TPS0
    REAL                                           :: cost,OWR,PW,OWS0,INJS0,TPCAP,Transportation_Cap,Min_cost
    CHARACTER(200)                                 :: RT
    TYPE(OG_WELLS)                                 :: Oilw,OilF,Gasw,GasF
    TYPE(SALT_WATER_DISPOSAL)                      :: SDW
    TYPE(TREATMENT_PLANT)                          :: T_plant,T_plantF
    REAL,DIMENSION(:,:),ALLOCATABLE                :: M_dist,M_cost,M_dist_2,M_cost_2,AI,AO,AT,AG,ATV
	REAL,DIMENSION(:),ALLOCATABLE                  :: Toil,Tgas,Tpw,ToilF,TgasF,TpwF,array,years,CountOW,CountGW,Excess,ExcessF
	REAL,DIMENSION(:),ALLOCATABLE                  :: aOD,bOD,cOD
	REAL,DIMENSION(:),ALLOCATABLE                  :: aOV,bOV,cOV
	REAL,DIMENSION(:),ALLOCATABLE                  :: aGV,bGV,cGV
	REAL,DIMENSION(:),ALLOCATABLE                  :: aOH,bOH,cOH
	REAL,DIMENSION(:),ALLOCATABLE                  :: aGH,bGH,cGH
    CHARACTER (LEN=12),ALLOCATABLE                 :: XIJ(:),XO(:),XT(:),XG(:),XTV(:)
	INTEGER,DIMENSION(:),ALLOCATABLE			   :: array1,TtimeOW,TtimeGW
    
    WRITE(*,*)"******* The simulation has started *******"


	!Read input values from input file text
    RT="./Input_file.txt"
    CALL read_input_file(ow,gw,iw,tini,tend,dt,Ncols,Nrows,Xcor,Ycor,delta,nan,RT,cost,OWS0&
	&,INJS0,NTP,TPCAP,RTPT,RTPO,TPS0)
    t=int((tend-tini)/dt)
	ny=t/12
	
    !Allocating Salt Disposal wells arrays
	!The prefix iw is being used for the old name of the SDW

    ALLOCATE(SDW%X(iw))
    ALLOCATE(SDW%Y(iw))
    ALLOCATE(SDW%TD(iw))
	ALLOCATE(SDW%TV(iw))
    ALLOCATE(SDW%ID(iw))
    ALLOCATE(SDW%RC(t+1,iw))
    ALLOCATE(SDW%INF(t+1,iw))
    ALLOCATE(SDW%OUTF(t+1,iw))
    ALLOCATE(SDW%STOCK(t+1,iw))    
    ALLOCATE(SDW%STATE(iw))   
    
    !Allocating oil wells arrays

    ALLOCATE(Oilw%In_Year(ow))
	ALLOCATE(Oilw%Dir(ow))
	ALLOCATE(Oilw%X(ow))
    ALLOCATE(Oilw%Y(ow))
    ALLOCATE(Oilw%ID(ow))
    ALLOCATE(Oilw%PW(t+1,ow))
    ALLOCATE(Oilw%OUTF(t+1,ow))
    ALLOCATE(Oilw%STOCK(t+1,ow))
	
	!Allocating gas wells arrays

    ALLOCATE(Gasw%In_Year(gw))
	ALLOCATE(Gasw%Dir(gw))
	ALLOCATE(Gasw%X(gw))
    ALLOCATE(Gasw%Y(gw))
    ALLOCATE(Gasw%ID(gw))
    ALLOCATE(Gasw%PW(t+1,gw))
    ALLOCATE(Gasw%OUTF(t+1,gw))
    ALLOCATE(Gasw%STOCK(t+1,gw))
	
    !Allocating Treatment Plant arrays
    ALLOCATE(T_plant%X(NTP))
    ALLOCATE(T_plant%Y(NTP))
    ALLOCATE(T_plant%ID(NTP))
    ALLOCATE(T_plant%INF(t+1,NTP))
    ALLOCATE(T_plant%OUTF(t+1,NTP))
    ALLOCATE(T_plant%STOCK(t+1,NTP))
	ALLOCATE(T_plant%RC(t+1,NTP))

    !Allocating other arrays
    ALLOCATE(M_dist(ow,iw))
	ALLOCATE(M_dist_2(ow,NTP))
    ALLOCATE(M_cost(ow,iw))
	ALLOCATE(M_cost_2(ow,NTP))
	ALLOCATE(Toil(t),Tgas(t),Tpw(t)) !Acumulated total volumens for each time step in historical part
	ALLOCATE(TtimeOW(t),TtimeGW(t)) !Total number of oil and gas wells for each time step (integer)
	ALLOCATE(years(ny))
	ALLOCATE(CountOW(ow),CountGW(gw)) !Total number of oil and gas wells for each time step (real)
	ALLOCATE(Excess(t),ExcessF(t))
	
	!Allocating exponential functions coefficients 
	
	ALLOCATE(aOD(ny),bOD(ny),cOD(ny))
	ALLOCATE(aOV(ny),bOV(ny),cOV(ny))
	ALLOCATE(aGV(ny),bGV(ny),cGV(ny))
	ALLOCATE(aOH(ny),bOH(ny),cOH(ny))
	ALLOCATE(aGH(ny),bGH(ny),cGH(ny))
	
	!Reading files section

	RT="../Input/SALT_WATER_DISPOSAL.csv"
    CALL READ_IW_FILE(iw,RT,SDW)
	RT="../Input/Oil_wells.csv"
    CALL READ_OW_FILE(RT,Oilw,ow)
	RT="../Input/Gas_wells.csv"
    CALL READ_OW_FILE(RT,Gasw,gw)
	RT="../Input/O-D.csv"
	CALL READ_OW_COEFF(ny,RT,years,aOD,bOD,cOD)
	RT="../Input/O-V.csv"
	CALL READ_OW_COEFF(ny,RT,years,aOV,bOV,cOV)
	RT="../Input/G-V.csv"
	CALL READ_OW_COEFF(ny,RT,years,aGV,bGV,cGV)
	RT="../Input/O-H.csv"
	CALL READ_OW_COEFF(ny,RT,years,aOH,bOH,cOH)
	RT="../Input/G-H.csv"
	CALL READ_OW_COEFF(ny,RT,years,aGH,bGH,cGH)
   
    !Treatment plant
	!Location created as the average of the x and y coordinates of the oil wells
    T_plant%X(1)=Oilw%X(1)
    T_plant%y(1)=Oilw%Y(1)
    DO i=2,ow
        T_plant%X(1)=T_plant%X(1)+Oilw%X(i)
        T_plant%y(1)=T_plant%y(1)+Oilw%Y(i)
    END DO
    T_plant%X(1)=T_plant%X(1)/ow
    T_plant%Y(1)=T_plant%Y(1)/ow

    !Calculating matrices of distance OW-IJW
	!This matrix is created to take decisions of PW disposal based in the less cost or distance
    DO i=1,ow
        DO j=1,iw
            M_dist(i,j)=SQRT(((Oilw%X(i)-SDW%X(j))**2+(Oilw%Y(i)-SDW%Y(j))**2))
            M_cost(i,j)=M_dist(i,j)*cost
        END DO
    END DO
	
	!*********************************************
	!********* CREATING FUTURE OIL WELLS *********
	!*********************************************
	
	!This commented block creates a set of future random located wells
	!the "array variable" has the number of new wells for each year of the future period
	
	!Array with future number wells
	! tf=300/12
	! ALLOCATE(array(tf))
	! ALLOCATE(array1(tf))
	! CALL CREATE_RANDOM_NUMBERS(tf,16.0,27.0,array)
	! array1=int(array)
	! owf=sum(array1)
	! WRITE(*,*)"total new oil wells",owf
	
	!Allocating future oil wells arrays
	owf=2678
	tf=t
    ALLOCATE(OilF%In_Year(owf))
	ALLOCATE(OilF%Dir(owf))
	ALLOCATE(OilF%X(owf))
    ALLOCATE(OilF%Y(owf))
    ALLOCATE(OilF%ID(owf))
    ALLOCATE(OilF%PW(tf+1,owf))
    ALLOCATE(OilF%OUTF(tf+1,owf))
    ALLOCATE(OilF%STOCK(tf+1,owf))
	RT="../Input/Future_Oil_wells.csv"
	CALL READ_OW_FILE(RT,OilF,owf)
	
	gwf=805
	ALLOCATE(GasF%In_Year(gwf))
	ALLOCATE(GasF%Dir(gwf))
	ALLOCATE(GasF%X(gwf))
    ALLOCATE(GasF%Y(gwf))
    ALLOCATE(GasF%ID(gwf))
    ALLOCATE(GasF%PW(tf+1,gwf))
    ALLOCATE(GasF%OUTF(tf+1,gwf))
    ALLOCATE(GasF%STOCK(tf+1,gwf))
	RT="../Input/Future_Gas_wells.csv"
	CALL READ_OW_FILE(RT,GasF,gwf)
	
	ALLOCATE(ToilF(tf),TgasF(tf),TpwF(tf))
	
	ALLOCATE(T_plantF%X(NTP))
    ALLOCATE(T_plantF%Y(NTP))
    ALLOCATE(T_plantF%ID(NTP))
    ALLOCATE(T_plantF%INF(t+1,NTP))
    ALLOCATE(T_plantF%OUTF(t+1,NTP))
    ALLOCATE(T_plantF%STOCK(t+1,NTP))
	ALLOCATE(T_plantF%RC(t+1,NTP))
	
	! CALL CREATE_RANDOM_LOCATION(Ncols,Nrows,Xcor,Ycor,owf,delta,OilF%X,OilF%Y)
	
	!*********************************************
	!********* INITIAL VALUES SD SIMULATION ******
	!*********************************************
	
	!Injection wells
    SDW%STOCK(0,:)=INJS0
    SDW%STATE(:)=.TRUE.
	SDW%TV(:)=SDW%TD(:)*100
	
	!OIL WELL 
	Oilw%STOCK(0,:)=OWS0
	Toil(:)=0
	!TREATMET PLANT 
	T_plant%TC=TPCAP
	Transportation_Cap=TPCAP*RTPT
	T_plant%STOCK(0,1)=TPS0
	!Total values arrays
	Toil(:)=0
	Tgas(:)=0
	Tpw(:)=0
	TtimeOW(:)=0
	TtimeGW(:)=0
	CountOW(:)=0
	CountGW(:)=0
	!Years counter
	m=1
	acm=0
	CountE=0
    DO k=0,t-1
		!This block accounts the years in the monthly iteration
		IF (acm<12) THEN
			acm=acm+1
		ELSE
			acm=1
			m=m+1
		END IF
		!OIL WATER RATIO
		OWR=0.229*m**(0.7088)
		!OIL WELLS SECTION
		!CountOW(:)=counter that starts in the year n, allows proper evaluation of the exponential functions
		DO i=1,ow
			IF ((k/12+1995).GE.(Oilw%In_Year(i))) THEN
				DO j=1,ny
					IF (Oilw%In_Year(i).EQ.years(j)) THEN
						n=j
					END IF
				END DO
				TtimeOW(k)=TtimeOW(k)+1
				CountOW(i)=CountOW(i)+1
				IF (Oilw%Dir(i).EQ."H") THEN
					Oilw%OUTF(k,i)=aOH(n)*exp(-bOH(n)*(CountOW(i)))+cOH(n)
				ELSE IF (Oilw%Dir(i).EQ."V") THEN
					Oilw%OUTF(k,i)=aOV(n)*exp(-bOV(n)*(CountOW(i)))+cOV(n)
				ELSE IF (Oilw%Dir(i).EQ."D") THEN
					Oilw%OUTF(k,i)=aOD(n)*exp(-bOD(n)*(CountOW(i)))+cOD(n)
				END IF
				IF (Oilw%OUTF(k,i)<0) THEN
					Oilw%OUTF(k,i)=0
				END IF
			ELSE 
				Oilw%OUTF(k,i)=0
			END IF
			Oilw%PW(k,i)=Oilw%OUTF(k,i)*OWR
			Oilw%STOCK(k+1,i)=Oilw%STOCK(k,i)-Oilw%OUTF(k,i)*dt
			Toil(k)=Oilw%OUTF(k,i)+Toil(k)
			Tpw(k)=Oilw%PW(k,i)+Tpw(k)
		END DO
		!GAS WELLS SECTION
		DO i=1,gw
			IF ((k/12+1995).GE.(Gasw%In_Year(i))) THEN
				DO j=1,ny
					IF (Gasw%In_Year(i).EQ.years(j)) THEN
						TtimeGW(k)=TtimeGW(k)+1
						n=j
						CountGW(i)=CountGW(i)+1
					END IF
				END DO
				IF (Gasw%Dir(i).EQ."H") THEN
					Gasw%OUTF(k,i)=aGH(n)*exp(-bGH(n)*(CountGW(i)))+cGH(n)
				ELSE IF (Gasw%Dir(i).EQ."V") THEN
					Gasw%OUTF(k,i)=aGV(n)*exp(-bGV(n)*(CountGW(i)))+cGV(n)
				END IF
				IF (Gasw%OUTF(k,i)<0) THEN
					Gasw%OUTF(k,i)=0
				END IF
			ELSE 
				Gasw%OUTF(k,i)=0
			END IF
			Gasw%PW(k,i)=Gasw%OUTF(k,i)*OWR
			Gasw%STOCK(k+1,i)=Gasw%STOCK(k,i)-Gasw%OUTF(k,i)*dt
			Tgas(k)=Gasw%OUTF(k,i)+Tgas(k)
		END DO
		
		IF ((m+1995).LE.(2011)) THEN
			Excess(k)=0
		ELSE
			CountE=CountE+1
			Excess(k)=((-0.38*CountE+85)/100)*Tpw(k)
		END IF
		
		!Uncomment this block to take decisions of PW disposal between SWD and TP
		!Calculating matrices of distance OW-TP
		! DO i=1,ow
			! M_dist_2(i,1)=SQRT(((Oilw%X(i)-T_plant%X(1))**2+(Oilw%Y(i)-T_plant%Y(1))**2))
			! M_cost_2(i,1)=M_dist_2(i,1)*cost
		! END DO
		! i=1
		! T_plant%INF(k,1)=0
		! DO WHILE (i<=ow)
			! Min_cost=MINVAL(M_cost_2(:,1))
			! IF ((Min_cost .eq. M_cost_2(i,1)).AND.(Transportation_Cap>=Oilw%PW(k,i)+T_plant%INF(k,1)))THEN
				! T_plant%INF(k,1)=Oilw%PW(k,i)+T_plant%INF(k,1)
				! M_cost_2(i,1)=10**6
				! i=1
			! END IF
			! i=i+1
		! END DO
		T_plant%INF(k,1)=Excess(k)
		T_plant%OUTF(k,1)=T_plant%STOCK(k,1)*RTPO
		T_plant%STOCK(k+1,1)=T_plant%STOCK(k,1)+T_plant%INF(k,1)*dt-T_plant%OUTF(k,1)*dt-T_plant%STOCK(k,1)*(1-RTPO)
		
		!uncomment this block to make the simulation of SDW 
		! !cost=minimum cost in each row
		! SDW%INF(k,:)=0

		! DO i=1,iw
			! SDW%RC(k,i)=SDW%TV(i)-SDW%STOCK(k,i)
			! !VERIFICAR MAS ADELANTE
			! IF (SDW%RC(k,i)<=0) THEN
				! SDW%STATE(i)=.FALSE.
				! M_cost(:,i)=10**6
			! ELSE
				! SDW%STATE(i)=.TRUE.
			! END IF
		! END DO
		! DO i=1,ow
			! IF (M_cost_2(i,1)/=10**6) THEN		
				! Min_cost=MINVAL(M_cost(i,:))
				! DO j=1,iw
					! IF (Min_cost .eq. M_cost(i,j)) THEN
						! aux=j
					! END IF
				! END DO
				! SDW%INF(k,aux)=SDW%INF(k,aux)+Oilw%PW(k,i)
			! END IF
		! END DO
		! DO i=1,iw
			! IF (SDW%STATE(i) .eqv. .FALSE.) THEN
				! SDW%INF(k,i)=0
				! SDW%STOCK(k+1,i)=SDW%STOCK(k,i)
			! ELSE
				! SDW%INF(k,i)=SDW%INF(k,i)
				! SDW%STOCK(k+1,i)=SDW%STOCK(k,i)+SDW%INF(k,i)*dt
			! END IF
		! END DO
	
    END DO
	
	!*****************************************************************
	!************************FUTURE SIMULATION************************
	!*****************************************************************
	acm=0
	n=25
	ToilF(:)=0
	TgasF(:)=0
	TpwF(:)=0
	TtimeOW(:)=0
	TtimeGW(:)=0
	CountOW(:)=0
	CountGW(:)=0
	T_plantF%STOCK(1,1)=T_plant%STOCK(t,1)
	DO k=0,tf-1
		!yearly loop
		IF (acm<12) THEN
			acm=acm+1
		ELSE
			acm=1
			m=m+1
		END IF
		!OIL WATER RATIO
		OWR=0.229*m**(0.7088)
		!OIL WELLS SECTION
		!CountOW(:)=counter that starts in the year n, allows proper evaluation of the exponential functions
		DO i=1,ow
			IF ((k/12+2021).GE.(OilF%In_Year(i))) THEN
				DO j=1,ny
					IF (OilF%In_Year(i).EQ.years(j)) THEN
						n=j
					END IF
				END DO
				TtimeOW(k)=TtimeOW(k)+1
				CountOW(i)=CountOW(i)+1
				IF (OilF%Dir(i).EQ."H") THEN
					OilF%OUTF(k,i)=aOH(n)*exp(-bOH(n)*(CountOW(i)))+cOH(n)
				ELSE IF (OilF%Dir(i).EQ."V") THEN
					OilF%OUTF(k,i)=aOV(n)*exp(-bOV(n)*(CountOW(i)))+cOV(n)
				ELSE IF (OilF%Dir(i).EQ."D") THEN
					OilF%OUTF(k,i)=aOD(n)*exp(-bOD(n)*(CountOW(i)))+cOD(n)
				END IF
				IF (OilF%OUTF(k,i)<0) THEN
					OilF%OUTF(k,i)=0
				END IF
			ELSE 
				OilF%OUTF(k,i)=0
			END IF
			OilF%PW(k,i)=OilF%OUTF(k,i)*OWR
			OilF%STOCK(k+1,i)=OilF%STOCK(k,i)-OilF%OUTF(k,i)*dt
			ToilF(k)=OilF%OUTF(k,i)+ToilF(k)
			TpwF(k)=OilF%PW(k,i)+TpwF(k)
		END DO
		!GAS WELLS SECTION
		DO i=1,gw
			IF ((k/12+2021).GE.(GasF%In_Year(i))) THEN
				DO j=1,ny
					IF (GasF%In_Year(i).EQ.years(j)) THEN
						n=j
					END IF
				END DO
				TtimeGW(k)=TtimeGW(k)+1
				CountGW(i)=CountGW(i)+1
				IF (GasF%Dir(i).EQ."H") THEN
					GasF%OUTF(k,i)=aGH(n)*exp(-bGH(n)*(CountGW(i)))+cGH(n)
				ELSE IF (GasF%Dir(i).EQ."V") THEN
					GasF%OUTF(k,i)=aGV(n)*exp(-bGV(n)*(CountGW(i)))+cGV(n)
				END IF
				IF (GasF%OUTF(k,i)<0) THEN
					GasF%OUTF(k,i)=0
				END IF
			ELSE 
				GasF%OUTF(k,i)=0
			END IF
			GasF%PW(k,i)=GasF%OUTF(k,i)*OWR
			GasF%STOCK(k+1,i)=GasF%STOCK(k,i)-GasF%OUTF(k,i)*dt
			TgasF(k)=GasF%OUTF(k,i)+TgasF(k)
		END DO
		CountE=CountE+1
		ExcessF(k)=((-0.38*CountE+85)/100)*TpwF(k)
		IF (ExcessF(k).LE.0) THEN
			ExcessF(k)=0
		END IF
		T_plantF%INF(k,1)=ExcessF(k)
		T_plantF%OUTF(k,1)=T_plantF%STOCK(k,1)*RTPO
		T_plantF%STOCK(k+1,1)=T_plantF%STOCK(k,1)+T_plantF%INF(k,1)*dt-T_plantF%OUTF(k,1)*dt-T_plantF%STOCK(k,1)*(1-RTPO)
		WRITE(*,*)TtimeOW(k),TtimeGW(k),k
	END DO
	
    !WRITING RESULTS SECTION

    ! ALLOCATE(XIJ(iw))
    ! ALLOCATE(AI(t+1,iw))
    ! !Capacity file
    ! RT="../Output/IJW_CAP.csv"
    ! DO i=1,iw
        ! XIJ(i)=SDW%ID(i)
        ! DO j=0,t
            ! AI(j,i)=SDW%STOCK(j,i)
        ! END DO
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(iw,t,AI,XIJ,RT)
    ! !Remaining capacity
    ! RT="../Output/IJW_RC.csv"
    ! DO i=1,iw
        ! XIJ(i)=SDW%ID(i)
        ! DO j=0,t
            ! AI(j,i)=SDW%RC(j,i)
        ! END DO
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(iw,t,AI,XIJ,RT)
    ! !Inflow
    ! RT="../Output/IJW_INF.csv"
    ! DO i=1,iw
        ! XIJ(i)=SDW%ID(i)
        ! DO j=0,t
            ! AI(j,i)=SDW%INF(j,i)
        ! END DO
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(iw,t,AI,XIJ,RT)
    
    ALLOCATE(XO(ow))
    ALLOCATE(AO(t+1,ow))

    !OIL WELLS
    RT="../Output/OW_STOCK.csv"
    DO i=1,ow
		XO(i)=Oilw%ID(i)
        DO j=0,t
            AO(j,i)=Oilw%STOCK(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(ow,t,AO,XO,RT)
    !Outflow
    RT="../Output/OW_OUTF.csv"
    DO i=1,ow
		XO(i)=Oilw%ID(i)
        DO j=0,t
            AO(j,i)=Oilw%OUTF(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(ow,t,AO,XO,RT)
    !Produce water
    RT="../Output/OW_PW.csv"
    DO i=1,ow
		XO(i)=Oilw%ID(i)
        DO j=0,t
            AO(j,i)=Oilw%PW(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(ow,t,AO,XO,RT)
	
	!GAS WELLS
	
	ALLOCATE(XG(gw))
    ALLOCATE(AG(t+1,gw))
    RT="../Output/GW_STOCK.csv"
    DO i=1,gw
		XG(i)=Gasw%ID(i)
        DO j=0,t
            AG(j,i)=Gasw%STOCK(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(gw,t,AG,XG,RT)
    !Outflow
    RT="../Output/GW_OUTF.csv"
    DO i=1,gw
		XG(i)=Gasw%ID(i)
        DO j=0,t
            AG(j,i)=Gasw%OUTF(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(gw,t,AG,XG,RT)
    !Produce water
    RT="../Output/GW_PW.csv"
    DO i=1,gw
		XG(i)=Gasw%ID(i)
        DO j=0,t
            AG(j,i)=Gasw%PW(j,i)
        END DO
    END DO
    CALL  WRITE_TXT_RESULTS(gw,t,AG,XG,RT)
	
	!TOTAL values
	ALLOCATE(XTV(10),ATV(t,10))
	XTV(1)="H_Total_oil "
	XTV(2)="H_Total_gas "
	XTV(3)="H_Total_PW  "
	XTV(4)="H_Excess"
	XTV(5)="H_OutFTPlant"
	XTV(6)="F_Total_Oil "
	XTV(7)="F_Total_gas "
	XTV(8)="F_Total_PW  "
	XTV(9)="F_Excess   "
	XTV(10)="F_OutFTPlant"
	RT="../Output/Total_values.csv"
	DO i=1,t
		ATV(i,1)=Toil(i)
		ATV(i,2)=Tgas(i)
		ATV(i,3)=Tpw(i)
		ATV(i,4)=Excess(i)
		ATV(i,5)=T_plant%OUTF(i,1)
		ATV(i,6)=ToilF(i)
		ATV(i,7)=TgasF(i)
		ATV(i,8)=TpwF(i)
		ATV(i,9)=ExcessF(i)
		ATV(i,10)=T_plantF%OUTF(i,1)
    END DO
	CALL  WRITE_TXT_RESULTS(10,t,ATV,XTV,RT)
	
	!Treatment plat section
	! ALLOCATE(XT(NTP))
	! ALLOCATE(AT(t+1,NTP))
	
	! XT(1)="Treatment_Pl"
	
    ! RT="../Output/TM_STOCK.csv"
    ! DO j=0,t
        ! AT(j,1)=T_plant%STOCK(j,1)
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(1,t,AT,XT,RT)
	
    ! !Outflow
    ! RT="../Output/TM_OUTF.csv"
    ! DO j=0,t
        ! AT(j,1)=T_plant%OUTF(j,1)
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(1,t,AT,XT,RT)
	
    ! !Remaining capacity
    ! RT="../Output/TM_RC.csv"
    ! DO j=0,t
        ! AT(j,1)=T_plant%RC(j,1)
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(1,t,AT,XT,RT)
	
	! !Remaining capacity
    ! RT="../Output/TM_INF.csv"
    ! DO j=0,t
        ! AT(j,1)=T_plant%INF(j,1)
    ! END DO
    ! CALL  WRITE_TXT_RESULTS(1,t,AT,XT,RT)

    WRITE(*,*)"******* The simulation has ended *******"

END PROGRAM