!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_IO
!! Contenent: Collections of Subroutines related to I/O and 
!!            Initialization
!! Author: PY Chuang
!!#####################################################################
MODULE IO
USE VAR
IMPLICIT NONE
 CONTAINS
!======================================================================


!======================================================================
!----------------------------------------------------------------------
! Initialize
!----------------------------------------------------------------------
SUBROUTINE initialize
IMPLICIT NONE


    OPEN(UNIT = 120, FILE = InputFileName)
    READ(UNIT = 120, NML = Initial_1)
    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( phn(FNph) )
    ALLOCATE( PoolL(Ne(2), Ne(3), PoolSize) )
    ALLOCATE( PoolR(Ne(2), Ne(3), PoolSize) )
    ALLOCATE( iNPoolL(Ne(2), Ne(3)) )
    ALLOCATE( iNPoolR(Ne(2), Ne(3)) )
    READ(UNIT = 120, NML = Initial_2)
    CLOSE( 120 )
    
    iMid = Ne(1) / 2

    ALLOCATE( NAdd(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( Tavg(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( Eavg(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( phId(FNph) )
    ALLOCATE( EmptyID(FNph) )
    NAdd = 0
    Tavg = 0
    Eavg = 0
    phId = -1
    EmptyID = -1
    
    ALLOCATE( EinjectL(Ne(1), Ne(2)) )
    ALLOCATE( EinjectR(Ne(1), Ne(2)) )
    ALLOCATE( qL(Ne(1), Ne(2)) )
    ALLOCATE( qC(Ne(1), Ne(2)) )
    ALLOCATE( qR(Ne(1), Ne(2)) )
    
    SELECTCASE( BCs(1) )
    CASE(3)
        WRITE(*, *) "This BC is not supported yet."
        WRITE(*, *) "Program will be shut down in 5 seconds."
        CALL SLEEP(5)
        STOP
    CASE DEFAULT
        EinjectL = 0
        EinjectR = 0
        qL = 0
        qC = 0
        qR = 0
        TBCL = 0
        TBCR = 0
        qfluxL = 0
        qfluxC = 0 
        qfluxR = 0
        VBCL = 0
        VBCR = 0
        EBCL = 0
        EBCR = 0
    END SELECT
    
    
    WRITE(*, *) "Case name: ", CaseName
    WRITE(*, *) "Domain size: ", L
    WRITE(*, *) "Numbers of elements: ", Ne
    WRITE(*, *) "Total number of phonons: ", RNph
    WRITE(*, *) "The size of phonon array: ", FNph
    WRITE(*, *) "Bundle: ", bundle
    WRITE(*, *) "dt: ", dt

END SUBROUTINE initialize

!======================================================================
END MODULE IO