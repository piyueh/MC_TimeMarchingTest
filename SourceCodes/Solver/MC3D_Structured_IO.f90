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

    iter = iter0
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
    WRITE(*, *) "TimeStep: ", TimeStep

END SUBROUTINE initialize


!======================================================================
!----------------------------------------------------------------------
! Output Transient Temperature
!----------------------------------------------------------------------
SUBROUTINE Output_Trans
IMPLICIt NONE

    WRITE(OutputFileName, "(I8.8, '_Transient.txt')") iter
    OPEN( UNIT = 150, FILE = OutputFileName )
    WRITE(150, *) time
    WRITE(150, *) ele(:, 1, 1)%T
    CLOSE( 150 )

END SUBROUTINE Output_Trans


!======================================================================
!----------------------------------------------------------------------
! Regular Output
!----------------------------------------------------------------------
SUBROUTINE Output
IMPLICIt NONE
REAL(KIND=8):: R1
INTEGER(KIND=4), SAVE:: ct = 0
INTEGER(KIND=4):: i, j, k, mt

    Eavg = Eavg + ele%E
    ct = ct + 1

    IF ( MOD( iter, nOutput ).eq.0 ) THEN
        Eavg = Eavg / DBLE( ct )

        DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)
            R1 = Eavg(i, j, k) / dV
            mt = ele(i, j, k)%Mat
            CALL Etable( mt, 1, R1, Tavg(i, j, k) )
        ENDDO; ENDDO; ENDDO

        WRITE(OutputFileName, "(I8.8, '.txt')") iter
        OPEN( UNIT = 150, FILE = OutputFileName )
        WRITE(150, *) ct
        WRITE(150, NML = DataOutput)
        CLOSE( 150 )

        Eavg = 0
        ct = 0

        CALL Output_Trans
    ENDIF

END SUBROUTINE Output


!======================================================================
!----------------------------------------------------------------------
! Show Informations on the Screen
!----------------------------------------------------------------------
SUBROUTINE ShowOnScreen( CPUTime )
IMPLICIT NONE
REAL(KIND=8):: CPUTime(4)

    WRITE(*, "(I7, 2X, I8, 2X, 3(F6.3, 2X))") iter, RNph, &
                                             CPUTime(2) - CPUTime(1), &
                                             CPUTime(3) - CPUTime(2), &
                                             CPUTime(4) - CPUTime(3)

END SUBROUTINE ShowOnScreen



!======================================================================
END MODULE IO