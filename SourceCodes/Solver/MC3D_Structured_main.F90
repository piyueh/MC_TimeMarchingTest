!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_main
!! Contenent: Main Program
!! Author: PY Chuang
!!#####################################################################
PROGRAM main
USE RNG
USE VAR_ALL
USE ADV
USE IO
USE ROUTINES
IMPLICIT NONE
REAL(KIND=8):: R1, t(5)
INTEGER(KIND=4):: i, NCores
TYPE(rng_t), ALLOCATABLE:: SeedMP(:)


    WRITE(*, '("Enter CASE Name: ")', ADVANCE = 'NO')
    !READ(*, *) CaseName
    CaseName = "test"
    InputFileName = casename(1:LEN_TRIM(casename))//'_initial.txt'


    WRITE(*, '("Enter BCs (BCx, BCy, BCz): ")', ADVANCE = 'NO')
    !READ(*, *) BCs
    BCs = (/ 1, 2, 2 /)


    WRITE(*, '("Enter the Numer of CPU Cores: ")', ADVANCE = 'NO')
    !READ(*, *) NCores
    NCores = 1
    ALLOCATE( SeedMP(NCores) )
    CALL RANDOM_SEED()
    DO i = 1, NCores
        CALL RANDOM_NUMBER( R1 )
        CALL RNG_SEED( SeedMP(i), INT( R1*1D8 ) )
    ENDDO

    WRITE(*, '("Enter the Maximum Iterations: ")', ADVANCE = 'NO')
    !READ(*, *) iterations
    iterations = 20000000


    WRITE(*, '("Enter How Many Steps It Outputs Once: ")', ADVANCE = 'NO')
    !READ(*, *) nOutput
    nOutput = 1000

    WRITE(*, '("Enter the Method for Time Marching: ")', ADVANCE = 'NO')
    !READ(*, *) WAY_FlightTime
    WAY_FlightTime = 1

    CALL Initialize_Ge( Ge_table, Ge_start, dU_Ge, N_Ge1, N_Ge2 )
    CALL Initialize_Si( Si_table, Si_start, dU_Si, N_Si1, N_Si2 )

    CALL initialize
    CALL Reorder_CellInfo
    CALL CreateDelete( SeedMP(1) )

    CALL Output_Trans

    DO iter = iter0 + 1, iterations

        CALL advance( NCores, SeedMP, t )
        time = time + TimeStep

        CALL Output

        CALL ShowOnScreen( t )

    ENDDO

    DEALLOCATE( SeedMP )

END PROGRAM main


!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Controller
!----------------------------------------------------------------------
SUBROUTINE advance( NCores, SeedMP, t )
USE RNG
USE HEAT
USE ADV
USE VAR_BC, ONLY: BCs
USE VAR_ph, ONLY: FNph, phn
USE VAR_Others, ONLY: TimeStep, WAY_FlightTime
USE ROUTINES, ONLY: Reorder_CellInfo
REAL(KIND=8):: dtRemain, t(5)
INTEGER(KIND=4):: NCores, iCPU, i
TYPE(rng_t):: SeedMP(NCores)

    iCPU = 1

    CALL CPU_TIME( t(1) )

    DO i = 1, FNph
        IF ( phn(i).Exist ) THEN
            dtRemain = TimeStep
            DO WHILE ( dtRemain.gt.0 )
                SELECTCASE( WAY_FlightTime )
                CASE(1)
                    CALL phn_adv_CT( phn(i), SeedMP(iCPU), dtRemain )
                CASE(2)
                    CALL phn_adv_RT( phn(i), SeedMP(iCPU), dtRemain )
                END SELECT
            ENDDO
        ENDIF
    ENDDO

    CALL CPU_TIME( t(2) )
    IF ( BCs(1).eq.3 ) CALL Heat_Control( SeedMP(iCPU) )

    CALL CPU_TIME( t(3) )
    CALL Reorder_CellInfo
    
    CALL CPU_TIME( t(4) )
    CALL CreateDelete( SeedMP(iCPU) )
    
    CALL CPU_TIME( t(5) )

END SUBROUTINE advance
