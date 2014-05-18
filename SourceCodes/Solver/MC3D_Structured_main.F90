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
REAL(KIND=8):: R1
INTEGER(KIND=4):: i, FLAG
TYPE(rng_t):: SeedMP

    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( R1 )
    CALL RNG_SEED( SeedMP, INT( R1*1D8 ) )
    
    WRITE(*, '("Enter CASE Name: ")', ADVANCE = 'NO')
    READ(*, *) CaseName
    InputFileName = casename(1:LEN_TRIM(casename))//'_initial.txt'

    WRITE(*, '("Enter BCs (BCx, BCy, BCz): ")', ADVANCE = 'NO')
    READ(*, *) BCs

    WRITE(*, '("Enter the Maximum Iterations: ")', ADVANCE = 'NO')
    READ(*, *) iterations
        
    WRITE(*, '("Enter How Many Steps It Outputs Once: ")', ADVANCE = 'NO')
    READ(*, *) nOutput

    WRITE(*, '("Enter the Method for Time Marching: ")')
    WRITE(*, '("(1: Constant Flight Time, 2: Random Flight Time) ")', ADVANCE = 'NO')
    READ(*, *) WAY_FlightTime

    
    CALL Initialize_Ge( Ge_table, Ge_start, dU_Ge, N_Ge1, N_Ge2 )
    CALL Initialize_Si( Si_table, Si_start, dU_Si, N_Si1, N_Si2 )

    
    CALL initialize
    CALL Reorder_CellInfo
    CALL CreateDelete( SeedMP )

    
    CALL Output_Trans
    
    
    CALL Flag_Gen( BCs(1), WAY_FlightTime, FLAG )
    
    
    CALL Iter_CTRL( FLAG )
    

END PROGRAM main

!======================================================================
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE Flag_Gen( BCx, WFT, FLAG )
IMPLICIT NONE
INTEGER(KIND=4), INTENT(IN):: BCx, WFT
INTEGER(KIND=4), INTENT(OUT):: FLAG

    SELECT CASE( BCx )
    CASE(1:2)
        SELECT CASE( WFT )
        CASE(1)
            FLAG = 1
        CASE(2)
            FLAG = 2
        CASE DEFAULT
            ! left for error msg.
        END SELECT
    CASE(3)
        SELECT CASE( WFT )
        CASE(1)
            FLAG = 3
        CASE(2)
            FLAG = 4
        CASE DEFAULT
            ! left for error msg.
        END SELECT
    CASE DEFAULT
        ! left for error msg.
    END SELECT 

END SUBROUTINE Flag_Gen