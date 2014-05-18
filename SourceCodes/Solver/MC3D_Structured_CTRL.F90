!======================================================================
!----------------------------------------------------------------------
! Iteration Controller
!----------------------------------------------------------------------
SUBROUTINE Iter_CTRL( FLAG, SeedMP )
IMPLICIT NONE
INTEGER(KIND=4), INTENT(IN):: FLAG
REAL(KIND=8):: t(5)
        
    SELECT CASE( FLAG )
    CASE(1)
        DO iter = iter0 + 1, iterations
            CALL ADV_CFFT
            CALL Output
            CALL ShowOnScreen( t )
        ENDDO
    CASE(2)
        DO iter = iter0 + 1, iterations
            CALL ADV_RFFT
            CALL Output
            CALL ShowOnScreen( t )
        ENDDO
    CASE(3)
        DO iter = iter0 + 1, iterations
            CALL ADV_CFFT_H
            CALL Output
            CALL ShowOnScreen( t )
        ENDDO
    CASE(4)
        DO iter = iter0 + 1, iterations
            CALL ADV_RFFT_H
            CALL Output
            CALL ShowOnScreen( t )
        ENDDO
    CASE DEFAULT
        ! left for error msg.
    END SELECT

END SUBROUTINE Iter_CTRL


!======================================================================
!----------------------------------------------------------------------
! ADV_CFFT: Controller for Constant Free Flight Time w/o Hear Injection
!----------------------------------------------------------------------
SUBROUTINE ADV_CFFT(  )
IMPLICIT NONE

    CALL CPU_TIME( t(1) )
    DO i = 1, FNph
        IF ( phn(i).Exist ) THEN
            dtRemain = TimeStep
            DO WHILE ( dtRemain.gt.0 )
                CALL phn_adv_CT( phn(i), SeedMP, dtRemain, 1 )
            ENDDO
        ENDIF
    ENDDO
    CALL CPU_TIME( t(2) )
    
END SUBROUTINE ADV_CFFT


!======================================================================
!----------------------------------------------------------------------
! ADV_RFFT: Controller for Random Free Flight Time w/o Hear Injection
!----------------------------------------------------------------------
SUBROUTINE ADV_RFFT(  )
IMPLICIT NONE

    CALL CPU_TIME( t(1) )
    DO i = 1, FNph
        IF ( phn(i).Exist ) THEN
            dtRemain = TimeStep
            DO WHILE ( dtRemain.gt.0 )
                CALL phn_adv_RT( phn(i), SeedMP, dtRemain, 1 )
            ENDDO
        ENDIF
    ENDDO
    CALL CPU_TIME( t(2) )
    
END SUBROUTINE ADV_RFFT







    CALL CPU_TIME( t(2) )

    CALL CPU_TIME( t(3) )
    CALL Reorder_CellInfo

    CALL CPU_TIME( t(4) )
    CALL CreateDelete( SeedMP(iCPU) )

    CALL CPU_TIME( t(5) )
    time = time + TimeStep



!======================================================================
!----------------------------------------------------------------------
! ADV_CFFT: Controller for Constant Free Flight Time w/ Hear Injection
!----------------------------------------------------------------------
SUBROUTINE ADV_CFFT
IMPLICIT NONE

    CALL CPU_TIME( t(1) )

    DO i = 1, FNph
        IF ( phn(i).Exist ) THEN
            dtRemain = TimeStep
            DO WHILE ( dtRemain.gt.0 )
                CALL phn_adv_CT( phn(i), SeedMP, dtRemain, 1 )
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
    time = time + TimeStep

END SUBROUTINE ADV_CFFT



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
                    CALL phn_adv_CT( phn(i), SeedMP(iCPU), dtRemain, 1 )
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
    time = time + TimeStep

END SUBROUTINE advance
