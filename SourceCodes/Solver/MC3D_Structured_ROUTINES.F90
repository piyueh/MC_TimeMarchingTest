!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_ROUTINE
!! Contenent: Collection of miscellaneous subroutines.
!! Author: PY Chuang
!!#####################################################################
MODULE ROUTINES
USE RNG
IMPLICIT NONE
CONTAINS
!======================================================================


!======================================================================
!----------------------------------------------------------------------
! Phonon reorder
!----------------------------------------------------------------------
SUBROUTINE Reorder_CellInfo
USE VAR_ALL
IMPLICIT NONE
REAL(KIND=8):: R1
INTEGER(KIND=4):: i, j, k, tmpI1, mt
TYPE(Phonon), ALLOCATABLE:: tmpPhn(:)


    NEmpty = 0
    ele%E = 0D0
    ele%Ntol = 0
    EmptyID = -1
    phId = -1

    RNph = COUNT( phn%Exist )
    RNHeatPh = COUNT( HeatPhn%Exist )

    IF ( (RNHeatPh + RNph).gt.FNph ) THEN

        mt = RNHeatPh + RNph

        ALLOCATE( tmpPhn(mt) )

        mt = 0
        DO tmpI1 = 1, FNph
            IF ( phn(tmpI1)%Exist ) THEN
                i = phn(tmpI1)%eID(1)
                j = phn(tmpI1)%eID(2)
                k = phn(tmpI1)%eID(3)

                ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
                ele(i, j, k)%E = ele(i, j, k)%E + phn(tmpI1)%E

                mt = mt + 1
                tmpPhn(mt) = phn(tmpI1)
            ENDIF
        ENDDO

        DO tmpI1 = 1, FNHeatPh
            IF ( HeatPhn(tmpI1)%Exist ) THEN
                i = HeatPhn(tmpI1)%eID(1)
                j = HeatPhn(tmpI1)%eID(2)
                k = HeatPhn(tmpI1)%eID(3)

                ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
                ele(i, j, k)%E = ele(i, j, k)%E + HeatPhn(tmpI1)%E

                mt = mt + 1
                tmpPhn(mt) = HeatPhn(tmpI1)
                HeatPhn(tmpI1)%Exist = .FALSE.
            ENDIF
        ENDDO

        DEALLOCATE( phn )
        FNph = INT( mt * 1.2 + 0.5 )
        ALLOCATE( phn(FNph) )
        phn(1:mt) = tmpPhn
        DEALLOCATE( tmpPhn )

        NEmpty = FNph - mt
        EmptyID(1:NEmpty) = (/ (i, i = mt+1, FNph) /)

    ELSE

        DO tmpI1 = 1, FNph
            IF ( phn(tmpI1)%Exist ) THEN
                i = phn(tmpI1)%eID(1)
                j = phn(tmpI1)%eID(2)
                k = phn(tmpI1)%eID(3)
                
!                 IF (  ANY( (phn(tmpI1)%xyz - ele(i, j, k)%BD(1,:)).lt.0 ).OR.ANY( (phn(tmpI1)%xyz - ele(i, j, k)%BD(2,:)).gt.0 ) ) THEN
!                     WRITE(*, *) phn(tmpI1)
!                     CALL Errors(5001)
!                 ENDIF
                
                ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
                ele(i, j, k)%E = ele(i, j, k)%E + phn(tmpI1)%E
            ELSE
                NEmpty = NEmpty + 1
                EmptyID(NEmpty) = tmpI1
            ENDIF
        ENDDO

        DO tmpI1 = 1, FNHeatPh
            IF ( HeatPhn(tmpI1)%Exist ) THEN
                i = HeatPhn(tmpI1)%eID(1)
                j = HeatPhn(tmpI1)%eID(2)
                k = HeatPhn(tmpI1)%eID(3)

                ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
                ele(i, j, k)%E = ele(i, j, k)%E + HeatPhn(tmpI1)%E

                phn(EmptyID(NEmpty)) = HeatPhn(tmpI1)
                NEmpty = NEmpty - 1
                HeatPhn(tmpI1)%Exist = .FALSE.
            ENDIF
        ENDDO

    ENDIF

    RNph = RNph + RNHeatPh

    ele%E = ele%E + ele%Ediff

!      IF ( ANY( ele%E.le.0D0 ) ) CALL Errors(3)
!      IF ( RNph.ne.SUM( ele%Ntol ) ) CALL Errors(100)
!      IF ( NEmpty.ne.(FNph - RNph) ) CALL Errors(4)

    tmpI1 = 0
    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)

        ele(i, j, k)%Nbg = tmpI1 + 1
        ele(i, j, k)%Ned = tmpI1 + ele(i, j, k)%Ntol
        tmpI1 = tmpI1 + ele(i, j, k)%Ntol

        R1 = ele(i, j, k)%E / dV
        mt = ele(i, j, k)%Mat

        CALL Etable( mt, 1, R1, ele(i, j, k)%T )
        CALL Etable( mt, 2, R1, ele(i, j, k)%ND )
        ele(i, j, k)%Eph = R1 / ele(i, j, k)%ND * bundle
        CALL Etable( mt, 4, R1, ele(i, j, k)%Vph )
        CALL Etable( mt, 5, R1, ele(i, j, k)%MFP )
        ele(i, j, k)%SCR = ele(i, j, k)%Vph / ele(i, j, k)%MFP

        IF ( ele(i, j, k)%SCR.ge.GammaT ) THEN
            WRITE(*, *) ele(i, j, k)%T
            WRITE(*, *) ele(i, j, k)%Vph
            WRITE(*, *) ele(i, j, k)%MFP
            WRITE(*, *) ele(i, j, k)%SCR, GammaT
            CALL Errors(9)
        ENDIF

    ENDDO; ENDDO; ENDDO

    ele%Ntol = 0
    DO tmpI1 = 1, FNph
        IF ( phn(tmpI1)%Exist ) THEN
            i = phn(tmpI1)%eID(1)
            j = phn(tmpI1)%eID(2)
            k = phn(tmpI1)%eID(3)
            ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
            phId(ele(i, j, k)%Nbg + ele(i, j, k)%Ntol - 1) = tmpI1
        ENDIF
    ENDDO
!      IF ( ele(Ne(1), Ne(2), Ne(3))%Ned.ne.RNph ) CALL Errors(5)
!      IF ( SUM( ele%Ntol ).ne.RNph ) CALL Errors(6)

END SUBROUTINE Reorder_CellInfo


!======================================================================
!----------------------------------------------------------------------
! Determine Pseudo Scattering Rate
!----------------------------------------------------------------------
SUBROUTINE Pseudo_ScatteringRate( T1, T2, T3, PseudoG )
USE VAR_MAT
IMPLICIT NONE
REAL(KINd=8), INTENT(IN):: T1, T2, T3
REAL(KIND=8), INTENT(OUT):: PseudoG
REAL(KIND=8):: T, E(2), V(2), MFP(2), G(2)
INTEGER(KIND=4):: i

    T = MAX( T1, T2, T3 ) + 1D2
    DO i = 1, 2
        CALL energy( i, T, E(i) )
        CALL Etable( i, 4, E(i), V(i) )
        CALL Etable( i, 5, E(i), MFP(i) )
        G(i) = V(i) / MFP(i)
    ENDDO
    PseudoG = MAXVAL( G )
    WRITE(*, *) "Use Temperature ", T, " K to generate pseudo gamma."
    WRITE(*, *) "The pseudo scattering rate is: ", PseudoG

END SUBROUTINE Pseudo_ScatteringRate


!======================================================================
!----------------------------------------------------------------------
! Determine Random Free Flight Time
!----------------------------------------------------------------------
SUBROUTINE FreeFlightTime( dt, Rseed )
USE VAR_Others, ONLY: GammaT
IMPLICIT NONE
REAL(KIND=8):: R, dt
TYPE(rng_t):: Rseed

    CALL RAN_NUM( Rseed, R )

    dt = - DLOG(R) / GammaT

END SUBROUTINE FreeFlightTime


!======================================================================
!----------------------------------------------------------------------
! Errors:
! code = 1: SUBROUTINE initialize
!      = 2: SUBROUTINE initialize
!      = 3: SUBROUTINE Reorder_CellInfo
!      = 4: SUBROUTINE Reorder_CellInfo
!      = 5: SUBROUTINE Reorder_CellInfo
!      = 6: SUBROUTINE Reorder_CellInfo
!      = 7: SUBROUTINE CreateDelete
!      = 8: SUBROUTINE CreateDelete
!----------------------------------------------------------------------
SUBROUTINE Errors( code )
IMPLICIT NONE
INTEGER(KIND=4), INTENT(IN):: code

    WRITE(*, *) "ERROR: ERROR CODE IS ", code, " !"
    WRITE(*, *) "PRESS ENTER TO STOP THE PROGRAM!"
    READ(*, *)
    STOP

END SUBROUTINE Errors


!======================================================================
!----------------------------------------------------------------------
! Return target propertie (out) under specific material (mat) and
! energy density (E).
! mat = 1: Ge,  2: Si
! idx = 1: temperature,  2: number density,  4: velocity,
!       5: MFP,  6: specific heat,  7: thermal conductivity
!----------------------------------------------------------------------
SUBROUTINE Etable( mat, idx, E, out )
USE VAR_MAT
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: E
REAL(KIND=8), INTENT(OUT):: out
INTEGER(KIND=4), INTENT(IN):: mat, idx
REAL(KIND=8):: out_a, out_b
INTEGER(KIND=4):: iU

    IF ( mat.eq.1 ) THEN
        iU = INT( (E - Ge_start) / dU_Ge ) + 1
        out_a = Ge_table(idx, iU)
        out_b = Ge_table(idx, iU+1)
        out = out_a + (out_b - out_a) * (E - Ge_table(3, iU)) / dU_Ge
    ELSE IF ( mat.eq.2 ) THEN
        iU = int( (E - Si_start) / dU_Si ) + 1
        out_a = Si_table(idx, iU)
        out_b = Si_table(idx, iU+1)
        out = out_a + (out_b - out_a) * (E - Si_table(3, iU)) / dU_Si
    END IF

END SUBROUTINE Etable


!======================================================================
!----------------------------------------------------------------------
! Return energy density (Eout) under specific material (nc) and
! specific temperature (T0)
!----------------------------------------------------------------------
SUBROUTINE energy( nc, T0, Eout )
USE VAR_MAT
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: T0
REAL(KIND=8), INTENT(OUT):: Eout
INTEGER(KIND=4), INTENT(IN):: nc
REAL(KIND=8):: E1, E2, T1, T2, Tout

    IF ( nc.eq.1 ) THEN

        E1 = Ge_table(3, 1)
        E2 = Ge_table(3, N_Ge2)
        T1 = Ge_table(1, 1) - T0
        T2 = Ge_table(1, N_Ge2) - T0
    ELSE IF ( nc.eq.2 ) THEN

        E1 = Si_table(3, 1)
        E2 = Si_table(3, N_Si2)
        T1 = Si_table(1, 1) - T0
        T2 = Si_table(1, N_Si2) - T0

    ENDIF

    Eout = (E1 + E2) / 2D0
    CALL ETable( nc, 1, Eout, Tout )
    Tout = Tout - T0

    DO WHILE ( ABS( Tout ).gt.1D-12 )

        IF ( (Tout * T1).lt.0D0 ) THEN
            E2 = Eout
            T2 = Tout
        ELSE
            E1 = Eout
            T1 = Tout
        ENDIF

        Eout = (E1 + E2) / 2D0
        CALL ETable( nc, 1, Eout, Tout )
        Tout = Tout - T0

    ENDDO

END SUBROUTINE energy

!======================================================================
END MODULE ROUTINES