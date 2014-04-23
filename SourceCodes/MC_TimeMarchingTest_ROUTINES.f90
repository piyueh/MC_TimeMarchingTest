!!#####################################################################
!! Module Name: ROUTINES
!! Purpose: Collection of miscellaneous subroutines.
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
SUBROUTINE reorder
USE VAR, ONLY: N, bg, ed, Nph, phn, Ne, phID
IMPLICIT NONE
INTEGER(KIND=4):: i, j, k, tmpI1

    N = 0
    DO i = 1, Nph
        IF ( phn(i)%E.gt.0 ) THEN
            N(phn(i)%eID(1), phn(i)%eID(2), phn(i)%eID(3)) = &
                     N(phn(i)%eID(1), phn(i)%eID(2), phn(i)%eID(3)) + 1
        ENDIF
    ENDDO
    
    tmpI1 = 0
    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)
                bg(i, j, k) = tmpI1 + 1
                ed(i, j, k) = tmpI1 + N(i, j, k)
                tmpI1 = tmpI1 + N(i, j, k)
            ENDDO
        ENDDO
    ENDDO
    
    N = 0
    DO i = 1, Nph
        IF ( phn(i)%E.gt.0 ) THEN
            N(phn(i)%eID(1), phn(i)%eID(2), phn(i)%eID(3)) = &
                     N(phn(i)%eID(1), phn(i)%eID(2), phn(i)%eID(3)) + 1
            phID(bg(i, j, k) + &
                N(phn(i)%eID(1), phn(i)%eID(2), phn(i)%eID(3)) - 1) = i
        ENDIF
    ENDDO
    
    Nph = SUM( N )
    
END SUBROUTINE reorder


!======================================================================
!----------------------------------------------------------------------
! Get informations of elements from phonons
!----------------------------------------------------------------------
SUBROUTINE cellInfo
USE VAR, ONLY: Nph, phn, Ne, ele, dV, bundle
IMPLICIT NONE
INTEGER(KIND=4):: i, j, k
INTEGER(KIND=4):: I1, I2, I3

    ele%E = 0D0
    
    DO i = 1, Nph
        I1 = phn(i)%eID(1)
        I2 = phn(i)%eID(2)
        I3 = phn(i)%eID(3)
        ele(I1, I2, I3)%E = ele(I1, I2, I3)%E + phn(i)%E
    ENDDO
    
    ele%E = (ele%E + ele%Ediff) / dV
    
    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)
            
                CALL Etable( ele(i, j, k)%Mat, 1, &
                             ele(i, j, k)%E, ele(i, j, k)%T )
                
                CALL Etable( ele(i, j, k)%Mat, 2, &
                             ele(i, j, k)%E, ele(i, j, k)%N )
                             
                ele(i, j, k)%Eph = ele(i, j, k)%E / ele(i, j, k)%N * &
                                   bundle
                                   
                CALL Etable( ele(i, j, k)%Mat, 4, &
                             ele(i, j, k)%E, ele(i, j, k)%Vph )
                
                CALL Etable( ele(i, j, k)%Mat, 5, &
                             ele(i, j, k)%E, ele(i, j, k)%MFP )
                             
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE cellInfo


!======================================================================
!----------------------------------------------------------------------
! Return target propertie (out) under specific material (mat) and
! energy density (E).
! mat = 1: Ge,  2: Si
! idx = 1: temperature,  2: number density,  4: velocity,
!       5: MFP,  6: specific heat,  7: thermal conductivity
!----------------------------------------------------------------------
SUBROUTINE Etable( mat, idx, E, out )
USE VAR, ONLY: Ge_table, Si_table, dU_Ge, dU_Si, Ge_start, Si_start
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
USE VAR, ONLY: Ge_table, Si_table, N_Si2, N_Ge2
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
!----------------------------------------------------------------------
! Errors
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
END MODULE ROUTINES