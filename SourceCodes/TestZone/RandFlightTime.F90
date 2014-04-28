!######################################################################
! Test of Random Free Flight Time
!######################################################################

!======================================================================
!======================================================================
MODULE VAR
USE MATERIALS
IMPLICIT NONE

REAL(KIND=8):: GammaT, dtP
REAL(KIND=8):: TBCL, TBCR

ENd MODULE VAR
!======================================================================
!======================================================================
PROGRAM main
USE VAR
IMPLICIT NONE
INTEGER(KIND=4):: i

    CALL RANDOM_SEED()
    CALL Initialize_Materials

    TBCL = 4D2
    TBCR = 2D2

    CALL Pseudo_ScatteringRate

!     DO i= 1, 10
!         CALL FreeFlightTime( dtP, GammaT )
!         WRITE(*, *) i, dtP
!     ENDDO

END PROGRAM main

!======================================================================
!======================================================================
SUBROUTINE Pseudo_ScatteringRate
USE VAR
IMPLICIT NONE
REAL(KIND=8):: T, E(2), V(2), MFP(2), G(2)
INTEGER(KIND=4):: i, j

    OPEN(UNIT = 120, FILE = "test.txt")
    DO j = 1, 8001

        WRITE(120, *) Si_table(1, j), Si_table(4, j) / Si_table(5, j)
!     T = MAX( TBCL, TBCR )
!         T = Ge_table(1, j)
!         DO i = 1, 2
!             CALL energy( i, T, E(i) )
!             CALL Etable( i, 4, E(i), V(i) )
!             CALL Etable( i, 5, E(i), MFP(i) )
!             G(i) = V(i) / MFP(i)
!         ENDDO
!         WRITE(*, *) G
    ENDDO
    CLOSE( 120 )
!     T = MIN( TBCL, TBCR )
!     DO i = 1, 2
!         CALL energy( i, T, E(i) )
!         CALL Etable( i, 4, E(i), V(i) )
!         CALL Etable( i, 5, E(i), MFP(i) )
!         G(i) = V(i) / MFP(i)
!     ENDDO
!     WRITE(*, *) G
!     GammaT = MAXVAL( G )

END SUBROUTINE Pseudo_ScatteringRate
!======================================================================
!======================================================================
SUBROUTINE FreeFlightTime( dt, GammaT )
IMPLICIT NONE
REAL(KIND=8):: R1, dt, GammaT

    CALL RANDOM_NUMBER( R1 )
    dt = - DLOG(R1) / GammaT

END SUBROUTINE FreeFlightTime