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
SUBROUTINE Reorder_CellInfo
USE VAR
IMPLICIT NONE
REAL(KIND=8):: R1
INTEGER(KIND=4):: i, j, k, tmpI1, mt

    NEmpty = 0
    ele%E = 0D0
    ele%Ntol = 0
    EmptyID = -1
    phId = -1
    i = 0
    
    DO tmpI1 = 1, FNph
        IF ( phn(tmpI1)%Exist ) THEN
            i = phn(tmpI1)%eID(1)
            j = phn(tmpI1)%eID(2)
            k = phn(tmpI1)%eID(3)
            
            ele(i, j, k)%Ntol = ele(i, j, k)%Ntol + 1
            
            ele(i, j, k)%E = ele(i, j, k)%E + phn(tmpI1)%E
        ELSE
            NEmpty = NEmpty + 1
            EmptyID(NEmpty) = tmpI1
        ENDIF
    ENDDO
    
    ele%E = ele%E + ele%Ediff
    IF ( ANY( ele%E.le.0D0 ) ) CALL Errors(3)
    
    RNph = SUM( ele%Ntol )
    IF ( NEmpty.ne.(FNph - RNph) ) CALL Errors(4)
    
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
    IF ( ele(Ne(1), Ne(2), Ne(3))%Ned.ne.RNph ) CALL Errors(5)
    IF ( SUM( ele%Ntol ).ne.RNph ) CALL Errors(6)
    
END SUBROUTINE Reorder_CellInfo


!======================================================================
!----------------------------------------------------------------------
! Get informations of elements from phonons
!----------------------------------------------------------------------
! SUBROUTINE cellInfo
! USE VAR, ONLY: Nph, phn, Ne, ele, dV, bundle
! IMPLICIT NONE
! INTEGER(KIND=4):: i, j, k
! INTEGER(KIND=4):: I1, I2, I3
! 
!     ele%E = 0D0
!     
!     DO i = 1, Nph
!         I1 = phn(i)%eID(1)
!         I2 = phn(i)%eID(2)
!         I3 = phn(i)%eID(3)
!         ele(I1, I2, I3)%E = ele(I1, I2, I3)%E + phn(i)%E
!     ENDDO
!     
!     ele%E = (ele%E + ele%Ediff) / dV
!     
!     DO k = 1, Ne(3)
!         DO j = 1, Ne(2)
!             DO i = 1, Ne(1)
!             
!                 CALL Etable( ele(i, j, k)%Mat, 1, &
!                              ele(i, j, k)%E, ele(i, j, k)%T )
!                 
!                 CALL Etable( ele(i, j, k)%Mat, 2, &
!                              ele(i, j, k)%E, ele(i, j, k)%N )
!                              
!                 ele(i, j, k)%Eph = ele(i, j, k)%E / ele(i, j, k)%N * &
!                                    bundle
!                                    
!                 CALL Etable( ele(i, j, k)%Mat, 4, &
!                              ele(i, j, k)%E, ele(i, j, k)%Vph )
!                 
!                 CALL Etable( ele(i, j, k)%Mat, 5, &
!                              ele(i, j, k)%E, ele(i, j, k)%MFP )
!                              
!             ENDDO
!         ENDDO
!     ENDDO
! 
! END SUBROUTINE cellInfo


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
END MODULE ROUTINES