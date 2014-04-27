MODULE A
IMPLICIT NONE
INTEGER, PARAMETER:: aa=1
    CONTAINS

    SUBROUTINE aaa
    IMPLICIT NONE
    
        WRITE(*, *) aa
    
    END SUBROUTINE aaa

END MODULE A

MODULE B
USE A
IMPLICIT NONE
INTEGER, PARAMETER:: bb=1
    CONTAINS
    SUBROUTINE bbb
    IMPLICIT NONE
    
        WRITE(*, *) bb
    
    END SUBROUTINE bbb

END MODULE B

PROGRAM main
USE B
IMPLICIT NONE

    CALL aaa

END PROGRAM main




! PROGRAM main
! IMPLICIT NONE
! REAL(KIND=8):: Ge_start, Si_start, dU_Ge, dU_Si
! INTEGER(KIND=4):: N_Ge1, N_Ge2, N_Si1, N_Si2, i
! REAL(KIND=8), ALLOCATABLE:: Ge_table(:, :), Si_table(:, :)
! 
!     OPEN( UNIT = 110, FILE = "Ge_real_table.txt" )
!     OPEN( UNIT = 120, FILE = "Si_real_table.txt" )
! 
!     READ(110,*) N_Ge1, N_Ge2
!     READ(120,*) N_Si1, N_Si2
! 
!     ALLOCATE( Ge_table(N_Ge1, N_Ge2) )
!     ALLOCATE( Si_table(N_Si1, N_Si2) )
!     READ(110,*) Ge_table
!     READ(120,*) Si_table
! 
!     CLOSE( 110 )
!     CLOSE( 120 )
! 
!     Ge_start = Ge_table(3, 1)
!     dU_Ge = Ge_table(3, 2) - Ge_table(3, 1)
!     Si_start = Si_table(3, 1)
!     dU_Si = Si_table(3, 2) - Si_table(3, 1)
!     
!     OPEN(120, FILE = "MC3D_Structured_MAT.f90")
!     WRITE(120, *) "MODULE MAT"
!     WRITE(120, *) "IMPLICIT NONE"
!     WRITE(120, *) "REAL(KIND=8):: Ge_start = ", Ge_start
!     WRITE(120, *) "REAL(KIND=8):: Si_start = ", Si_start
!     WRITE(120, *) "REAL(KIND=8):: dU_Ge = ", dU_Ge
!     WRITE(120, *) "REAL(KIND=8):: dU_Si = ", dU_Si
!     WRITE(120, *) "INTEGER(KIND=4):: N_Ge1 = 7, N_Ge2 = 8001"
!     WRITE(120, *) "INTEGER(KIND=4):: N_Si1 = 7, N_Si2 = 8001"
!     WRITE(120, *)
!     WRITE(120, FMT = *) "REAL(KIND=8):: Ge_table(7, 8001)"
!     WRITE(120, FMT = *) "REAL(KIND=8):: Si_table(7, 8001)"
!     WRITE(120, *) "CONTAINS"
!     WRITE(120, *)
!     WRITE(120, *) "SUBROUTINE Initialize_Materials"
!     WRITE(120, *) "IMPLICIT NONE"
!     WRITE(120, *) 
!     DO i = 1, 8001
!         WRITE(UNIT = 120, FMT = "('Ge_table(:,', I4, ') = (/ ')", ADVANCE='NO')  i
!         WRITE(120, FMT = 100) Ge_table(:, i)
!     ENDDO
!     WRITE(120, *)
!     DO i = 1, 8001
!         WRITE(UNIT = 120, FMT = "('Si_table(:,', I4, ') = (/ ')", ADVANCE='NO')  i
!         WRITE(120, FMT = 100) Si_table(:, i)
!     ENDDO
!     WRITE(120, *)
!     WRITE(120, *) "END SUBROUTINE Initialize_Materials"
!     WRITE(120, *)
!     WRITE(120, *) "END MODULE MAT"
!     CLOSE( 120 )
! 
!     100 FORMAT(3(F, ","), "&", /, 3(F, ","), F, " /)")
! 
! END PROGRAM main
! PROGRAM main
! IMPLICIT NONE
! INTEGER(KIND=4):: N
! INTEGER(KIND=4), ALLOCATABLE:: A(:)
! NAMELIST /test/ N, A
! NAMELIST /t1/ N
! NAMELIST /t2/ A
!     
!     OPEN( UNIT = 120, FILE = "Namelist.txt" )
!     READ(120, NML = t1)
!     ALLOCATE( A(N) )
!     READ(120, NML = t2)
!     CLOSE( 120 )
!     
!     WRITE(*, *) N
!     WRITE(*, *) A
! !     OPEN( UNIT = 120, FILE = "N_A.txt" )
! !     WRITE(120, *) N
! !     WRITE(120, *) A
! !     CLOSE( 120 )
!     
! !     B = (/ 1, 2, 3 /)
! !     SELECTCASE(B(1))
! !     CASE(1)
! !         WRITE(*, *) '1'
! !     CASE(2)
! !         WRITE(*, *) '2'
! !     END SELECT
! 
! !     DO k = 1, 3; DO j = 1, 3; DO i = 1, 3
! !         A(i, j, k) = i * 100 + j * 10 + k
! !     ENDDO; ENDDO; ENDDO
! !     
! !     DO k = 1, 3; DO j = 1, 3; DO i = 1, 3
! !         WRITE(*, *) A(a, j, k)
! !     ENDDO; ENDDO; ENDDO
! 
! END PROGRAM main

! PROGRAM main
! IMPLICIT NONE
! INTEGER(KIND=4):: i, phn(15), s
! 
! TYPE ElePh
!     INTEGER(KIND=4):: Idx
!     INTEGER(KIND=4):: Value
!     TYPE(ElePh), POINTER:: prevPh => NULL()
!     TYPE(ElePh), POINTER:: nextPh => NULL()
! END TYPE ElePh
! 
! TYPE ElePhbged
!     TYPE(ElePh), POINTER:: bgPh => NULL()
!     TYPE(ElePh), POINTER:: edPh => NULL()
! END TYPE ElePhbged
! TYPE(ElePhbged):: ePh(3)
! TYPE(ElePh), POINTER:: aPh
! 
!     phn = (/ 1, 2, 2, 2, 1, 1, 3, 3, 2, 3, 1, 2, 3, 2, 3 /)
!     ! 1:  1, 5, 6, 11   => 4
!     ! 2:  2, 3, 4, 9, 12, 14  =>6
!     ! 3:  7, 8, 10, 13, 15  => 5
!     
!     DO i = 1, 3 
!         ALLOCATE( ePh(i)%bgPh )
!         ePh(i)%bgPh%Idx = 1
!         ePh(i)%edPh => ePh(i)%bgPh
!     ENDDO
!     
!     DO i = 1, 15
!         ePh(phn(i))%edPh%Value = i
!         ALLOCATE( ePh(phn(i))%edPh%nextPh )
!         ePh(phn(i))%edPh%nextPh%prevPh => ePh(phn(i))%edPh
!         ePh(phn(i))%edPh => ePh(phn(i))%edPh%nextPh
!         ePh(phn(i))%edPh%Idx = ePh(phn(i))%edPh%prevPh%Idx + 1
!     ENDDO
!     
!     DO i = 1, 3
!         ePh(i)%edPh => ePh(i)%edPh%prevPh
!         DEALLOCATE( ePh(i)%edPh%nextPh )
!     ENDDO
!     
!     DO i = 1, 3
!         WRITE(*, *) "Element ", i
!         aPh => ePh(i)%bgPh
!         DO s = ePh(i)%bgPh%Idx, ePh(i)%edPh%Idx
!             WRITE(*, *) aPh%Idx, aPh%Value
!             IF ( ASSOCIATED( aPh%nextPh ) ) aPh => aPh%nextPh
!         ENDDO
!         WRITE(*, *)
!     ENDDO
! 
! END PROGRAM main