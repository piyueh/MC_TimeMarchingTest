!!#####################################################################
!! Module Name: INIT
!! Purpose: Collection of initialization subroutines
!! Author: PY Chuang
!!#####################################################################
MODULE INIT
USE RNG
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
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: Npercell
INTEGER(KIND=4):: tmpI1
INTEGER(KIND=4), ALLOCATABLE:: N(:, :, :)

    dt = 1D-1
    iter0 = 0

    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( tmpR1 )

    ALLOCATE( SeedMP(NCores) )
    DO i = 1, NCores
        CALL RNG_SEED( SeedMP(i), INT( tmpR1*1D8 ) )
    ENDDO

    L = (/ 1D3, 1D3, 1D1 /)
    Ne = (/ 5D1, 5D1, 1D0 /)
    dL = L / DBLE( Ne )

    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( N(Ne(1), Ne(2), Ne(3)) )
    ele%T = 330D0
    ele%Mat = 1
    ele%Ediff = 0D0

    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)

                tmpI1 = ele(i, j, k)%Mat
                ele(i, j, k)%BD(:, 1) = DBLE( (/ i-1, i /) ) * dL(1)
                ele(i, j, k)%BD(:, 2) = DBLE( (/ j-1, j /) ) * dL(2)
                ele(i, j, k)%BD(:, 3) = DBLE( (/ k-1, k /) ) * dL(3)
                CALL energy( tmpI1, ele(i, j, k)%T, tmpR1 )
                ele(i, j, k)%E = tmpR1
                CALL Etable( tmpI1, 2, tmpR1, ele(i, j, k)%N )
                ele(i, j, k)%Eph = ele(i, j, k)%E / ele(i, j, k)%N
                CALL Etable( tmpI1, 4, tmpR1, ele(i, j, k)%Vph )
                CALL Etable( tmpI1, 5, tmpR1, ele(i, j, k)%MFP )

                N(i, j, k) = INT( ele(i, j, k)%N * dV + 0.5D0 )

            ENDDO
        ENDDO
    ENDDO

    Npercell = 400

    DEALLOCATE( N )


END SUBROUTINE initialize


!======================================================================
!----------------------------------------------------------------------
! Read material tables.
!----------------------------------------------------------------------
SUBROUTINE ReadTables
IMPLICIT NONE

    OPEN( UNIT = 110, FILE = "../Material_Tables/Ge_real_table.txt" )
    OPEN( UNIT = 120, FILE = "../Material_Tables/Si_real_table.txt" )

    READ(110,*) N_Ge1, N_Ge2
    READ(120,*) N_Si1, N_Si2

    ALLOCATE( Ge_table(N_Ge1, N_Ge2) )
    ALLOCATE( Si_table(N_Si1, N_Si2) )
    READ(110,*) Ge_table
    READ(120,*) Si_table

    CLOSE( 110 )
    CLOSE( 120 )

    Ge_start = Ge_table(3, 1)
    dU_Ge = Ge_table(3, 2) - Ge_table(3, 1)
    Si_start = Si_table(3, 1)
    dU_Si = Si_table(3, 2) - Si_table(3, 1)

END SUBROUTINE ReadTables


!======================================================================
END MODULE INIT