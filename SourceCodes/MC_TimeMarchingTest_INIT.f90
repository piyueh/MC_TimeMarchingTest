!!#####################################################################
!! Module Name: INIT
!! Purpose: Collection of initialization subroutines
!! Author: PY Chuang
!!#####################################################################
MODULE INIT
USE RNG
USE VAR
USE ROUTINES
IMPLICIT NONE
CONTAINS
!======================================================================


!======================================================================
!----------------------------------------------------------------------
! Initialize
!----------------------------------------------------------------------
SUBROUTINE initialize
IMPLICIT NONE
REAL(KIND=8):: tmpR1
INTEGER(KIND=4):: Npercell
INTEGER(KIND=4):: tmpI1, Loc(3), i, j, k
REAL(KIND=8), ALLOCATABLE:: rannum(:, :)

    dt = 1D-1
    iter0 = 0

    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( tmpR1 )

    ALLOCATE( SeedMP(NCores) )
    DO i = 1, NCores
        CALL RNG_SEED( SeedMP(i), INT( tmpR1*1D8 ) )
    ENDDO

    L = (/ 1D3, 1D3, 1D1 /)
    Ne = (/ 50, 50, 1 /)
    dL = L / DBLE( Ne )
    dV = dL(1) * dL(2) *dL(3)

    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( N(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( bg(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( ed(Ne(1), Ne(2), Ne(3)) )
    ele%T = 330D0
    ele%Mat = 1
    ele%Ediff = 0D0
    ele%E = 0D0

    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)

                tmpI1 = ele(i, j, k)%Mat
                ele(i, j, k)%BD(:, 1) = DBLE( (/ i-1, i /) ) * dL(1)
                ele(i, j, k)%BD(:, 2) = DBLE( (/ j-1, j /) ) * dL(2)
                ele(i, j, k)%BD(:, 3) = DBLE( (/ k-1, k /) ) * dL(3)
                CALL energy( tmpI1, ele(i, j, k)%T, tmpR1 )
                CALL Etable( tmpI1, 2, tmpR1, ele(i, j, k)%N )
                ele(i, j, k)%Eph = tmpR1 / ele(i, j, k)%N
                CALL Etable( tmpI1, 4, tmpR1, ele(i, j, k)%Vph )
                CALL Etable( tmpI1, 5, tmpR1, ele(i, j, k)%MFP )

                N(i, j, k) = INT( ele(i, j, k)%N * dV + 0.5D0 )

            ENDDO
        ENDDO
    ENDDO

    Npercell = 400
    Loc = MINLOC( N )
    bundle = DBLE( N(Loc(1), Loc(2), Loc(3)) ) / DBLE( Npercell )
    ele%Eph = ele%Eph * bundle
    N = INT( DBLE( N ) / bundle + 0.5D0 )
    Nph = SUM( N )
    
    ALLOCATE( phn(Nph) )
    ALLOCATE( phID(Nph) )
    
    tmpI1 = 0
    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)
                
                ALLOCATE( rannum(N(i, j, k), 5) )
                CALL RAN_NUM( SeedMP(1), rannum )
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(1) = &
                           rannum(:, 1) * dL(1) + ele(i, j, k)%BD(1, 1)
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(2) = &
                           rannum(:, 2) * dL(2) + ele(i, j, k)%BD(1, 2)
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(3) = &
                           rannum(:, 3) * dL(3) + ele(i, j, k)%BD(1, 3)
                           
                phn(tmpI1+1:tmpI1+N(i, j, k))%V = ele(i, j, k)%Vph
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%Vxyz(1) = &
                                    phn(tmpI1+1:tmpI1+N(i, j, k))%V * &
                                    DCOS( rannum(:, 4) * M_PI )
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%Vxyz(2) = &
                                    phn(tmpI1+1:tmpI1+N(i, j, k))%V * &
                                    DSIN( rannum(:, 4) * M_PI ) * &
                                    DCOS( rannum(:, 5) * M_PI * 2D0 )
                                    
                phn(tmpI1+1:tmpI1+N(i, j, k))%Vxyz(3) = &
                                    phn(tmpI1+1:tmpI1+N(i, j, k))%V * &
                                    DSIN( rannum(:, 4) * M_PI ) * &
                                    DSIN( rannum(:, 5) * M_PI * 2D0 )
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%E = ele(i, j, k)%Eph
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%Org(1) = &
                                      phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(1)
                phn(tmpI1+1:tmpI1+N(i, j, k))%Org(2) = &
                                      phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(2)
                phn(tmpI1+1:tmpI1+N(i, j, k))%Org(3) = &
                                      phn(tmpI1+1:tmpI1+N(i, j, k))%xyz(3)
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%Dis(1) = 0D0
                phn(tmpI1+1:tmpI1+N(i, j, k))%Dis(2) = 0D0
                phn(tmpI1+1:tmpI1+N(i, j, k))%Dis(3) = 0D0
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%Mat = ele(i, j, k)%Mat
                
                phn(tmpI1+1:tmpI1+N(i, j, k))%eID(1) = i
                phn(tmpI1+1:tmpI1+N(i, j, k))%eID(2) = j
                phn(tmpI1+1:tmpI1+N(i, j, k))%eID(3) = k
                
                DEALLOCATE( rannum )
                
                CALL energy( ele(i, j, k)%Mat, ele(i, j, k)%T, tmpR1 )
                ele(i, j, k)%E = SUM( phn(tmpI1+1:tmpI1+N(i, j, k))%E )
                ele(i, j, k)%Ediff = tmpR1 * dV - ele(i, j, k)%E
                
                bg(i, j, k) = tmpI1 + 1
                ed(i, j, k) = tmpI1 + N(i, j, k)
                
                tmpI1 = tmpI1 + N(i, j, k)
            ENDDO
        ENDDO
    ENDDO
    
    FORALL( i = 1:10 ) phID(i) = i
    
    IF ( tmpI1.ne.Nph ) CALL Errors( 1 )
    IF ( ed(Ne(1), Ne(2), Ne(3)).ne.Nph ) CALL Errors( 2 )
    
    WRITE(*, *) "Total number of phonons: ", Nph

END SUBROUTINE initialize


!======================================================================
!----------------------------------------------------------------------
! Read material tables.
!----------------------------------------------------------------------
SUBROUTINE ReadTables
IMPLICIT NONE

    OPEN( UNIT = 110, FILE = "Material_Tables/Ge_real_table.txt" )
    OPEN( UNIT = 120, FILE = "Material_Tables/Si_real_table.txt" )

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