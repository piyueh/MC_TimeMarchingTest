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
INTEGER(KIND=4):: Npercell, bgN, edN
INTEGER(KIND=4):: tmpI1, Loc(3), i, j, k
REAL(KIND=8), ALLOCATABLE:: R(:, :)

    time = 0D0
    iter0 = 0

    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( tmpR1 )

    ALLOCATE( SeedMP(NCores) )
    DO i = 1, NCores
        CALL RNG_SEED( SeedMP(i), INT( tmpR1*1D8 ) )
    ENDDO

    L = (/ 1D3, 1D1, 1D1 /)
    Ne = (/ 100, 4, 4 /)
    dL = L / DBLE( Ne )
    dV = dL(1) * dL(2) *dL(3)

    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( NAdd(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( Tavg(Ne(1), Ne(2), Ne(3)) )
    ALLOCATE( Eavg(Ne(1), Ne(2), Ne(3)) )
    
    ele(1:50, :, :)%T = 200D0
    ele(51:100, :, :)%T = 400D0
    ele%Mat = 2
    ele%Ediff = 0D0
    ele%E = 0D0

    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)

        tmpI1 = ele(i, j, k)%Mat
        ele(i, j, k)%BD(:, 1) = DBLE( (/ i-1, i /) ) * dL(1)
        ele(i, j, k)%BD(:, 2) = DBLE( (/ j-1, j /) ) * dL(2)
        ele(i, j, k)%BD(:, 3) = DBLE( (/ k-1, k /) ) * dL(3)
        CALL energy( tmpI1, ele(i, j, k)%T, tmpR1 )
        CALL Etable( tmpI1, 2, tmpR1, ele(i, j, k)%ND )
        ele(i, j, k)%Eph = tmpR1 / ele(i, j, k)%ND
        CALL Etable( tmpI1, 4, tmpR1, ele(i, j, k)%Vph )
        CALL Etable( tmpI1, 5, tmpR1, ele(i, j, k)%MFP )

        ele(i, j, k)%Ntol = INT( ele(i, j, k)%ND * dV + 0.5D0 )

    ENDDO; ENDDO; ENDDO

    Npercell = 8000 / 16
    Loc = MINLOC( ele%Ntol )
    bundle = DBLE( ele(Loc(1), Loc(2), Loc(3))%Ntol ) / DBLE( Npercell )
    ele%Eph = ele%Eph * bundle
    ele%Ntol = INT( DBLE( ele%Ntol ) / bundle + 0.5D0 )
    
    !dt = MINVAL( dL ) / MAXVAL( ele%Vph ) / 2D0
    dt = 1D0
    
    RNph = SUM( ele%Ntol )
    FNph = RNph * 1.2D0
    
    ALLOCATE( phn(FNph) )
    ALLOCATE( phId(FNph) )
    ALLOCATE( EmptyID(FNph) )
    phn%Exist = .FALSE.
    EmptyID = -1
    phId = -1
    
    tmpI1 = 0
    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)
     
        ele(i, j, k)%Nbg = tmpI1 + 1
        ele(i, j, k)%Ned = tmpI1 + ele(i, j, k)%Ntol
        bgN = ele(i, j, k)%Nbg
        edN = ele(i, j, k)%Ned
     
        ALLOCATE( R(ele(i, j, k)%Ntol, 5) )
        CALL RAN_NUM( SeedMP(1), R )
                
        phn(bgN:edN)%xyz(1) = R(:, 1) * dL(1) + ele(i, j, k)%BD(1, 1)
        phn(bgN:edN)%xyz(2) = R(:, 2) * dL(2) + ele(i, j, k)%BD(1, 2)
        phn(bgN:edN)%xyz(3) = R(:, 3) * dL(3) + ele(i, j, k)%BD(1, 3)
                           
        phn(bgN:edN)%V = ele(i, j, k)%Vph
        phn(bgN:edN)%Vxyz(1) = phn(bgN:edN)%V * DCOS( R(:, 4) * M_PI )
        phn(bgN:edN)%Vxyz(2) = phn(bgN:edN)%V * &
                  DSIN( R(:, 4) * M_PI ) * DCOS( R(:, 5) * M_PI * 2D0 )
        phn(bgN:edN)%Vxyz(3) = phn(bgN:edN)%V * &
                  DSIN( R(:, 4) * M_PI ) * DSIN( R(:, 5) * M_PI * 2D0 )
                
        phn(bgN:edN)%E = ele(i, j, k)%Eph
                
        phn(bgN:edN)%Org(1) = phn(bgN:edN)%xyz(1)
        phn(bgN:edN)%Org(2) = phn(bgN:edN)%xyz(2)
        phn(bgN:edN)%Org(3) = phn(bgN:edN)%xyz(3)
                
        phn(bgN:edN)%Dis(1) = 0D0
        phn(bgN:edN)%Dis(2) = 0D0
        phn(bgN:edN)%Dis(3) = 0D0
                
        phn(bgN:edN)%Mat = ele(i, j, k)%Mat
                
        phn(bgN:edN)%eID(1) = i
        phn(bgN:edN)%eID(2) = j
        phn(bgN:edN)%eID(3) = k
        
        phn(bgN:edN)%Exist = .TRUE.
        
        DEALLOCATE( R )
                
        CALL energy( ele(i, j, k)%Mat, ele(i, j, k)%T, tmpR1 )
        ele(i, j, k)%E = SUM( phn(bgN:edN)%E )
        ele(i, j, k)%Ediff = tmpR1 * dV - ele(i, j, k)%E
                
        tmpI1 = tmpI1 + ele(i, j, k)%Ntol
        
    ENDDO; ENDDO; ENDDO
    
    FORALL( i=1:RNph ) phId(i) = i
    
    IF ( tmpI1.ne.RNph ) CALL Errors( 1 )
    IF ( ele(Ne(1), Ne(2), Ne(3))%Ned.ne.RNph ) CALL Errors( 2 )
    
    WRITE(*, *) "Total number of phonons: ", RNph
    WRITE(*, *) "The size of phonon array: ", FNph
    WRITE(*, *) "Bundle: ", bundle
    WRITE(*, *) "The Energy of Phonon Bundle: ", MINVAL( ele%Eph )
    WRITE(*, *) "dt: ", dt

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