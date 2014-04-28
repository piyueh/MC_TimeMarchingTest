!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Pre-processing
!! File Name: MC3D_Structured_PreProcessing
!! Contenent: Main Program of Pre-processing
!! Author: PY Chuang
!!#####################################################################


!======================================================================
!----------------------------------------------------------------------
! Main Program
!----------------------------------------------------------------------
PROGRAM main
USE VAR
IMPLICIT NONE
INTEGER(KIND=4):: Npercell


    CALL RANDOM_SEED()

    time = 0D0
    iter0 = 0

    WRITE(*, '("Enter CASE Name: ")', ADVANCE = 'NO')
    !READ(*, *) CaseName
    CaseName = "Adiabatic_Transien_Ge"
    InputFileName = casename(1:LEN_TRIM(casename))//'_initial.txt'

    WRITE(*, '("Enter Domain, Lx, Ly, Lz: ")', ADVANCE = 'NO')
    READ(*, *) L

    CALL Initialize_Materials

    CALL Initialize_Mesh( NperCell )

    CALL TimeStep_setup

    CALL Initialize_Phonon

    CALL Initialize_Pools

    WRITE(*, *) "Case name: ", CaseName
    WRITE(*, *) "Domain size: ", L
    WRITE(*, *) "Numbers of elements: ", Ne
    WRITE(*, *) "NperCell: ", NperCell
    WRITE(*, *) "Total number of phonons: ", RNph
    WRITE(*, *) "The size of phonon array: ", FNph
    WRITE(*, *) "Bundle: ", bundle
    WRITE(*, *) "TimeStep: ", TimeStep

    CALL Output

END PROGRAM main


!======================================================================
!----------------------------------------------------------------------
! Initialize the mesh
!----------------------------------------------------------------------
SUBROUTINE Initialize_Mesh( NperCell )
USE VAR
IMPLICIT NONE
REAL(KIND=8):: tmpR1
INTEGER(KIND=4):: i, j, k, tmpI1, Loc(3), NperCell

    WRITE(*, '("Enter Element Numbers, Nx, Ny, Nz: ")', ADVANCE = 'NO')
    READ(*, *) Ne

    WRITE(*, '("Enter NperCell: ")', ADVANCE = 'NO')
    READ(*, *) NperCell

    dL = L / DBLE( Ne )
    dV = dL(1) * dL(2) * dL(3)

    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )

    CALL Initial_Temperature

    CALL Initial_Material

    ele%Ediff = 0D0
    ele%E = 0D0
    ele%Nbg = -1
    ele%Ned = -1

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
        ele(i, j, k)%SCR = ele(i, j, k)%Vph / ele(i, j, k)%MFP

        ele(i, j, k)%Ntol = INT( ele(i, j, k)%ND * dV + 0.5D0 )

    ENDDO; ENDDO; ENDDO

    Loc = MINLOC( ele%Ntol )
    bundle = DBLE( ele(Loc(1), Loc(2), Loc(3))%Ntol ) / DBLE( NperCell )
    ele%Eph = ele%Eph * bundle
    ele%Ntol = INT( DBLE( ele%Ntol ) / bundle + 0.5D0 )
    ele%E = ele%Eph * ele%Ntol

    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)
        CALL energy( tmpI1, ele(i, j, k)%T, tmpR1 )
        ele(i, j, k)%Ediff = tmpR1 * dV - ele(i, j, k)%E
    ENDDO; ENDDO; ENDDO

    !==========================================================
    CONTAINS
    !----------------------------------------------------------
    ! Set up initial temperature field
    !----------------------------------------------------------
    SUBROUTINE Initial_Temperature
    IMPLICIT NONE

        ele(1:Ne(1)/2, :, :)%T = 2D2
        ele(Ne(1)/2:Ne(1), :, :)%T = 4D2

    END SUBROUTINE Initial_Temperature

    !----------------------------------------------------------
    ! Set up initial material distribution
    !----------------------------------------------------------
    SUBROUTINE Initial_Material
    IMPLICIT NONE

        ele%Mat = 1

    END SUBROUTINE Initial_Material
    !==========================================================

END SUBROUTINE Initialize_Mesh


!======================================================================
!----------------------------------------------------------------------
! Time step setup
!----------------------------------------------------------------------
SUBROUTINE TimeStep_setup
USE VAR
IMPLICIT NONE
!CHARACTER:: YN

    !TimeStep = MINVAL( dL ) / MAXVAL( ele%Vph ) / 2D0
    !
    !WRITE(*, '("Current TimeStep is: ", F)', ADVANCE = 'NO') TimeStep
    !WRITE(*, '("Modify TimeStep Manually? (Y/N): ")', ADVANCE = 'NO')
    !READ(*, *) YN
    !IF ( YN.eq.'Y' ) THEN
    !    WRITE(*, '("Enter the New TimeStep: ")', ADVANCE = 'NO')
    !    READ(*, *) TimeStep
    !ENDIF
    TimeStep = 1D0

END SUBROUTINE TimeStep_setup


!======================================================================
!----------------------------------------------------------------------
! Initialize the phonons
!----------------------------------------------------------------------
SUBROUTINE Initialize_Phonon
USE VAR
IMPLICIT NONE
INTEGER(KIND=4):: tmpI1, bgN, edN, i, j, k
REAL(KIND=8), ALLOCATABLE:: R(:, :)

    RNph = SUM( ele%Ntol )
    FNph = RNph * 1.2D0

    ALLOCATE( phn(FNph) )
    phn%Exist = .FALSE.

    tmpI1 = 0

    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)

        bgN = tmpI1 + 1
        edN = tmpI1 + ele(i, j, k)%Ntol

        ALLOCATE( R(ele(i, j, k)%Ntol, 5) )
        CALL RANDOM_NUMBER( R )

        !--------------------------------------------------------------
        ! R(:, 4) = COS(theta) => SIN(theta) = (1-R(:, 4)**2)**0.5
        ! R(:, 5) = phi
        ! Vx = V * COS(theta)
        ! Vy = V * SIN(theta) * COS(phi)
        ! Vz = V * SIN(theta) * SIN(phi)
        !--------------------------------------------------------------
        R(:, 4) = 2D0 * R(:, 4) - 1D0
        R(:, 5) = R(:, 5) * M_PI * 2D0


        phn(bgN:edN)%xyz(1) = R(:, 1) * dL(1) + ele(i, j, k)%BD(1, 1)
        phn(bgN:edN)%xyz(2) = R(:, 2) * dL(2) + ele(i, j, k)%BD(1, 2)
        phn(bgN:edN)%xyz(3) = R(:, 3) * dL(3) + ele(i, j, k)%BD(1, 3)

        phn(bgN:edN)%V = ele(i, j, k)%Vph

        phn(bgN:edN)%Vxyz(1) = phn(bgN:edN)%V * R(:, 4)
        phn(bgN:edN)%Vxyz(2) = phn(bgN:edN)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DCOS( R(:, 5) )
        phn(bgN:edN)%Vxyz(3) = phn(bgN:edN)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DSIN( R(:, 5) )

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

        tmpI1 = tmpI1 + ele(i, j, k)%Ntol

    ENDDO; ENDDO; ENDDO

    IF ( tmpI1.ne.RNph ) WRITE(*, *) "Errors( 1 )"


END SUBROUTINE Initialize_Phonon


!======================================================================
!----------------------------------------------------------------------
! Initialize Pools
!----------------------------------------------------------------------
SUBROUTINE Initialize_Pools
USE VAR
IMPLICIT NONE

    PoolSize = 5 * MAXVAL( ele%Ntol )

    ALLOCATE( PoolL(Ne(2), Ne(3), PoolSize) )
    ALLOCATE( PoolR(Ne(2), Ne(3), PoolSize) )
    ALLOCATE( iNPoolL(Ne(2), Ne(3)) )
    ALLOCATE( iNPoolR(Ne(2), Ne(3)) )

    PoolL%y = -1D0
    PoolL%z = -1D0
    PoolL%direction(1) = 0D0
    PoolL%direction(2) = 0D0
    PoolL%direction(3) = 0D0
    PoolL%dtRemain = 0D0
    PoolL%Mat = -1

    PoolR%y = -1D0
    PoolR%z = -1D0
    PoolR%direction(1) = 0D0
    PoolR%direction(2) = 0D0
    PoolR%direction(3) = 0D0
    PoolR%dtRemain = 0D0
    PoolR%Mat = -1

    iNPoolL = 0
    iNpoolR = 0

END SUBROUTINE Initialize_Pools



!======================================================================
!----------------------------------------------------------------------
! Output
!----------------------------------------------------------------------
SUBROUTINE output
USE VAR
IMPLICIT NONE

    WRITE(*, *)
    WRITE(*, *) "The output file name is: ", InputFileName
    WRITE(*, *) "Start to Output..."
    OPEN(UNIT = 120, FILE = InputFileName)
    WRITE(UNIT = 120, NML = Initial_1)
    WRITE(UNIT = 120, NML = Initial_2)
    CLOSE( 120 )
    WRITE(*, *) "Finished."

END SUBROUTINE output
!======================================================================
!======================================================================

