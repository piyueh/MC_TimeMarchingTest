!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_ADV
!! Contenent: Collection of subroutines related to phonon advancement.
!! Author: PY Chuang
!!#####################################################################
MODULE ADV
USE RNG
IMPLICIT NONE
CONTAINS


!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Random Free Flight Time
!----------------------------------------------------------------------
SUBROUTINE phn_adv_RT( phm, RSeed, dtRemain )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele
USE VAR_BC, ONLY: BCs
USE VAR_Others, ONLY: GammaT
USE ROUTINES, ONLY: FreeFlightTime
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(INOUT):: dtRemain
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: dt(3), R
INTEGER(KIND=4):: hit, Loc(1)
INTEGER(KIND=4):: I1, I2, I3
LOGICAL:: TF

    dt(1) = dtRemain
    CALL FreeFlightTime( dt(2), Rseed )
    CALL hit_ElementSurface( phm, dt(3), hit )

    Loc = MINLOC( dt )

    SELECTCASE( Loc(1) )
    CASE(1)
        phm%xyz = phm%xyz + dtRemain * phm%Vxyz
        dtRemain = 0D0
    CASE(2)
        I1 = phm%eID(1)
        I2 = phm%eID(2)
        I3 = phm%eID(3)
        phm%xyz = phm%xyz + dt(2) * phm%Vxyz
        dtRemain = dtRemain - dt(2)

        CALL RAN_NUM( Rseed, R )
        R = R * GammaT
        IF ( R.le.ele(I1, I2, I3)%SCR ) &
                    CALL intrinsic_scattering( phm, Rseed, I1, I2, I3 )
    CASE(3)
        I1 = phm%eID(1)
        I2 = phm%eID(2)
        I3 = phm%eID(3)
        phm%xyz = phm%xyz + dt(3) * phm%Vxyz
        phm%xyz(IABS( hit )) = &
                 ele(I1, I2, I3)%BD((ISIGN(1, hit)+1)/2+1, IABS( hit ))
        dtRemain = dtRemain - dt(3)

        CALL touchDomainBC( phm, hit, TF )

        SELECTCASE( TF )
        CASE(.TRUE.)
            SELECTCASE( BCs(IABS( hit )) )
            CASE(1)
                CALL Diffused( phm, hit, Rseed )
            CASE(2)
                CALL Perodic( phm, hit )
            CASE(3)
                CALL OutDoamin( phm, dtRemain, hit )
            END SELECT
        CASE(.FALSE.)
            phm%eID(IABS( hit )) = phm%eID(IABS( hit )) + ISIGN(1, hit)
        ENDSELECT
    ENDSELECT

END SUBROUTINE phn_adv_RT


!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Constant Free Flight Time
!----------------------------------------------------------------------
SUBROUTINE phn_adv_CT( phm, RSeed, dtRemain )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele
USE VAR_BC, ONLY: BCs
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(INOUT):: dtRemain
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: dtUsed
INTEGER(KIND=4):: hit
LOGICAL:: TF

    CALL hit_ElementSurface( phm, dtUsed, hit )

    IF ( dtUsed.gt.dtRemain ) THEN

        dtUsed = dtRemain
        dtRemain = 0D0
        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        CALL Adjust_scattering( phm, dtUsed, Rseed, TF )

    ELSEIF ( dtUsed.le.dtRemain ) THEN

        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        phm%xyz(IABS( hit )) = ele(phm%eID(1), &
                                   phm%eID(2), &
                                   phm%eID(3)) &
                                   %BD((ISIGN(1, hit)+1)/2+1, &
                                       IABS( hit ))
        dtRemain = dtRemain - dtUsed

        CALL Adjust_scattering( phm, dtUsed, Rseed, TF )
        ! TF =
        !   true: scattering occured. The phonon may or may not
        !         continue going through the element surface.
        !  false: scattering didn't occured. The phonon must
        !         continue going through the element surface.

        CALL still_GoThrough( TF, hit, phm%Vxyz(IABS( hit )) )
        ! TF =
        !   true: the phonon continues going through element surface
        !  false: otherwise

        IF ( TF ) THEN
            CALL touchDomainBC( phm, hit, TF )

            SELECTCASE( TF )
            CASE(.TRUE.)
                SELECTCASE( BCs(IABS( hit )) )
                CASE(1)
                    CALL Diffused( phm, hit, Rseed )
                CASE(2)
                    CALL Perodic( phm, hit )
                CASE(3)
                    CALL OutDoamin( phm, dtRemain, hit )
                END SELECT
            CASE(.FALSE.)
                phm%eID(IABS( hit )) = phm%eID(IABS( hit )) + &
                                       ISIGN(1, hit)
            ENDSELECT
        ENDIF

    ENDIF

END SUBROUTINE phn_adv_CT


!======================================================================
!----------------------------------------------------------------------
! Adjust the element surface the phonon will hit first.  And calculate
! the time it needs.
!----------------------------------------------------------------------
SUBROUTINE hit_ElementSurface( phm, dtUsed, hit )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele
IMPLICIT NONE
TYPE(Phonon), INTENT(IN):: phm
REAL(KIND=8), INTENT(OUT):: dtUsed
INTEGER(KIND=4), INTENT(OUT):: hit
REAL(KIND=8):: ds(3)
INTEGER(KIND=4):: i, tmpI1, tmpI2, tmpI3, Loc(1), face(3)

    tmpI1 = phm%eID(1)
    tmpI2 = phm%eID(2)
    tmpI3 = phm%eID(3)

    DO i = 1, 3
        IF ( phm%Vxyz(i).gt.0 ) THEN
            ds(i) = (ele(tmpI1, tmpI2, tmpI3)%BD(2, i) - phm%xyz(i)) /&
                    phm%Vxyz(i)
            face(i) = i
        ELSEIF ( phm%Vxyz(i).lt.0 ) THEN
            ds(i) = (ele(tmpI1, tmpI2, tmpI3)%BD(1, i) - phm%xyz(i)) /&
                    phm%Vxyz(i)
            face(i) = -i
        ELSE
            ds(i) = 1D8
            face(i) = 0
        ENDIF
    ENDDO

    Loc = MINLOC( ds )
    dtUsed = ds(Loc(1))
    hit = face(Loc(1))

END SUBROUTINE hit_ElementSurface


!======================================================================
!----------------------------------------------------------------------
! Adjust Whether Intrinsic Scattering Occurs
!----------------------------------------------------------------------
SUBROUTINE Adjust_scattering( phm, dt, Rseed, TF )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8), INTENT(IN):: dt
REAL(KIND=8):: prob, R1
INTEGER(KIND=4):: tmpI1, tmpI2, tmpI3
LOGICAL:: TF

    tmpI1 = phm%eID(1)
    tmpI2 = phm%eID(2)
    tmpI3 = phm%eID(3)

    prob = 1D0 - DEXP( -dt * phm%V / ele(tmpI1, tmpI2, tmpI3)%MFP )

    CALL RAN_NUM( Rseed, R1 )

    IF ( R1.le.prob ) THEN
        TF = .TRUE.
        CALL intrinsic_scattering( phm, Rseed, tmpI1, tmpI2, tmpI3 )
    ELSE
        TF = .FALSE.
    ENDIF


END SUBROUTINE Adjust_scattering
!======================================================================
!----------------------------------------------------------------------
! Phonon Intrinsic Scattering
!----------------------------------------------------------------------
SUBROUTINE intrinsic_scattering( phm, Rseed, I1, I2, I3 )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele
USE VAR_Others, ONLY: M_PI
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: R1, R2
INTEGER(KIND=4):: I1, I2, I3

    CALL RAN_NUM( Rseed, R1 )
    CALL RAN_NUM( Rseed, R2 )

    R1 = 2D0 * R1 - 1D0
    R2 = R2 * M_PI * 2D0
    !--------------------------------------------------------------
    ! R1 = COS(theta) => SIN(theta) = (1-R1**2)**0.5
    ! R2 = phi
    ! Vx = V * COS(theta)
    ! Vy = V * SIN(theta) * COS(phi)
    ! Vz = V * SIN(theta) * SIN(phi)
    !--------------------------------------------------------------

    ele(I1, I2, I3)%Ediff = ele(I1, I2, I3)%Ediff + &
                                        phm%E - ele(I1, I2, I3)%Eph
    phm%V = ele(I1, I2, I3)%Vph
    phm%Vxyz(1) = phm%V * R1
    phm%Vxyz(2:3) = phm%V * DSQRT(1D0 - R1**2) * &
                                        (/ DCOS( R2 ), DSIN( R2 ) /)
    phm%E = ele(I1, I2, I3)%Eph
    phm%Mat = ele(I1, I2, I3)%Mat

END SUBROUTINE intrinsic_scattering


!======================================================================
!----------------------------------------------------------------------
! Adjust whether the phonon still goes through the element surface.
!----------------------------------------------------------------------
SUBROUTINE still_GoThrough( TF, hit, V )
IMPLICIT NONE
LOGICAL, INTENT(INOUT):: TF
INTEGER(KIND=4), INTENT(IN):: hit
REAL(KIND=8), INTENT(IN):: V

    IF ( TF ) THEN
        IF ( (V * hit).le.0 ) TF = .FALSE.
    ELSE
        TF = .TRUE.
    ENDIF

END SUBROUTINE still_GoThrough


!======================================================================
!----------------------------------------------------------------------
! Out Domain
!----------------------------------------------------------------------
SUBROUTINE touchDomainBC( phm, hit, TF )
USE VAR_TYPES
USE VAR_SPACES, ONLY: Ne
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit
LOGICAL, INTENT(OUT):: TF
INTEGER(KIND=4):: tmpI1

    TF = .FALSE.

    tmpI1 = phm%eID(IABS( hit )) + ISIGN(1, hit)

    IF ( (tmpI1.gt.Ne(IABS( hit ))) .OR. (tmpI1.lt.1) ) TF = .TRUE.

END SUBROUTINE touchDomainBC


!======================================================================
!----------------------------------------------------------------------
! Diffused Response on domain boundaries
!----------------------------------------------------------------------
SUBROUTINE Diffused( phm, hit, Rseed )
USE VAR_TYPES
USE VAR_Others, ONLY: M_PI
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit
TYPE(rng_t), INTENT(INOUT):: Rseed
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: tmpI1

    !------------------------------------------------------------------
    ! V = (+ or -) (r1**0.5) * n +
    !     [(1 - r1)**0.5] * COS(2 * pi * r2) * t1 +
    !     [(1 - r1)**0.5] * SIN(2 * pi * r2) * t1
    !
    ! 1. r1 & r2 are uniform random numbers in [0, 1].
    ! 2. n is normal vector of element face
    ! 3. t1 & t2 are two unit vectors located on the element surface
    !    and perpendicular to each other
    ! 4. '+' or '-' depend on reflection of transmission occurs
    ! 5. For example, in cubic elements, if hit = 1,
    !    then the n = ( 1, 0, 0); t1 = (0, 1, 0); t2 = (0, 0, 1)
    !    and use '-' for reflection. Therefore:
    !
    !    (Vx, Vy, Vz) = -r1^0.5 * (1, 0, 0) +
    !                   (1-r1)^0.5 * cos(2pi*r2) * (0, 1, 0) +
    !                   (1-r1)^0.5 * sin(2pi*r2) * (0, 0, 1)
    !--------------------------------------------------------------

    CALL RAN_NUM( Rseed, tmpR1 )
    CALL RAN_NUM( Rseed, tmpR2 )

    tmpI1 = ISIGN( 1, hit )
    SELECTCASE( tmpI1)
    CASE(1)
        tmpR1 = - DSQRT( tmpR1 )
    CASE(-1)
        tmpR1 = DSQRT( tmpR1 )
    ENDSELECT
    tmpR2 = tmpR2 * M_PI * 2D0

    !------------------------------------------------------------------
    ! tmpR1 = (+ or -) r1^0.5
    ! tmpR2 = 2 * pi * r2
    !------------------------------------------------------------------

    tmpI1 = IABS( hit )
    SELECTCASE( tmpI1 )
    CASE(1)
        phm%Vxyz(1) = phm%V * tmpR1
        phm%Vxyz(2) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
        phm%Vxyz(3) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
    CASE(2)
        phm%Vxyz(1) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
        phm%Vxyz(2) = phm%V * tmpR1
        phm%Vxyz(3) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
    CASE(3)
        phm%Vxyz(3) = phm%V * tmpR1
        phm%Vxyz(2) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
        phm%Vxyz(1) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
    ENDSELECT

END SUBROUTINE Diffused


!======================================================================
!----------------------------------------------------------------------
! Perodic response on boundary
!----------------------------------------------------------------------
SUBROUTINE Perodic( phm, hit )
USE VAR_TYPES
USE VAR_SPACES, ONLY: L, Ne
USE ROUTINES, ONLY: Errors
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit
INTEGER(KIND=4):: tmpI1, tmpI2

    tmpI1 = IABS( hit )
    tmpI2 = ISIGN( 1, hit )

    SELECTCASE( tmpI2 )
    CASE(1)
        IF ( phm%eID(tmpI1).ne.Ne(tmpI1) ) CALL Errors(9)
        phm%xyz(tmpI1) = 0D0
        phm%eID(tmpI1) = 1
    CASE(-1)
        IF ( phm%eID(tmpI1).ne.1 ) CALL Errors(10)
        phm%xyz(tmpI1) = L(tmpI1)
        phm%eID(tmpI1) = Ne(tmpI1)
    END SELECT


END SUBROUTINE Perodic


!======================================================================
!----------------------------------------------------------------------
! Create or delete phonons in elements for energy conservation
!----------------------------------------------------------------------
SUBROUTINE CreateDelete( Rseed )
USE VAR_TYPES
USE VAR_Others, ONLY: M_PI
USE VAR_SPACES, ONLY: NAdd, Ne, ele, dL
USE VAR_ph, ONLY: phID, phn, NEmpty, EmptyID
USE ROUTINES
IMPLICIT NONE
TYPE(rng_t), INTENT(INOUT):: Rseed
REAL(KIND=8):: R1, rannum
REAL(KIND=8), ALLOCATABLE:: R(:, :)
INTEGER(KIND=4):: i, j, k, s1, s2, m
INTEGER(KIND=4):: NAddTol

    !------------------------------------------------------------------
    ! Must be executed after calling SUBROUTINE Reorder_CellInfo
    ! Theredore, ele%E = (summation of energy of phonons inside this
    !                    element) + ele%Ediff
    !------------------------------------------------------------------

    NAdd = 0

    DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)

        R1 = ele(i, j, k)%Eph * 0.5D0

        IF ( ele(i, j, k)%Ediff.gt.R1 ) THEN

            NAdd(i, j, k) = INT( ele(i, j, k)%Ediff / &
                                             ele(i, j, k)%Eph + 0.5D0 )

            ele(i, j, k)%Ediff = ele(i, j, k)%Ediff - &
                                       NAdd(i, j, k) * ele(i, j, k)%Eph

        ELSEIF ( ele(i, j, k)%Ediff.lt.(-R1) ) THEN

            DO WHILE ( ele(i, j, k)%Ediff.lt.(-R1) )

                CALL RAN_NUM( Rseed, rannum )

                m = MIN( INT( rannum * ele(i, j, k)%Ntol + 0.5D0 ) + &
                                   ele(i, j, k)%Nbg, ele(i, j, k)%Ned )

                IF ( phn(phId(m))%Exist ) THEN

                    ele(i, j, k)%Ediff = ele(i, j, k)%Ediff + &
                                                         phn(phId(m))%E

                    phn(phID(m))%Exist = .FALSE.
                    NEmpty = NEmpty + 1
                    EmptyID(NEmpty) = phId(m)

                ENDIF

            ENDDO

        ENDIF

    ENDDO; ENDDO; ENDDO

    NAddTol = SUM( NAdd )

    IF ( NAddTol.le.NEmpty ) THEN
        s1 = 1
        s2 = 0
        DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)
            IF ( NAdd(i, j, k).gt.0 ) THEN

                s2 = s1 + NAdd(i, j, k) - 1

                ALLOCATE( R(NAdd(i, j, k), 5) )

                R(:, 4) = 2D0 * R(:, 4) - 1D0
                R(:, 5) = R(:, 5) * M_PI * 2D0

                phn(EmptyID(s1:s2))%xyz(1) = R(:, 1) * dL(1) + &
                                             ele(i, j, k)%BD(1, 1)
                phn(EmptyID(s1:s2))%xyz(2) = R(:, 2) * dL(2) + &
                                             ele(i, j, k)%BD(1, 2)
                phn(EmptyID(s1:s2))%xyz(3) = R(:, 3) * dL(3) + &
                                             ele(i, j, k)%BD(1, 3)

                phn(EmptyID(s1:s2))%V = ele(i, j, k)%Vph

                phn(EmptyID(s1:s2))%Vxyz(1) = &
                                        phn(EmptyID(s1:s2))%V * R(:, 4)
                phn(EmptyID(s1:s2))%Vxyz(2) = phn(EmptyID(s1:s2))%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DCOS( R(:, 5) )
                phn(EmptyID(s1:s2))%Vxyz(3) = phn(EmptyID(s1:s2))%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DSIN( R(:, 5) )

                phn(EmptyID(s1:s2))%E = ele(i, j, k)%Eph

                phn(EmptyID(s1:s2))%Org(1) = phn(EmptyID(s1:s2))%xyz(1)
                phn(EmptyID(s1:s2))%Org(2) = phn(EmptyID(s1:s2))%xyz(2)
                phn(EmptyID(s1:s2))%Org(3) = phn(EmptyID(s1:s2))%xyz(3)

                phn(EmptyID(s1:s2))%Dis(1) = 0D0
                phn(EmptyID(s1:s2))%Dis(2) = 0D0
                phn(EmptyID(s1:s2))%Dis(3) = 0D0

                phn(EmptyID(s1:s2))%Mat = ele(i, j, k)%Mat

                phn(EmptyID(s1:s2))%eID(1) = i
                phn(EmptyID(s1:s2))%eID(2) = j
                phn(EmptyID(s1:s2))%eID(3) = k

                phn(EmptyID(s1:s2))%Exist = .TRUE.

                DEALLOCATE( R )

                s1 = s2 + 1

            ENDIF
        ENDDO; ENDDO; ENDDO

        IF ( s2.gt.NEmpty ) CALL Errors(7)

    ELSE
        CALL Errors(8)
    ENDIF


END SUBROUTINE CreateDelete


!======================================================================
!----------------------------------------------------------------------
! Phonon encounter heat control BC
!----------------------------------------------------------------------
SUBROUTINE OutDoamin( phm, dtRemain, hit )
USE VAR_TYPES
USE VAR_BC, ONLY: PoolSize, PoolL, PoolR, iNPoolL, iNPoolR, qL, qR
IMPLICIT NONE
REAL(KIND=8), INTENT(INOUT):: dtRemain
INTEGER(KIND=4), INTENT(IN):: hit
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4):: I1, I2, I3

    I2 = phm%eID(2)
    I3 = phm%eID(3)

    SELECTCASE( hit )
    CASE(1)
        IF ( iNPoolL(I2, I3).eq.PoolSize ) iNPoolL(I2, I3) = 0
        iNPoolL(I2, I3) = iNPoolL(I2, I3) + 1
        I1 = iNPoolL(I2, I3)
        PoolL(I1, I2, I3)%y = phm%xyz(2)
        PoolL(I1, I2, I3)%z = phm%xyz(3)
        PoolL(I1, I2, I3)%direction(1) = phm%Vxyz(1) / phm%V
        PoolL(I1, I2, I3)%direction(2) = phm%Vxyz(2) / phm%V
        PoolL(I1, I2, I3)%direction(3) = phm%Vxyz(3) / phm%V
        PoolL(I1, I2, I3)%dtRemain = dtRemain
        PoolL(I1, I2, I3)%Mat = phm%Mat
        qL(I2, I3) = qL(I2, I3) + phm%E
    CASE(-1)
        IF ( iNPoolR(I2, I3).eq.PoolSize ) iNPoolR(I2, I3) = 0
        iNPoolR(I2, I3) = iNPoolR(I2, I3) + 1
        I1 = iNPoolR(I2, I3)
        PoolR(I1, I2, I3)%y = phm%xyz(2)
        PoolR(I1, I2, I3)%z = phm%xyz(3)
        PoolR(I1, I2, I3)%direction(1) = phm%Vxyz(1) / phm%V
        PoolR(I1, I2, I3)%direction(2) = phm%Vxyz(2) / phm%V
        PoolR(I1, I2, I3)%direction(3) = phm%Vxyz(3) / phm%V
        PoolR(I1, I2, I3)%dtRemain = dtRemain
        PoolR(I1, I2, I3)%Mat = phm%Mat
        qR(I2, I3) = qR(I2, I3) + phm%E
    END SELECT

    phm%Exist = .FALSE.
    dtRemain = 0D0

END SUBROUTINE OutDoamin

!======================================================================
END MODULE ADV