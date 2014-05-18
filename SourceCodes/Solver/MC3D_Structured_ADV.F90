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
USE VAR_SPACES, ONLY: ele, dV
USE VAR_BC, ONLY: BCs
USE VAR_Others, ONLY: GammaT
USE ROUTINES, ONLY: FreeFlightTime
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(INOUT):: dtRemain
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: dt(3), R
INTEGER(KIND=4):: hit1, hit2, Loc(1)
INTEGER(KIND=4):: OldID(3), I1, I2, I3
LOGICAL:: TF

    dt(1) = dtRemain
    CALL FreeFlightTime( dt(2), Rseed )
    CALL hit_ElementSurface( phm, dt(3), hit1, hit2 )

    Loc = MINLOC( dt )

    SELECTCASE( Loc(1) )
    CASE(1)
        phm%xyz = phm%xyz + dtRemain * phm%Vxyz
        dtRemain = 0D0
    CASE(2)
        phm%xyz = phm%xyz + dt(2) * phm%Vxyz
        dtRemain = dtRemain - dt(2)

        CALL RAN_NUM( Rseed, R )
        R = R * GammaT
        IF ( R.le.ele(phm%eID(1), phm%eID(2), phm%eID(3))%SCR ) &
                    CALL intrinsic_scattering( phm, Rseed, &
                                   phm%eID(1), phm%eID(2), phm%eID(3) )
    CASE(3)
    
        OldID = phm%eID
        phm%xyz = phm%xyz + dt(3) * phm%Vxyz
        phm%xyz(hit2) = ele(OldID(1), OldID(2), OldID(3)) &
                                                %BD((hit1+1)/2+1, hit2)
        phm%eID(hit2) = phm%eID(hit2) + hit1
        dtRemain = dtRemain - dt(3)

        CALL touchDomainBC( phm, hit2, TF )

        SELECTCASE( TF )
        CASE(.TRUE.)
            SELECTCASE( BCs(hit2) )
            CASE(1)
                
                phm%eID(hit2) = OldID(hit2)
                I1 = phm%eID(1)
                I2 = phm%eID(2)
                I3 = phm%eID(3)
                ele(I1, I2, I3)%Ediff = ele(I1, I2, I3)%Ediff + &
                                        phm%E - &
                                        ele(I1, I2, I3)%Eph
                phm%V = ele(I1, I2, I3)%Vph
                CALL Diffused( phm, hit1, hit2, Rseed, -1 )
                phm%E = ele(I1, I2, I3)%Eph
                phm%Mat = ele(I1, I2, I3)%Mat
            CASE(2)
                CALL Perodic( phm, hit1, hit2 )
            CASE(3)
                CALL OutDoamin( phm, dtRemain, hit1 )
            END SELECT
        CASE(.FALSE.)
            CALL Mat_Interface( phm, hit1, hit2, OldID,           &
                                   ele(OldID(1), OldID(2), OldID(3)), &
                             ele(phm%eID(1), phm%eID(2), phm%eID(3)), &
                                                            dV, RSeed )
        ENDSELECT
    ENDSELECT

END SUBROUTINE phn_adv_RT


!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Constant Free Flight Time
!----------------------------------------------------------------------
SUBROUTINE phn_adv_CT( phm, RSeed, dtRemain, flag )
USE VAR_TYPES
USE VAR_SPACES, ONLY: ele, dV
USE VAR_BC, ONLY: BCs
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(INOUT):: dtRemain
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: dtUsed
INTEGER(KIND=4):: hit1, hit2, OldID(3), I1, I2, I3, flag
LOGICAL:: TF
TYPE(Phonon):: OldPhn

    OldPhn = phm
    
    CALL hit_ElementSurface( phm, dtUsed, hit1, hit2 )

    IF ( dtUsed.gt.dtRemain ) THEN

        phm%xyz = phm%xyz + dtRemain * phm%Vxyz
        CALL Adjust_scattering( phm, dtRemain, Rseed, TF )
        dtRemain = 0D0

    ELSEIF ( dtUsed.le.dtRemain ) THEN
        
        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        phm%xyz(hit2) = ele(phm%eID(1), phm%eID(2), phm%eID(3)) &
                                                %BD((hit1+1)/2+1, hit2)
        dtRemain = dtRemain - dtUsed

        CALL Adjust_scattering( phm, dtUsed, Rseed, TF )
        ! TF =
        !   true: scattering occured. The phonon may or may not
        !         continue going through the element surface.
        !  false: scattering didn't occured. The phonon must
        !         continue going through the element surface.

        CALL still_GoThrough( TF, hit1, phm%Vxyz(hit2) )
        ! TF =
        !   true: the phonon continues going through element surface
        !  false: otherwise

        IF ( TF ) THEN
            
            OldID = phm%eID
            phm%eID(hit2) = phm%eID(hit2) + hit1
            CALL touchDomainBC( phm, hit2, TF )

            SELECTCASE( TF )
            CASE(.TRUE.)
                SELECTCASE( BCs(hit2) )
                CASE(1)
                    phm%eID(hit2) = OldID(hit2)
                    I1 = phm%eID(1)
                    I2 = phm%eID(2)
                    I3 = phm%eID(3)
                    ele(I1, I2, I3)%Ediff = ele(I1, I2, I3)%Ediff + &
                                            phm%E - &
                                            ele(I1, I2, I3)%Eph
                    phm%V = ele(I1, I2, I3)%Vph
                    CALL Diffused( phm, hit1, hit2, Rseed, -1 )
                    phm%E = ele(I1, I2, I3)%Eph
                    phm%Mat = ele(I1, I2, I3)%Mat
                CASE(2)
                    CALL Perodic( phm, hit1, hit2 )
                CASE(3)
                    ! Only possible when hit2 = 1
                    SELECTCASE(flag)
                    CASE(1)
                        CALL OutDoamin( phm, dtRemain, hit1 )
                    CASE(2)
                        phm = OldPhn
                        dtRemain = dtRemain + dtUsed
                    END SELECT
                END SELECT
            CASE(.FALSE.)
                CALL Mat_Interface( phm, hit1, hit2, OldID,           &
                                   ele(OldID(1), OldID(2), OldID(3)), &
                             ele(phm%eID(1), phm%eID(2), phm%eID(3)), &
                                                            dV, RSeed )
            ENDSELECT
        ENDIF

    ENDIF

END SUBROUTINE phn_adv_CT


!======================================================================
!----------------------------------------------------------------------
! Adjust the element surface the phonon will hit first.  And calculate
! the time it needs.
!----------------------------------------------------------------------
SUBROUTINE hit_ElementSurface( phm, elm, dtUsed, hit1, hit2 )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(IN):: phm
TYPE(Element), INTENT(IN):: elm
REAL(KIND=8), INTENT(OUT):: dtUsed
INTEGER(KIND=4), INTENT(OUT):: hit1, hit2
REAL(KIND=8):: ds(3)
INTEGER(KIND=4):: i, Loc(1), dir(3)

    DO i = 1, 3
        IF ( phm%Vxyz(i).gt.0 ) THEN
            ds(i) = (elm%BD(2, i) - phm%xyz(i)) / phm%Vxyz(i)
            dir(i) = 1
        ELSEIF ( phm%Vxyz(i).lt.0 ) THEN
            ds(i) = (elm%BD(1, i) - phm%xyz(i)) / phm%Vxyz(i)
            dir(i) = -1
        ELSE
            ds(i) = 1D8
            dir(i) = 0
        ENDIF
    ENDDO

    Loc = MINLOC( ds )
    hit2 = Loc(1)
    dtUsed = ds(hit2)
    hit1 = dir(hit2)
    
END SUBROUTINE hit_ElementSurface


!======================================================================
!----------------------------------------------------------------------
! Adjust Whether Intrinsic Scattering Occurs
!----------------------------------------------------------------------
SUBROUTINE Adjust_scattering( phm, elm, dt, Rseed, TF )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(Element), INTENT(IN):: elm
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8), INTENT(IN):: dt
REAL(KIND=8):: prob, R1
LOGICAL:: TF

    prob = 1D0 - DEXP( -dt * phm%V / elm%MFP )

    CALL RAN_NUM( Rseed, R1 )

    IF ( R1.le.prob ) THEN
        TF = .TRUE.
        CALL intrinsic_scattering( phm, elm, Rseed )
    ELSE
        TF = .FALSE.
    ENDIF


END SUBROUTINE Adjust_scattering
!======================================================================
!----------------------------------------------------------------------
! Phonon Intrinsic Scattering
!----------------------------------------------------------------------
SUBROUTINE intrinsic_scattering( phm, elm, Rseed )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(Element), INTENT(IN):: elm
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8), PARAMETER:: M_PI = 3.14159265359D0
REAL(KIND=8):: R1, R2

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

    elm%Ediff = elm%Ediff + phm%E - elm%Eph
    phm%V = elm%Vph
    phm%Vxyz(1) = phm%V * R1
    phm%Vxyz(2:3) = phm%V * DSQRT(1D0 - R1**2) * &
                    (/ DCOS( R2 ), DSIN( R2 ) /)
    phm%E = elm%Eph
    phm%Mat = elm%Mat

END SUBROUTINE intrinsic_scattering


!======================================================================
!----------------------------------------------------------------------
! Adjust whether the phonon still goes through the element surface.
!----------------------------------------------------------------------
SUBROUTINE still_GoThrough( TF, dir, V )
IMPLICIT NONE
LOGICAL, INTENT(INOUT):: TF
INTEGER(KIND=4), INTENT(IN):: dir
REAL(KIND=8), INTENT(IN):: V

    IF ( TF ) THEN
        IF ( (V * dir).le.0 ) TF = .FALSE.
    ELSE
        TF = .TRUE.
    ENDIF

END SUBROUTINE still_GoThrough


!======================================================================
!----------------------------------------------------------------------
! Out Domain
!----------------------------------------------------------------------
SUBROUTINE touchDomainBC( phm, Ne, hit2, TF )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: Ne(3)
INTEGER(KIND=4), INTENT(IN):: hit2
LOGICAL, INTENT(OUT):: TF
INTEGER(KIND=4):: I1

    I1 = phm%eID(hit2)

    IF ( (I1.gt.Ne(hit2)) .OR. (I1.lt.1) ) THEN
        TF = .TRUE.
    ELSE
        TF = .FALSE.
    ENDIF

END SUBROUTINE touchDomainBC


!======================================================================
!----------------------------------------------------------------------
! Diffused Response on domain boundaries
!----------------------------------------------------------------------
SUBROUTINE Diffused( phm, hit1, hit2, Rseed, D_Type )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit1, hit2, D_Type
TYPE(rng_t), INTENT(INOUT):: Rseed
REAL(KIND=8), PARAMETER:: M_PI = 3.14159265359D0
REAL(KIND=8):: tmpR1, tmpR2

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

    tmpR1 = hit1 * D_Type * DSQRT( tmpR1 )
    tmpR2 = tmpR2 * M_PI * 2D0

    !------------------------------------------------------------------
    ! tmpR1 = (+ or -) r1^0.5
    ! tmpR2 = 2 * pi * r2
    !------------------------------------------------------------------

    SELECTCASE( hit2 )
    CASE(1)
        phm%Vxyz(1) = phm%V * tmpR1
        phm%Vxyz(2) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
        phm%Vxyz(3) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
    CASE(2)
        phm%Vxyz(1) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
        phm%Vxyz(2) = phm%V * tmpR1
        phm%Vxyz(3) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
    CASE(3)
        phm%Vxyz(1) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DCOS( tmpR2 )
        phm%Vxyz(2) = phm%V * DSQRT( 1D0 - tmpR1**2 ) * DSIN( tmpR2 )
        phm%Vxyz(3) = phm%V * tmpR1
    ENDSELECT

END SUBROUTINE Diffused


!======================================================================
!----------------------------------------------------------------------
! Perodic response on boundary
!----------------------------------------------------------------------
SUBROUTINE Perodic( phm, L, Ne, hit1, hit2 )
USE VAR_TYPES
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(IN):: L(3)
INTEGER(KIND=4), INTENT(IN):: Ne(3)
INTEGER(KIND=4), INTENT(IN):: hit1, hit2

    SELECTCASE( hit1 )
    CASE(1)
        phm%xyz(hit2) = 0D0
        phm%eID(hit2) = 1
    CASE(-1)
        phm%xyz(hit2) = L(hit2)
        phm%eID(hit2) = Ne(hit2)
    END SELECT

END SUBROUTINE Perodic


!======================================================================
!----------------------------------------------------------------------
! Phonon encounter heat control BC
!----------------------------------------------------------------------
SUBROUTINE OutDoamin( phm, dtRemain, dir )
USE ROUTINES
USE VAR_TYPES
USE VAR_BC, ONLY: PoolSize, PoolL, PoolR, iNPoolL, iNPoolR, qL, qR
IMPLICIT NONE
REAL(KIND=8), INTENT(INOUT):: dtRemain
INTEGER(KIND=4), INTENT(IN):: dir
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4):: I1, I2, I3

    I2 = phm%eID(2)
    I3 = phm%eID(3)

    SELECTCASE( dir )
    CASE(1)
        IF ( phm%Vxyz(1).le.0D0 ) THEN
            WRITE(*, *) dir
            WRITE(*, *) phm
            CALL Errors(511)
        ENDIF
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
        qR(I2, I3) = qR(I2, I3) - phm%E
    CASE(-1)
        IF ( phm%Vxyz(1).ge.0D0 ) THEN
            WRITE(*, *) dir
            WRITE(*, *) phm
            CALL Errors(512)
        ENDIF
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
        qL(I2, I3) = qL(I2, I3) - phm%E
    CASE DEFAULT
        CALL Errors(513)
    END SELECT

    phm%Exist = .FALSE.
    dtRemain = 0D0

END SUBROUTINE OutDoamin


!======================================================================
!----------------------------------------------------------------------
! Phonon hit the material interface
!----------------------------------------------------------------------
SUBROUTINE Mat_Interface( phm, hit1, hit2, OldID, elm1, elm2, dV, RSeed )
USE VAR_TYPES
USE ROUTINES, ONLY: energy, Etable
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: dV
INTEGER(KIND=4), INTENT(IN):: hit1, hit2, OldID(3)
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(Element), INTENT(IN):: elm1, elm2
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: R1
REAL(KIND=8):: tau12, U2, Vg2, U1, Vg1
INTEGER(KIND=4):: mt1, mt2
    
    mt1 = elm1%Mat
    mt2 = elm2%Mat

    IF ( mt1.ne.mt2 ) THEN
        
               
        U1 = elm1%E / dV
        Vg1 = elm1%Vph
        CALL energy( mt2, elm1%T, U2 )
        CALL Etable( mt2, 4, U2, Vg2 )
        tau12 = U2 * Vg2 / (U1 * Vg1 + U2 * Vg2)
        
        CALL RAN_NUM( RSeed, R1 )
        
        IF ( R1.le.tau12 ) THEN ! Diffused Transmission
            phm%V = elm2%Vph
            CALL Diffused( phm, hit1, hit2, Rseed, 1 )
        ELSE ! Diffused Reflection
            phm%eID(hit2) = OldID(hit2)
            phm%V = elm1%Vph
            CALL Diffused( phm, hit1, hit2, Rseed, -1 )
        ENDIF 
        
    ENDIF

END SUBROUTINE Mat_Interface


!======================================================================
END MODULE ADV