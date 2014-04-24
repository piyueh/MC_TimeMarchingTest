!!#####################################################################
!! Module Name: ADV
!! Purpose: Collection of subroutines related to phonon advancement.
!! Author: PY Chuang
!!#####################################################################
MODULE ADV
USE RNG
IMPLICIT NONE
CONTAINS
!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Controller
!----------------------------------------------------------------------
SUBROUTINE advance
USE VAR
USE ROUTINES, ONLY: Reorder_CellInfo
REAL(KIND=8):: dtRemain
INTEGER(KIND=4):: iCPU, i
    
    iCPU = 1
    
    DO i = 1, FNph
        IF ( phn(i).Exist ) THEN
            dtRemain = dt
            DO WHILE ( dtRemain.gt.0 )
                CALL phn_advance( phn(i), SeedMP(iCPU), dtRemain )
            ENDDO
        ENDIF
    ENDDO
    
    CALL Reorder_CellInfo
    CALL CreateDelete( SeedMP(iCPU) )

END SUBROUTINE advance


!======================================================================
!----------------------------------------------------------------------
! Phonon Advance Subroutine - Single Time
!----------------------------------------------------------------------
SUBROUTINE phn_advance( phm, RSeed, dtRemain )
USE VAR, ONLY: Phonon, Element, ele, BCs
USE RNG
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
        CALL intrinsic_scattering( phm, dtUsed, Rseed, TF )

    ELSEIF ( dtUsed.le.dtRemain ) THEN

        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        phm%xyz(IABS( hit )) = ele(phm%eID(1), &
                                   phm%eID(2), &
                                   phm%eID(3)) &
                                   %BD((ISIGN(1, hit)+1)/2+1, &
                                       IABS( hit ))
        dtRemain = dtRemain - dtUsed

        CALL intrinsic_scattering( phm, dtUsed, Rseed, TF )
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
            CALL outDomain( phm, hit, TF )

            SELECTCASE( TF )
            CASE(.TRUE.)
                SELECTCASE( BCs(IABS( hit )) )
                CASE(1)
                    CALL Diffused( phm, hit, Rseed )
                CASE(2)
                    CALL Perodic( phm, hit )
                END SELECT
            CASE(.FALSE.)
                phm%eID(IABS( hit )) = phm%eID(IABS( hit )) + &
                                       ISIGN(1, hit)
            ENDSELECT
        ENDIF

    ENDIF

END SUBROUTINE phn_advance


!======================================================================
!----------------------------------------------------------------------
! Adjust the element surface the phonon will hit first.  And calculate
! the time it needs.
!----------------------------------------------------------------------
SUBROUTINE hit_ElementSurface( phm, dtUsed, hit )
USE VAR, ONLY: Phonon, Element, ele
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
! Phonon Intrinsic Scattering
!----------------------------------------------------------------------
SUBROUTINE intrinsic_scattering( phm, dt, Rseed, TF )
USE VAR, ONLY: Phonon, Element, ele, M_PI
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8), INTENT(IN):: dt
REAL(KIND=8):: rannum(3), prob
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: tmpI1, tmpI2, tmpI3
LOGICAL:: TF

    TF = .FALSE.

    tmpI1 = phm%eID(1)
    tmpI2 = phm%eID(2)
    tmpI3 = phm%eID(3)

    prob = 1D0 - DEXP( -dt * phm%V / ele(tmpI1, tmpI2, tmpI3)%MFP )
    CALL RAN_NUM( Rseed, rannum )

    IF ( rannum(1).le.prob ) THEN

        tmpR1 = rannum(2) * M_PI
        tmpR2 = rannum(3) * M_PI * 2D0

        ele(tmpI1, tmpI2, tmpI3)%Ediff = &
                                   ele(tmpI1, tmpI2, tmpI3)%Ediff + &
                                   phm%E - ele(tmpI1, tmpI2, tmpI3)%Eph
        phm%V = ele(tmpI1, tmpI2, tmpI3)%Vph
        phm%Vxyz(1) = phm%V * DCOS( tmpR1 )
        phm%Vxyz(2:3) = phm%V * DSIN( tmpR1 ) * &
                                     (/ DCOS( tmpR2 ), DSIN( tmpR2 ) /)
        phm%E = ele(tmpI1, tmpI2, tmpI3)%Eph
        phm%Mat = ele(tmpI1, tmpI2, tmpI3)%Mat

        TF = .TRUE.

    ENDIF

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
SUBROUTINE outDomain( phm, hit, TF )
USE VAR, ONLY: Phonon, Ne
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit
LOGICAL, INTENT(OUT):: TF
INTEGER(KIND=4):: tmpI1

    TF = .FALSE.

    tmpI1 = phm%eID(IABS( hit )) + ISIGN(1, hit)

    IF ( (tmpI1.gt.Ne(IABS( hit ))) .OR. (tmpI1.lt.1) ) TF = .TRUE.

END SUBROUTINE outDomain


!======================================================================
!----------------------------------------------------------------------
! Diffused Response on domain boundaries
!----------------------------------------------------------------------
SUBROUTINE Diffused( phm, hit, Rseed )
USE VAR, ONLY: Phonon, M_PI
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
INTEGER(KIND=4), INTENT(IN):: hit
TYPE(rng_t), INTENT(INOUT):: Rseed
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: tmpI1

    CALL RAN_NUM( Rseed, tmpR1 )
    CALL RAN_NUM( Rseed, tmpR2 )

    tmpI1 = ISIGN( 1, hit )
    SELECTCASE( tmpI1)
    CASE(1)
        tmpR1 = (tmpR1 + 1D0) * M_PI / 2D0
    CASE(-1)
        tmpR1 = tmpR1 * M_PI / 2D0
    ENDSELECT
    tmpR2 = tmpR2 * M_PI * 2D0

    tmpI1 = IABS( hit )
    SELECTCASE( tmpI1 )
    CASE(1)
        phm%Vxyz(1) = phm%V * DCOS( tmpR1 )
        phm%Vxyz(2) = phm%V * DSIN( tmpR1 ) * DCOS( tmpR2 )
        phm%Vxyz(3) = phm%V * DSIN( tmpR1 ) * DSIN( tmpR2 )
    CASE(2)
        phm%Vxyz(1) = phm%V * DSIN( tmpR1 ) * DCOS( tmpR2 )
        phm%Vxyz(2) = phm%V * DCOS( tmpR1 )
        phm%Vxyz(3) = phm%V * DSIN( tmpR1 ) * DSIN( tmpR2 )
    CASE(3)
        phm%Vxyz(3) = phm%V * DCOS( tmpR1 )
        phm%Vxyz(2) = phm%V * DSIN( tmpR1 ) * DCOS( tmpR2 )
        phm%Vxyz(1) = phm%V * DSIN( tmpR1 ) * DSIN( tmpR2 )
    ENDSELECT

END SUBROUTINE Diffused


!======================================================================
!----------------------------------------------------------------------
! Perodic response on boundary
!----------------------------------------------------------------------
SUBROUTINE Perodic( phm, hit )
USE VAR, ONLY: Phonon, L, Ne
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
        IF ( phm%eID(tmpI1).ne.Ne(tmpI1) ) CALL Errors(10)
        phm%xyz(tmpI1) = L(tmpI1)
        phm%eID(tmpI1) = Ne(tmpI1)
    END SELECT
    

END SUBROUTINE Perodic


!======================================================================
!----------------------------------------------------------------------
! Create or delete phonons in elements for energy conservation
!----------------------------------------------------------------------
SUBROUTINE CreateDelete( Rseed )
USE ROUTINES
USE VAR
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
                
                phn(EmptyID(s1:s2))%xyz(1) = R(:, 1) * dL(1) + &
                                             ele(i, j, k)%BD(1, 1)
                phn(EmptyID(s1:s2))%xyz(2) = R(:, 2) * dL(2) + &
                                             ele(i, j, k)%BD(1, 2)
                phn(EmptyID(s1:s2))%xyz(3) = R(:, 3) * dL(3) + &
                                             ele(i, j, k)%BD(1, 3)
                
                phn(EmptyID(s1:s2))%V = ele(i, j, k)%Vph
                phn(EmptyID(s1:s2))%Vxyz(1) = phn(EmptyID(s1:s2))%V * &
                                              DCOS( R(:, 4) * M_PI )
                phn(EmptyID(s1:s2))%Vxyz(2) = phn(EmptyID(s1:s2))%V * &
                  DSIN( R(:, 4) * M_PI ) * DCOS( R(:, 5) * M_PI * 2D0 )
                phn(EmptyID(s1:s2))%Vxyz(3) = phn(EmptyID(s1:s2))%V * &
                  DSIN( R(:, 4) * M_PI ) * DSIN( R(:, 5) * M_PI * 2D0 )
                
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
END MODULE ADV