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

                CALL RAN_NUM( RSeed, R )

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

                phn(EmptyID(s1:s2))%Dis(1) = -10D0
                phn(EmptyID(s1:s2))%Dis(2) = -10D0
                phn(EmptyID(s1:s2))%Dis(3) = -100D0

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