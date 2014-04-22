!!#####################################################################
!! Module Name: rng
!! Purpose: Random number generator (RNG)
!! Author: PY Chuang
!!#####################################################################
MODULE RNG
IMPLICIT NONE
PRIVATE
PUBLIC:: rng_t, RNG_SEED, RAN_NUM
!----------------------------------------------------------------------
INTEGER(KIND=4), PARAMETER:: ns = 4
INTEGER(KIND=4), PARAMETER, DIMENSION(ns):: default_seed = &
                            (/521288629,362436069,16163801,1131199299/)
!----------------------------------------------------------------------
! rng_t:
!       A data type for storing the state of the random number
!       generator. "state" is a seed vector.
!----------------------------------------------------------------------
TYPE:: rng_t
    INTEGER(KIND=4), DIMENSION(ns):: state
END TYPE rng_t
!----------------------------------------------------------------------
! RAN_NUM:
!       An overloading subroutine.  It returns random numbers.
!----------------------------------------------------------------------
INTERFACE RAN_NUM
    MODULE PROCEDURE RAN_NUM_ONE_DOUBLE
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY_2D
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY_3D
END INTERFACE RAN_NUM
!----------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------

!======================================================================
!======================================================================
    SUBROUTINE RNG_SEED( self, seed )
    IMPLICIT NONE
    TYPE(rng_t), INTENT(INOUT):: self
    INTEGER(KIND=4), INTENT(IN):: seed
    !------------------------------------------------------------------
    ! This subroutine seeds the RNG using a double integer and a
    ! default seed vector.
    !
    ! self: A seed vector. It will be set to a default seed vector
    !       first and then be modified by a double integer seed. The
    !       subroutine will return the modified seed vector.
    ! seed: A double integer seed.
    !------------------------------------------------------------------

        self%state = default_seed
        self%state(1) = seed
        self%state(2:ns) = default_seed(2:ns)

    END SUBROUTINE RNG_SEED
!======================================================================
!======================================================================
    SUBROUTINE RAN_NUM_ONE_DOUBLE( self, u )
    IMPLICIT NONE
    TYPE(rng_t), INTENT(INOUT):: self
    REAL(KIND=8):: u
    INTEGER(KIND=4):: imz
    !------------------------------------------------------------------
    ! This subroutine randomly draws an uniform real number on [0,1].
    !
    ! self: An input seed vector.  The subroutine will modified the
    !       seed to next state after generating the random number and
    !       then return it.
    ! u: An output.  Random number.
    !------------------------------------------------------------------

        imz = self%state(1) - self%state(3)

        IF ( imz < 0 ) imz = imz + 2147483579

        self%state(1) = self%state(2)
        self%state(2) = self%state(3)
        self%state(3) = imz
        self%state(4) = 69069 * self%state(4) + 1013904243

        imz = imz + self%state(4)
        u = 0.5d0 + 0.23283064d-9 * DBLE( imz )

    END SUBROUTINE RAN_NUM_ONE_DOUBLE
!======================================================================
!======================================================================
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY( self, A )
    IMPLICIT NONE
    TYPE(rng_t),INTENT(INOUT):: self
    INTEGER(KIND=4):: i, N(1)
    REAL(KIND=8):: A(:)

        N = SIZE( A )
        DO i = 1, N(1)
            CALL RAN_NUM_ONE_DOUBLE( self, A(i) )
        ENDDO

    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY
!======================================================================
!======================================================================
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY_2D( self, A )
    IMPLICIT NONE
    TYPE(rng_t),INTENT(INOUT):: self
    INTEGER(KIND=4):: i, j, N(2)
    REAL(KIND=8):: A(:, :)

        N = SHAPE( A )
        DO j = 1, N(2)
            DO i = 1, N(1)
                CALL RAN_NUM_ONE_DOUBLE( self, A(i, j) )
            ENDDO
        ENDDO

    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY_2D
!======================================================================
!======================================================================
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY_3D( self, A )
    IMPLICIT NONE
    TYPE(rng_t),INTENT(INOUT):: self
    INTEGER(KIND=4):: i, j, k, N(3)
    REAL(KIND=8):: A(:, :, :)

        N = SHAPE( A )
        DO k = 1, N(3)
            DO j = 1, N(2)
                DO i = 1, N(1)
                    CALL RAN_NUM_ONE_DOUBLE( self, A(i, j, k) )
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY_3D
!======================================================================
!======================================================================
END MODULE RNG
