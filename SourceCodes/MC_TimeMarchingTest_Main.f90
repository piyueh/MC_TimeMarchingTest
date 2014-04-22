!!#####################################################################
!! Module Name: Main
!! Purpose: Main program
!! Author: PY Chuang
!!#####################################################################
PROGRAM main
USE RNG
USE variables
IMPLICIT NONE
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: i, icpu

    NCores = 1
    dt = 1D-1
    iter0 = 0
    iterations = 20000000

    ALLOCATE( SeedMP(NCores) )
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( tmpR1 )

    DO i = 1, NCores
        CALL RNG_SEED( SeedMP(i), INT( tmpR1*1D8 ) )
    ENDDO

    icpu = 1
    tmpR2 = 0D0
    i = 0
    DO iter = iter0 + 1, iterations
        CALL RAN_NUM( SeedMP(icpu), tmpR1 )
        tmpR2 = tmpR2 + tmpR1
        i = i + 1
    ENDDO

    WRITE(*, *) i, tmpR2 / DBLE( i )

END PROGRAM main
