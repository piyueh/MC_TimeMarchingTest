!!#####################################################################
!! Module Name: Main
!! Purpose: Main program
!! Author: PY Chuang
!!#####################################################################
PROGRAM main
USE RNG
USE VAR
USE ADV
USE INIT
IMPLICIT NONE
INTEGER(KIND=4):: i
INTEGER(KIND=4):: iCPU


    CALL ReadTables
    NCores = 1
    iterations = 20000000

    iCPU = 1

    CALL initialize


    DO iter = iter0 + 1, iterations
        DO i = 1, Nph
            CALL advance( phn(i), SeedMP(iCPU), dt )
        ENDDO
        WRITE(*, *) iter
    ENDDO

END PROGRAM main
