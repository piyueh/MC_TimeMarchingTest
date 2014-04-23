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


    CALL ReadTables
    NCores = 1
    iterations = 20000000

    CALL initialize
    CALL cellInfo

    DO iter = iter0 + 1, iterations
        CALL advance
        WRITE(*, *) iter, SUM( ele%Ediff ) / DBLE( PRODUCT( Ne ) )
        WRITE(*, *) "   ", MAXLOC( ele%Ediff ), MAXVAL( ele%Ediff )
        WRITE(*, *) "   ", MINLOC( ele%Ediff ), MINVAL( ele%Ediff )
    ENDDO

END PROGRAM main
