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
REAL(KIND=8):: R1, t(4)
INTEGER(KIND=4):: mt, i, j, k

    CaseName = "Test"

    BCs = (/ 1, 2, 2 /)

    NCores = 1
    iterations = 20000000

    nOutput = 400

    CALL ReadTables

    CALL initialize
    CALL Reorder_CellInfo
    CALL CreateDelete( SeedMP(1) )

    WRITE(OutputFileName, "(I8.8, '_Transient.txt')") 0
    OPEN( UNIT = 150, FILE = OutputFileName )
    WRITE(150, *) time
    WRITE(150, *) ele(:, 1, 1)%T
    CLOSE( 150 )
    
    Tavg = 0D0
    Eavg = 0D0
    ct = 0
    DO iter = iter0 + 1, iterations

        CALL advance( t )
        time = time + dt
        Eavg = Eavg + ele%E
        ct = ct + 1

        IF ( MOD( iter, nOutput ).eq.0 ) THEN
            Eavg = Eavg / DBLE( ct )
            DO k = 1, Ne(3); DO j = 1, Ne(2); DO i = 1, Ne(1)
                R1 = Eavg(i, j, k) / dV
                mt = ele(i, j, k)%Mat
                CALL Etable( mt, 1, R1, Tavg(i, j, k) )
            ENDDO; ENDDO; ENDDO

            WRITE(OutputFileName, "(I8.8, '.txt')") iter
            OPEN( UNIT = 150, FILE = OutputFileName )
            WRITE(150, NML = Output)
            CLOSE( 150 )

            ct = 0
            Eavg = 0

            WRITE(OutputFileName, "(I8.8, '_Transient.txt')") iter
            OPEN( UNIT = 150, FILE = OutputFileName )
            WRITE(150, *) time
            WRITE(150, *) ele(:, 1, 1)%T
            CLOSE( 150 )
        ENDIF

        WRITE(*, "(I7, 2X, I8, 2X, 3(F6.3, 2X))") iter, RNph, t(2) - t(1), t(3) - t(2), t(4) - t(3)

    ENDDO

END PROGRAM main
