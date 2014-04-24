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
REAL(KIND=8):: R1
INTEGER(KIND=4):: mt, i, j, k

    CaseName = "Test"

    BCs = (/ 1, 2, 2 /)
    
    NCores = 1
    iterations = 20000000
    
    nOutput = 100

    CALL ReadTables
    
    CALL initialize
    CALL Reorder_CellInfo
    CALL CreateDelete( SeedMP(1) )

    Tavg = 0D0
    Eavg = 0D0
    ct = 0
    DO iter = iter0 + 1, iterations
        
        CALL advance
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
        ENDIF
        
        WRITE(*, *) iter, RNph, SUM( ele%E ) / DBLE( PRODUCT( Ne ) ), &
                    SUM( ele%T ) / DBLE( PRODUCT( Ne ) )
        
    ENDDO

END PROGRAM main
