!!#####################################################################
!! Module Name: Main
!! Purpose: Main program
!! Author: PY Chuang
!!#####################################################################
PROGRAM main
USE RNG
USE VAR
USE ADV
USE IO
USE ROUTINES
IMPLICIT NONE
REAL(KIND=8):: R1, t(4)
INTEGER(KIND=4):: mt, i, j, k
INTEGER(KIND=4):: NCores
TYPE(rng_t), ALLOCATABLE:: SeedMP(:)


    WRITE(*, '("Enter CASE Name: ")', ADVANCE = 'NO')
    !READ(*, *) CaseName
    CaseName = "Adiabatic_Transien_Ge"
    InputFileName = casename(1:LEN_TRIM(casename))//'_initial.txt'
    
    
    !WRITE(*, '("Enter BCs (BCx, BCy, BCz): ")', ADVANCE = 'NO')
    !READ(*, *) BCs
    BCs = (/ 1, 2, 2 /)
    
    
    !WRITE(*, '("Enter the Numer of CPU Cores: ")', ADVANCE = 'NO')
    !READ(*, *) NCores
    NCores = 1
    ALLOCATE( SeedMP(NCores) )
    CALL RANDOM_SEED()
    DO i = 1, NCores
        CALL RANDOM_NUMBER( R1 )
        CALL RNG_SEED( SeedMP(i), INT( R1*1D8 ) )
    ENDDO
    
    !WRITE(*, '(Enter the Maximum Iterations: )', ADVANCE = 'NO')
    !READ(*, *) iterations
    iterations = 20000000
    
    
    WRITE(*, '("Enter How Many Steps It Outputs Once: ")', ADVANCE = 'NO')
    !READ(*, *) nOutput
    nOutput = 100

    CALL Initialize_Materials

    CALL initialize
    
    CALL Reorder_CellInfo

    CALL CreateDelete( SeedMP(1) )
    
    
    WRITE(OutputFileName, "(I8.8, '_Transient.txt')") 0
    OPEN( UNIT = 150, FILE = OutputFileName )
    WRITE(150, *) time
    WRITE(150, *) ele(:, 1, 1)%T
    CLOSE( 150 )
    
    ct = 0
    DO iter = iter0 + 1, iterations

        CALL advance( NCores, SeedMP, t )
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
            WRITE(150, NML = DataOutput)
            CLOSE( 150 )

            ct = 0
            Eavg = 0
            Tavg = 0

            WRITE(OutputFileName, "(I8.8, '_Transient.txt')") iter
            OPEN( UNIT = 150, FILE = OutputFileName )
            WRITE(150, *) time
            WRITE(150, *) ele(:, 1, 1)%T
            CLOSE( 150 )
        ENDIF

        WRITE(*, "(I7, 2X, I8, 2X, 3(F6.3, 2X))") iter, RNph, &
                                                  t(2) - t(1), &
                                                  t(3) - t(2), &
                                                  t(4) - t(3)

    ENDDO
    
    DEALLOCATE( SeedMP )

END PROGRAM main
