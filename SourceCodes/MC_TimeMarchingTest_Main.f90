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
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: i, j, k
INTEGER(KIND=4):: iCPU
REAL(KIND=8), ALLOCATABLE:: rannum(:)


    CALL ReadTables
    NCores = 1
    iterations = 20000000

    iCPU = 1

    CALL initialize



    Nph = 10000
    ALLOCATE( phn(Nph) )
    ALLOCATE( rannum(5) )
    DO i = 1, Nph
        CALL RAN_NUM( SeedMP(iCPU), rannum )
        tmpR1 = rannum(4) * M_PI
        tmpR2 = rannum(5) * M_PI * 2D0

        phn(i)%xyz = rannum(1:3) * L
        phn(i)%V = 1.30822382162595D0
        phn(i)%Vxyz(1) = phn(i)%V * DCOS( tmpR1 )
        phn(i)%Vxyz(2:3) = phn(i)%V * DSIN( tmpR1 ) * &
                          (/ DCOS( tmpR2 ), DSIN( tmpR2 ) /)
        phn(i)%E = 1535.5D0 / 141.493541620377D0
        phn(i)%Org = phn(i)%xyz
        phn(i)%Dis = 0D0
        phn(i)%Mat = 1
        phn(i)%eID = INT( phn(i)%xyz / dL ) + (/ 1, 1, 1 /)
    ENDDO
    DEALLOCATE( rannum )


     DO iter = iter0 + 1, iterations
        DO i = 1, Nph
            CALL advance( phn(i), SeedMP(iCPU), dt, i )
        ENDDO
        WRITE(*, *) iter
     ENDDO

END PROGRAM main
