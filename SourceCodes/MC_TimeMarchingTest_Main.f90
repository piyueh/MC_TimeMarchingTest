!!#####################################################################
!! Module Name: Main
!! Purpose: Main program
!! Author: PY Chuang
!!#####################################################################
PROGRAM main
USE RNG
USE VAR
USE ADV
IMPLICIT NONE
REAL(KIND=8):: tmpR1, tmpR2
INTEGER(KIND=4):: i, j, k
INTEGER(KIND=4):: iCPU
REAL(KIND=8), ALLOCATABLE:: rannum(:)

    NCores = 1
    dt = 1D-1
    iter0 = 0
    iterations = 20000000

    iCPU = 1

    ALLOCATE( SeedMP(NCores) )
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER( tmpR1 )

    DO i = 1, NCores
        CALL RNG_SEED( SeedMP(i), INT( tmpR1*1D8 ) )
    ENDDO

    L = (/ 1D3, 1D3, 1D1 /)
    Ne = (/ 5D1, 5D1, 1D0 /)
    dL = L / DBLE( Ne )

    ALLOCATE( ele(Ne(1), Ne(2), Ne(3)) )
    DO k = 1, Ne(3)
        DO j = 1, Ne(2)
            DO i = 1, Ne(1)
                ele(i, j, k)%BD(:, 1) = DBLE( (/ i-1, i /) ) * dL(1)
                ele(i, j, k)%BD(:, 2) = DBLE( (/ j-1, j /) ) * dL(2)
                ele(i, j, k)%BD(:, 3) = DBLE( (/ k-1, k /) ) * dL(3)
                ele(i, j, k)%T = 330.034003189333D0
                ele(i, j, k)%E = 1535.50000000000D0
                ele(i, j, k)%N = 141.493541620377D0
                ele(i, j, k)%Eph = ele(i, j, k)%E / ele(i, j, k)%N
                ele(i, j, k)%Ediff = 0D0
                ele(i, j, k)%Vph = 1.30822382162595D0
                ele(i, j, k)%MFP = 133.386434452556D0
            ENDDO
        ENDDO
    ENDDO

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
