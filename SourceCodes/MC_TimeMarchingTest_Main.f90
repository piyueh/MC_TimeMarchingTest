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
INTEGER(KIND=4):: i, j, k, ct
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
                ele(i, j, k)%T = 0D0
                ele(i, j, k)%E = 0D0
                ele(i, j, k)%Ediff = 0D0
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
    
    
!     DO iter = iter0 + 1, iterations
        DO i = 1, Nph
            CALL phn_advance( phn(i),SeedMP(iCPU), dt )
        ENDDO
!     ENDDO

END PROGRAM main


!----------------------------------------------------------------------
! Phonon Advance Subroutine
!----------------------------------------------------------------------
SUBROUTINE phn_advance( phm, RSeed, dt )
USE variables, ONLY: Phonon, Element, ele
USE RNG
IMPLICIT NONE
TYPE(Phonon), INTENT(INOUT):: phm
REAL(KIND=8), INTENT(IN):: dt
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: xyz(3), ds(3), dtUsed, dtRemain
INTEGER(KIND=4):: i, tmpI1, tmpI2, tmpI3, Loc(1)
    
    tmpI1 = phm%eID(1)
    tmpI2 = phm%eID(2)
    tmpI3 = phm%eID(3)
    
    dtRemain = dt
    
    DO i = 1, 3
        IF ( phm%Vxyz(i).gt.0 ) THEN
            ds(i) = (ele(tmpI1, tmpI2, tmpI3)%BD(2, i) - phm%xyz(i)) /&
                    phm%Vxyz(i)
        ELSEIF ( phm%Vxyz(i).lt.0 ) THEN
            ds(i) = (ele(tmpI1, tmpI2, tmpI3)%BD(1, i) - phm%xyz(i)) /&
                    phm%Vxyz(i)
        ELSE
            ds(i) = 1D8
        ENDIF
    ENDDO
    
    Loc = MINLOC( ds )
    dtUsed = ds(Loc(1))
    
    !##################################################################
    ! Stop at here.  04/23/2014
    !##################################################################
    
    IF ( dtUsed.gt.dt ) THEN
        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        dtRemain = 0D0
        CALL intrinsic_scattering
    ELSEIF ( dtUsed.lt.dt ) THEN
        dtRemain = dtRemain - dtUsed
        phm%xyz = phm%xyz + dtUsed * phm%Vxyz
        CALL intrinsic_scattering
    ELSE
    ENDIF
    

END SUBROUTINE phn_advance