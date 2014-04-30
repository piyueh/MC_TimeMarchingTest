!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_HEAT
!! Contenent: Collections of Subroutines Related to Heat Flux Control
!! Author: PY Chuang
!!#####################################################################
MODULE HEAT
IMPLICIT NONE
CONTAINS


!======================================================================
!----------------------------------------------------------------------
! Initialize heat control boundary conditions.
!----------------------------------------------------------------------
SUBROUTINE init_heat_BC
USE VAR_BC
USE VAR_SPACES, ONLY: Ne, ele, dA_heat
USE VAR_Others, ONLY: bundles, TimeStep
IMPLICIT NONE

    WRITE(*, '("Enter Boundary Temperatures (TL, TR): )', ADVANCE = 'NO')
    READ(*, *) TBCL, TBCR

    CALL init_inj_ph_prop( TBCL, VphBCL, EphBCL, bundles)
    CALL init_inj_ph_prop( TBCR, VphBCR, EphBCR, bundles)
    CALL calc_inj_heat( TBCL, Ne(2), Ne(3), ele(1, :, :)%Mat,  &
                        HinjectL, dA_heat, TimeStep )
    CALL calc_inj_heat( TBCR, Ne(2), Ne(3), ele(Ne(1), :, :)%Mat,  &
                        HinjectR, dA_heat, TimeStep )

END SUBROUTINE init_heat_BC

!======================================================================
!----------------------------------------------------------------------
! Get the informations of phonons injected from the boundary at
! specific temperature and material.
!----------------------------------------------------------------------
SUBROUTINE init_inj_ph_prop( TBC, VphBC, EphBC, w )
USE ROUTINES, ONLY: energy, Etable
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: TBC, w
REAL(KIND=8), INTENT(OUT):: VphBC(2), EphBC(2)
REAL(KIND=8):: U
INTEGER(KIND=4):: i

    DO i = 1, 2
        CALL energy( i, TBC, U )
        CALL Etable( i, 4, U, VphBC(i) )
        CALL Etable( i, 2, U, EphBC(i) )
        EphBC(i) = w * U / EphBC(i)
    ENDDO

END SUBROUTINE init_inj_ph_prop


!======================================================================
!----------------------------------------------------------------------
! Calculate how much energy the boundary elements should obtain at each
! time step.
!----------------------------------------------------------------------
SUBROUTINE calc_inj_heat( T, Ny, Nz, Mat, Hinject, A, dt )
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: T, A, dt
REAL(KIND=8), INTENT(OUT):: Hinject(Ny, Nz)
INTEGER(KIND=4), INTENT(IN):: Ny, Nx, Mat(Ny, Nz)
REAL(KIND=8):: U, V
INTEGER(KIND=8):: j, k

    DO k = 1, Nz; DO j = 1, Ny
        CALL energy( Mat(j, k), T, U )
        CALL Etable( Mat(j, k), 4, U, V)
        Hinject(j, k) = U * V * A * dt * 0.25
    ENDDO; ENDDO

END SUBROUTINE calc_inj_heat


!======================================================================
!----------------------------------------------------------------------
! Deal with heat control
!----------------------------------------------------------------------
SUBROUTINE Heat_Control
IMPLICIT NONE

    HeatPhn%Exist = .FALSE.
    RNHeatPh = 0
    dtHeat = 0D0

    DO k = 1, Ne(3); DO j = 1, Ne(2)
        H = HinjectL(j, k)
        SELECTCASE( PoolL(j, k)%Mat)
        CASE(-1)
            N = INT( H / EphBCL(j, k) + 0.5 )
            ele(1, j, k)%Ediff = ele(1, j, k)%Ediff + H - &
                                                       N * EphBCL(j, k)

            bg = RNHeatPh + 1
            ed = RNHeatPh + N

            ALLOCATE( R(N, 5) )
            CALL RAN_NUM( R )

            SELECTCASE( WAY_FlightTime )
            CASE(1)
                dtHeat(bg:ed) = R(:, 1) * TimeStep
            CASE(2)
                dtHeat(bg:ed) = - DLOG( R(:, 1) ) / GammaT
            END SELECT

            R(:, 4) = DSQRT( R(:, 4) )
            R(:, 5) = R(:, 5) * M_PI * 2D0

            HeatPhn(bg:ed)%xyz(1) = 0D0
            HeatPhn(bg:ed)%xyz(2) = R(:, 2) * dL(2)
            HeatPhn(bg:ed)%xyz(3) = R(:, 3) * dL(3)

            HeatPhn(bg:ed)%V = VphVCL(ele(1, j, k)%Mat)
            HeatPhn(bg:ed)%Vxyz(1) = HeatPhn(bg:ed)%V * R(:, 4)
            HeatPhn(bg:ed)%Vxyz(2) = HeatPhn(bg:ed)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DCOS( R(:, 5) )
            HeatPhn(bg:ed)%Vxyz(3) = HeatPhn(bg:ed)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DSIN( R(:, 5) )

            HeatPhn(bg:ed)%E = EphBCL(ele(1, j, k)%Mat)

            HeatPhn(bg:ed)%Org(1) = HeatPhn(bg:ed)%xyz(1)
            HeatPhn(bg:ed)%Org(2) = HeatPhn(bg:ed)%xyz(2)
            HeatPhn(bg:ed)%Org(3) = HeatPhn(bg:ed)%xyz(3)

            HeatPhn(bg:ed)%Dis(1) = 0D0
            HeatPhn(bg:ed)%Dis(2) = 0D0
            HeatPhn(bg:ed)%Dis(3) = 0D0

            HeatPhn(bg:ed)%MAt = ele(1, j, k)%Mat

            HeatPhn(bg:ed)%eID(1) = 1
            HeatPhn(bg:ed)%eID(2) = j
            HeatPhn(bg:ed)%eID(3) = k

            HeatPhn(bg:ed)%Exist = .TRUE.

            DEALLOCATE( R )

            RNHeatPh = ed

        CASE DEFAULT
            N = 0
            DO WHILE(  )
                CALL RAN_NUM( Rseed, R1 )
                m = MIN(INT( R1 * PoolSize + 0.5 ), 1)
                N = N + 1
                H = H - PoolL(j, k, m)
            ENDDO
        END SELECT
    ENDDO; ENDDO

END SUBROUTINE Heat_Control


!======================================================================
!----------------------------------------------------------------------
! Random Injection
!----------------------------------------------------------------------
! SUBROUTINE Random_Injection( Pool )
! IMPLICIT NONE



END SUBROUTINE Random_Injection

!======================================================================
!----------------------------------------------------------------------
! Inject phonons in pools
!----------------------------------------------------------------------

END MODULE HEAT