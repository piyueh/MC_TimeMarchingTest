!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver
!! File Name: MC3D_Structured_HEAT
!! Contenent: Collections of Subroutines Related to Heat Flux Control
!! Author: PY Chuang
!!#####################################################################
MODULE HEAT
USE RNG
IMPLICIT NONE
CONTAINS


!======================================================================
!----------------------------------------------------------------------
! Initialize heat control boundary conditions.
!----------------------------------------------------------------------
SUBROUTINE init_heat_BC
USE VAR_BC
USE VAR_ph, ONLY: bundle
USE VAR_SPACES, ONLY: Ne, ele, dA_heat
USE VAR_Others, ONLY: TimeStep
IMPLICIT NONE

    WRITE(*, '("Enter Boundary Temperatures (TL, TR): ")', ADVANCE = 'NO')
    !READ(*, *) TBCL, TBCR
    TBCL = 350D0
    TBCR = 310D0

    CALL init_inj_ph_prop( TBCL, VphBCL, EphBCL, bundle)
    CALL init_inj_ph_prop( TBCR, VphBCR, EphBCR, bundle)
    CALL calc_inj_heat( TBCL, Ne(2), Ne(3), ele(1, :, :)%Mat,  &
                        HinjectL, dA_heat, TimeStep )
    CALL calc_inj_heat( TBCR, Ne(2), Ne(3), ele(Ne(1), :, :)%Mat,  &
                        HinjectR, dA_heat, TimeStep )
    
    CALL calc_approx_NHeatPh( TBCL, TBCR, dA_Heat * Ne(2) * Ne(3), &
                              TimeStep, bundle, FNHeatPh )

    ALLOCATE( HeatPhn(FNHeatPh) )
    ALLOCATE( dtHeat(FNHeatPh) )
    RNHeatPh = 0
    HeatPhn%Exist = .FALSE.
                              
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
USE ROUTINES, ONLY: energy, Etable
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: T, A, dt
REAL(KIND=8), INTENT(OUT):: Hinject(Ny, Nz)
INTEGER(KIND=4), INTENT(IN):: Ny, Nz, Mat(Ny, Nz)
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
! Calculate approximated number of injected phonons
!----------------------------------------------------------------------
SUBROUTINE calc_approx_NHeatPh( T1, T2, A, dt, w, N )
USE ROUTINES, ONLY: energy, Etable
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: T1, T2, A, dt, w
INTEGER(KIND=4), INTENT(OUT):: N
REAL(KIND=8):: U1(2), N1(2), V1(2), Eph1(2)
REAL(KIND=8):: U2(2), N2(2), V2(2), Eph2(2)
INTEGER(KIND=4):: i

    DO i = 1, 2
        CALL energy( i, T1, U1(i) )
        CALL Etable( i, 2, U1(i), N1(i) )
        CALL Etable( i, 4, U1(i), V1(i) )
        Eph1(i) = U1(i) / N1(i) * w
    ENDDO
    
    DO i = 1, 2
        CALL energy( i, T2, U2(i) )
        CALL Etable( i, 2, U2(i), N2(i) )
        CALL Etable( i, 4, U2(i), V2(i) )
        Eph2(i) = U2(i) / N2(i) * w
    ENDDO
    
    N = MAX(MAXVAL(U1), MAXVAL(U2)) * MAX(MAXVAL(V1), MAXVAL(V2)) * &
        A * dt * 0.25 / MIN(MINVAL(Eph1), MINVAL(Eph2))
        
    N = N * 1.2

END SUBROUTINE calc_approx_NHeatPh


!======================================================================
!----------------------------------------------------------------------
! Deal with heat control
!----------------------------------------------------------------------
SUBROUTINE Heat_Control( RSeed )
USE ADV
USE VAR_ALL
USE ROUTINES
IMPLICIT NONE
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: H, R1
INTEGER(KIND=4):: i, j, k, m, N, bg, ed
REAL(KIND=8), ALLOCATABLE:: R(:, :)

    dtHeat = 0D0

    DO k = 1, Ne(3); DO j = 1, Ne(2)
        H = HinjectL(j, k)
        SELECTCASE( PoolL(PoolSize, j, k)%Mat)
        CASE(-1)
            N = INT( H / EphBCL(ele(1, j, k)%Mat) + 0.5 )
            IF (N.gt.0) THEN
            ele(1, j, k)%Ediff = ele(1, j, k)%Ediff + H - &
                                                       N * EphBCL(ele(1, j, k)%Mat)

            bg = RNHeatPh + 1
            ed = RNHeatPh + N

            ALLOCATE( R(N, 5) )
            CALL RAN_NUM( RSeed, R )

            R(:, 4) = DSQRT( R(:, 4) )
            R(:, 5) = R(:, 5) * M_PI * 2D0

            HeatPhn(bg:ed)%xyz(1) = 0D0
            HeatPhn(bg:ed)%xyz(2) = R(:, 2) * dL(2)
            HeatPhn(bg:ed)%xyz(3) = R(:, 3) * dL(3)

            HeatPhn(bg:ed)%V = VphBCL(ele(1, j, k)%Mat)
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

            HeatPhn(bg:ed)%Mat = ele(1, j, k)%Mat

            HeatPhn(bg:ed)%eID(1) = 1
            HeatPhn(bg:ed)%eID(2) = j
            HeatPhn(bg:ed)%eID(3) = k

            HeatPhn(bg:ed)%Exist = .TRUE.
            
            SELECTCASE( WAY_FlightTime )
            CASE(1)
                dtHeat(bg:ed) = R(:, 1) * TimeStep
                DO i = bg, ed
                    CALL phn_adv_CT( HeatPhn(i), RSeed, dtHeat(i) )
                ENDDO
            CASE(2)
                dtHeat(bg:ed) = - DLOG( R(:, 1) ) / GammaT
                DO i = bg, ed
                    CALL phn_adv_RT( HeatPhn(i), RSeed, dtHeat(i) )
                ENDDO
            END SELECT
            
            DEALLOCATE( R )

            RNHeatPh = ed
            ENDIF
        CASE DEFAULT
            WRITE(*, *) "L", j, k
            DO WHILE( H.gt.0D0 )
                CALL RAN_NUM( Rseed, R1 )
                m = MAX(INT( R1 * PoolSize + 0.5 ), 1)
                RNHeatPh = RNHeatPh + 1
                
                HeatPhn(RNHeatPh)%xyz(1) = 0D0
                HeatPhn(RNHeatPh)%xyz(2) = PoolL(m, j, k)%y
                HeatPhn(RNHeatPh)%xyz(3) = PoolL(m, j, k)%z
                
                HeatPhn(RNHeatPh)%V = VphBCL(ele(1, j, k)%Mat)
                HeatPhn(RNHeatPh)%Vxyz(1) = HeatPhn(RNHeatPh)%V * &
                                            PoolL(m, j, k)%direction(1)
                HeatPhn(RNHeatPh)%Vxyz(2) = HeatPhn(RNHeatPh)%V * &
                                            PoolL(m, j, k)%direction(2)
                HeatPhn(RNHeatPh)%Vxyz(3) = HeatPhn(RNHeatPh)%V * &
                                            PoolL(m, j, k)%direction(3)
                
                HeatPhn(RNHeatPh)%E = EphBCL(PoolL(m, j, k)%Mat)
                
                HeatPhn(RNHeatPh)%Org(1) = HeatPhn(RNHeatPh)%xyz(1)
                HeatPhn(RNHeatPh)%Org(2) = HeatPhn(RNHeatPh)%xyz(2)
                HeatPhn(RNHeatPh)%Org(3) = HeatPhn(RNHeatPh)%xyz(3)
                
                HeatPhn(RNHeatPh)%Dis(1) = 0D0
                HeatPhn(RNHeatPh)%Dis(2) = 0D0
                HeatPhn(RNHeatPh)%Dis(3) = 0D0
                
                HeatPhn(RNHeatPh)%Mat = PoolL(m, j, k)%Mat
                
                HeatPhn(RNHeatPh)%eID(1) = 1
                HeatPhn(RNHeatPh)%eID(2) = j
                HeatPhn(RNHeatPh)%eID(3) = k
                
                HeatPhn(RNHeatPh)%Exist = .TRUE.
                
                SELECTCASE(WAY_FlightTime)
                CASE(1)
                    dtHeat(RNHeatPh) = PoolL(m, j, k)%dtRemain
                    CALL phn_adv_CT( HeatPhn(RNHeatPh), RSeed, dtHeat(RNHeatPh) )
                CASE(2)
                    CALL RAN_NUM( Rseed, R1 )
                    dtHeat(RNHeatPh) = - DLOG( R1 ) / GammaT
                    CALL phn_adv_RT( HeatPhn(RNHeatPh), RSeed, dtHeat(RNHeatPh) )
                END SELECT
                
                H = H - HeatPhn(RNHeatPh)%E
                
            ENDDO
            
            ele(1, j, k)%Ediff = ele(1, j, k)%Ediff + H
            
        END SELECT
    ENDDO; ENDDO
    
    
    DO k = 1, Ne(3); DO j = 1, Ne(2)
        H = HinjectR(j, k)
        SELECTCASE( PoolR(PoolSize, j, k)%Mat )
        CASE(-1)
            N = INT( H / EphBCR(ele(Ne(1), j, k)%Mat) + 0.5 )
            IF (N.gt.0) THEN
            ele(1, j, k)%Ediff = ele(1, j, k)%Ediff + H - &
                                                       N * EphBCR(ele(Ne(1), j, k)%Mat)

            bg = RNHeatPh + 1
            ed = RNHeatPh + N

            ALLOCATE( R(N, 5) )
            CALL RAN_NUM( RSeed, R )

            R(:, 4) = - DSQRT( R(:, 4) )
            R(:, 5) = R(:, 5) * M_PI * 2D0

            HeatPhn(bg:ed)%xyz(1) = L(1)
            HeatPhn(bg:ed)%xyz(2) = R(:, 2) * dL(2)
            HeatPhn(bg:ed)%xyz(3) = R(:, 3) * dL(3)

            HeatPhn(bg:ed)%V = VphBCR(ele(Ne(1), j, k)%Mat)
            HeatPhn(bg:ed)%Vxyz(1) = HeatPhn(bg:ed)%V * R(:, 4)
            HeatPhn(bg:ed)%Vxyz(2) = HeatPhn(bg:ed)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DCOS( R(:, 5) )
            HeatPhn(bg:ed)%Vxyz(3) = HeatPhn(bg:ed)%V * &
                              DSQRT(1D0 - R(:, 4)**2) * DSIN( R(:, 5) )

            HeatPhn(bg:ed)%E = EphBCR(ele(Ne(1), j, k)%Mat)

            HeatPhn(bg:ed)%Org(1) = HeatPhn(bg:ed)%xyz(1)
            HeatPhn(bg:ed)%Org(2) = HeatPhn(bg:ed)%xyz(2)
            HeatPhn(bg:ed)%Org(3) = HeatPhn(bg:ed)%xyz(3)

            HeatPhn(bg:ed)%Dis(1) = 0D0
            HeatPhn(bg:ed)%Dis(2) = 0D0
            HeatPhn(bg:ed)%Dis(3) = 0D0

            HeatPhn(bg:ed)%Mat = ele(Ne(1), j, k)%Mat

            HeatPhn(bg:ed)%eID(1) = Ne(1)
            HeatPhn(bg:ed)%eID(2) = j
            HeatPhn(bg:ed)%eID(3) = k

            HeatPhn(bg:ed)%Exist = .TRUE.

            SELECTCASE( WAY_FlightTime )
            CASE(1)
                dtHeat(bg:ed) = R(:, 1) * TimeStep
                DO i = bg, ed
                    CALL phn_adv_CT( HeatPhn(i), RSeed, dtHeat(i) )
                ENDDO
            CASE(2)
                dtHeat(bg:ed) = - DLOG( R(:, 1) ) / GammaT
                DO i = bg, ed
                    CALL phn_adv_RT( HeatPhn(i), RSeed, dtHeat(i) )
                ENDDO
            END SELECT
            
            DEALLOCATE( R )

            RNHeatPh = ed
            ENDIF
        CASE DEFAULT
            WRITE(*, *) "R", j, k
            DO WHILE( H.gt.0D0 )
                CALL RAN_NUM( Rseed, R1 )
                m = MAX(INT( R1 * PoolSize + 0.5 ), 1)
                RNHeatPh = RNHeatPh + 1
                
                HeatPhn(RNHeatPh)%xyz(1) = L(1)
                HeatPhn(RNHeatPh)%xyz(2) = PoolR(m, j, k)%y
                HeatPhn(RNHeatPh)%xyz(3) = PoolR(m, j, k)%z
                
                HeatPhn(RNHeatPh)%V = VphBCR(ele(Ne(1), j, k)%Mat)
                HeatPhn(RNHeatPh)%Vxyz(1) = HeatPhn(RNHeatPh)%V * &
                                            PoolR(m, j, k)%direction(1)
                HeatPhn(RNHeatPh)%Vxyz(2) = HeatPhn(RNHeatPh)%V * &
                                            PoolR(m, j, k)%direction(2)
                HeatPhn(RNHeatPh)%Vxyz(3) = HeatPhn(RNHeatPh)%V * &
                                            PoolR(m, j, k)%direction(3)
                
                HeatPhn(RNHeatPh)%E = EphBCR(PoolL(m, j, k)%Mat)
                
                HeatPhn(RNHeatPh)%Org(1) = HeatPhn(RNHeatPh)%xyz(1)
                HeatPhn(RNHeatPh)%Org(2) = HeatPhn(RNHeatPh)%xyz(2)
                HeatPhn(RNHeatPh)%Org(3) = HeatPhn(RNHeatPh)%xyz(3)
                
                HeatPhn(RNHeatPh)%Dis(1) = 0D0
                HeatPhn(RNHeatPh)%Dis(2) = 0D0
                HeatPhn(RNHeatPh)%Dis(3) = 0D0
                
                HeatPhn(RNHeatPh)%Mat = PoolR(m, j, k)%Mat
                
                HeatPhn(RNHeatPh)%eID(1) = Ne(1)
                HeatPhn(RNHeatPh)%eID(2) = j
                HeatPhn(RNHeatPh)%eID(3) = k
                
                HeatPhn(RNHeatPh)%Exist = .TRUE.
                
                SELECTCASE(WAY_FlightTime)
                CASE(1)
                    dtHeat(RNHeatPh) = PoolR(m, j, k)%dtRemain
                    CALL phn_adv_CT( HeatPhn(RNHeatPh), RSeed, dtHeat(RNHeatPh) )
                CASE(2)
                    CALL RAN_NUM( Rseed, R1 )
                    dtHeat(RNHeatPh) = - DLOG( R1 ) / GammaT
                    CALL phn_adv_RT( HeatPhn(RNHeatPh), RSeed, dtHeat(RNHeatPh) )
                END SELECT
            
                H = H - HeatPhn(RNHeatPh)%E
            
            ENDDO
            
            ele(Ne(1), j, k)%Ediff = ele(Ne(1), j, k)%Ediff + H
            
        END SELECT
    ENDDO; ENDDO
    
END SUBROUTINE Heat_Control


!======================================================================
!----------------------------------------------------------------------
! Random Injection
!----------------------------------------------------------------------
! SUBROUTINE Random_Injection( Pool )
! IMPLICIT NONE



!END SUBROUTINE Random_Injection

!======================================================================
!----------------------------------------------------------------------
! Inject phonons in pools
!----------------------------------------------------------------------

END MODULE HEAT