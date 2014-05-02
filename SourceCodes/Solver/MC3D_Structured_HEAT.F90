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
! Make the temperature of boundary element constant.
!----------------------------------------------------------------------
SUBROUTINE Fixed_Boundary_T
USE VAR_ALL
USE ROUTINES
IMPLICIT NONE
REAL(KIND=8):: U
INTEGER(KIND=4):: i, j, k, mt
    
    i = 1
    DO k = 1, Ne(3); DO j = 1, Ne(2)
        phn(phID(ele(i, j, k)%Nbg:ele(i, j, k)%Ned))%Exist = .FALSE.
        mt = ele(i, j, k)%Mat
        CALL energy( mt, TBCL, U )
        ele(i, j, k)%T = TBCL
        
        CALL Etable( mt, 2, U, ele(i, j, k)%ND )
        ele(i, j, k)%Eph = U / ele(i, j, k)%ND * bundle
        ele(i, j, k)%Ediff = U * dV
        ele(i, j, k)%E = U * dV
        CALL Etable( mt, 4, U, ele(i, j, k)%Vph )
        CALL Etable( mt, 5, U, ele(i, j, k)%MFP )
        ele(i, j, k)%SCR = ele(i, j, k)%Vph / ele(i, j, k)%MFP
    ENDDO; ENDDO
    
    i = Ne(1)
    DO k = 1, Ne(3); DO j = 1, Ne(2)
        phn(phID(ele(i, j, k)%Nbg:ele(i, j, k)%Ned))%Exist = .FALSE.
        mt = ele(i, j, k)%Mat
        CALL energy( mt, TBCR, U )
        ele(i, j, k)%T = TBCR
        
        CALL Etable( mt, 2, U, ele(i, j, k)%ND )
        ele(i, j, k)%Eph = U / ele(i, j, k)%ND * bundle
        ele(i, j, k)%Ediff = U * dV
        ele(i, j, k)%E = U * dV
        CALL Etable( mt, 4, U, ele(i, j, k)%Vph )
        CALL Etable( mt, 5, U, ele(i, j, k)%MFP )
        ele(i, j, k)%SCR = ele(i, j, k)%Vph / ele(i, j, k)%MFP
    ENDDO; ENDDO

END SUBROUTINE Fixed_Boundary_T


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
    TBCR = 350D0

    CALL init_inj_ph_prop( TBCL, VphBCL, EphBCL, bundle)
    CALL init_inj_ph_prop( TBCR, VphBCR, EphBCR, bundle)
    CALL calc_inj_heat( TBCL, dEHeatL, dA_heat, TimeStep )
    CALL calc_inj_heat( TBCR, dEHeatR, dA_heat, TimeStep )

    HinjectL = 0D0
    HinjectR = 0D0

    CALL calc_approx_NHeatPh( TBCL, TBCR, dA_Heat * Ne(2) * Ne(3), &
                              TimeStep, bundle, FNHeatPh )

    ALLOCATE( HeatPhn(FNHeatPh) )
    ALLOCATE( dtHeat(FNHeatPh) )
    RNHeatPh = 0
    HeatPhn%Exist = .FALSE.

    WRITE(*, *) "At Left Boundary: "
    WRITE(*, *) "   TBCL = ", TBCL
    WRITE(*, *) "   VphBCL = ", VphBCL
    WRITE(*, *) "   EphBCL = ", EphBCL
    WRITE(*, *) "   dEHeatL = ", dEHeatL
    WRITE(*, *) "At Right Boundary: "
    WRITE(*, *) "   TBCR = ", TBCR
    WRITE(*, *) "   VphBCR = ", VphBCR
    WRITE(*, *) "   EphBCR = ", EphBCR
    WRITE(*, *) "   dEHeatR = ", dEHeatR
    WRITE(*, *) "RNHeatPh = ", RNHeatPh
    WRITE(*, *) "FNHeatPh = ", FNHeatPh

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
SUBROUTINE calc_inj_heat( T, Hinj, A, dt )
USE ROUTINES, ONLY: energy, Etable
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: T, A, dt
REAL(KIND=8), INTENT(OUT):: Hinj(2)
REAL(KIND=8):: U, V
INTEGER(KIND=4):: i

    DO i = 1, 2
        CALL energy( i, T, U )
        CALL Etable( i, 4, U, V)
        Hinj(i) = U * V * A * dt * 0.25
    ENDDO

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
REAL(KIND=8):: R1, x, Vph, Eph
INTEGER(KIND=4):: i, j, k, mt, eID(3)
    
    PoolL(PoolSize, :, :)%Mat = -1
    PoolR(PoolSize, :, :)%Mat = -1
    
    RNHeatPh = 0

    CALL RAN_NUM( RSeed, dtHeat )
    dtHeat = dtHeat * TimeStep

    x = 0D0
    i = 1
    DO k = 1, Ne(3); DO j = 1, Ne(2)
        mt = ele(i, j, k)%Mat
        HinjectL(j, k) = dEHeatL(mt) + HinjectL(j, k)
        Vph = VphBCL(mt)
        eID = (/ i, j, k /)
        SELECTCASE( PoolL(PoolSize, j, k)%Mat )
        CASE(-1)
            Eph = EphBCL(mt)

            CALL Random_Injection_Control( RNHeatPh, FNHeatPh, &
                                           HeatPhn, HinjectL(j, k), &
                                           x, dL, Vph, Eph, mt, eID, &
                                           1, qL(j, k), RSeed )
        CASE DEFAULT
            R1 = EphBCL(mt) * 0.5

            CALL Perodic_Injection_Control( RNHeatph, FNHeatPh, &
                                            HeatPhn, PoolSize, &
                                            PoolL(:, j, k), &
                                            dtHeat, HinjectL(j, k), &
                                            x, Vph, EphBCL, R1, eID, &
                                            qL(j, k), RSeed )
        END SELECT
    ENDDO; ENDDO

    x = L(1)
    i = Ne(1)
    DO k = 1, Ne(3); DO j = 1, Ne(2)
        mt = ele(i, j, k)%Mat
        HinjectR(j, k) = dEHeatR(mt) + HinjectR(j, k)
        Vph = VphBCR(mt)
        eID = (/ i, j, k /)
        SELECTCASE( PoolR(PoolSize, j, k)%Mat )
        CASE(-1)
            Eph = EphBCR(mt)

            CALL Random_Injection_Control( RNHeatPh, FNHeatPh, &
                                           HeatPhn, HinjectR(j, k), &
                                           x, dL, Vph, Eph, mt, eID, &
                                           -1, qR(j, k), RSeed )
        CASE DEFAULT
            R1 = EphBCR(mt) * 0.5

            CALL Perodic_Injection_Control( RNHeatph, FNHeatPh, &
                                            HeatPhn, PoolSize, &
                                            PoolR(:, j, k), &
                                            dtHeat, HinjectR(j, k), &
                                            x, Vph, EphBCR, R1, eID, &
                                            qR(j, k), RSeed )
        END SELECT
    ENDDO; ENDDO

    DO i = 1, FNHeatph
        IF ( HeatPhn(i).Exist ) THEN
            DO WHILE ( dtHeat(i).gt.0 )
                SELECTCASE( WAY_FlightTime )
                CASE(1)
                    CALL phn_adv_CT( HeatPhn(i), RSeed, dtHeat(i) )
                CASE(2)
                    CALL phn_adv_RT( HeatPhn(i), RSeed, dtHeat(i) )
                END SELECT
            ENDDO
        ENDIF
    ENDDO

END SUBROUTINE Heat_Control


!======================================================================
!----------------------------------------------------------------------
! Perodic Injection Controller
!----------------------------------------------------------------------
SUBROUTINE Perodic_Injection_Control( NHph, S, HPhn, PLSize, PL, dt, &
                                      Hinj, x, V, EphBC, tmp, eID, &
                                      q, RSeed )
USE VAR_TYPES
IMPLICIT NONE
INTEGER(KIND=4), INTENT(IN):: S, PLSize, eID(3)
INTEGER(KIND=4), INTENT(INOUT):: NHph
REAL(KIND=8), INTENT(IN):: x, V, EphBC(2), tmp
REAL(KIND=8), INTENT(INOUT):: dt(S), Hinj, q
TYPE(Phonon), INTENT(INOUT):: HPhn(S)
TYPE(phnPool), INTENT(IN):: PL(PLSize)
TYPE(rng_t), INTENT(INOUT):: RSeed
REAL(KIND=8):: R1, E
INTEGER(KIND=4):: m


    DO WHILE( Hinj.gt.tmp )
        CALL RAN_NUM( RSeed, R1 )
        m = MAX( INT( R1 * PLSize + 0.5 ), 1 )
        NHph = NHph + 1
        E = EphBC(PL(m)%Mat)
        CALL S_Perodic_Injection( HPhn(NHph), dt(NHph), PL(m), &
                                  x, V, E, eID )
        Hinj = Hinj - E
        q = q + E
    ENDDO

END SUBROUTINE Perodic_Injection_Control


!======================================================================
!----------------------------------------------------------------------
! Random Injection Controller
!----------------------------------------------------------------------
SUBROUTINE Random_Injection_Control( NHPh, S, HPhn, Hinj, x, dL, &
                                     V, E, mt, eID, dir, q, RSeed )
USE VAR_TYPES
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: E, x, dL(3), V
REAL(KIND=8), INTENT(INOUT):: Hinj, q
INTEGER(KIND=4), INTENT(IN):: S, mt, eID(3), dir
INTEGER(KIND=4), INTENT(INOUT):: NHPH
TYPE(Phonon), INTENT(INOUT):: HPhn(S)
TYPE(rng_t), INTENT(INOUT):: RSeed
INTEGER(KIND=4):: N, bg, ed, i
!-------------------------------------------------
! dir = 1, Phonons are injected from left boundary
!      -1, ..........................right........
!-------------------------------------------------

    N = INT( Hinj / E + 0.5 )

    IF (N.gt.0) THEN
        Hinj = Hinj - N * E
        q = q + N * E

        bg = NHPh + 1
        ed = NHPh + N

        DO i = bg, ed
            CALL S_Random_Injection( HPhn(i), x, dL, V, &
                                             E, mt, eID, dir, RSeed )
        ENDDO

        NHPh = ed
    ENDIF

END SUBROUTINE Random_Injection_Control


!======================================================================
!----------------------------------------------------------------------
! Random Injection of a single phonon
!----------------------------------------------------------------------
SUBROUTINE S_Random_Injection( HPhm, x, dL, V, Eph, mt, &
                                                      eID, dir, RSeed )
USE VAR_TYPES
USE VAR_Others, ONLY: M_PI
IMPLICIT NONE
TYPE(rng_t), INTENT(INOUT):: RSeed
TYPE(Phonon), INTENT(OUT):: HPhm
REAL(KIND=8), INTENT(IN):: x, dL(3), V, Eph
INTEGER(KIND=4), INTENT(IN):: dir, mt, eID(3)
REAL(KIND=8):: R(4)
!-------------------------------------------------
! dir = 1, Phonons are injected from left boundary
!      -1, ..........................right........
!-------------------------------------------------

    CALL RAN_NUM( RSeed, R )

    R(3) = DBLE( dir ) * DSQRT( R(3) )
    R(4) = R(4) * M_PI * 2D0

    HPhm%xyz(1) = x
    HPhm%xyz(2) = (DBLE( eID(2) - 1 ) + R(1)) * dL(2)
    HPhm%xyz(3) = (DBLE( eID(3) - 1 ) + R(2)) * dL(3)

    HPhm%V = V
    HPhm%Vxyz(1) = HPhm%V * R(3)
    HPhm%Vxyz(2) = HPhm%V * DSQRT(1D0 - R(3)**2) * DCOS( R(4) )
    HPhm%Vxyz(3) = HPhm%V * DSQRT(1D0 - R(3)**2) * DSIN( R(4) )

    HPhm%E = Eph

    HPhm%Org(1) = HPhm%xyz(1)
    HPhm%Org(2) = HPhm%xyz(2)
    HPhm%Org(3) = HPhm%xyz(3)

    HPhm%Dis(1) = 0D0
    HPhm%Dis(2) = 0D0
    HPhm%Dis(3) = 0D0

    HPhm%Mat = mt

    HPhm%eID(1) = eID(1)
    HPhm%eID(2) = eID(2)
    HPhm%eID(3) = eID(3)

    HPhm%Exist = .TRUE.

END SUBROUTINE S_Random_Injection


!======================================================================
!----------------------------------------------------------------------
! Random Injection of a single phonon
!----------------------------------------------------------------------
SUBROUTINE S_Perodic_Injection( HPhm, dt, HPl, &
                                x, Vph, Eph, eID )
USE VAR_TYPES
IMPLICIT NONE
REAL(KIND=8), INTENT(IN):: x, Vph, Eph
REAL(KIND=8), INTENT(INOUT):: dt
INTEGER(KIND=4), INTENT(IN):: eID(3)
TYPE(Phonon), INTENT(OUT):: HPhm
TYPE(phnPool), INTENT(IN):: HPl

    HPhm%xyz(1) = x
    HPhm%xyz(2) = HPl%y
    HPhm%xyz(3) = HPl%z

    HPhm%V = Vph
    HPhm%Vxyz(1) = HPhm%V * HPl%direction(1)
    HPhm%Vxyz(2) = HPhm%V * HPl%direction(2)
    HPhm%Vxyz(3) = HPhm%V * HPl%direction(3)

    HPhm%E = Eph

    HPhm%Org(1) = HPhm%xyz(1)
    HPhm%Org(2) = HPhm%xyz(2)
    HPhm%Org(3) = HPhm%xyz(3)

    HPhm%Dis(1) = 0D0
    HPhm%Dis(2) = 0D0
    HPhm%Dis(3) = 0D0

    HPhm%Mat = HPl%Mat

    HPhm%eID(1) = eID(1)
    HPhm%eID(2) = eID(2)
    HPhm%eID(3) = eID(3)

    HPhm%Exist = .TRUE.

    dt = Hpl%dtRemain

END SUBROUTINE S_Perodic_Injection



!======================================================================
!----------------------------------------------------------------------
! Inject phonons in pools
!----------------------------------------------------------------------

END MODULE HEAT