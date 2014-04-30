!!#####################################################################
!! Projrct: 3D Monte-Carlo Simulator Using Structured Grids
!! Program: Solver & PreProcessing
!! File Name: MC3D_Structured_VAR
!! Contenent: Collections of Public Variable Modules
!! Author: PY Chuang
!!#####################################################################


!======================================================================
!----------------------------------------------------------------------
! All other variables
!----------------------------------------------------------------------
MODULE VAR_Others
IMPLICIT NONE

    !------------------------------------------------------------------
    ! Constant parameters
    !------------------------------------------------------------------
    REAL(KIND=8), PARAMETER:: M_PI = 3.14159265359D0
    REAL(KIND=8), PARAMETER:: kB = 8.617D-2
    REAL(KIND=8), PARAMETER:: barh = 6.5822D-1 ! meV/K, meV.ps

    !------------------------------------------------------------------
    ! GammaT: Constant Pseudo Scattering Rate
    ! WAY_FlightTime: 1: Constant Free Flight Time
    !                 2: Random Free Flight Time
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: GammaT
    INTEGER(KIND=4), SAVE:: WAY_FlightTime

    !------------------------------------------------------------------
    ! Variable related computation
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: TimeStep, time
    INTEGER(KIND=4), SAVE:: iter, iter0, iterations

END MODULE VAR_Others


!======================================================================
!----------------------------------------------------------------------
! User defined variable types
!----------------------------------------------------------------------
MODULE VAR_TYPES
IMPLICIT NONE

    TYPE:: Phonon
        REAL(KIND=8):: xyz(3)
        REAL(KIND=8):: V, Vxyz(3)
        REAL(KIND=8):: E
        REAL(KIND=8):: Org(3), Dis(3)
        INTEGER(KIND=4):: Mat
        INTEGER(KIND=4):: eID(3)
        LOGICAL:: Exist
    END TYPE Phonon

    TYPE:: Element
        REAL(KIND=8):: BD(2, 3)
        REAL(KIND=8):: T, E, ND, Eph, Ediff, Vph, MFP, SCR
        INTEGER(KIND=4):: Mat
        INTEGER(KIND=4):: Ntol, Nbg, Ned
    END TYPE Element

    TYPE:: phnPool
        REAL(KIND=8):: y, z
        REAL(KIND=8):: direction(3)
        REAL(KIND=8):: dtRemain
        INTEGER(KIND=4):: Mat
    END TYPE phnPool

END MODULE VAR_TYPES


!======================================================================
!----------------------------------------------------------------------
! Materials
!----------------------------------------------------------------------
MODULE VAR_MAT
IMPLICIT NONE

    REAL(KIND=8), SAVE:: Ge_start, Si_start
    REAL(KIND=8), SAVE:: dU_Ge, dU_Si
    INTEGER(KIND=4), SAVE:: N_Ge1 = 7, N_Ge2 = 8001
    INTEGER(KIND=4), SAVE:: N_Si1 = 7, N_Si2 = 8001
    REAL(KIND=8), SAVE:: Ge_table(7, 8001)
    REAL(KIND=8), SAVE:: Si_table(7, 8001)

END MODULE VAR_MAT


!======================================================================
!----------------------------------------------------------------------
! Domains and Elements
!----------------------------------------------------------------------
MODULE VAR_SPACES
USE VAR_TYPES
IMPLICIT NONE

    !------------------------------------------------------------------
    ! Computational model
    !       L: Dimensions of the computational model
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: L(3)

    !------------------------------------------------------------------
    ! Element
    !       dL: Dimensions of cubic elements
    !       dV: Volume of cubic elementa
    !       Ne: Numbers of cubic elements in 3 directions (x, y, z)
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: dL(3), dV
    INTEGER(KIND=4), SAVE:: Ne(3)
    INTEGER(KIND=4), SAVE:: NEmpty
    INTEGER(KIND=8), SAVE:: iMid
    TYPE(Element), ALLOCATABLE, SAVE:: ele(:, :, :)
    INTEGER(KIND=4), ALLOCATABLE, SAVE:: NAdd(:, :, :)
    REAL(KIND=8), ALLOCATABLE, SAVE:: Tavg(:, :, :), Eavg(:, :, :)

END MODULE VAR_SPACES


!======================================================================
!----------------------------------------------------------------------
! Phonons
!----------------------------------------------------------------------
MODULE VAR_ph
USE VAR_TYPES
IMPLICIT NONE

    !------------------------------------------------------------------
    ! Nph: Number of phonons in computational model.
    ! Nprop: Number of properties carried with a phonon.
    !        Nprop =
    !        1: x, 2: y, 3: z, 4: vx, 5: vy, 6:vz, 7: velocity,
    !        8: energy, 9: energy_material, 10: current element ID
    ! phn: Phonon property array.  Dimension: Nprop x Nph.
    ! phID:
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: bundle
    INTEGER(KIND=4), SAVE:: RNph, FNph
    TYPE(Phonon), ALLOCATABLE, SAVE:: phn(:)
    INTEGER(KIND=4), ALLOCATABLE, SAVE:: phId(:), EmptyID(:)

END MODULE VAR_ph


!======================================================================
!----------------------------------------------------------------------
! BCs and Heat Control
!----------------------------------------------------------------------
MODULE VAR_BC
USE VAR_TYPES
IMPLICIT NONE

    !------------------------------------------------------------------
    ! Boundary Conditions
    !       BC = 1: Adiabatic
    !            2: Periodic
    !            3: Heat control
    !------------------------------------------------------------------
    INTEGER(KIND=4), SAVE:: BCs(3)

    !------------------------------------------------------------------
    ! Variables related to heat control
    ! TBCL, TBCR: Temperature on boundaries
    ! QBCL, QBCR: Net heat through boundaries
    ! qfluxL, qfluxC, qfluxR: Net heat flux through boundaries
    ! VBCL, VBCR, EBCL, EBCR: The property of injected phonons of each
    !                         material
    ! EinjectL, EinjectR: The energy which should be injected from
    !                     boundary elements
    ! qL, qC, qR: The energy carried by phonons transfer through
    !             boundaries and central section elements
    !------------------------------------------------------------------
    REAL(KIND=8), SAVE:: TBCL, TBCR
    REAL(KIND=8), SAVE:: QBCL, QBCR
    REAL(KIND=8), SAVE:: qfluxL, qfluxC, qfluxR
    REAL(KIND=8), SAVE:: VphBCL(2), VphBCR(2), EphBCL(2), EphBCR(2)
    REAL(KIND=8), ALLOCATABLE, SAVE:: EinjectL(:, :), EinjectR(:, :)
    REAL(KIND=8), ALLOCATABLE, SAVE:: qL(:, :), qC(:, :), qR(:, :)

    !------------------------------------------------------------------
    ! Variables with regard to random heat injection method
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Variables with regard to periodic heat injection method
    !------------------------------------------------------------------
    INTEGER(KIND=4), SAVE:: PoolSize
    TYPE(phnPool), ALLOCATABLE, SAVE:: PoolL(:, :, :), PoolR(:, :, :)
    INTEGER(KIND=4), ALLOCATABLE, SAVE:: iNPoolL(:, :), iNPoolR(:, :)

END MODULE VAR_BC


!======================================================================
!----------------------------------------------------------------------
! Variable related to output
!----------------------------------------------------------------------
MODULE VAR_Output
USE VAR_TYPES
USE VAR_Others
USE VAR_SPACES
USE VAR_ph
USE VAR_BC
IMPLICIT NONE

    INTEGER(KIND=4), SAVE:: nOutput
    CHARACTER(LEN=72), SAVE:: CaseName, OutputFileName
    CHARACTER(LEN=72), SAVE:: InputFileName, RestartFileName
    NAMELIST /DataOutput/ iter, time, L, dL, Ne, Tavg, Eavg
    NAMELIST /Initial_1/ TimeStep, time, iter0, L, Ne, dL, dV, bundle, &
                         FNph, RNph, PoolSize
    NAMELIST /Initial_2/ ele, phn, PoolL, PoolR, iNPoolL, iNPoolR

END MODULE VAR_Output


!======================================================================
!----------------------------------------------------------------------
! Use all variables
!----------------------------------------------------------------------
MODULE VAR_ALL
USE VAR_TYPES
USE VAR_MAT
USE VAR_Output
USE VAR_SPACES
USE VAR_ph
USE VAR_BC
USE VAR_Others
IMPLICIT NONE
END MODULE VAR_ALL