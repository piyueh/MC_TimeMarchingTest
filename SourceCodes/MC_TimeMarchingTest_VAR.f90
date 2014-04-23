!!#####################################################################
!! Module Name: VAR
!! Purpose: Global variable collection.
!! Author: PY Chuang
!!#####################################################################
MODULE VAR
USE RNG
IMPLICIT NONE
!----------------------------------------------------------------------
! User defined variable types
!----------------------------------------------------------------------
TYPE:: Phonon
    REAL(KIND=8):: xyz(3)
    REAL(KIND=8):: V, Vxyz(3)
    REAL(KIND=8):: E
    REAL(KIND=8):: Org(3), Dis(3)
    INTEGER(KIND=4):: Mat
    INTEGER(KIND=4):: eID(3)
END TYPE Phonon

TYPE:: Element
    REAL(KIND=8):: BD(2, 3)
    REAL(KIND=8):: T, E, N, Eph, Ediff, Vph, MFP
    INTEGER(KIND=4):: Mat
END TYPE Element


!----------------------------------------------------------------------
! Constant parameters
!----------------------------------------------------------------------
REAL(KIND=8), PARAMETER:: M_PI = 3.14159265359D0


!----------------------------------------------------------------------
! Computational model
!       L: Dimensions of the computational model
!----------------------------------------------------------------------
REAL(KIND=8):: L(3)


!----------------------------------------------------------------------
! Material tables
!----------------------------------------------------------------------
REAL(KIND=8):: Ge_start, Si_start
REAL(KIND=8):: dU_Ge, dU_Si
INTEGER(KIND=4):: N_Ge1, N_Ge2, N_Si1, N_Si2
REAL(KIND=8), ALLOCATABLE:: Ge_table(:, :), Si_table(:, :)


!----------------------------------------------------------------------
! Phonons
!       Nph: Number of phonons in computational model.
!       Nprop: Number of properties carried with a phonon.
!              Nprop =
!              1: x, 2: y, 3: z, 4: vx, 5: vy, 6:vz, 7: velocity,
!              8: energy, 9: energy_material, 10: current element ID
!       phn: Phonon property array.  Dimension: Nprop x Nph.
!       phID:
!----------------------------------------------------------------------
REAL(KIND=8):: bundle
INTEGER(KIND=4):: Nph
TYPE(Phonon), ALLOCATABLE:: phn(:)
REAL(KIND=8), ALLOCATABLE:: phID(:)


!----------------------------------------------------------------------
! Element
!       dL: Dimensions of cubic elements
!       dV: Volume of cubic elementa
!       Ne: Numbers of cubic elements in 3 directions (x, y, z)
!----------------------------------------------------------------------
REAL(KIND=8):: dL(3), dV
INTEGER(KIND=4):: Ne(3)
TYPE(Element), ALLOCATABLE:: ele(:, :, :)
INTEGER(KIND=4), ALLOCATABLE:: N(:, :, :), bg(:, :, :), ed(:, :, :)


!----------------------------------------------------------------------
! Variable related to random number generation
!----------------------------------------------------------------------
TYPE(rng_t), ALLOCATABLE:: SeedMP(:)


!----------------------------------------------------------------------
! Variable related computation
!----------------------------------------------------------------------
REAL(KIND=8):: dt
INTEGER(KIND=4):: NCores
INTEGER(KIND=4):: iter, iter0, iterations


!----------------------------------------------------------------------
!----------------------------------------------------------------------
END MODULE VAR
