!!#####################################################################
!! Module Name: variables
!! Purpose: Global variable collection.
!! Author: PY Chuang
!!#####################################################################
MODULE variables
USE RNG
IMPLICIT NONE


!----------------------------------------------------------------------
! Computational model
!       L: Dimensions of the computational model
!----------------------------------------------------------------------
REAL(KIND=8):: L(3)



!----------------------------------------------------------------------
! Phonons
!       Nph: Number of phonons in computational model.
!       Nprop: Number of properties carried with a phonon.
!              Nprop =
!              1: x, 2: y, 3: z, 4: vx, 5: vy, 6:vz, 7: velocity,
!              8: energy, 9: energy_material, 10: current element ID
!              11: displacement
!       phn: Phonon property array.  Dimension: Nprop x Nph.
!       phID:
!----------------------------------------------------------------------
INTEGER(KIND=4):: Nph
INTEGER(KIND=4):: Nprop
REAL(KIND=8), ALLOCATABLE:: phn(:, :)
REAL(KIND=8), ALLOCATABLE:: phID(:)


!----------------------------------------------------------------------
! Element
!       dL: Dimensions of cubic elements
!       Ne: Numbers of cubic elements in 3 directions (x, y, z)
!       Te: Current temperature of elements
!       Ee: Current total energy of phonons located in the elements
!       EeDiff: Current energy difference of elements
!       eID: The relationship between scalar ID and 3D index (i, j, k)
!       IDe: The relationship between 3D index (i, j, k) and scalar ID
!----------------------------------------------------------------------
REAL(KIND=8):: dL(3)
INTEGER(KIND=4):: Ne(3)
REAL(KIND=8), ALLOCATABLE:: Te(:), Ee(:), EeDiff(:)
INTEGER(KIND=4), ALLOCATABLE:: eID(:,:), IDe(:, :, :)


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
END MODULE variables
