SUBROUTINE ReadData1( i, ct, L, dL, Ne, qfluxL, qfluxC, qfluxR )
IMPLICIT NONE
REAL(KIND=8), INTENT(OUT) :: L(3), dL(3), qfluxL, qfluxC, qfluxR
INTEGER(KIND=4), INTENT(OUT):: ct, Ne(3)
INTEGER(KIND=4), INTENT(IN):: i
REAL(KIND=8):: time
INTEGER(KIND=4):: iter
CHARACTER(LEN=72):: FileName
NAMELIST /DataOutput_1/ iter, time, L, dL, Ne, qfluxL, qfluxC, qfluxR

    WRITE(FileName, '(I8.8, ".txt")') i
    OPEN( UNIT = 150, FILE = FileName )
    READ( UNIT = 150, FMT = *) ct
    READ( UNIT = 150, NML = DataOutput_1 )
    CLOSE( 150 )

END SUBROUTINE ReadData1

!======================================================================

SUBROUTINE ReadData2( i, Nx, Ny, Nz, qL, qC, qR, Tavg, Eavg)
IMPLICIT NONE
INTEGER(KIND=4), INTENT(IN):: i
INTEGER(KIND=4), INTENT(IN):: Nx, Ny, Nz
REAL(KIND=8), INTENT(OUT):: Tavg(Nx, Ny, Nz)
REAL(KIND=8), INTENT(OUT):: Eavg(Nx, Ny, Nz)
REAL(KIND=8), INTENT(OUT):: qL(Ny, Nz)
REAL(KIND=8), INTENT(OUT):: qC(Ny, Nz)
REAL(KIND=8), INTENT(OUT):: qR(Ny, Nz)
CHARACTER(LEN=72):: FileName
NAMELIST /DataOutput_2/ qL, qC, qR, Tavg, Eavg
    
    WRITE(FileName, '(I8.8, ".txt")') i
    OPEN( UNIT = 150, FILE = FileName )
    READ( UNIT = 150, NML = DataOutput_2 )
    CLOSE( 150 )

END SUBROUTINE ReadData2
