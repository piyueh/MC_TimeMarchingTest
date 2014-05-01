PROGRAM main
IMPLICIT NONE
TYPE test
    REAL(KIND=8):: E
    LOGICAL:: Exist
END TYPE test
TYPE(test):: TT(5000000)
REAL(KIND=8):: R(2), time(2)
INTEGER(KIND=4):: i, tmp(5000000), N, tmp2(5000000)
    
    tmp2 = (/ (i, i=1,5000000) /)
    
    CALL RANDOM_SEED()
    DO i = 1, 5000000
        CALL RANDOM_NUMBER( R )
        TT(i)%E = R(1) * 100D0
        TT(i)%Exist = INT( R(2) * 1.9 + 0.5 + 0.5 )
    ENDDO
    
    tmp = -1
    
    CALL CPU_TIME( time(1) )
    N = COUNT( TT%Exist )
    CALL CPU_TIME( time(2) )
    WRITE(*, *) N, time(2) - time(1)
    
    CALL CPU_TIME( time(1) )
    N = 0
    DO i = 1, 5000000
        IF (TT(i)%Exist) THEN
            N = N + 1
        ENDIF
    ENDDO
    CALL CPU_TIME( time(2) )
    WRITE(*, *) N, time(2) - time(1)

END PROGRAM main