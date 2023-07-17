
SUBROUTINE ZERO(X,N)
    use global
    REAL X(N)
    DO 10 I=1,N
        X(I) = 0.0
10 CONTINUE
END
