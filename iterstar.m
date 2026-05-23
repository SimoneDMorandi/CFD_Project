function ITERSTAR(A,B,C)
        
        USE PIGAMMA_PAR
        
        INTEGER KIP
        REAL*8 TRY1,TRY2,A,B,C
        
        TRY1 = .5*A/GB
        KIP = 0
1       CONTINUE
        KIP = KIP+1
        IF(KIP.GE.500) THEN
        WRITE(*,*)'ITERATION FAILS IN ITERSTAR !!!!!'
        STOP
        ENDIF
C       TRY2 = B-A*DLOG(TRY1)
        TRY2 = EXP((B-TRY1)/A)
        IF(ABS(TRY2-TRY1).LT.1.E-5) GO TO 2
        TRY1 = TRY2
        GO TO 1
2       CONTINUE
        C = TRY2
        RETURN
end
