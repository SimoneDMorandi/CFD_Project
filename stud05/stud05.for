C       PROGRAM STUD05 (ESPERIMENTI 1D e QUASI 1D).

C       (gas perfetto (gamma=1.4), variabili P,U,H)

        INCLUDE 'COMST05.INC'
        COMMON /CONVERG/ RTRMS,RTMAX,NRT
144     FORMAT(5X,I4,4X,F5.2,4X,F7.4,4X,F7.4,4X,F7.4)
        GAMMA=1.4
        PI=4.*ATAN(1.)
        CONV=180./PI
        GA=GAMMA/(GAMMA-1)
        GB=1./(GAMMA-1.)
        GC=(GAMMA+1.)/(GAMMA-1.)
        GD=(GAMMA-1.)/2.
        GE=(GAMMA+1.)/2.
        GF=SQRT(GAMMA)
        GG=2./(GAMMA-1.)
        GH=(GAMMA+1.)/(2.*GAMMA)
        GI=(GAMMA-1.)/(2.*GAMMA)
        GJ=(GAMMA-1.)/GAMMA

        WRITE(*,*)'PROGRAMMA STUD05 (FDS - 2nd Order)'
        K=0
        TIME=.0

        CALL INIT

        WRITE(*,*)'IORD (   1 = first  order )'
        WRITE(*,*)'     (  21 = second order (ENO 1st level))'
        WRITE(*,*)'     (  22 = second order (ENO 2nd level))'
        READ(*,*)IORD

        NMID = NA/2+2

        IF(ITEST.EQ.4) THEN
c       WRITE(*,*)'dare il numero J di identificazione '
c       WRITE(*,*)'del file per DIA05.FOR              '
c       READ(*,*)J
c       WRITE(10+J,*)NA,IORD,KA
        ENDIF

        CALL OUTPUT

1       K=K+1
        IF(ITEST.EQ.1)                    WRITE(*,*)'K = ',K
        IF(ITEST.EQ.2.AND.MOD(K,50).EQ.0) WRITE(*,*)'K = ',K
        IF(ITEST.EQ.3.AND.MOD(K,50).EQ.0) WRITE(*,*)'K = ',K
        IF(ITEST.EQ.4)                    WRITE(*,*)'K = ',K
        IF(ITEST.EQ.5.AND.MOD(K,10).EQ.0) WRITE(*,*)'K = ',K
        IF(ITEST.EQ.6.AND.MOD(K,50).EQ.0) WRITE(*,*)'K = ',K
        IF(ITEST.EQ.7.AND.MOD(K,50).EQ.0) WRITE(*,*)'K = ',K
        IF(ITEST.EQ.8.AND.MOD(K,50).EQ.0) WRITE(*,*)'K = ',K

        IF((ITEST.EQ.2.OR.ITEST.EQ.3).AND.MOD(K,100).EQ.0) THEN
        WRITE(*,*)'=',K,' NRT RTRMS = ',NRT,RTRMS
        ENDIF

        DT=1.E5
        DO 2 N=2,NCM
        ALAM=(ABS(U(N))+A(N))
        DTD=DX/ALAM
        IF(DTD.LT.DT) DT=DTD
2       CONTINUE

        DT=DT*STAB
        IF(ITEST.EQ.1) DT = DTODX * DX
        IF(ITEST.EQ.4) DT = DTITEST4
        TIME=TIME+DT

        IF(IORD.EQ.21) CALL ENO1
        IF(IORD.EQ.22) CALL ENO2

        CALL SPLIT
        CALL MARCH

        if(itest.eq.1) then
        sum1=0.0
        sum2=0.0
        sum3=0.0
        do 3333 n=2,ncm
        sum1=sum1+w1(n)
        sum2=sum2+w2(n)
        sum3=sum3+w3(n)
3333    continue
        sum1=sum1/na
        sum2=sum2/na+(p(ncm)-p(2))*time
        sum3=sum3/na
        write(*,*) k,sum1,sum2,sum3
        endif

        IF(ITEST.EQ.4) THEN
        IF(K.EQ.1)          WRITE(*,*)'K,TIME P U TTOT at NMID  '
        IF(MOD(K,100).EQ.0) WRITE(*,*)'K,TIME P U TTOT at NMID  '
        KPRINT=4*(NA+1)/40
        IF(MOD(K,KPRINT).EQ.0)
     1          WRITE(6,144) K,TIME,P(NMID),U(NMID),TTOT(NMID)
        WRITE(10+J,*)K,TIME,P(NMID),U(NMID),TTOT(NMID)
        ENDIF

        IF(MOD(K,KOUT).EQ.0) THEN

        IF(ITEST.EQ.2) THEN
        NCMM=NCM-1
        DO 200 N=2,NCMM
        DIFF=(AMACH(N+1)-1.)*(AMACH(N)-1.)
        IF(DIFF.GT.0) GO TO 200
        XSH = X(N)+
     1   (1.-AMACH(N))/(AMACH(N+1)-AMACH(N))*(X(N+1)-X(N))
        XSH=XSH/C
200     CONTINUE
        WRITE(*,*)'the shock is located at XSH/C = ',XSH
        ENDIF
        CALL OUTPUT
        ENDIF
        IF(K.LT.KA) GO TO 1
        WRITE(*,*)'KA,KOUT (old values = ',KA,KOUT,'  )'
        READ(*,*)KA,KOUT
        WRITE(*,*)'old STAB was ',STAB,'  give the new one'
        READ(*,*)STAB
        WRITE(*,*)'old IORD was ',IORD,'  give the new one'
        READ(*,*)IORD

        IF(ITEST.EQ.2) THEN
        WRITE(*,*)'if you like, give the new PEX  '
        WRITE(*,*)'PEXINF,PEXSUP PEXPP=', PEXINF,PEXSUP,PEXPP
        WRITE(*,*)'the old PEX was ',PEX
        READ(*,*)PEX
        ENDIF

        IF(ITEST.EQ.3) THEN
        WRITE(*,*)'if you like, give the new PEX '
        WRITE(*,*)'PEXINF,PEXSUP PEXPP=', PEXINF,PEXSUP,PEXPP
        WRITE(*,*)'the old PEX was ',PEX
        READ(*,*)PEX
        ENDIF

        GO TO 1
        END

C....................................................................

        SUBROUTINE INIT

        INCLUDE 'COMST05.INC'
        COMMON / ENO / PXENO(NMAX),UXENO(NMAX),HXENO(NMAX),
     1           PPXENO(NMAX),HHXENO(NMAX),
     2           PPT(NMAX),UT(NMAX),HHT(NMAX),
     3           PPTT(NMAX),UTT(NMAX),HHTT(NMAX)

        WRITE(*,*)'esempio numerico: ITEST = 1 (shock tube)'
        WRITE(*,*)'        ITEST = 2 (supersonic nozzle)'
        WRITE(*,*)'        ITEST = 3 (transonic nozzle)'
        WRITE(*,*)'        ITEST = 4 (discharge from a reservoir)'
        WRITE(*,*)'        ITEST = 5 (expansion at end of a duct)'
        WRITE(*,*)'        ITEST = 6 (opposite shocks interaction)'
        WRITE(*,*)'        ITEST = 7 (same-fam.shocks interaction)'
        WRITE(*,*)'        ITEST = 8 (shock contact-surface inter.)'
        READ(*,*)ITEST
        IF(ITEST.EQ.1) GO TO 1111
        IF(ITEST.EQ.2) GO TO 2222
        IF(ITEST.EQ.3) GO TO 3333
        IF(ITEST.EQ.4) GO TO 4444
        IF(ITEST.EQ.5) GO TO 5555
        IF(ITEST.EQ.6) GO TO 6666
        IF(ITEST.EQ.7) GO TO 7777
        IF(ITEST.EQ.8) GO TO 8888

1111    CONTINUE
        DTODX = 0.411
        NA    = 100
        KA    =  35
        RATRHO=   8.
        RATP  =  10.
        WRITE(*,*)'Test case di Sod (Dt/Dx=0.411 , NA=100 , KA=35)'
        WRITE(*,*)'for double number of points (NA=200), give:    '
        WRITE(*,*)'                  Dt/Dx=0.411 , NA=200 , KA=70)'
        WRITE(*,*)'type : DTODX NA KA '
        READ(*,*)DTODX,NA,KA
        KOUT=KA
        WRITE(*,*)'ratrho , ratp (Sod Test RATRHO=8. RATP=10.)'
        READ(*,*)RATRHO,RATP
        WRITE(*,*)'STAB'
        READ(*,*)STAB
        PMIN   =1.
        RHOMIN =1.
        TMIN   =PMIN/RHOMIN
        PMAX   =PMIN*RATP
        RHOMAX =RHOMIN*RATRHO
        TMAX   =PMAX/RHOMAX
        AMIN   =GF*SQRT(TMIN)
        AMAX   =GF*SQRT(TMAX)
        SMIN   =ALOG(PMIN)-GAMMA*ALOG(RHOMIN)
        SMAX   =ALOG(PMAX)-GAMMA*ALOG(RHOMAX)
        EMIN   =GB*PMIN
        EMAX   =GB*PMAX
        NC=NA+2
        NCM=NC-1
        NHIGH=NA/2+2
        B=.0
        C=1.0
        DX=(C-B)/NA
        DO 110 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)    =1.0
        AMEDIA(N)  =1.0
        DIFFAREA(N)=0.0
C       WRITE(*,*)N,X(N),AREA(N)
110     CONTINUE
        DO 120 N=2,NCM
        A(N)   = AMIN
        P(N)   = PMIN
        RHO(N) = RHOMIN
        T(N)   = TMIN
        S(N)   = SMIN
        E(N)   = EMIN
        IF(N.LT.NHIGH) THEN
        A(N)   = AMAX
        P(N)   = PMAX
        RHO(N) = RHOMAX
        T(N)   = TMAX
        S(N)   = SMAX
        E(N)   = EMAX
        ENDIF

        U(N)   = .0
        AMACH(N)=.0
        H(N)  =3.5*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
120     CONTINUE
        RETURN
2222    CONTINUE
        WRITE(*,*)'NA,STAB'
        READ(*,*)NA,STAB
        WRITE(*,*)'DIVERGENCE OF THE NOZZLE (AEND/AINI)'
        READ(*,*)DIVERG
        WRITE(*,*)'ACHIN'
        READ(*,*)ACHIN
        WRITE(*,*)'KA,KOUT'
        READ(*,*)KA,KOUT
        B = 0.
        C = 1.
        NC=NA+2
        NCM=NC-1
        DX=(C-B)/NA
        DO 210 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.+(DIVERG-1.)*X(N)/(C-B)
C       WRITE(*,*)N,X(N),AREA(N)
210     CONTINUE
C
        TIN    = 1./(1.+.2*ACHIN**2)
        AREAIN = 1.
        AREAEN = AREAIN*DIVERG
        AMEDIA(1)=AREAIN
        COST=.2*(2./2.4)**6
        ATHR = AREAIN*SQRT(TIN**5*(1.-TIN)/COST)
        DO 220 N=2,NCM
        AMEDIA(N)  =.5*(AREA(N+1)+AREA(N))
        DIFFAREA(N)= AMEDIA(N)-AMEDIA(N-1)

        PTOTS= 1.
        RATAREA=AREA(N)/ATHR
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        T(N)=TEMP
        P(N)=T(N)**3.5*PTOTS
        AMACH(N)=SQRT(5.*(1./T(N)-1.))
        A(N)=SQRT(1.4*T(N))
        U(N)=SQRT(7.*(1.-T(N)))
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
220     CONTINUE
        RATAREA=AREAEN/ATHR
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TDUM=TEMP
        PDUM=TDUM**3.5*PTOTS
        AMACHDUM=SQRT(5.*(1./TDUM-1.))
        AM2=AMACHDUM**2
        CALL SHOCK(AM2,RATP,RATR,RATPT)
        PEXINF=PDUM
        PEXSUP=PDUM*RATP
        WRITE(*,*)'PRESSIONI DI USCITA (PEXINF,PEXSUP)'
        WRITE(*,*)PEXINF,PEXSUP

        CALL SHOCK (ACHIN*ACHIN,AA,BB,CC)
        CALL SUB (COST,TT,AREAEN/ATHR,CC)
        PEXPP = TT**3.5*CC
        WRITE(*,*) 'PRESSIONE USCITA MAX =',PEXPP

        WRITE(*,*)'PRESSIONE DI USCITA'
        READ(*,*)PEX
        XTRY=-1.E4
        IF(PEX.GT.PEXSUP) GO TO 230
        IF(PEX.LE.PEXSUP) GO TO 240
230     CONTINUE
        KIP=0
        XTRY=0
        DXTRY=.2
260     CONTINUE
        KIP=KIP+1
        if(kip.gt.1000) WRITE(*,*)'pressione a valle troppo alta'
        if(kip.gt.1000) stop
        XTRY=XTRY+DXTRY
        AREADUM=1.+(DIVERG-1.)*XTRY/(C-B)
        RATAREA=AREADUM/ATHR
        PTOTS=1.
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TDUM= TEMP
        TU  = TDUM
        AMADUM=SQRT(5.*(1./TDUM-1.))
        AM2=AMADUM**2
        CALL SHOCK(AM2,RATP,RATR,RATPT)
        PTOTS=RATPT
        RATAREA=AREAEN/ATHR
        CALL SUB(COST,TTRY,RATAREA,PTOTS)
        PTRY=PTOTS*(TTRY/1.)**3.5
C       WRITE(*,*)'KIP,XTRY,PTRY,PTOTS',KIP,XTRY,PTRY,PTOTS
        IF(ABS(PTRY-PEX)/PEX.LT.1.E-5) GO TO 250
        IF(PTRY.GT.PEX) GO TO 260
        XTRY=XTRY-DXTRY
        DXTRY=DXTRY/2.
        GO TO 260
250     CONTINUE
        DO 270 N=2,NCM
        IF(X(N).LT.XTRY) GO TO 270
        RATAREA=AREA(N)/ATHR
        CALL SUB(COST,TEMP,RATAREA,PTOTS)
        T(N) = TEMP
        P(N) = PTOTS*(TEMP/1.)**3.5
        AMACH(N)=SQRT(5.*(1./T(N)-1.))
        A(N)=SQRT(1.4*T(N))
        U(N)=SQRT(7.*(1.-T(N)))
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
270     CONTINUE

        P(NC) = PEX
        S(NC) = S(NCM)
        T(NC) = T(NCM)*(PEX/P(NCM))**(1./3.5)
        AMACH(NC)=SQRT(5.*(1./T(NC)-1.))
        A(NC)=SQRT(1.4*T(NC))
        U(NC)=SQRT(7.*(1.-T(NC)))
        RHO(NC)=P(NC)/T(NC)
        E(NC)  =GB*P(NC)+.5*RHO(NC)*U(NC)**2
        H(NC)  =GA*T(NC)
        PTOT(NC)  =P(NC)*(1.+GD*AMACH(NC)**2)**GA
        TTOT(NC)  =T(NC)*(1.+GD*AMACH(NC)**2)
        FLOW(NC)= RHO(NC)*U(NC)*AMEDIA(NCM)
        FLHT(NC)= U(NC)*(P(NC)+E(NC))*AMEDIA(N)

        IF(XTRY.GE.0) THEN
        PU=TU**3.5
        RHOU=TU**2.5
        SU=.0
        UU=SQRT((1.-TU)*7.)
        AU=SQRT(1.4*TU)
        AMU=UU/AU
        TTOTU=TU+1./7.*UU**2
        PTOTU=PU*(1.+.2*AMU**2)**3.5
        AMU2=AMU**2
        CALL SHOCK(AMU2,RATP,RATR,RATPT)
        PD=PU*RATP
        RHOD=RHOU*RATR
        TD=PD/RHOD
        SD=-.4*ALOG(RATPT)
        UD=SQRT((1.-TD)*7.)
        AD=SQRT(1.4*TD)
        AMD=UD/AD
        TTOTD=TD+1./7.*UD**2
        PTOTD=PD*(1.+.2*AMD**2)**3.5
        XSH=XTRY
        WRITE(*,*)'PROPRIETA DELLO URTO'
        WRITE(*,*)'XSH,PU,AMU,SU,PTOTU'
        WRITE(*,*) XSH,PU,AMU,SU,PTOTU
        WRITE(*,*)'XSH,PD,AMD,SD,PTOTD'
        WRITE(*,*) XSH,PD,AMD,SD,PTOTD
        ENDIF
240     CONTINUE
        RETURN

3333    CONTINUE
        WRITE(*,*)'NA,STAB'
        READ(*,*)NA,STAB
        WRITE(*,*)'KA,KOUT'
        READ(*,*)KA,KOUT
C
C       AREA=CCC1*(X(N)+X000)+CCC2/(X(N)+X000)
C
        CCC1=2.5
        CCC2=0.3
        X000=.1
C
        B = 0.
        C = 1.
        NC=NA+2
        NCM=NC-1
        DX=(C-B)/NA
        DO 310 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=CCC1*(X(N)+X000)+CCC2/(X(N)+X000)
C       WRITE(*,*)N,X(N),AREA(N)
310     CONTINUE

        AREAIN = CCC1*(0.00+X000)+CCC2/(0.00+X000)
        AREAEN = CCC1*(1.00+X000)+CCC2/(1.00+X000)
        XTHR   = SQRT(CCC2/CCC1)-X000
        ATHR   = CCC1*(XTHR+X000)+CCC2/(XTHR+X000)

        WRITE(*,*)'Convergent-divergent Nozzle'
        WRITE(*,*)'AREA=2.5*(X+.1)+0.3/(X+.1)'
        WRITE(*,*)'ATHR at XTRH  ',ATHR,XTHR
        WRITE(*,*)'AREAIN AREAEN ',AREAIN,AREAEN

        PTOTS= 1.
        COST=.2*(2./2.4)**6
        AMEDIA(1)  = AREAIN
        DO 320 N=2,NCM
        AMEDIA(N)  =.5*(AREA(N+1)+AREA(N))
        IF(N.EQ.NCM) AMEDIA(N)=AREAEN
        DIFFAREA(N)= AMEDIA(N)-AMEDIA(N-1)
320     CONTINUE

        PTOTS= 1.
        RATAREA=AREAEN/ATHR
        CALL SUB (COST,TEMP,RATAREA,PTOTS)
        PASUB=TEMP**3.5
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        PASUP=TEMP**3.5
        AMEN2=5.*((1./PASUP)**(1./3.5)-1.)
        PINT=PASUP*(7.*AMEN2-1.)/6.
        WRITE(*,*)'PRESSIONI DI ADATTAMENTO'
        WRITE(*,*)'PASUB,PINT,PASUP'
        WRITE(*,*)PASUB,PINT,PASUP
        WRITE(*,*)'PRESSIONE DI USCITA'
        READ(*,*)PEX
        PEXINF = PASUP
        PEXSUP = PINT
        PEXPP  = PASUB
        IF(PEX.GT.PASUB) GO TO 311
        IF(PEX.LE.PASUB) GO TO 312
311     CONTINUE
        TEMP=PEX**(1./3.5)
        FLOWDUM=AREAEN*SQRT(7.)*TEMP**2.5*SQRT(1.-TEMP)
        AREF=FLOWDUM/((2./2.4)**2.5*(2.8/2.4)**.5)
        DO 305 N=2,NCM
        RATAREA=AREA(N)/AREF
        CALL SUB (COST,TEMP,RATAREA,PTOTS)
        T(N) = TEMP
        P(N) = PTOTS*(TEMP/1.)**3.5
        AMACH(N)=SQRT(5.*(1./T(N)-1.))
        A(N)=SQRT(1.4*T(N))
        U(N)=SQRT(7.*(1.-T(N)))
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
305     CONTINUE
        GO TO 350
312     CONTINUE
        DO 306 N=2,NCM
        RATAREA=AREA(N)/ATHR
        IF(X(N).LT.XTHR) CALL SUB  (COST,TEMP,RATAREA,PTOTS)
        IF(X(N).GE.XTHR) CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        T(N) = TEMP
        P(N) = PTOTS*(TEMP/1.)**3.5
        AMACH(N)=SQRT(5.*(1./T(N)-1.))
        A(N)=SQRT(1.4*T(N))
        U(N)=SQRT(7.*(1.-T(N)))
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
306     CONTINUE
        IF(PEX.LE.PINT) GO TO 350

        KIP=0
        XTRY=XTHR
        DXTRY=.1
360     CONTINUE
        KIP=KIP+1
        if(kip.gt.1000) WRITE(*,*)'pressione a valle troppo alta'
        if(kip.gt.1000) stop
        XTRY=XTRY+DXTRY
        AREADUM=CCC1*(XTRY+X000)+CCC2/(XTRY+X000)
        RATAREA=AREADUM/ATHR
        PTOTS=1.
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TDUM= TEMP
        TU  = TDUM
        AMADUM=SQRT(5.*(1./TDUM-1.))
        AM2=AMADUM**2
        CALL SHOCK(AM2,RATP,RATR,RATPT)
        PTOTS=RATPT
        RATAREA=AREAEN/ATHR
        CALL SUB(COST,TTRY,RATAREA,PTOTS)
        PTRY=PTOTS*(TTRY/1.)**3.5
C       WRITE(*,*)'KIP,XTRY,PTRY,PTOTS',KIP,XTRY,PTRY,PTOTS
        IF(ABS(PTRY-PEX)/PEX.LT.1.E-5) GO TO 361
        IF(PTRY.GT.PEX) GO TO 360
        XTRY=XTRY-DXTRY
        DXTRY=DXTRY/2.
        GO TO 360
361     CONTINUE
        DO 362 N=2,NCM
        IF(X(N).LT.XTRY) GO TO 362
        RATAREA=AREA(N)/ATHR
        CALL SUB(COST,TEMP,RATAREA,PTOTS)
        T(N) = TEMP
        P(N) = PTOTS*(TEMP/1.)**3.5
        AMACH(N)=SQRT(5.*(1./T(N)-1.))
        A(N)=SQRT(1.4*T(N))
        U(N)=SQRT(7.*(1.-T(N)))
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
362     CONTINUE
        IF(XTRY.GE.XTHR.AND.XTRY.LT.1.0) THEN
        PU=TU**3.5
        RHOU=TU**2.5
        SU=.0
        UU=SQRT((1.-TU)*7.)
        AU=SQRT(1.4*TU)
        AMU=UU/AU
        TTOTU=TU+1./7.*UU**2
        PTOTU=PU*(1.+.2*AMU**2)**3.5
        AMU2=AMU**2
        CALL SHOCK(AMU2,RATP,RATR,RATPT)
        PD=PU*RATP
        RHOD=RHOU*RATR
        TD=PD/RHOD
        SD=-.4*ALOG(RATPT)
        UD=SQRT((1.-TD)*7.)
        AD=SQRT(1.4*TD)
        AMD=UD/AD
        TTOTD=TD+1./7.*UD**2
        PTOTD=PD*(1.+.2*AMD**2)**3.5
        XSH=XTRY
        WRITE(*,*)'PROPRIETA DELLO URTO'
        WRITE(*,*)'XSH,PU,AMU,SU,PTOTU'
        WRITE(*,*) XSH,PU,AMU,SU,PTOTU
        WRITE(*,*)'XSH,PD,AMD,SD,PTOTD'
        WRITE(*,*) XSH,PD,AMD,SD,PTOTD
        ENDIF

350     CONTINUE
        RETURN

4444    CONTINUE

        WRITE(*,*)'NA (dispari) '
        READ(*,*)NA
        WRITE(*,*)'exit pressure (PEX)'
        READ(*,*)PEX

        WRITE(*,*)'Il valore di DT risultera costante'
        WRITE(*,*)'per NA=(40-1) est DT=.0125'
        WRITE(*,*)'per diverso NA est DT = 0.0125*40/(NA+1)'
        DTITEST4 = 0.0125*40/(NA+1)

        WRITE(*,*)'KA,KOUT'
        READ(*,*)KA,KOUT
        STAB=1.E6
        B = 0.
        C = 1.
        NC   = NA+2
        NCM  = NC-1
        DX=(C-B)/NA
        DO 410 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.
        AMEDIA(N)=1.
        DIFFAREA(N)=0.0

410     CONTINUE

        DO 420 N=2,NC
        P(N)=1.0
        IF(N.EQ.NC) P(N)=PEX
        T(N)=P(N)**3.5
        A(N)=SQRT(1.4*T(N))
        U(N)=0.0
        AMACH(N)=0.0
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
420     CONTINUE
        RETURN

5555    CONTINUE

        WRITE(*,*)'NA,STAB'
        READ(*,*)NA,STAB
        WRITE(*,*)'exit pressure (PEX)'
        READ(*,*)PEX
        WRITE(*,*)'KA,KOUT'
        READ(*,*)KA,KOUT
        B = 0.
        C = 1.
        NC   = NA+2
        NCM  = NC-1
        DX=(C-B)/NA
        DO 510 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.
        AMEDIA(N)=1.
        DIFFAREA(N)=0.0

510     CONTINUE

        DO 520 N=2,NC
        P(N)=1.0
        IF(N.EQ.NC) P(N)=PEX
        T(N)=P(N)**3.5
        A(N)=SQRT(1.4*T(N))
        U(N)=0.0
        AMACH(N)=0.0
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
520     CONTINUE

        RETURN

6666    CONTINUE

        WRITE(*,*)'NA,STAB (200 .8)'
        READ(*,*)NA,STAB

        WRITE(*,*)'KA,KOUT (100, clean, 150, interact., 300)'
        READ(*,*)KA,KOUT
        B = 0.
        C = 1.
        NC   = NA+2
        NCM  = NC-1
        DX=(C-B)/NA
        DO 61 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.
        AMEDIA(N)=1.
        DIFFAREA(N)=0.0
61      CONTINUE

        XSH10= 0.10
        XSH20= 0.90
        RATP1= 4.5
        RATP2= 2.5

        WRITE(*,*)'location of shocks XSH10  XSH20  (0.10 0.90)'
        READ(*,*)XSH10,XSH20
        WRITE(*,*)'strength of shocks RATP1  RATP2  (4.5  2.5 )'
        READ(*,*)RATP1,RATP2
        P1  =  RATP1
        AM1 =  -SQRT((6.0*RATP1+1.0)/7.0)
        RHO1=  6.0*AM1**2/(5.0+AM1**2)
        VSH1=  -AM1*GF
        U1  =  VSH1*(1-1./RHO1)
        P2  =  RATP2
        AM2 =  SQRT((6.0*RATP2+1.0)/7.0)
        RHO2=  6.0*AM2**2/(5.0+AM2**2)
        VSH2=  -AM2*GF
        U2  =  VSH2*(1-1./RHO2)

        DO 62 N=2,NC
        IF(X(N).LT.XSH10) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        T(N)  = P(N)/RHO(N)
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GE.XSH10.AND.X(N).LE.XSH20) THEN
        P(N)   =1.0
        RHO(N) =1.0
        U(N)   =0.0
        T(N)   =P(N)/RHO(N)
        A(N)=SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GT.XSH20) THEN
        P(N)  = P2
        RHO(N)= RHO2
        U(N)  = U2
        T(N)  = P(N)/RHO(N)
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        AMACH(N)=U(N)/A(N)
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
62      CONTINUE

        RETURN

7777    CONTINUE

        WRITE(*,*)'NA,STAB (200 .8)'
        READ(*,*)NA,STAB

        WRITE(*,*)'KA,KOUT (90, clean, 200, interact., 280)'
        READ(*,*)KA,KOUT
        B = 0.
        C = 1.
        NC   = NA+2
        NCM  = NC-1
        DX=(C-B)/NA
        DO 71 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.
        AMEDIA(N)=1.
        DIFFAREA(N)=0.0
71      CONTINUE

        XSH10= 0.35
        XSH20= 0.05
        RATP1= 2.5
        RATP2= 4.0

        WRITE(*,*)'location of shocks XSH10  XSH20  (0.35 0.05)'
        READ(*,*)XSH10,XSH20
        WRITE(*,*)'strength of shocks RATP1  RATP2  (2.5  4.0 )'
        READ(*,*)RATP1,RATP2
        P1  =  RATP1
        AM1 =  -SQRT((6.0*RATP1+1.0)/7.0)
        RHO1=  6.0*AM1**2/(5.0+AM1**2)
        VSH1=  -GF*AM1
        U1  =  VSH1*(1-1./RHO1)
        P2  =  P1*RATP2
        AM2 =  -SQRT((6.0*RATP2+1.0)/7.0)
        RHO2=  RHO1*6.0*AM2**2/(5.0+AM2**2)
        VSH2=  U1-SQRT(1.4*P1/RHO1)*AM2
        U2  =  VSH2*(1-RHO1/RHO2)+U1*RHO1/RHO2

        DO 72 N=2,NC
        IF(X(N).LT.XSH20) THEN
        P(N)  = P2
        RHO(N)= RHO2
        U(N)  = U2
        T(N)  = P(N)/RHO(N)
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GE.XSH20.AND.X(N).LE.XSH10) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        T(N)  = P(N)/RHO(N)
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GT.XSH10) THEN
        P(N)   =1.0
        RHO(N) =1.0
        U(N)   =0.0
        T(N)   =P(N)/RHO(N)
        A(N)=SQRT(1.4*T(N))
        ENDIF

        AMACH(N)=U(N)/A(N)
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
72      CONTINUE

        RETURN

8888    CONTINUE

        NA  = 200
        STAB= 0.8
        WRITE(*,*)'NA,STAB (200 .8)'
        READ(*,*)NA,STAB

        B = 0.
        C = 1.
        NC   = NA+2
        NCM  = NC-1
        DX=(C-B)/NA
        DO 81 N=1,NC
        X(N)=B+0.5*DX+DX*(N-2)
        AREA(N)=1.
        AMEDIA(N)=1.
        DIFFAREA(N)=0.0
81      CONTINUE

        XSH10= 0.05
        RATP1= 4.50
        XCSF2= 0.30
        RATT2= 2.50

        WRITE(*,*)'shock location XSH10 and strength RATP1 (0.05 4.50)'
        READ(*,*)XSH10,RATP1
        WRITE(*,*)'c.s.  location XCSF2 and strength RATT2 '
        WRITE(*,*)'                        (0.30 2.50/0.40)'
        READ(*,*)XCSF2,RATT2
        WRITE(*,*)'KA,KOUT (50, clean, 85, interact., 280/350)'
        READ(*,*)KA,KOUT


        P1  =  RATP1
        AM1 =  -SQRT((6.0*RATP1+1.0)/7.0)
        RHO1=  6.0*AM1**2/(5.0+AM1**2)
        VSH1=  -AM1*GF
        U1  =  VSH1*(1-1./RHO1)
        P2  =  1.0
        RHO2=  1.0/RATT2
        VCS2=  0.0
        U2  =  0.0

        DO 82 N=2,NC
        IF(X(N).LT.XSH10) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        T(N)  = P(N)/RHO(N)
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GE.XSH10.AND.X(N).LE.XCSF2) THEN
        P(N)   =1.0
        RHO(N) =1.0
        U(N)   =0.0
        T(N)   =P(N)/RHO(N)
        A(N)=SQRT(1.4*T(N))
        ENDIF

        IF(X(N).GT.XCSF2) THEN
        P(N)  = P2
        T(N)  = RATT2
        RHO(N)= RHO2
        U(N)  = U2
        A(N)  = SQRT(1.4*T(N))
        ENDIF

        AMACH(N)=U(N)/A(N)
        RHO(N)=P(N)/T(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        PXENO(N)=0.0
        UXENO(N)=0.0
        HXENO(N)=0.0
        PPXENO(N)=0.0
        HHXENO(N)=0.0
82      CONTINUE

        RETURN

        END

C....................................................................

        SUBROUTINE MARCH

        INCLUDE 'COMST05.INC'
        COMMON /CONVERG/ RTRMS,RTMAX,NRT

        DTODX=DT/DX

        RTRMS=0.
        RTMAX=0.
        NRT  =0.
        NIN  = 2
        DO 1 N=NIN,NCM
        ENU=DTODX
        W1OLD=W1(N)
        W1(N)=W1(N)
     1       -ENU/AREA(N)*(PHI1(N)*AMEDIA(N)-PHI1(N-1)*AMEDIA(N-1))
        W2(N)=W2(N)
     1       -ENU/AREA(N)*(PHI2(N)*AMEDIA(N)-PHI2(N-1)*AMEDIA(N-1))
     2       -ENU/AREA(N)*(-P(N) * DIFFAREA(N))
        W3(N)=W3(N)
     1       -ENU/AREA(N)*(PHI3(N)*AMEDIA(N)-PHI3(N-1)*AMEDIA(N-1))

        RT    = (W1(N)-W1OLD)/RHO(2)
        RTRMS = RTRMS + RT**2
        IF (ABS(RT).GE.ABS(RTMAX)) THEN
        RTMAX = ABS(RT)
        NRT = N
        ENDIF
1       CONTINUE

        RTRMS = SQRT(RTRMS / FLOAT(NCM-4))

        DO 2 N=NIN,NCM
        RHO(N)   =W1(N)
        U(N)     =W2(N)/W1(N)
        E(N)     =W3(N)
        P(N)     =(E(N)-.5*RHO(N)*U(N)*U(N))/GB

        IF (P(N).LT.0.) THEN
        WRITE(*,*) 'PRESSIONE NEGATIVA', N,P(N),P(N-1)
        ENDIF

        T(N)     =P(N)/RHO(N)
        A(N)     =GF*SQRT(T(N))
        H(N)     =GA*T(N)
        AMACH(N) =U(N)/A(N)
        S(N)     =ALOG(P(N))-GAMMA*ALOG(RHO(N))
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        HTOT(N)  =H(N)+.5*U(N)**2
        F1(N)    =RHO(N)*U(N)
        F2(N)    =P(N)+RHO(N)*U(N)*U(N)
        F3(N)    =U(N)*(P(N)+E(N))
        FLOW(N)  =RHO(N)*U(N)*AREA(N)
        FLHT(N)  =U(N)*(P(N)+E(N))*AREA(N)
2       CONTINUE
        RETURN
        END

C....................................................................

        SUBROUTINE ENO1

        INCLUDE 'COMST05.INC'
        COMMON / ENO / PXENO(NMAX),UXENO(NMAX),HXENO(NMAX),
     1           PPXENO(NMAX),HHXENO(NMAX),
     2           PPT(NMAX),UT(NMAX),HHT(NMAX),
     3           PPTT(NMAX),UTT(NMAX),HHTT(NMAX)

        NCMM =NC-2
        NCMMM=NC-3

        DO 1 N=3,NCMM

        NM1=N-1
        N00=N
        NP1=N+1

C       PPXP=ALOG(P(NP1))-ALOG(P(N00))
C       PPXM=ALOG(P(N00))-ALOG(P(NM1))
C       CALL MINMOD(PPXP,PPXM,PPXDUM)

C       UXP=U(NP1)-U(N00)
C       UXM=U(N00)-U(NM1)
C       CALL MINMOD(UXP,UXM,UXDUM)

C       HHXP=ALOG(H(NP1))-ALOG(H(N00))
C       HHXM=ALOG(H(N00))-ALOG(H(NM1))
C       CALL MINMOD(HHXP,HHXM,HHXDUM)

        GOA    = GAMMA/A(N00)
        R1XP   = (ALOG(P(NP1))-GOA*U(NP1))-(ALOG(P(N00))-GOA*U(N00))
        R1XM   = (ALOG(P(N00))-GOA*U(N00))-(ALOG(P(NM1))-GOA*U(NM1))
        CALL MINMOD(R1XP,R1XM,R1XDUM)

        R2XP   = (ALOG(H(NP1))-GJ*ALOG(P(NP1)))-
     1           (ALOG(H(N00))-GJ*ALOG(P(N00)))
        R2XM   = (ALOG(H(N00))-GJ*ALOG(P(N00)))-
     1           (ALOG(H(NM1))-GJ*ALOG(P(NM1)))
        CALL MINMOD(R2XP,R2XM,R2XDUM)

        R3XP   = (ALOG(P(NP1))+GOA*U(NP1))-(ALOG(P(N00))+GOA*U(N00))
        R3XM   = (ALOG(P(N00))+GOA*U(N00))-(ALOG(P(NM1))+GOA*U(NM1))
        CALL MINMOD(R3XP,R3XM,R3XDUM)

        PPXDUM  = .5*(R3XDUM+R1XDUM)
        UXDUM   = .5*(R3XDUM-R1XDUM)/GOA
        HHXDUM  = R2XDUM+GJ*PPXDUM

        PPXENO(N) =  PPXDUM
        UXENO(N)  =   UXDUM
        HHXENO(N) =  HHXDUM

1       CONTINUE

        N = 2
        PPXENO(N)= ALOG(P(3))-ALOG(P(2))
        UXENO(N) = U(3)-U(2)
        HHXENO(N)= ALOG(H(3))-ALOG(H(2))
        N = NCM
        PPXENO(N)= ALOG(P(NCM))-ALOG(P(NCMM))
        UXENO(N) = U(NCM)-U(NCMM)
        HHXENO(N)= ALOG(H(NCM))-ALOG(H(NCMM))

        DO 2 N=2,NCM
        PPT(N) = -(U(N)*PPXENO(N)+GAMMA*UXENO(N)+
     1             GAMMA*U(N)*(DIFFAREA(N)/AREA(N)))
        UT(N)  = -(U(N)*UXENO(N)+T(N)*PPXENO(N))
        HHT(N) =   PPT(N)/GA - U(N)*(HHXENO(N)-PPXENO(N)/GA)
2       CONTINUE

        DO 3 N=2,NCM
        PPTT(N) = 0.0
        UTT(N)  = 0.0
        HHTT(N) = 0.0
3       CONTINUE

        RETURN
        END

C..................................................................
        SUBROUTINE ENO2

        INCLUDE 'COMST05.INC'
        COMMON / ENO / PXENO(NMAX),UXENO(NMAX),HXENO(NMAX),
     1           PPXENO(NMAX),HHXENO(NMAX),
     2           PPT(NMAX),UT(NMAX),HHT(NMAX),
     3           PPTT(NMAX),UTT(NMAX),HHTT(NMAX)

        NCMM =NC-2
        NCMMM=NC-3
        DO 1 N=4,NCMMM
        NM2=N-2
        NM1=N-1
        N00=N
        NP1=N+1
        NP2=N+2

C       PPNM2 = ALOG(P(NM2))
C       PPNM1 = ALOG(P(NM1))
C       PPN00 = ALOG(P(N00))
C       PPNP1 = ALOG(P(NP1))
C       PPNP2 = ALOG(P(NP2))
C       CALL DECONV(PPNM2,PPNM1,PPN00,PPNP1,PPNP2,PPXDUM)

C       UNM2 = U(NM2)
C       UNM1 = U(NM1)
C       UN00 = U(N00)
C       UNP1 = U(NP1)
C       UNP2 = U(NP2)
C       CALL DECONV(UNM2,UNM1,UN00,UNP1,UNP2,UXDUM)

C       HHNM2 = ALOG(H(NM2))
C       HHNM1 = ALOG(H(NM1))
C       HHN00 = ALOG(H(N00))
C       HHNP1 = ALOG(H(NP1))
C       HHNP2 = ALOG(H(NP2))
C       CALL DECONV(HHNM2,HHNM1,HHN00,HHNP1,HHNP2,HHXDUM)

        GOA    = GAMMA/A(N00)
        R1NM2 = ALOG(P(NM2))-GOA*U(NM2)
        R1NM1 = ALOG(P(NM1))-GOA*U(NM1)
        R1N00 = ALOG(P(N00))-GOA*U(N00)
        R1NP1 = ALOG(P(NP1))-GOA*U(NP1)
        R1NP2 = ALOG(P(NP2))-GOA*U(NP2)
        CALL DECONV(R1NM2,R1NM1,R1N00,R1NP1,R1NP2,R1XDUM)

        R2NM2 = ALOG(H(NM2))-GJ*ALOG(P(NM2))
        R2NM1 = ALOG(H(NM1))-GJ*ALOG(P(NM1))
        R2N00 = ALOG(H(N00))-GJ*ALOG(P(N00))
        R2NP1 = ALOG(H(NP1))-GJ*ALOG(P(NP1))
        R2NP2 = ALOG(H(NP2))-GJ*ALOG(P(NP2))
        CALL DECONV(R2NM2,R2NM1,R2N00,R2NP1,R2NP2,R2XDUM)

        R3NM2 = ALOG(P(NM2))+GOA*U(NM2)
        R3NM1 = ALOG(P(NM1))+GOA*U(NM1)
        R3N00 = ALOG(P(N00))+GOA*U(N00)
        R3NP1 = ALOG(P(NP1))+GOA*U(NP1)
        R3NP2 = ALOG(P(NP2))+GOA*U(NP2)
        CALL DECONV(R3NM2,R3NM1,R3N00,R3NP1,R3NP2,R3XDUM)

        PPXDUM  = .5*(R3XDUM+R1XDUM)
        UXDUM   = .5*(R3XDUM-R1XDUM)/GOA
        HHXDUM  = R2XDUM+GJ*PPXDUM

        PPXENO(N) = PPXDUM
        UXENO(N)  =  UXDUM
        HHXENO(N) = HHXDUM

1       CONTINUE

        N=2
        PPXENO(N) = ALOG(P(3))-ALOG(P(2))
        UXENO(N)  = U(3)-U(2)
        HHXENO(N) = ALOG(H(3))-ALOG(H(2))

        N = 3
        PPXD1 = ALOG(P(3)) -ALOG(P(2))
        PPXD2 = ALOG(P(4)) -ALOG(P(3))
        CALL MINMOD(PPXD1,PPXD2,PPXDUM)

        UXD1 = U(3) -U(2)
        UXD2 = U(4) -U(3)
        CALL MINMOD(UXD1,UXD2,UXDUM)

        HHXD1 = ALOG(H(3)) -ALOG(H(2))
        HHXD2 = ALOG(H(4)) -ALOG(H(3))
        CALL MINMOD(HHXD1,HHXD2,HHXDUM)

        PPXENO(N) = PPXDUM
        UXENO(N)  =  UXDUM
        HHXENO(N) = HHXDUM


        N = NCMM
        PPXD1 = ALOG(P(NCMM)) -ALOG(P(NCMMM))
        PPXD2 = ALOG(P(NCM )) -ALOG(P(NCMM ))
        CALL MINMOD(PPXD1,PPXD2,PPXDUM)

        UXD1 = U(NCMM) -U(NCMMM)
        UXD2 = U(NCM ) -U(NCMM )
        CALL MINMOD(UXD1,UXD2,UXDUM)

        HHXD1 = ALOG(H(NCMM)) -ALOG(H(NCMMM))
        HHXD2 = ALOG(H(NCM )) -ALOG(H(NCMM ))
        CALL MINMOD(HHXD1,HHXD2,HHXDUM)

        PPXENO(N) = PPXDUM
        UXENO(N)  =  UXDUM
        HHXENO(N) = HHXDUM

        N = NCM
        PPXENO(N)= ALOG(P(NCM))-ALOG(P(NCMM))
        UXENO(N) = U(NCM)-U(NCMM)
        HHXENO(N)= ALOG(H(NCM))-ALOG(H(NCMM))

        DO 2 N=2,NCM
        PPT(N) = -(U(N)*PPXENO(N)+GAMMA*UXENO(N)+
     1             GAMMA*U(N)*(DIFFAREA(N)/AREA(N)))
        UT(N)  = -(U(N)*UXENO(N)+T(N)*PPXENO(N))
        HHT(N) =   PPT(N)/GA - U(N)*(HHXENO(N)-PPXENO(N)/GA)
2       CONTINUE

        DO 3 N=2,NCM
        IF(N.GE.3.AND.N.LE.NCM-1) THEN
        PPTXP = (PPT(N+1)-PPT(N))/DX
        UTXP  = (UT (N+1)-UT (N))/DX
        HHTXP = (HHT(N+1)-HHT(N))/DX
        PPTXM = (PPT(N)-PPT(N-1))/DX
        UTXM  = (UT (N)-UT (N-1))/DX
        HHTXM = (HHT(N)-HHT(N-1))/DX
        CALL MINMOD(PPTXP,PPTXM,PPTX)
        CALL MINMOD(UTXP ,UTXM ,UTX )
        CALL MINMOD(PPTXP,PPTXM,PPTX)
        ENDIF
        IF(N.EQ.2) THEN
        PPTX = (PPT(N+1)-PPT(N))/DX
        UTX  = (UT (N+1)-UT (N))/DX
        HHTX = (HHT(N+1)-HHT(N))/DX
        ENDIF
        IF(N.EQ.NCM) THEN
        PPTX = (PPT(N)-PPT(N-1))/DX
        UTX  = (UT (N)-UT (N-1))/DX
        HHTX = (HHT(N)-HHT(N-1))/DX
        ENDIF

        TT = HHT(N)*H(N)/GA
        PPTT(N) = -(U(N)*PPTX+UT(N)*PPXENO(N)+GAMMA*UTX+
     1              GAMMA*UT(N)*(DIFFAREA(N)/AREA(N)))
        UTT(N)  = -(U(N)*UTX+UT(N)*UXENO(N)+TT*PPXENO(N)+T(N)*PPTX)
        HHTT(N) = -(U(N)*HHTX+UT(N)*HHXENO(N)-
     1             (PPTT(N)+U(N)*PPTX+UT(N)*PPXENO(N))/GA)

3       CONTINUE

        RETURN
        END

C..................................................................

        SUBROUTINE DECONV(RM2,RM1,R00,RP1,RP2,RX)

        CURV00 = ABS(RP1-2.*R00+RM1)
        CURVP1 = ABS(RP2-2.*RP1+R00)
        IF(CURV00.LT.CURVP1) THEN
        RXPLUS = (RP1-RM1)/2.
        ELSE
        RXPLUS =-(RP2-4.*RP1+3.*R00)/2.
        ENDIF

        CURV00 = ABS(RP1-2.*R00+RM1)
        CURVM1 = ABS(R00-2.*RM1+RM2)
        IF(CURV00.LT.CURVM1) THEN
        RXMINUS = (RP1-RM1)/2.
        ELSE
        RXMINUS = (3.*R00-4.*RM1+RM2)/2.
        ENDIF

        CALL MINMOD(RXPLUS,RXMINUS,RX)

        RETURN
        END

C....................................................................

        SUBROUTINE MINMOD(A,B,C)

        C = 0.0
        IF(A*B.GE..0) C=SIGN(1.,A)*AMIN1(ABS(A),ABS(B))
        RETURN
        END
C..................................................................


        SUBROUTINE SPLIT

        INCLUDE 'COMST05.INC'
        COMMON / ENO / PXENO(NMAX),UXENO(NMAX),HXENO(NMAX),
     1           PPXENO(NMAX),HHXENO(NMAX),
     2           PPT(NMAX),UT(NMAX),HHT(NMAX),
     3           PPTT(NMAX),UTT(NMAX),HHTT(NMAX)

        ENUO2 = 0.5*DT/DX

        NCMM = NCM-1
        DO 1 N=2,NCMM
        NM=N
        NP=NM+1

        PA  =  P(NM)
        PB  =  P(NP)
        UA  =  U(NM)
        UB  =  U(NP)
        HA  =  H(NM)
        HB  =  H(NP)

        IF(N.EQ.2) THEN
        P002   =  PA
        U002   =  UA
        H002   =  HA
        RHO002 =  P002/H002*GA
        A002   =  SQRT(GAMMA*P002/RHO002)
        ENDIF

        IF(N.EQ.NCMM) THEN
        PNCM   =  PB
        UNCM   =  UB
        HNCM   =  HB
        RHONCM =  PNCM/HNCM*GA
        ANCM   =  SQRT(GAMMA*PNCM/RHONCM)
        ENDIF

        IF(IORD.NE.1) THEN
        PPA = ALOG(PA)
        PPB = ALOG(PB)
        HHA = ALOG(HA)
        HHB = ALOG(HB)
        PPA  =  PPA + .5* PPXENO(N)
        PPB  =  PPB - .5* PPXENO(N+1)
        UA   =  UA  + .5* UXENO(N)
        UB   =  UB  - .5* UXENO(N+1)
        HHA  =  HHA + .5* HHXENO(N)
        HHB  =  HHB - .5* HHXENO(N+1)

        PPA = PPA + ENUO2*(PPT(NM)+PPTT(NM)*DT/3.)
        UA  = UA  + ENUO2*(UT(NM) +UTT(NM) *DT/3.)
        HHA = HHA + ENUO2*(HHT(NM)+HHTT(NM)*DT/3.)

        IF(N.EQ.2) THEN
        PP002  =  ALOG(P(NM))
        U002   =  U(NM)
        HH002  =  ALOG(H(NM))
        PP002  =  PP002 - .5* PPXENO(N)
        U002   =  U002  - .5* UXENO(N)
        HH002  =  HH002 - .5* HHXENO(N)
        PP002  =  PP002 + ENUO2*(PPT(N)+PPTT(N)*DT/3.)
        U002   =  U002  + ENUO2*(UT(N) +UTT (N)*DT/3.)
        HH002  =  HH002 + ENUO2*(HHT(N)+HHTT(N)*DT/3.)
        P002   =  EXP(PP002)
        H002   =  EXP(HH002)
        RHO002 =  P002/H002*GA
        A002   =  SQRT(GAMMA*P002/RHO002)
        ENDIF

        PPB = PPB + ENUO2*(PPT(NP)+PPTT(NP)*DT/3.)
        UB  = UB  + ENUO2*(UT(NP) +UTT(NP) *DT/3.)
        HHB = HHB + ENUO2*(HHT(NP)+HHTT(NP)*DT/3.)

        IF(N.EQ.NCMM) THEN
        PPNCM  =  ALOG(P(NP))
        UNCM   =  U(NP)
        HHNCM  =  ALOG(H(NP))
        PPNCM  =  PPNCM    + .5* PPXENO(N+1)
        UNCM   =  UNCM     + .5* UXENO(N+1)
        HHNCM  =  HHNCM    + .5* HHXENO(N+1)
        PPNCM  =  PPNCM    + ENUO2*(PPT(N+1)+PPTT(N+1)*DT/3.)
        UNCM   =  UNCM     + ENUO2*(UT(N+1) +UTT (N+1)*DT/3.)
        HHNCM  =  HHNCM    + ENUO2*(HHT(N+1)+HHTT(N+1)*DT/3.)
        PNCM   =  EXP(PPNCM)
        HNCM   =  EXP(HHNCM)
        RHONCM =  PNCM/HNCM*GA
        ANCM   =  SQRT(GAMMA*PNCM/RHONCM)
        ENDIF

        PA = EXP(PPA)
        PB = EXP(PPB)
        HA = EXP(HHA)
        HB = EXP(HHB)
        ENDIF

        RHOA=  PA/HA*GA
        RHOB=  PB/HB*GA
        AA  =  SQRT(GAMMA*PA/RHOA)
        AB  =  SQRT(GAMMA*PB/RHOB)

        ICALC=0

        R3A=PA+(RHOA*AA)*UA
        R2A=HA-PA/RHOA
        R2B=HB-PB/RHOB
        R1B=PB-(RHOB*AB)*UB
        UC = (R3A-R1B)/((RHOA*AA)+(RHOB*AB))
        UD = UC
        PC = R3A-(RHOA*AA)*UC
        PD = PC
        HC = PC/RHOA+R2A
        HD = PD/RHOB+R2B

        IF((PC.LT.0.0).OR.(PD.LT.0.0).
     1                 OR.(HC.LT.0.0).OR.(HD.LT.0.0)) ICALC=1

        IF(ICALC.EQ.1) THEN
        PPA = ALOG(PA)
        PPB = ALOG(PB)
        HHA = ALOG(HA)
        HHB = ALOG(HB)
        R3A=PPA+GAMMA/AA*UA
        R2A=HHA-PPA/GA
        R2B=HHB-PPB/GA
        R1B=PPB-GAMMA/AB*UB
        UC = (R3A-R1B)/(GAMMA/AA+GAMMA/AB)
        UD = UC
        PPC= R3A-GAMMA/AA*UC
        PPD= PPC
        HHC= PPC/GA+R2A
        HHD= PPD/GA+R2B
        PC = EXP(PPC)
        PD = EXP(PPD)
        HC = EXP(HHC)
        HD = EXP(HHD)
        WRITE(*,*)'WARNING ICALC=1 AT K,N = ',K,N
        ENDIF

        AC = SQRT(2.*GD*HC)
        AD = SQRT(2.*GD*HD)
        ALA=UA-AA
        ALC=UC-AC
        ALD=UD+AD
        ALB=UB+AB
        ALX=UC
        CALL DECOD(PA,UA,HA,F1A,F2A,F3A)
        CALL DECOD(PB,UB,HB,F1B,F2B,F3B)
        CALL DECOD(PC,UC,HC,F1C,F2C,F3C)
        CALL DECOD(PD,UD,HD,F1D,F2D,F3D)
C       ...................................................
        IF((ALA.LE..0.AND.ALC.LE..0).AND.
     1               (ALD.GT..0.AND.ALB.GT..0)) GO TO 10
        IF((ALA.GT..0.AND.ALC.GT..0).AND.
     1               (ALD.GT..0.AND.ALB.GT..0)) GO TO 11
        IF((ALA.LE..0.AND.ALC.LE..0).AND.
     1               (ALD.LE..0.AND.ALB.LE..0)) GO TO 12
        IF((ALA.GT..0.AND.ALC.LE..0).AND.
     1               (ALD.GT..0.AND.ALB.GT..0)) GO TO 13
        IF((ALA.LE..0.AND.ALC.GT..0).AND.
     1               (ALD.GT..0.AND.ALB.GT..0)) GO TO 14
        IF((ALD.GT..0.AND.ALB.LE..0).AND.
     1               (ALC.LT..0.AND.ALA.LT..0)) GO TO 15
        IF((ALD.LE..0.AND.ALB.GT..0).AND.
     1               (ALC.LT..0.AND.ALA.LT..0)) GO TO 16
        IF((ALA.GT..0.AND.ALC.LE..0).AND.
     1               (ALD.GT..0.AND.ALB.LT..0)) GO TO 17

        WRITE(*,*)'attenzione scelta delle caratteristiche ciucca'
        WRITE(*,*)'ala,alc,alx,ald,alb',ALA,ALC,ALX,ALD,ALB
        STOP
10      CONTINUE
        IF(ALX.GT..0) GO TO 20
        DF1R=F1B-F1D
        DF2R=F2B-F2D
        DF3R=F3B-F3D
        DF1L=F1D-F1A
        DF2L=F2D-F2A
        DF3L=F3D-F3A
        GO TO 2
20      CONTINUE
        DF1R=F1B-F1C
        DF2R=F2B-F2C
        DF3R=F3B-F3C
        DF1L=F1C-F1A
        DF2L=F2C-F2A
        DF3L=F3C-F3A
        GO TO 2
11      CONTINUE
        DF1R=F1B-F1A
        DF2R=F2B-F2A
        DF3R=F3B-F3A
        DF1L=.0
        DF2L=.0
        DF3L=.0
        GO TO 2
12      CONTINUE
        DF1R=.0
        DF2R=.0
        DF3R=.0
        DF1L=F1B-F1A
        DF2L=F2B-F2A
        DF3L=F3B-F3A
        GO TO 2
13      CONTINUE
        GDGD=2.*GD
        IF(ICALC.EQ.0) THEN
        SQRHST=.5*(-SQRT(GDGD)*AA+SQRT(GDGD*AA*AA+4.*(R3A/RHOA+R2A)))
        HST = SQRHST**2
        AST = SQRT(GDGD)*SQRHST
        UST = AST
        PST = RHOA*(-R2A+HST)
        ENDIF
        IF(ICALC.EQ.1) THEN
        ADUM = 2.*AA/GDGD**1.5
        BDUM = AA/GAMMA/SQRT(GDGD)*(R3A+GA*R2A)
        CALL ITERSTAR (ADUM,BDUM,CDUM)
        AST = SQRT(GDGD)*CDUM
        UST = AST
        HST = GB*AST**2
        PST= EXP((ALOG(HST)-R2A)*GA)
        ENDIF
        CALL DECOD(PST,UST,HST,F1ST,F2ST,F3ST)
        IF(ALX.LE.0) GO TO 21
        DF1R=F1B-F1C+F1ST-F1A
        DF2R=F2B-F2C+F2ST-F2A
        DF3R=F3B-F3C+F3ST-F3A
        DF1L=F1C-F1ST
        DF2L=F2C-F2ST
        DF3L=F3C-F3ST
        GO TO 2
21      CONTINUE
        DF1R=F1B-F1D+F1ST-F1A
        DF2R=F2B-F2D+F2ST-F2A
        DF3R=F3B-F3D+F3ST-F3A
        DF1L=F1D-F1ST
        DF2L=F2D-F2ST
        DF3L=F3D-F3ST
        GO TO 2
14      CONTINUE
        GDGD=2.*GD
        SQRHST=.5*(-SQRT(GDGD)*AA+SQRT(GDGD*AA*AA+4.*(R3A/RHOA+R2A)))
        HST = SQRHST**2
        AST = SQRT(GDGD)*SQRHST
        UST = AST
        PST = RHOA*(-R2A+HST)
        CALL DECOD(PST,UST,HST,F1ST,F2ST,F3ST)
        DF1R=F1B-F1ST
        DF2R=F2B-F2ST
        DF3R=F3B-F3ST
        DF1L=F1ST-F1A
        DF2L=F2ST-F2A
        DF3L=F3ST-F3A
        GO TO 2
15      CONTINUE
        GDGD=2.*GD
        SQRHST=.5*(-SQRT(GDGD)*AB+SQRT(GDGD*AB*AB+4.*(R1B/RHOB+R2B)))
        HST = SQRHST**2
        AST = SQRT(GDGD)*SQRHST
        UST = -AST
        PST = RHOB*(-R2B+HST)
        CALL DECOD(PST,UST,HST,F1ST,F2ST,F3ST)
        IF(ALX.LE.0) GO TO 22
        DF1R=F1ST-F1C
        DF2R=F2ST-F2C
        DF3R=F3ST-F3C
        DF1L=F1B-F1ST+F1C-F1A
        DF2L=F2B-F2ST+F2C-F2A
        DF3L=F3B-F3ST+F3C-F3A
        GO TO 2
22      CONTINUE
        DF1R=F1ST-F1D
        DF2R=F2ST-F2D
        DF3R=F3ST-F3D
        DF1L=F1B-F1ST+F1D-F1A
        DF2L=F2B-F2ST+F2D-F2A
        DF3L=F3B-F3ST+F3D-F3A
        GO TO 2
16      CONTINUE
        GDGD=2.*GD
        SQRHST=.5*(-SQRT(GDGD)*AB+SQRT(GDGD*AB*AB+4.*(R1B/RHOB+R2B)))
        HST = SQRHST**2
        AST = SQRT(GDGD)*SQRHST
        UST = -AST
        PST = RHOB*(-R2B+HST)
        CALL DECOD(PST,UST,HST,F1ST,F2ST,F3ST)
        DF1R=F1B-F1ST
        DF2R=F2B-F2ST
        DF3R=F3B-F3ST
        DF1L=F1ST-F1A
        DF2L=F2ST-F2A
        DF3L=F3ST-F3A
        GO TO 2
17      CONTINUE
        GDGD=2.*GD
        SQRHST=.5*(-SQRT(GDGD)*AA+SQRT(GDGD*AA*AA+4.*(R3A/RHOA+R2A)))
        HST1 = SQRHST**2
        AST1 = SQRT(GDGD)*SQRHST
        UST1 = AST1
        PST1 = RHOA*(-R2A+HST1)
        CALL DECOD(PST1,UST1,HST1,F1ST1,F2ST1,F3ST1)
        SQRHST=.5*(-SQRT(GDGD)*AB+SQRT(GDGD*AB*AB+4.*(R1B/RHOB+R2B)))
        HST2 = SQRHST**2
        AST2 = SQRT(GDGD)*SQRHST
        UST2 = -AST2
        PST2 = RHOB*(-R2B+HST2)
        CALL DECOD(PST2,UST2,HST2,F1ST2,F2ST2,F3ST2)

        IF(ALX.LE.0) GO TO 171
        DF1L=F1C-F1ST1+F1B-F1ST2
        DF2L=F2C-F2ST1+F2B-F2ST2
        DF3L=F3C-F3ST1+F3B-F3ST2
        GO TO 2
171     CONTINUE
        DF1L=F1D-F1ST+F1B-F1ST2
        DF2L=F2D-F2ST+F2B-F2ST2
        DF3L=F3D-F3ST+F3B-F3ST2
        GO TO 2

2       CONTINUE
C       PHI1(N)=F1B-DF1R
C       PHI2(N)=F2B-DF2R
C       PHI3(N)=F3B-DF3R
        PHI1(N)=F1A+DF1L
        PHI2(N)=F2A+DF2L
        PHI3(N)=F3A+DF3L
1       CONTINUE

        IF(ITEST.EQ.1) THEN
        R1DUM  = P002-RHO002*A002*U002
        R2DUM  = H002-P002/RHO002
        PIN    = R1DUM
        UIN    = 0.0
        HIN    = R2DUM+PIN/RHO002
        CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))
        R3DUM  = PNCM+RHONCM*ANCM*UNCM
        R2DUM  = HNCM-PNCM/RHONCM
        PEX    = R3DUM
        UEX    = 0.0
        HEX    = R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF

        IF(ITEST.EQ.2) THEN
        TIN    = 1./(1.+.2*ACHIN**2)
        PIN    = TIN**3.5
        AIN    = SQRT(1.4*TIN)
        UIN    = ACHIN*AIN
        HIN    = 3.5*TIN
        CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))

        AMA2DUM=(UNCM/ANCM)**2
        IF(AMA2DUM.LT.1.) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        IF(AMA2DUM.GE.1.) THEN
        CALL SHOCK(AMA2DUM,RATP,RATR,RATPT)
        PRH    = PNCM*RATP
        IF(PRH.LE.PEX) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        IF(UEX.LT.0) UEX=.01
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ELSE
        CALL DECOD(P(NCM),U(NCM),H(NCM),PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        ENDIF

        ENDIF

        IF(ITEST.EQ.3) THEN
        R1DUM  = GG*A002-U002
        AIN    = (R1DUM+SQRT(GAMMA*GC-R1DUM**2*GD))/GC
        UIN    = GG*AIN-R1DUM
        SIN    = 0.0
        TIN    = 1./(1.+.2*(UIN/AIN)**2)
        PIN    = TIN**GA
        HIN    = GA*TIN
        CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))

        AMA2DUM=(UNCM/ANCM)**2
        IF(AMA2DUM.LT.1.) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        IF(AMA2DUM.GE.1.) THEN
        CALL SHOCK(AMA2DUM,RATP,RATR,RATPT)
        PRH    = PNCM*RATP
        IF(PRH.LE.PEX) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        IF(UEX.LT.0) UEX=.01
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ELSE
        CALL DECOD(P(NCM),U(NCM),H(NCM),PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        ENDIF

        ENDIF

        IF(ITEST.EQ.4) THEN
        R1DUM  = GG*A002-U002
        AIN    = (R1DUM+SQRT(GAMMA*GC-R1DUM**2*GD))/GC
        UIN    = GG*AIN-R1DUM
        SIN    = 0.0
        TIN    = 1./(1.+.2*(UIN/AIN)**2)
        PIN    = TIN**GA
        HIN    = GA*TIN
        CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF

        IF(ITEST.EQ.5) THEN

C....... simple wave boundary condition at the inlet.......
        R3IN   = GG*GF
        R1DUM  = GG*A002-U002
        AIN    = 0.5*GD*(R3IN+R1DUM)
        UIN    = 0.5   *(R3IN-R1DUM)
        SIN    = 0.0
        TIN    = AIN*AIN/GAMMA
        PIN    = TIN**GA
        HIN    = GA*TIN
        CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))

C....... prescribed total temperature and entropy at the inlet..
C       R1DUM  = GG*A002-U002
C       AIN    = (R1DUM+SQRT(GAMMA*GC-R1DUM**2*GD))/GC
C       UIN    = GG*AIN-R1DUM
C       SIN    = 0.0
C       TIN    = 1./(1.+.2*(UIN/AIN)**2)
C       PIN    = TIN**GA
C       HIN    = GA*TIN
C       CALL DECOD(PIN,UIN,HIN,PHI1(1),PHI2(1),PHI3(1))
C......................................................
        AMA2DUM=(UNCM/ANCM)**2
        IF(AMA2DUM.LT.1.) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        GDGD=2.*GD
        SQRHST=.5*(-SQRT(GDGD)*ANCM+SQRT(GDGD*ANCM**2+
     1         4.*(R3DUM/RHONCM+R2DUM)))
        HSTEX = SQRHST**2
        ASTEX = SQRT(GDGD)*SQRHST
        USTEX = ASTEX
        PSTEX = RHONCM*(-R2DUM+HSTEX)
        IF(PEX.GE.PSTEX) THEN
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ELSE
        CALL DECOD(PSTEX,USTEX,HSTEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        ENDIF
        IF(AMA2DUM.GE.1.) THEN
        CALL SHOCK(AMA2DUM,RATP,RATR,RATPT)
        PRH    = PNCM*RATP
        IF(PRH.LE.PEX) THEN
        R3DUM = PNCM+RHONCM*ANCM*UNCM
        R2DUM = HNCM-PNCM/RHONCM
        UEX=(R3DUM-PEX)/(RHONCM*ANCM)
        IF(UEX.LT.0) UEX=.01
        HEX=R2DUM+PEX/RHONCM
        CALL DECOD(PEX,UEX,HEX,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ELSE
        CALL DECOD(P(NCM),U(NCM),H(NCM),PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF
        ENDIF

        ENDIF

        IF(ITEST.EQ.6) THEN
        AIN    = SQRT(GAMMA*P1/RHO1)
        UIN    = U1
        ATIN1S = AIN**2+GD*UIN**2
        R1DUM  = GG*A002-U002
        AIN1   = (R1DUM+SQRT(ATIN1S*GC-R1DUM**2*GD))/GC
        UIN1   = GG*AIN1-R1DUM
        PIN1   = P1*(AIN1/AIN)**GI
        TIN1   = AIN1**2/GAMMA
        HIN1   = GA*TIN1
        CALL DECOD(PIN1,UIN1,HIN1,PHI1(1),PHI2(1),PHI3(1))
        AIN    = SQRT(GAMMA*P2/RHO2)
        UIN    = U2
        ATIN2S = AIN**2+GD*UIN**2
        R3DUM  = GG*ANCM+UNCM
        AIN2   = (R3DUM+SQRT(ATIN2S*GC-R3DUM**2*GD))/GC
        UIN2   = -(GG*AIN2-R3DUM)
        PIN2   = P2*(AIN2/AIN)**GI
        TIN2   = AIN2**2/GAMMA
        HIN2   = GA*TIN2
        CALL DECOD(PIN2,UIN2,HIN2,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF

        IF(ITEST.EQ.7) THEN
        AIN    = SQRT(GAMMA*P2/RHO2)
        UIN    = U2
        ATIN1S = AIN**2+GD*UIN**2
        R1DUM  = GG*A002-U002
        AIN1   = (R1DUM+SQRT(ATIN1S*GC-R1DUM**2*GD))/GC
        UIN1   = GG*AIN1-R1DUM
        PIN1   = P2*(AIN1/AIN)**GI
        TIN1   = AIN1**2/GAMMA
        HIN1   = GA*TIN1
        CALL DECOD(PIN1,UIN1,HIN1,PHI1(1),PHI2(1),PHI3(1))
        R1DUM  = GG*GF
        R3DUM  = GG*ANCM+UNCM
        AIN2   = 0.5*GD*(R3DUM+R1DUM)
        UIN2   = 0.5*   (R3DUM-R1DUM)
        PIN2   = (AIN2/GF)**GI
        TIN2   = AIN2**2/GAMMA
        HIN2   = GA*TIN2
        CALL DECOD(PIN2,UIN2,HIN2,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF

        IF(ITEST.EQ.8) THEN
        AIN    = SQRT(GAMMA*P1/RHO1)
        UIN    = U1
        ATIN1S = AIN**2+GD*UIN**2
        R1DUM  = GG*A002-U002
        AIN1   = (R1DUM+SQRT(ATIN1S*GC-R1DUM**2*GD))/GC
        UIN1   = GG*AIN1-R1DUM
        PIN1   = P1*(AIN1/AIN)**GI
        TIN1   = AIN1**2/GAMMA
        HIN1   = GA*TIN1
        CALL DECOD(PIN1,UIN1,HIN1,PHI1(1),PHI2(1),PHI3(1))
        R1DUM  = GG*A(NC)
        R3DUM  = GG*ANCM+UNCM
        AIN2   = 0.5*GD*(R3DUM+R1DUM)
        UIN2   = 0.5*   (R3DUM-R1DUM)
        PIN2   = (AIN2/A(NC))**GI
        TIN2   = AIN2**2/GAMMA
        HIN2   = GA*TIN2
        CALL DECOD(PIN2,UIN2,HIN2,PHI1(NCM),PHI2(NCM),PHI3(NCM))
        ENDIF

        RETURN
        END

C....................................................................

        SUBROUTINE DECOD(P,U,H,F1,F2,F3)
        GAMMA=1.4
        GB=1./(GAMMA-1.)
        GA=GAMMA/(GAMMA-1.)
        RHO=P/H*GA
        E=P*GB+.5*RHO*U**2
        F1=RHO*U
        F2=P+RHO*U**2
        F3=U*(P+E)
        RETURN
        END

C....................................................................

        SUBROUTINE ITERSTAR(A,B,C)
        GAMMA=1.4
        GB=1./(GAMMA-1.)
        TRY1 = .5*A/GB
        KIP = 0
1       CONTINUE
        KIP = KIP+1
        IF(KIP.GE.500) THEN
        WRITE(*,*)'ITERATION FAILS IN ITERSTAR !!!!!'
        STOP
        ENDIF
C       TRY2 = B-A*ALOG(TRY1)
        TRY2 = EXP((B-TRY1)/A)
        IF(ABS(TRY2-TRY1).LT.1.E-5) GO TO 2
        TRY1 = TRY2
        GO TO 1
2       CONTINUE
        C = TRY2
        RETURN
        END

C....................................................................

        SUBROUTINE SUB(COST,TEMP,RATAREA,PTOTS)
        KIP=0
        TR1=1.
1       KIP=KIP+1
        TR2=1.-(COST/PTOTS**2/TR1**5/RATAREA**2)
        IF(ABS(TR2-TR1).LT.1.E-5) GO TO 2
        IF(KIP.GT.1000) WRITE(*,*)'ITERAZIONE CIUCCA IN SUB'
        IF(KIP.GT.1000) STOP
        TR1=TR2
        GO TO 1
2       CONTINUE
        TEMP=TR2
        RETURN
        END

C....................................................................

        SUBROUTINE SUPER(COST,TEMP,RATAREA,PTOTS)
        KIP=0
        TR1=.7
1       KIP=KIP+1
        TR2=(COST/PTOTS**2/(1.-TR1)/RATAREA**2)**.2
C       WRITE(*,*)KIP,TR1,TR2
        IF(ABS(TR2-TR1).LT.1.E-5) GO TO 2
        IF(KIP.GT.1000) WRITE(*,*)'ITERAZIONE CIUCCA IN SUPER'
        IF(KIP.GT.1000) STOP
        TR1 =TR2
        GO TO 1
2       CONTINUE
        TEMP=TR2
        RETURN
        END

C....................................................................

        SUBROUTINE SHOCK(AM2,RP,RR,RPT)
        RP=(7.*AM2-1.)/6.
        RR=6.*AM2/(AM2+5.)
        RPT=RR**3.5/RP**2.5
        RETURN
        END

C....................................................................

        SUBROUTINE OUTPUT

        INCLUDE 'COMST05.INC'

        COMMON PTH(NMAX),TTH(NMAX),UTH(NMAX),STH(NMAX),RHOTH(NMAX),
     1         ATH(NMAX),ETH(NMAX),AMACHTH(NMAX),PTOTTH(NMAX),
     2         TTOTTH(NMAX),FLOWTH(NMAX),FLHTTH(NMAX),HTH(NMAX),
     3         HTOTTH(NMAX),R1TH(NMAX),R2TH(NMAX),R3TH(NMAX),
     4         XSHTH,FTHSHU,FTHSHD,PU,PD,RHOU,RHOD,TU,TD,SU,SD,
     5         UU,UD,AMU,AMD,HTOTU,HTOTD,FLOWU,FLOWD,FLHTU,FLHTD,
     6         R1U,R1D,R2U,R2D,R3U,R3D,X1,X2,X3,X4,
     7  P1L,P1R,U1L,U1R,S1L,S1R,A1L,A1R,T1L,T1R,RHO1L,RHO1R,E1L,E1R,
     8  P2L,P2R,U2L,U2R,S2L,S2R,A2L,A2R,T2L,T2R,RHO2L,RHO2R,E2L,E2R,
     9  P3L,P3R,U3L,U3R,S3L,S3R,A3L,A3R,T3L,T3R,RHO3L,RHO3R,E3L,E3R,
     9  P4L,P4R,U4L,U4R,S4L,S4R,A4L,A4R,T4L,T4R,RHO4L,RHO4R,E4L,E4R
        COMMON/OUT/LINE

1001    FORMAT(A1)

        CALL THESOL

        AMSH=1.
        DO 200 N=1,NA
        DIFF=(AMACH(N+1)-AMSH)*(AMACH(N)-AMSH)
        IF(DIFF.GT.0) GO TO 200
        XSH = X(N)+
     1       (AMSH-AMACH(N))/(AMACH(N+1)-AMACH(N))*(X(N+1)-X(N))
        XSH=XSH/C
200     CONTINUE
        IF(ITEST.EQ.2.OR.ITEST.EQ.3) THEN
        WRITE(*,*)'the length of the nozzle is  C = ',C
        IF(XSHTH.GT..0) THEN
        WRITE(*,*)'the shock is located at XSH/C = ',XSH
        WRITE(*,*)'its theoretical location is   = ',XSHTH
        ENDIF
        ENDIF
                
10      CONTINUE

        WRITE(*,*)'INPUT:         1: GNUPLOT '
        WRITE(*,*)'               2: STAMPA  '
        WRITE(*,*)'               3: RETURN  '
        WRITE(*,*)'               4: STOP    '
        WRITE(*,*)'               5: CLEAN   '
        READ (*,*)INPUT

	IF(INPUT.EQ.3) RETURN
        IF(INPUT.EQ.4) STOP
        IF(INPUT.EQ.5) CALL CLEAN
        WRITE(*,*)'Input LINE for the diagram'
	WRITE(*,*)'P(1)   RHO(2) T(3)   U(4)   MACHN(5) '
	WRITE(*,*)'S(6)   HT(7)  FL(8)  FLHT(9)         '
	WRITE(*,*)'R1(10)  R2(11)  R3(12)                '
	READ(*,*)LINE

        IF(ITEST.EQ.1) THEN
        IF(LINE.EQ.1) THEN
        FTH1L = P1L
        FTH1R = P1R
        FTH2L = P2L
        FTH2R = P2R
        FTH3L = P3L
        FTH3R = P3R
        FTH4L = P4L
        FTH4R = P4R
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHO1L
        FTH1R = RHO1R
        FTH2L = RHO2L
        FTH2R = RHO2R
        FTH3L = RHO3L
        FTH3R = RHO3R
        FTH4L = RHO4L
        FTH4R = RHO4R
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = T1L
        FTH1R = T1R
        FTH2L = T2L
        FTH2R = T2R
        FTH3L = T3L
        FTH3R = T3R
        FTH4L = T4L
        FTH4R = T4R
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = U1L
        FTH1R = U1R
        FTH2L = U2L
        FTH2R = U2R
        FTH3L = U3L
        FTH3R = U3R
        FTH4L = U4L
        FTH4R = U4R
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = U1L/SQRT(GAMMA*T1L)
        FTH1R = U1R/SQRT(GAMMA*T1R)
        FTH2L = U2L/SQRT(GAMMA*T2L)
        FTH2R = U2R/SQRT(GAMMA*T2R)
        FTH3L = U3L/SQRT(GAMMA*T3L)
        FTH3R = U3R/SQRT(GAMMA*T3R)
        FTH4L = U4L/SQRT(GAMMA*T4L)
        FTH4R = U4R/SQRT(GAMMA*T4R)
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = S1L
        FTH1R = S1R
        FTH2L = S2L
        FTH2R = S2R
        FTH3L = S3L
        FTH3R = S3R
        FTH4L = S4L
        FTH4R = S4R
        ENDIF
        DO 1100 N=2,NCM
        IF(LINE.EQ.1) DUMDUM(N)  = PTH(N)
        IF(LINE.EQ.2) DUMDUM(N)  = RHOTH(N)
        IF(LINE.EQ.3) DUMDUM(N)  = TTH(N)
        IF(LINE.EQ.4) DUMDUM(N)  = UTH(N)
        IF(LINE.EQ.5) DUMDUM(N)  = AMACHTH(N)
        IF(LINE.EQ.6) DUMDUM(N)  = STH(N)
1100    CONTINUE
        IK=1
        DO 1101 N=2,NCM
        IF(X(N-1).LT.X1.AND.X(N).GE.X1)THEN
        IK=IK+1
        XTH(IK)     = X1
        FPLOTTH(IK) = FTH1L
        IK=IK+1
        XTH(IK)     = X1
        FPLOTTH(IK) = FTH1R
        ENDIF          
        IF(X(N-1).LT.X2.AND.X(N).GE.X2)THEN
        IK=IK+1
        XTH(IK)     = X2
        FPLOTTH(IK) = FTH2L
        IK=IK+1
        XTH(IK)     = X2
        FPLOTTH(IK) = FTH2R
        ENDIF          
        IF(X(N-1).LT.X3.AND.X(N).GE.X3)THEN
        IK=IK+1
        XTH(IK)     = X3
        FPLOTTH(IK) = FTH3L
        IK=IK+1
        XTH(IK)     = X3
        FPLOTTH(IK) = FTH3R
        ENDIF          
        IF(X(N-1).LT.X4.AND.X(N).GE.X4)THEN
        IK=IK+1
        XTH(IK)     = X4
        FPLOTTH(IK) = FTH4L
        IK=IK+1
        XTH(IK)     = X4
        FPLOTTH(IK) = FTH4R
        ENDIF          
        IK=IK+1
        XTH(IK)     = X(N)
        FPLOTTH(IK) = DUMDUM(N)
1101    CONTINUE
      
        OPEN(10)
        DO 1103 N=2,NCM+8
        WRITE( 10,*)XTH(N),FPLOTTH(N)
1103    CONTINUE
        WRITE(10,*)
        CLOSE(10)

        DO 1102 N=2,NCM
        IF(LINE.EQ.1) FPLOT(N)  = P(N)
        IF(LINE.EQ.2) FPLOT(N)  = RHO(N)
        IF(LINE.EQ.3) FPLOT(N)  = T(N)
        IF(LINE.EQ.4) FPLOT(N)  = U(N)
        IF(LINE.EQ.5) FPLOT(N)  = AMACH(N)
        IF(LINE.EQ.6) FPLOT(N)  = S(N)
1102    CONTINUE
        
        ENDIF

        IF(ITEST.EQ.2.OR.ITEST.EQ.3) THEN
        IF(LINE.EQ.1) THEN
        FTH1L = PU
        FTH1R = PD
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHOU
        FTH1R = RHOD
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = TU
        FTH1R = TD
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = UU
        FTH1R = UD
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = UU/SQRT(GAMMA*TU)
        FTH1R = UD/SQRT(GAMMA*TD)
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = SU
        FTH1R = SD
        ENDIF
        DO 2200 N=2,NCM
        IF(LINE.EQ.1) DUMDUM(N)  = PTH(N)
        IF(LINE.EQ.2) DUMDUM(N)  = RHOTH(N)
        IF(LINE.EQ.3) DUMDUM(N)  = TTH(N)
        IF(LINE.EQ.4) DUMDUM(N)  = UTH(N)
        IF(LINE.EQ.5) DUMDUM(N)  = AMACHTH(N)
        IF(LINE.EQ.6) DUMDUM(N)  = STH(N)
2200    CONTINUE
        IK=1
        DO 2201 N=2,NCM
        IF(X(N-1).LT.XSHTH.AND.X(N).GE.XSHTH)THEN
        IK=IK+1
        XTH(IK)     = XSHTH
        FPLOTTH(IK) = FTH1L
        IK=IK+1
        XTH(IK)     = XSHTH
        FPLOTTH(IK) = FTH1R
        ENDIF          
        IK=IK+1
        XTH(IK)     = X(N)
        FPLOTTH(IK) = DUMDUM(N)
2201    CONTINUE

        OPEN(10)
        DO 2203 N=2,NCM+2
        IF(XSHTH.GE.0.0) THEN
        WRITE(10,*)XTH(N),FPLOTTH(N)
        ELSE
        IF(N.LE.NCM) WRITE(10,*)XTH(N),FPLOTTH(N) 
        IF(N.GT.NCM) WRITE(10,*)XTH(NCM),FPLOTTH(NCM) 
        ENDIF
2203    CONTINUE
        WRITE(10,*)
        CLOSE(10)

        DO 2202 N=2,NCM
        IF(LINE.EQ.1) FPLOT(N)  = P(N)
        IF(LINE.EQ.2) FPLOT(N)  = RHO(N)
        IF(LINE.EQ.3) FPLOT(N)  = T(N)
        IF(LINE.EQ.4) FPLOT(N)  = U(N)
        IF(LINE.EQ.5) FPLOT(N)  = AMACH(N)
        IF(LINE.EQ.6) FPLOT(N)  = S(N)
2202    CONTINUE
        ENDIF

        IF(ITEST.EQ.5) THEN
        IF(LINE.EQ.1) THEN
        FTH1L = P1L
        FTH1R = P1R
        FTH2L = P2L
        FTH2R = P2R
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHO1L
        FTH1R = RHO1R
        FTH2L = RHO2L
        FTH2R = RHO2R
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = T1L
        FTH1R = T1R
        FTH2L = T2L
        FTH2R = T2R
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = U1L
        FTH1R = U1R
        FTH2L = U2L
        FTH2R = U2R
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = U1L/A1L
        FTH1R = U1R/A1R
        FTH2L = U2L/A2L
        FTH2R = U2R/A2R
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = S1L
        FTH1R = S1R
        FTH2L = S2L
        FTH2R = S2R
        ENDIF
        DO 5500 N=2,NCM
        IF(LINE.EQ.1) DUMDUM(N)  = PTH(N)
        IF(LINE.EQ.2) DUMDUM(N)  = RHOTH(N)
        IF(LINE.EQ.3) DUMDUM(N)  = TTH(N)
        IF(LINE.EQ.4) DUMDUM(N)  = UTH(N)
        IF(LINE.EQ.5) DUMDUM(N)  = AMACHTH(N)
        IF(LINE.EQ.6) DUMDUM(N)  = STH(N)
5500    CONTINUE
      
        IK=1
        DO 5501 N=2,NCM
        IF(X(N-1).LT.X1.AND.X(N).GE.X1)THEN
        IK=IK+1
        XTH(IK)     = X1
        FPLOTTH(IK) = FTH1L
        IK=IK+1
        XTH(IK)     = X1
        FPLOTTH(IK) = FTH1R
        ENDIF          
        IF(X(N-1).LT.X2.AND.X(N).GE.X2)THEN
        IK=IK+1
        XTH(IK)     = X2
        FPLOTTH(IK) = FTH2L
        IK=IK+1
        XTH(IK)     = X2
        FPLOTTH(IK) = FTH2R
        ENDIF          
        IK=IK+1
        XTH(IK)     = X(N)
        FPLOTTH(IK) = DUMDUM(N)
5501    CONTINUE


        OPEN(10)
        DO 5503 N=2,NCM+4
        WRITE( 10,*)XTH(N),FPLOTTH(N)
5503    CONTINUE
        WRITE(10,*)
        CLOSE(10)
        ENDIF
        
        IF(ITEST.EQ.6) THEN
        T1  = P1/RHO1
        S1  = 1.4*ALOG(T1)-.4*ALOG(P1)
        T2  = P2/RHO2
        S2  = 1.4*ALOG(T2)-.4*ALOG(P2)
        T3  = P3/RHO3
        S3  = 1.4*ALOG(T3)-.4*ALOG(P3)
        T4  = P4/RHO4
        S4  = 1.4*ALOG(T4)-.4*ALOG(P4)
        XSHOINT = (XSH10-XSH20*VSH1/VSH2)/(1.-VSH1/VSH2)
        TIMEINT = (XSHOINT-XSH10)/VSH1
        WRITE(*,*)'XSHOINT,TIMEINT',XSHOINT,TIMEINT
        IF(TIME.LE.TIMEINT) THEN
        XSH1 = XSH10+VSH1*TIME
        XSH2 = XSH20+VSH2*TIME
        ENDIF
        IF(TIME.GT.TIMEINT) THEN
        XSH3 = XSHOINT+VSH3*(TIME-TIMEINT)
        XSH4 = XSHOINT+VSH4*(TIME-TIMEINT)
        XSH5 = XSHOINT+VSH5*(TIME-TIMEINT)
        ENDIF
        IF(LINE.EQ.1) THEN
        FTH1L = P1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = P2
        FTH3L = P4
        FTH3R = P2
        FTH4L = P1
        FTH4R = P3
        FTH5L = P3
        FTH5R = P4
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHO1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = RHO2
        FTH3L = RHO4
        FTH3R = RHO2
        FTH4L = RHO1
        FTH4R = RHO3
        FTH5L = RHO3
        FTH5R = RHO4
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = T1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = T2
        FTH3L = T4
        FTH3R = T2
        FTH4L = T1
        FTH4R = T3
        FTH5L = T3
        FTH5R = T4
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = U1
        FTH1R = 0.0
        FTH2L = 0.0
        FTH2R = U2
        FTH3L = U4
        FTH3R = U2
        FTH4L = U1
        FTH4R = U3
        FTH5L = U3
        FTH5R = U4
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = U1/SQRT(GAMMA*T1)
        FTH1R = 1.0/SQRT(GAMMA)
        FTH2L = 1.0/SQRT(GAMMA)
        FTH2R = U2/SQRT(GAMMA*T2)
        FTH3L = U4/SQRT(GAMMA*T4)
        FTH3R = U2/SQRT(GAMMA*T2)
        FTH4L = U1/SQRT(GAMMA*T1)
        FTH4R = U3/SQRT(GAMMA*T3)
        FTH5L = U3/SQRT(GAMMA*T3)
        FTH5R = U4/SQRT(GAMMA*T4)
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = S1
        FTH1R = 0.0
        FTH2L = 0.0
        FTH2R = S2
        FTH3L = S4
        FTH3R = S2
        FTH4L = S1
        FTH4R = S3
        FTH5L = S3
        FTH5R = S4
        ENDIF
        IF(TIME.LE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH1L
        XTH(2)     = XSH1
        FPLOTTH(2) = FTH1L
        XTH(3)     = XSH1
        FPLOTTH(3) = FTH1R
        XTH(4)     = XSH2
        FPLOTTH(4) = FTH2L
        XTH(5)     = XSH2
        FPLOTTH(5) = FTH2R
        XTH(6)     = 1.0
        FPLOTTH(6) = FTH2R
        XTH(7)     = XTH(6)
        FPLOTTH(7) = FPLOTTH(6)
        XTH(8)     = XTH(7)
        FPLOTTH(8) = FPLOTTH(7)
        ENDIF
        IF(TIME.GE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH4L
        XTH(2)     = XSH4
        FPLOTTH(2) = FTH4L
        XTH(3)     = XSH4
        FPLOTTH(3) = FTH4R
        XTH(4)     = XSH5
        FPLOTTH(4) = FTH5L
        XTH(5)     = XSH5
        FPLOTTH(5) = FTH5R
        XTH(6)     = XSH3
        FPLOTTH(6) = FTH3L
        XTH(7)     = XSH3
        FPLOTTH(7) = FTH3R
        XTH(8)     = 1.0
        FPLOTTH(8) = FTH3R
        ENDIF
    
        OPEN(10)
        DO 6603 N=1,8
        WRITE( 10,*)XTH(N),FPLOTTH(N)
6603    CONTINUE
        WRITE(10,*)
        CLOSE(10)
        ENDIF

        IF(ITEST.EQ.7) THEN
        T1  = P1/RHO1
        S1  = 1.4*ALOG(T1)-.4*ALOG(P1)
        T2  = P2/RHO2
        S2  = 1.4*ALOG(T2)-.4*ALOG(P2)
        T3  = P3/RHO3
        S3  = 1.4*ALOG(T3)-.4*ALOG(P3)
        T4  = P4/RHO4
        S4  = 1.4*ALOG(T4)-.4*ALOG(P4)
        XSHOINT = (XSH10-XSH20*VSH1/VSH2)/(1.-VSH1/VSH2)
        TIMEINT = (XSHOINT-XSH10)/VSH1
        WRITE(*,*)'XSHOINT,TIMEINT',XSHOINT,TIMEINT
        IF(TIME.LE.TIMEINT) THEN
        XSH1 = XSH10+VSH1*TIME
        XSH2 = XSH20+VSH2*TIME
        ENDIF
        IF(TIME.GT.TIMEINT) THEN
        XSH3 = XSHOINT+VSH3 *(TIME-TIMEINT)
        XSH4L= XSHOINT+VSH4L*(TIME-TIMEINT)
        XSH4R= XSHOINT+VSH4R*(TIME-TIMEINT)
        XSH5 = XSHOINT+VSH5 *(TIME-TIMEINT)
        ENDIF
        IF(LINE.EQ.1) THEN
        FTH1L = P1
        FTH1R = 1.0
        FTH2L = P2
        FTH2R = P1
        FTH3L = P4
        FTH3R = 1.0
        FTH4L = P2
        FTH4R = P3
        FTH5L = P3
        FTH5R = P4
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHO1
        FTH1R = 1.0
        FTH2L = RHO2
        FTH2R = RHO1
        FTH3L = RHO4
        FTH3R = 1.0
        FTH4L = RHO2
        FTH4R = RHO3
        FTH5L = RHO3
        FTH5R = RHO4
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = T1
        FTH1R = 1.0
        FTH2L = T2
        FTH2R = T1
        FTH3L = T4
        FTH3R = 1.0
        FTH4L = T2
        FTH4R = T3
        FTH5L = T3
        FTH5R = T4
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = U1
        FTH1R = 0.0
        FTH2L = U2
        FTH2R = U1
        FTH3L = U4
        FTH3R = 0.0
        FTH4L = U2
        FTH4R = U3
        FTH5L = U3
        FTH5R = U4
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = U1/SQRT(GAMMA*T1)
        FTH1R = 0.0/SQRT(GAMMA)
        FTH2L = U2/SQRT(GAMMA*T2)
        FTH2R = U1/SQRT(GAMMA*T1)
        FTH3L = U4/SQRT(GAMMA*T4)
        FTH3R = 0.0/SQRT(GAMMA)
        FTH4L = U2/SQRT(GAMMA*T2)
        FTH4R = U3/SQRT(GAMMA*T3)
        FTH5L = U3/SQRT(GAMMA*T3)
        FTH5R = U4/SQRT(GAMMA*T4)
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = S1
        FTH1R = 0.0
        FTH2L = S2
        FTH2R = S1
        FTH3L = S4
        FTH3R = 0.0
        FTH4L = S2
        FTH4R = S3
        FTH5L = S3
        FTH5R = S4
        ENDIF

        IF(TIME.LE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH2L
        XTH(2)     = XSH2
        FPLOTTH(2) = FTH2L
        XTH(3)     = XSH2
        FPLOTTH(3) = FTH2R
        XTH(4)     = XSH1
        FPLOTTH(4) = FTH1L
        XTH(5)     = XSH1
        FPLOTTH(5) = FTH1R
        XTH(6)     = 1.0
        FPLOTTH(6) = FTH1R
        XTH(7)     = XTH(6)
        FPLOTTH(7) = FPLOTTH(6)
        XTH(8)     = XTH(7)
        FPLOTTH(8) = FPLOTTH(7)
        ENDIF

        IF(TIME.GE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH4L
        XTH(2)     = XSH4L
        FPLOTTH(2) = FTH4L
        XTH(3)     = XSH4R
        FPLOTTH(3) = FTH4R
        XTH(4)     = XSH5
        FPLOTTH(4) = FTH5L
        XTH(5)     = XSH5
        FPLOTTH(5) = FTH5R
        XTH(6)     = XSH3
        FPLOTTH(6) = FTH3L
        XTH(7)     = XSH3
        FPLOTTH(7) = FTH3R
        XTH(8)     = 1.0
        FPLOTTH(8) = FTH3R
        ENDIF        

        OPEN(10)
        DO 7703 N=1,8
        WRITE( 10,*)XTH(N),FPLOTTH(N)
7703    CONTINUE
        WRITE(10,*)
        CLOSE(10)
        ENDIF

        IF(ITEST.EQ.8) THEN
        T1  = P1/RHO1
        S1  = 1.4*ALOG(T1)-.4*ALOG(P1)
        T2  = P2/RHO2
        S2  = 1.4*ALOG(T2)-.4*ALOG(P2)
        T3  = P3/RHO3
        S3  = 1.4*ALOG(T3)-.4*ALOG(P3)
        T4  = P4/RHO4
        S4  = 1.4*ALOG(T4)-.4*ALOG(P4)

        XSHOINT = XCSF2
        TIMEINT = (XSHOINT-XSH10)/VSH1
        WRITE(*,*)'XSHOINT,TIMEINT',XSHOINT,TIMEINT

        IF(TIME.LE.TIMEINT) THEN
        XSH1 = XSH10+VSH1*TIME
        XSH2 = XCSF2
        ENDIF

        IF(TIME.GT.TIMEINT) THEN
        IF(P2/RHO2.GT.1) THEN
        XSH3L= XSHOINT+VSH3L*(TIME-TIMEINT)
        XSH3R= XSHOINT+VSH3R*(TIME-TIMEINT)
        XSH4L= XSHOINT+VSH4 *(TIME-TIMEINT)
        XSH4R= XSHOINT+VSH4 *(TIME-TIMEINT)
        XSH5 = XSHOINT+VSH5 *(TIME-TIMEINT)
        ENDIF
        IF(P2/RHO2.LE.1) THEN
        XSH3L= XSHOINT+VSH3 *(TIME-TIMEINT)
        XSH3R= XSHOINT+VSH3 *(TIME-TIMEINT)
        XSH4L= XSHOINT+VSH4 *(TIME-TIMEINT)
        XSH4R= XSHOINT+VSH4 *(TIME-TIMEINT)
        XSH5 = XSHOINT+VSH5 *(TIME-TIMEINT)
        ENDIF
        ENDIF

        IF(LINE.EQ.1) THEN
        FTH1L = P1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = P2
        FTH3L = P1
        FTH3R = P3
        FTH4L = P4
        FTH4R = P2
        FTH5L = P3
        FTH5R = P4
        ENDIF
        IF(LINE.EQ.2) THEN
        FTH1L = RHO1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = RHO2
        FTH3L = RHO1
        FTH3R = RHO3
        FTH4L = RHO4
        FTH4R = RHO2
        FTH5L = RHO3
        FTH5R = RHO4
        ENDIF
        IF(LINE.EQ.3) THEN
        FTH1L = T1
        FTH1R = 1.0
        FTH2L = 1.0
        FTH2R = T2
        FTH3L = T1
        FTH3R = T3
        FTH4L = T4
        FTH4R = T2
        FTH5L = T3
        FTH5R = T4
        ENDIF
        IF(LINE.EQ.4) THEN
        FTH1L = U1
        FTH1R = 0.0
        FTH2L = 0.0
        FTH2R = U2
        FTH3L = U1
        FTH3R = U3
        FTH4L = U4
        FTH4R = U2
        FTH5L = U3
        FTH5R = U4
        ENDIF
        IF(LINE.EQ.5) THEN
        FTH1L = U1/SQRT(GAMMA*T1)
        FTH1R = 0.0/SQRT(GAMMA)
        FTH2L = 0.0/SQRT(GAMMA)
        FTH2R = U2/SQRT(GAMMA*T2)
        FTH3L = U1/SQRT(GAMMA*T1)
        FTH3R = U3/SQRT(GAMMA*T3)
        FTH4L = U4/SQRT(GAMMA*T4)
        FTH4R = U2/SQRT(GAMMA*T2)
        FTH5L = U3/SQRT(GAMMA*T3)
        FTH5R = U4/SQRT(GAMMA*T4)
        ENDIF
        IF(LINE.EQ.6) THEN
        FTH1L = S1
        FTH1R = 0.0
        FTH2L = 0.0
        FTH2R = S2
        FTH3L = S1
        FTH3R = S3
        FTH4L = S4
        FTH4R = S2
        FTH5L = S3
        FTH5R = S4
        ENDIF

        IF(TIME.LE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH1L
        XTH(2)     = XSH1
        FPLOTTH(2) = FTH1L
        XTH(3)     = XSH1
        FPLOTTH(3) = FTH1R
        XTH(4)     = XSH2
        FPLOTTH(4) = FTH2L
        XTH(5)     = XSH2
        FPLOTTH(5) = FTH2R
        XTH(6)     = 1.0
        FPLOTTH(6) = FTH2R
        XTH(7)     = XTH(6)
        FPLOTTH(6) = FPLOTTH(6)
        XTH(8)     = XTH(7)
        FPLOTTH(8) = FPLOTTH(7)
        ENDIF

        IF(TIME.GE.TIMEINT) THEN
        XTH(1)     = 0.0
        FPLOTTH(1) = FTH1L
        XTH(2)     = XSH3L
        FPLOTTH(2) = FTH3L
        XTH(3)     = XSH3R
        FPLOTTH(3) = FTH3R
        XTH(4)     = XSH5
        FPLOTTH(4) = FTH5L
        XTH(5)     = XSH5
        FPLOTTH(5) = FTH5R
        XTH(6)     = XSH4L
        FPLOTTH(6) = FTH4L
        XTH(7)     = XSH4R
        FPLOTTH(7) = FTH4R
        XTH(8)     = 1.0
        FPLOTTH(8) = FTH4R
        ENDIF

        OPEN(10)
        DO 8803 N=1,8
        WRITE( 10,*)XTH(N),FPLOTTH(N)
8803    CONTINUE
        WRITE(10,*)
        CLOSE(10)
        ENDIF
C...............................................................

        IF(ITEST.GE.4) THEN
        DO 4402 N=2,NCM
        IF(LINE.EQ.1) FPLOT(N)  = P(N)
        IF(LINE.EQ.2) FPLOT(N)  = RHO(N)
        IF(LINE.EQ.3) FPLOT(N)  = T(N)
        IF(LINE.EQ.4) FPLOT(N)  = U(N)
        IF(LINE.EQ.5) FPLOT(N)  = AMACH(N)
        IF(LINE.EQ.6) FPLOT(N)  = S(N)
4402    CONTINUE
        ENDIF

        
c????????????
         
        IF(INPUT.EQ.1) CALL GNUPLOT
        IF(INPUT.EQ.2) CALL STAMPA

        GO TO 10
	END
C.....................................................................

        SUBROUTINE GNUPLOT
        INCLUDE 'COMST05.INC'
        
c        OPEN(UNIT=10,ACCESS='APPEND')
C        OPEN(UNIT=11,ACCESS='APPEND')

        OPEN(11)
        DO 1 N=2,NCM
C        write( *,*)x(n),fplot(n)
        WRITE(11,*)X(N),FPLOT(N)
1       CONTINUE
        WRITE(11,*)
        CLOSE(11)

        RETURN
        END

C___________________________________________________________

        SUBROUTINE STAMPA
C
        INCLUDE 'COMST05.INC'
        COMMON PTH(NMAX),TTH(NMAX),UTH(NMAX),STH(NMAX),RHOTH(NMAX),
     1         ATH(NMAX),ETH(NMAX),AMACHTH(NMAX),PTOTTH(NMAX),
     2         TTOTTH(NMAX),FLOWTH(NMAX),FLHTTH(NMAX),HTH(NMAX),
     3         HTOTTH(NMAX),R1TH(NMAX),R2TH(NMAX),R3TH(NMAX),
     4         XSHTH,FTHSHU,FTHSHD,PU,PD,RHOU,RHOD,TU,TD,SU,SD,
     5         UU,UD,AMU,AMD,HTOTU,HTOTD,FLOWU,FLOWD,FLHTU,FLHTD,
     6         R1U,R1D,R2U,R2D,R3U,R3D,X1,X2,X3,X4,
     7  P1L,P1R,U1L,U1R,S1L,S1R,A1L,A1R,T1L,T1R,RHO1L,RHO1R,E1L,E1R,
     8  P2L,P2R,U2L,U2R,S2L,S2R,A2L,A2R,T2L,T2R,RHO2L,RHO2R,E2L,E2R,
     9  P3L,P3R,U3L,U3R,S3L,S3R,A3L,A3R,T3L,T3R,RHO3L,RHO3R,E3L,E3R,
     9  P4L,P4R,U4L,U4R,S4L,S4R,A4L,A4R,T4L,T4R,RHO4L,RHO4R,E4L,E4R
C
100     FORMAT(5X,'STEP K =  ',I5,'  AND  TIME =   ',F8.4)
101     FORMAT(10X,'PRESSURE')
102     FORMAT(10X,'DENSITY')
103     FORMAT(10X,'TEMPERATURE')
104     FORMAT(10X,'VELOCITY')
105     FORMAT(10X,'MACH NUMBER')
106     FORMAT(10X,'ENTROPY')
107     FORMAT(10X,'TOTAL ENTHALPY')
108     FORMAT(10X,'MASS FLOW')
109     FORMAT(10X,'TOTAL ENTHALPHY FLOW')
110     FORMAT(10X,'characteristic variable R1')
111     FORMAT(10X,'characteristic variable R2')
112     FORMAT(10X,'characteristic variable R3')
150     FORMAT(1X,10F7.3)
151     FORMAT(1X,10F8.0)

        INFORM=0
        WRITE(*,*)'information (WRITE(*,*) 1 to get them)'
        READ(*,*)INFORM
        IF(INFORM.EQ.1) THEN
        WRITE(*,*)'KA and NA =',KA,NA,'       STAB =',STAB
        IF(IORD.EQ.1)   WRITE(*,*)'IORD ( 1 = first order)'
        IF(IORD.EQ.21)  WRITE(*,*)'IORD (21 = ENO 1st level)'
        IF(IORD.EQ.22)  WRITE(*,*)'IORD (22 = ENO 2nd level)'
        WRITE(*,*)'nozzle DIVERGence and length =',DIVERG,C
        WRITE(*,*)'inlet Mach number    = ',AMACH(2)
        WRITE(*,*)'exit static pressure =',PEX
        ENDIF

99      CONTINUE
        WRITE(6,100) K,TIME
        WRITE(*,*)'Input LINE for printing'
        WRITE(*,*)'P(1)   RHO(2) T(3)   U(4)   MACHN(5) '
        WRITE(*,*)'S(6)   HT(7)  FL(8)  FLHT(9)         '
        WRITE(*,*)'R1(10)  R2(11)  R3(12)                '
        WRITE(*,*)'with theoretical values add 1000'
        READ(*,*)LINE
        IF(LINE.EQ.0) RETURN
        IMULT=0
        WRITE(*,*)'give IMULT to multiply by 10**IMULT (/=0)'
        READ(*,*)IMULT
        IF(LINE.EQ.1.OR.LINE.EQ.1001) GO TO 1
        IF(LINE.EQ.2.OR.LINE.EQ.1002 ) GO TO 2
        IF(LINE.EQ.3.OR.LINE.EQ.1003 ) GO TO 3
        IF(LINE.EQ.4.OR.LINE.EQ.1004 ) GO TO 4
        IF(LINE.EQ.5.OR.LINE.EQ.1005 ) GO TO 5
        IF(LINE.EQ.6.OR.LINE.EQ.1006 ) GO TO 6
        IF(LINE.EQ.7.OR.LINE.EQ.1007 ) GO TO 7
        IF(LINE.EQ.8.OR.LINE.EQ.1008 ) GO TO 8
        IF(LINE.EQ.9.OR.LINE.EQ.1009 ) GO TO 9
        IF(LINE.EQ.10.OR.LINE.EQ.1010 ) GO TO 10
        IF(LINE.EQ.11.OR.LINE.EQ.1011 ) GO TO 11
        IF(LINE.EQ.12.OR.LINE.EQ.1012 ) GO TO 12
1       WRITE(6,101)
        WRITE(6,150)(10.**IMULT* P(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* PTH(N),N=2,NCM)
        GO TO 99
2       WRITE(6,102)
        WRITE(6,150)(10.**IMULT* RHO(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* RHOTH(N),N=2,NCM)
        GO TO 99
3       WRITE(6,103)
        WRITE(6,150)(10.**IMULT* T(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* TTH(N),N=2,NCM)
        GO TO 99
4       WRITE(6,104)
        WRITE(6,150)(10.**IMULT* U(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* UTH(N),N=2,NCM)
        GO TO 99
5       WRITE(6,105)
        WRITE(6,150)(10.**IMULT* AMACH(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* AMACHTH(N),N=2,NCM)
        GO TO 99
6       WRITE(6,106)
        WRITE(6,150)(10.**IMULT* S(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* STH(N),N=2,NCM)
        GO TO 99
7       WRITE(6,107)
        WRITE(6,150)(10.**IMULT* HTOT(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* HTOTTH(N),N=2,NCM)
        GO TO 99
8       WRITE(6,108)
        WRITE(6,150)(10.**IMULT* FLOW(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* FLOWTH(N),N=2,NCM)
        GO TO 99
9       WRITE(6,109)
        WRITE(6,150)(10.**IMULT* FLHT(N),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* FLHTTH(N),N=2,NCM)
        GO TO 99
10      WRITE(6,110)
C       WRITE(6,150)(10.**IMULT*(P(N)-RHO(N)*A(N)*U(N)),N=2,NCM)
        WRITE(6,150)(10.**IMULT*(5.*A(N)-U(N)),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* R1TH(N),N=2,NCM)
        GO TO 99
11      WRITE(6,111)
        WRITE(6,150)(10.**IMULT*(H(N)-P(N)/RHO(N)),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* R2TH(N),N=2,NCM)
        GO TO 99
12      WRITE(6,112)
        WRITE(6,150)(10.**IMULT*(5.*A(N)+U(N)),N=2,NCM)
        WRITE(*,*)'theoretical values'
        IF(LINE.GT.1000) WRITE(6,150)(10.**IMULT* R3TH(N),N=2,NCM)
        GO TO 99
        END

C       ..................................................

        SUBROUTINE THESOL

        INCLUDE 'COMST05.INC'

        COMMON PTH(NMAX),TTH(NMAX),UTH(NMAX),STH(NMAX),RHOTH(NMAX),
     1         ATH(NMAX),ETH(NMAX),AMACHTH(NMAX),PTOTTH(NMAX),
     2         TTOTTH(NMAX),FLOWTH(NMAX),FLHTTH(NMAX),HTH(NMAX),
     3         HTOTTH(NMAX),R1TH(NMAX),R2TH(NMAX),R3TH(NMAX),
     4         XSHTH,FTHSHU,FTHSHD,PU,PD,RHOU,RHOD,TU,TD,SU,SD,
     5         UU,UD,AMU,AMD,HTOTU,HTOTD,FLOWU,FLOWD,FLHTU,FLHTD,
     6         R1U,R1D,R2U,R2D,R3U,R3D,X1,X2,X3,X4,
     7  P1L,P1R,U1L,U1R,S1L,S1R,A1L,A1R,T1L,T1R,RHO1L,RHO1R,E1L,E1R,
     8  P2L,P2R,U2L,U2R,S2L,S2R,A2L,A2R,T2L,T2R,RHO2L,RHO2R,E2L,E2R,
     9  P3L,P3R,U3L,U3R,S3L,S3R,A3L,A3R,T3L,T3R,RHO3L,RHO3R,E3L,E3R,
     9  P4L,P4R,U4L,U4R,S4L,S4R,A4L,A4R,T4L,T4R,RHO4L,RHO4R,E4L,E4R

        IF(ITEST.EQ.1) GO TO 1111
        IF(ITEST.EQ.2) GO TO 2222
        IF(ITEST.EQ.3) GO TO 3333
        IF(ITEST.EQ.4) RETURN
        IF(ITEST.EQ.5) GO TO 5555
        IF(ITEST.EQ.6) GO TO 6666
        IF(ITEST.EQ.7) GO TO 7777
        IF(ITEST.EQ.8) GO TO 8888

1111    CONTINUE
        IF(K.GT.100) RETURN
        AA=A(2)
        PA=P(2)
        UA=U(2)
        SA=S(2)
        AB=A(NC-1)
        PB=P(NC-1)
        UB=U(NC-1)
        SB=S(NC-1)
        UC=0.0
        DELU=.1
        KIP=0
11      CONTINUE
        KIP=KIP+1
        UCO= UC
        UC = UC+DELU
        UD = UC
        AC = AA-GD*UC
        PC = PA*(AC/AA)**(1./GI)
        UOA= UD/AB
        AMS= (GE*UOA+SQRT(GE*GE*UOA*UOA+4.))/2.
        PD = PB*(AMS*AMS-GI)/GH
        DIFF=PC-PD
        IF(ABS(DIFF).LT.1.E-5) GO TO 12
        IF(kip.gt.100) then
        write(*,*)'kip maggiore di 100 in iterazione THESOL'
        go to 12
        endif
        IF(PC.GT.PD) GO TO 11
        UC = UCO
        DELU=DELU/2.
        GO TO 11
12      CONTINUE
        SC = SA
        UD = UC
        ANUM=(AMS**2-GI)*(GI*AMS**2+1./GAMMA)
        ADEN=GH*GH*AMS**2
        AD = AB*SQRT(ANUM/ADEN)
        SD = 2.*GAMMA*ALOG(AD/AB)-(GAMMA-1.)*ALOG(PD/PB)
        ALAMA=UA-AA
        ALAMC=UC-AC
        ALAMX=UC
        VSH  =AMS*AB
C       TIMESOD=.411*KA/FLOAT(NA)
        TIMESOD=TIME
        X1=0.5+ALAMA*TIMESOD
        X2=0.5+ALAMC*TIMESOD
        X3=0.5+ALAMX*TIMESOD
        X4=0.5+VSH  *TIMESOD
        DO 13 N=2,NC
        IF(X(N).LT.X1) THEN
        PTH(N)=PA
        UTH(N)=UA
        STH(N)=SA
        ATH(N)=AA
        ENDIF
        IF(X(N).GE.X1.AND.X(N).LT.X2) THEN
        ALAM=(X(N)-0.5)/TIMESOD
        RIE =GG*AA
        ATH(N)=(RIE-ALAM)/GC
        UTH(N)=RIE-GG*ATH(N)
        STH(N)=SA
        PTH(N)=PA*(ATH(N)/AA)**(1./GI)
        ENDIF
        IF(X(N).GE.X2.AND.X(N).LT.X3) THEN
        PTH(N)=PC
        UTH(N)=UC
        STH(N)=SC
        ATH(N)=AC
        ENDIF
        IF(X(N).GE.X3.AND.X(N).LT.X4) THEN
        PTH(N)=PD
        UTH(N)=UD
        STH(N)=SD
        ATH(N)=AD
        ENDIF
        IF(X(N).GE.X4) THEN
        PTH(N)=PB
        UTH(N)=UB
        STH(N)=SB
        ATH(N)=AB
        ENDIF
        TTH(N)=ATH(N)**2/GAMMA
        RHOTH(N)=PTH(N)/TTH(N)
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        AMACHTH(N) = UTH(N)/ATH(N)
13      CONTINUE

        P1L=PA
        P1R=PA
        U1L=UA
        U1R=UA
        S1L=SA
        S1R=SA
        A1L=AA
        A1R=AA
        T1L=A1L**2/GAMMA
        T1R=A1R**2/GAMMA
        RHO1L=P1L/T1L
        RHO1R=P1R/T1R
        E1L=GB*P1L+.5*RHO1L*U1L**2
        E1R=GB*P1R+.5*RHO1R*U1R**2

        P2L=PC
        P2R=PC
        U2L=UC
        U2R=UC
        S2L=SC
        S2R=SC
        A2L=AC
        A2R=AC
        T2L=A2L**2/GAMMA
        T2R=A2R**2/GAMMA
        RHO2L=P2L/T2L
        RHO2R=P2R/T2R
        E2L=GB*P2L+.5*RHO2L*U2L**2
        E2R=GB*P2R+.5*RHO2R*U2R**2

        P3L=PC
        P3R=PD
        U3L=UC
        U3R=UD
        S3L=SC
        S3R=SD
        A3L=AC
        A3R=AD
        T3L=A3L**2/GAMMA
        T3R=A3R**2/GAMMA
        RHO3L=P3L/T3L
        RHO3R=P3R/T3R
        E3L=GB*P3L+.5*RHO3L*U3L**2
        E3R=GB*P3R+.5*RHO3R*U3R**2

        P4L=PD
        P4R=PB
        U4L=UD
        U4R=UB
        S4L=SD
        S4R=SB
        A4L=AD
        A4R=AB
        T4L=A4L**2/GAMMA
        T4R=A4R**2/GAMMA
        RHO4L=P4L/T4L
        RHO4R=P4R/T4R
        E4L=GB*P4L+.5*RHO4L*U4L**2
        E4R=GB*P4R+.5*RHO4R*U4R**2
        RETURN
2222    CONTINUE
        XSHTH  = -1.E6
        TIN    = 1./(1.+.2*ACHIN**2)
        AREAIN = 1.
        AREAEN = AREAIN*DIVERG
        AMEDIA(1) = AREAIN
        COST=.2*(2./2.4)**6
        ATHR = AREAIN*SQRT(TIN**5*(1.-TIN)/COST)
        DO 20 N=2,NCM
        AMEDIA(N)  =.5*(AREA(N+1)+AREA(N))
        IF(N.NE.1) DIFFAREA(N)= AMEDIA(N)-AMEDIA(N-1)

        PTOTS= 1.
        RATAREA=AREA(N)/ATHR
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TTH(N)=TEMP
        PTH(N)=TTH(N)**3.5*PTOTS
        AMACHTH(N)=SQRT(5.*(1./TTH(N)-1.))
        ATH(N)=SQRT(1.4*TTH(N))
        UTH(N)=SQRT(7.*(1.-TTH(N)))
        RHOTH(N)=PTH(N)/TTH(N)
        STH(N) = 1.4*ALOG(TTH(N))-.4*ALOG(PTH(N))
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        HTH(N)  =GA*TTH(N)
        PTOTTH(N)  =PTH(N)*(1.+GD*AMACHTH(N)**2)**GA
        TTOTTH(N)  =TTH(N)*(1.+GD*AMACHTH(N)**2)
        FLOWTH(N)= RHOTH(N)*UTH(N)*AREA(N)
        FLHTTH(N)= UTH(N)*(PTH(N)+ETH(N))*AREA(N)
        R1TH(N)  = PTH(N)-RHOTH(N)*ATH(N)*UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = PTH(N)+RHOTH(N)*ATH(N)*UTH(N)
        R1TH(N)  = 5.*ATH(N)-UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = 5.*ATH(N)+UTH(N)
20      CONTINUE

        XTRY = -1.E4
        XSHTH= -333.333
        IF(PEX.GT.PEXSUP) GO TO 30
        IF(PEX.LE.PEXSUP) GO TO 40
30      CONTINUE
        KIP=0
        XTRY=0
        DXTRY=.2
60      CONTINUE
        KIP=KIP+1
        IF(KIP.GT.1000) WRITE(*,*)'pressione a valle troppo alta'
        IF(KIP.GT.1000) STOP
        XTRY=XTRY+DXTRY
        AREADUM=1.+(DIVERG-1.)*XTRY/(C-B)
        RATAREA=AREADUM/ATHR
        PTOTS=1.
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TDUM= TEMP
        TU  = TDUM
        AMADUM=SQRT(5.*(1./TDUM-1.))
        AM2=AMADUM**2
        CALL SHOCK(AM2,RATP,RATR,RATPT)
        PTOTS=RATPT
        RATAREA=AREAEN/ATHR
        CALL SUB(COST,TTRY,RATAREA,PTOTS)
        PTRY=PTOTS*(TTRY/1.)**3.5
        IF(ABS(PTRY-PEX)/PEX.LT.1.E-5) GO TO 50
        IF(PTRY.GT.PEX) GO TO 60
        XTRY=XTRY-DXTRY
        DXTRY=DXTRY/2.
        GO TO 60
50      CONTINUE
        DO 70 N=2,NCM
        IF(X(N).LT.XTRY) GO TO 70
        RATAREA=AREA(N)/ATHR
        CALL SUB(COST,TEMP,RATAREA,PTOTS)
        TTH(N) = TEMP
        PTH(N) = PTOTS*(TEMP/1.)**3.5
        AMACHTH(N)=SQRT(5.*(1./TTH(N)-1.))
        ATH(N)=SQRT(1.4*TTH(N))
        UTH(N)=SQRT(7.*(1.-TTH(N)))
        RHOTH(N)=PTH(N)/TTH(N)
        STH(N) = 1.4*ALOG(TTH(N))-.4*ALOG(PTH(N))
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        HTH(N)  =GA*TTH(N)
        PTOTTH(N)  =PTH(N)*(1.+GD*AMACHTH(N)**2)**GA
        TTOTTH(N)  =TTH(N)*(1.+GD*AMACHTH(N)**2)
        FLOWTH(N)= RHOTH(N)*UTH(N)*AREA(N)
        FLHTTH(N)= UTH(N)*(PTH(N)+ETH(N))*AREA(N)
        R1TH(N)  = PTH(N)-RHOTH(N)*ATH(N)*UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = PTH(N)+RHOTH(N)*ATH(N)*UTH(N)
        R1TH(N)  = 5.*ATH(N)-UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = 5.*ATH(N)+UTH(N)
70      CONTINUE
        XSHTH=XTRY
        IF(XTRY.GE.0) THEN
        PU=TU**3.5
        RHOU=TU**2.5
        TU=PU/RHOU
        HU=3.5*TU
        SU=.0
        UU=SQRT((1.-TU)*7.)
        AU=SQRT(1.4*TU)
        AMU=UU/AU
        TTOTU=TU+1./7.*UU**2
        HTOTU=3.5*TTOTU
        PTOTU=PU*(1.+.2*AMU**2)**3.5
        FLOWU=FLOWTH(2)
        FLHTU=FLHTTH(2)
        R1U  = PU-RHOU*AU*UU
        R2U  = HU-PU/RHOU
        R3U  = PU+RHOU*AU*UU
        AMU2=AMU**2
        CALL SHOCK(AMU2,RATP,RATR,RATPT)
        PD=PU*RATP
        RHOD=RHOU*RATR
        TD=PD/RHOD
        HD=3.5*TD
        SD=-.4*ALOG(RATPT)
        UD=SQRT((1.-TD)*7.)
        AD=SQRT(1.4*TD)
        AMD=UD/AD
        TTOTD=TD+1./7.*UD**2
        HTOTD=3.5*TTOTD
        PTOTD=PD*(1.+.2*AMD**2)**3.5
        FLOWD=FLOWTH(2)
        FLHTD=FLHTTH(2)
        R1D  = PD-RHOD*AD*UD
        R2D  = HD-PD/RHOD
        R3D  = PD+RHOD*AD*UD
40      CONTINUE
        ENDIF
        RETURN

3333    CONTINUE

        XSHTH = -1.E6
        CCC1=2.5
        CCC2=0.3
        X000=.1
        AREAIN = CCC1*(0.00+X000)+CCC2/(0.00+X000)
        AREAEN = CCC1*(1.00+X000)+CCC2/(1.00+X000)
        XTHR   = SQRT(CCC2/CCC1)-X000
        ATHR   = CCC1*(XTHR+X000)+CCC2/(XTHR+X000)
        PTOTS= 1.
        COST=.2*(2./2.4)**6
        AMEDIA( 1 )=AREAIN
        DO 320 N=2,NCM
        AMEDIA(N)  =.5*(AREA(N+1)+AREA(N))
        IF(N.EQ.NCM) AMEDIA(N)=AREAEN
        DIFFAREA(N)= AMEDIA(N)-AMEDIA(N-1)
320     CONTINUE

        PTOTS= 1.
        RATAREA=AREAEN/ATHR
        CALL SUB (COST,TEMP,RATAREA,PTOTS)
        PASUB=TEMP**3.5
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        PASUP=TEMP**3.5
        AMEN2=5.*((1./PASUP)**(1./3.5)-1.)
        PINT=PASUP*(7.*AMEN2-1.)/6.
        XTRY=-1.E4
        IF(PEX.GT.PASUB) GO TO 311
        IF(PEX.LE.PASUB) GO TO 312
311     CONTINUE
        TEMP=PEX**(1./3.5)
        FLOWDUM=AREAEN*SQRT(7.)*TEMP**2.5*SQRT(1.-TEMP)
        AREF=FLOWDUM/((2./2.4)**2.5*(2.8/2.4)**.5)
        DO 305 N=2,NCM
        RATAREA=AREA(N)/AREF
        CALL SUB (COST,TEMP,RATAREA,PTOTS)
        TTH(N) = TEMP
        PTH(N) = PTOTS*(TEMP/1.)**3.5
        AMACHTH(N)=SQRT(5.*(1./TTH(N)-1.))
        ATH(N)=SQRT(1.4*TTH(N))
        UTH(N)=SQRT(7.*(1.-TTH(N)))
        RHOTH(N)=PTH(N)/TTH(N)
        STH(N) = 1.4*ALOG(TTH(N))-.4*ALOG(PTH(N))
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        HTH(N)  =GA*TTH(N)
        PTOTTH(N)  =PTH(N)*(1.+GD*AMACHTH(N)**2)**GA
        TTOTTH(N)  =TTH(N)*(1.+GD*AMACHTH(N)**2)
        FLOWTH(N)= RHOTH(N)*UTH(N)*AREA(N)
        FLHTTH(N)= UTH(N)*(PTH(N)+ETH(N))*AREA(N)
        R1TH(N)  = PTH(N)-RHOTH(N)*ATH(N)*UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = PTH(N)+RHOTH(N)*ATH(N)*UTH(N)
        R1TH(N)  = 5.*ATH(N)-UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = 5.*ATH(N)+UTH(N)
305     CONTINUE
        GO TO 350
312     CONTINUE
        DO 306 N=2,NCM
        RATAREA=AREA(N)/ATHR
        IF(X(N).LT.XTHR) CALL SUB  (COST,TEMP,RATAREA,PTOTS)
        IF(X(N).GE.XTHR) CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TTH(N)=TEMP
        PTH(N)=TTH(N)**3.5*PTOTS
        AMACHTH(N)=SQRT(5.*(1./TTH(N)-1.))
        ATH(N)=SQRT(1.4*TTH(N))
        UTH(N)=SQRT(7.*(1.-TTH(N)))
        RHOTH(N)=PTH(N)/TTH(N)
        STH(N) = 1.4*ALOG(TTH(N))-.4*ALOG(PTH(N))
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        HTH(N)  =GA*TTH(N)
        PTOTTH(N)  =PTH(N)*(1.+GD*AMACHTH(N)**2)**GA
        TTOTTH(N)  =TTH(N)*(1.+GD*AMACHTH(N)**2)
        FLOWTH(N)= RHOTH(N)*UTH(N)*AREA(N)
        FLHTTH(N)= UTH(N)*(PTH(N)+ETH(N))*AREA(N)
        R1TH(N)  = PTH(N)-RHOTH(N)*ATH(N)*UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = PTH(N)+RHOTH(N)*ATH(N)*UTH(N)
        R1TH(N)  = 5.*ATH(N)-UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = 5.*ATH(N)+UTH(N)
306     CONTINUE
        IF(PEX.LE.PINT) GO TO 350

        KIP=0
        XTRY=XTHR
        DXTRY=.1
360     CONTINUE
        KIP=KIP+1
        if(kip.gt.1000) WRITE(*,*)'pressione a valle troppo alta'
        if(kip.gt.1000) stop
        XTRY=XTRY+DXTRY
        AREADUM=CCC1*(XTRY+X000)+CCC2/(XTRY+X000)
        RATAREA=AREADUM/ATHR
        PTOTS=1.
        CALL SUPER(COST,TEMP,RATAREA,PTOTS)
        TDUM= TEMP
        TU  = TDUM
        AMADUM=SQRT(5.*(1./TDUM-1.))
        AM2=AMADUM**2
        CALL SHOCK(AM2,RATP,RATR,RATPT)
        PTOTS=RATPT
        RATAREA=AREAEN/ATHR
        CALL SUB(COST,TTRY,RATAREA,PTOTS)
        PTRY=PTOTS*(TTRY/1.)**3.5
C       WRITE(*,*)'KIP,XTRY,PTRY,PTOTS',KIP,XTRY,PTRY,PTOTS
        IF(ABS(PTRY-PEX)/PEX.LT.1.E-5) GO TO 361
        IF(PTRY.GT.PEX) GO TO 360
        XTRY=XTRY-DXTRY
        DXTRY=DXTRY/2.
        GO TO 360
361     CONTINUE
        DO 362 N=2,NCM
        IF(X(N).LT.XTRY) GO TO 362
        RATAREA=AREA(N)/ATHR
        CALL SUB(COST,TEMP,RATAREA,PTOTS)
        TTH(N) = TEMP
        PTH(N) = PTOTS*(TEMP/1.)**3.5
        AMACHTH(N)=SQRT(5.*(1./TTH(N)-1.))
        ATH(N)=SQRT(1.4*TTH(N))
        UTH(N)=SQRT(7.*(1.-TTH(N)))
        RHOTH(N)=PTH(N)/TTH(N)
        STH(N) = 1.4*ALOG(TTH(N))-.4*ALOG(PTH(N))
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
        HTH(N)  =GA*TTH(N)
        PTOTTH(N)  =PTH(N)*(1.+GD*AMACHTH(N)**2)**GA
        TTOTTH(N)  =TTH(N)*(1.+GD*AMACHTH(N)**2)
        FLOWTH(N)= RHOTH(N)*UTH(N)*AREA(N)
        FLHTTH(N)= UTH(N)*(PTH(N)+ETH(N))*AREA(N)
        R1TH(N)  = PTH(N)-RHOTH(N)*ATH(N)*UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = PTH(N)+RHOTH(N)*ATH(N)*UTH(N)
        R1TH(N)  = 5.*ATH(N)-UTH(N)
        R2TH(N)  = HTH(N)-PTH(N)/RHOTH(N)
        R3TH(N)  = 5.*ATH(N)+UTH(N)
362     CONTINUE
        IF(XTRY.GE.XTHR.AND.XTRY.LT.1.0) THEN
        XSHTH=XTRY
        PU=TU**3.5
        RHOU=TU**2.5
        TU=PU/RHOU
        HU=3.5*TU
        SU=.0
        UU=SQRT((1.-TU)*7.)
        AU=SQRT(1.4*TU)
        AMU=UU/AU
        TTOTU=TU+1./7.*UU**2
        HTOTU=3.5*TTOTU
        PTOTU=PU*(1.+.2*AMU**2)**3.5
        FLOWU=FLOWTH(2)
        FLHTU=FLHTTH(2)
        R1U  = PU-RHOU*AU*UU
        R2U  = HU-PU/RHOU
        R3U  = PU+RHOU*AU*UU
        AMU2=AMU**2
        CALL SHOCK(AMU2,RATP,RATR,RATPT)
        PD=PU*RATP
        RHOD=RHOU*RATR
        TD=PD/RHOD
        HD=3.5*TD
        SD=-.4*ALOG(RATPT)
        UD=SQRT((1.-TD)*7.)
        AD=SQRT(1.4*TD)
        AMD=UD/AD
        TTOTD=TD+1./7.*UD**2
        HTOTD=3.5*TTOTD
        PTOTD=PD*(1.+.2*AMD**2)**3.5
        FLOWD=FLOWTH(2)
        FLHTD=FLHTTH(2)
        R1D  = PD-RHOD*AD*UD
        R2D  = HD-PD/RHOD
        R3D  = PD+RHOD*AD*UD
        ENDIF

350     CONTINUE
        RETURN

5555    CONTINUE

        RATCRI=(1./GE)**(1./GI)
        AA=GF
        UA=0.0
        PA=1.0
        SA=0.0
        AC=GF*PEX**GI
        IF(PEX.LE.RATCRI) AC=GF/GE
        UC=GG*(GF-AC)
        PC=PA*(AC/AA)**(1./GI)
        SC=0.0
        ALAMA=UA-AA
        ALAMC=UC-AC
        X1   = 1.0 + ALAMA*TIME
        X2   = 1.0 + ALAMC*TIME

        WRITE(*,*)'X1 X2   ',X1,X2
        WRITE(*,*)'PA PC   ',PA,PC

        DO 51 N=2,NC

        IF(X(N).LE.X1) THEN
        ATH(N) = AA
        UTH(N) = 0.0
        STH(N) = 0.0
        PTH(N) = 1.0
        ENDIF

        IF(X(N).GE.X1.AND.X(N).LE.X2) THEN
        ALAM   = (X(N)-1.0)/TIME
        RIE    = GG*AA
        ATH(N) = (RIE-ALAM)/GC
        UTH(N) = RIE-GG*ATH(N)
        AMACHTH(N)=UTH(N)/ATH(N)
        STH(N) = 0.0
        PTH(N) = (ATH(N)/AA)**(1./GI)
        ENDIF

        IF(X(N).GE.X2) THEN
        ATH(N) = AC
        UTH(N) = UC
        AMACHTH(N)=UTH(N)/ATH(N)
        STH(N) = 0.0
        PTH(N) = (AC/AA)**(1./GI)
        ENDIF

        TTH(N)  =ATH(N)**2/GAMMA
        RHOTH(N)=PTH(N)/TTH(N)
        ETH(N)  =GB*PTH(N)+.5*RHOTH(N)*UTH(N)**2
51      CONTINUE

        P1L=PA
        P1R=PA
        U1L=UA
        U1R=UA
        S1L=SA
        S1R=SA
        A1L=AA
        A1R=AA
        T1L=A1L**2/GAMMA
        T1R=A1R**2/GAMMA
        RHO1L=P1L/T1L
        RHO1R=P1R/T1R
        E1L=GB*P1L+.5*RHO1L*U1L**2
        E1R=GB*P1R+.5*RHO1R*U1R**2

        P2L=PC
        P2R=PC
        U2L=UC
        U2R=UC
        S2L=SC
        S2R=SC
        A2L=AC
        A2R=AC
        T2L=A2L**2/GAMMA
        T2R=A2R**2/GAMMA
        RHO2L=P2L/T2L
        RHO2R=P2R/T2R
        E2L=GB*P2L+.5*RHO2L*U2L**2
        E2R=GB*P2R+.5*RHO2R*U2R**2

        RETURN

6666    CONTINUE
        A1   = SQRT(1.4*P1/RHO1)
        A2   = SQRT(1.4*P2/RHO2)
        PDUM = AMAX1(P1,P2)
        DELP = 0.5
        P3   = PDUM
        P4   = PDUM
        KIP  = 0
61      CONTINUE
        KIP = KIP+1
        IF(KIP.EQ.500) THEN
        WRITE(*,*)'ATTENZIONE KIP=500 !!!!'
        STOP
        ENDIF
        AM3 =-SQRT((6.*P4/P2+1.)/7.)
        WW2  = AM3*A2
        WW4  = WW2*(5.+AM3**2)/(6.*AM3**2)
        VSH3= U2-WW2
        RHO4=RHO2*WW2/WW4
        U4  = VSH3+WW4
        AM4 = SQRT((6.*P3/P1+1.)/7.)
        WW1  = AM4*A1
        WW3  = WW1*(5.+AM4**2)/(6.*AM4**2)
        VSH4= U1-WW1
        RHO3=RHO1*WW1/WW3
        U3  = VSH4+WW3
        VSH5= U3
        DIFF = U3-U4
        IF(ABS(DIFF).LT.1.E-5) GO TO 62
        IF(DIFF.LE.0.0) THEN
        P3  = P3-DELP
        DELP=DELP/2.
        ENDIF
        P3  = P3+DELP
        P4  = P3
        GO TO 61
62      CONTINUE
        RETURN

7777    CONTINUE

        A1   = SQRT(1.4*P1/RHO1)
        A2   = SQRT(1.4*P2/RHO2)
        PDUM = P1
        DELP = 0.5
        P3   = PDUM
        P4   = PDUM
        KIP  = 0
71      CONTINUE
        KIP = KIP+1
        IF(KIP.EQ.500) THEN
        WRITE(*,*)'ATTENZIONE KIP=500 !!!!'
        STOP
        ENDIF
        AM3 =-SQRT((6.*P4/1.0+1.)/7.)
        WW0  = AM3*GF
        WW4  = WW0*(5.+AM3**2)/(6.*AM3**2)
        VSH3= 0.0-WW0
        RHO4=1.0*WW0/WW4
        U4  = VSH3+WW4

        A3  = A2*(P3/P2)**GI
        U3  = U2+GG*(A2-A3)
        RHO3= RHO2*(P3/P2)**(1./GAMMA)
        VSH4L= U2-A2
        VSH4R= U3-A3
        VSH5 = U3

        DIFF = U3-U4
        IF(ABS(DIFF).LT.1.E-5) GO TO 72
        IF(DIFF.LE.0.0) THEN
        P3  = P3-DELP
        DELP=DELP/2.
        ENDIF
        P3  = P3+DELP
        P4  = P3
        GO TO 71
72      CONTINUE
        RETURN

8888    CONTINUE
        IF(P2/RHO2.GT.1) THEN
        A1   = SQRT(1.4*P1/RHO1)
        A2   = SQRT(1.4*P2/RHO2)
        PDUM = P1
        DELP = 0.5
        P3   = PDUM
        P4   = PDUM
        KIP  = 0
81      CONTINUE
        KIP = KIP+1
        IF(KIP.EQ.500) THEN
        WRITE(*,*)'ATTENZIONE KIP=500 !!!!'
        STOP
        ENDIF
        AM4 =-SQRT((6.*P4/1.0+1.)/7.)
        WW0  = AM4*GF*SQRT(P2/RHO2)
        WW4  = WW0*(5.+AM4**2)/(6.*AM4**2)
        VSH4= 0.0-WW0
        RHO4= RHO2*WW0/WW4
        U4  = VSH4+WW4

        A3  = A1*(P3/P1)**GI
        U3  = U1+GG*(A1-A3)
        RHO3= RHO1*(P3/P1)**(1./GAMMA)
        VSH3L= U1-A1
        VSH3R= U3-A3
        VSH5 = U3

        DIFF = U3-U4
        IF(ABS(DIFF).LT.1.E-5) GO TO 82
        IF(DIFF.GT.0.0) THEN
        P3  = P3+DELP
        DELP=DELP/2.
        ENDIF
        P3  = P3-DELP
        P4  = P3
        GO TO 81
82      CONTINUE
        RETURN
        ENDIF

        IF(P2/RHO2.LE.1) THEN
        A1   = SQRT(1.4*P1/RHO1)
        A2   = SQRT(1.4*P2/RHO2)
        PDUM = P1
        DELP = 0.5
        P3   = PDUM
        P4   = PDUM
        KIP  = 0
83      CONTINUE
        KIP = KIP+1
        IF(KIP.EQ.500) THEN
        WRITE(*,*)'ATTENZIONE KIP=500 !!!!'
        STOP
        ENDIF
        AM4 =-SQRT((6.*P4/1.0+1.)/7.)
        WW0  = AM4*GF*SQRT(P2/RHO2)
        WW4  = WW0*(5.+AM4**2)/(6.*AM4**2)
        VSH4= 0.0-WW0
        RHO4= RHO2*WW0/WW4
        U4  = VSH4+WW4

        AM3 = SQRT((6.*P3/P1+1.)/7.)
        WW1  = AM3*GF*SQRT(P1/RHO1)
        WW3  = WW1*(5.+AM3**2)/(6.*AM3**2)
        VSH3= U1-WW1
        RHO3= RHO1*WW1/WW3
        U3  = VSH3+WW3

        VSH5 = U4

        DIFF = U3-U4
        IF(ABS(DIFF).LT.1.E-5) GO TO 84
        IF(DIFF.LT.0.0) THEN
        P3  = P3-DELP
        DELP=DELP/2.
        ENDIF
        P3  = P3+DELP
        P4  = P3
        GO TO 83
84      CONTINUE
        RETURN
        ENDIF

        END

C................................................................

        SUBROUTINE CLEAN

        INCLUDE 'COMST05.INC'

        IF(ITEST.EQ.6) GO TO 6666
        IF(ITEST.EQ.7) GO TO 7777
        IF(ITEST.EQ.8) GO TO 8888
        RETURN
6666    CONTINUE

        NMID =NA/2+2
        NDIFF=10000
        DIFFX=0.0
        DO 1 N=3,NMID
        NP=N+1
        NM=N-1
        DIFF=ABS(P(NP)-P(NM))
        IF(DIFF.GT.DIFFX) THEN
        NDIFF=N
        DIFFX=DIFF
        ENDIF
1       CONTINUE

        DO 2 N=2,NMID

        IF(N.LT.(NDIFF-6)) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        ENDIF

        IF(N.GT.(NDIFF+6)) THEN
        P(N)  = 1.0
        RHO(N)= 1.0
        U(N)  = 0.0
        ENDIF

        IF(N.LE.(NDIFF+6).AND.N.GE.(NDIFF-6)) GO TO 2

        T(N)    = P(N)/RHO(N)
        A(N)    = SQRT(1.4*T(N))
        AMACH(N)= U(N)/A(N)
        S(N)    = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)    = GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)    = GA*T(N)
        PTOT(N) = P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N) = T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
2       CONTINUE

        NMID =NA/2+2
        NDIFF=10000
        DIFFX=0.0
        DO 3 N=NMID,NCM-1
        NP=N+1
        NM=N-1
        DIFF=ABS(P(NP)-P(NM))
        IF(DIFF.GT.DIFFX) THEN
        NDIFF=N
        DIFFX=DIFF
        ENDIF
3       CONTINUE

        DO 4 N=NMID,NCM

        IF(N.GT.(NDIFF+6)) THEN
        P(N)  = P2
        RHO(N)= RHO2
        U(N)  = U2
        ENDIF

        IF(N.LT.(NDIFF-6)) THEN
        P(N)  = 1.0
        RHO(N)= 1.0
        U(N)  = 0.0
        ENDIF

        IF(N.LE.(NDIFF+6).AND.N.GE.(NDIFF-6)) GO TO 4

        T(N)    = P(N)/RHO(N)
        A(N)=SQRT(1.4*T(N))
        AMACH(N)=U(N)/A(N)
        S(N) = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)  =GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)  =GA*T(N)
        PTOT(N)  =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N)  =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
4       CONTINUE
        RETURN

7777    CONTINUE

        NDIFF2 = 10000
        NDIFF1 = 10000
        DIFFX2 = 0.0
        DIFFX1 = 0.0
        DO 5 N=3,NCM-1
        NP=N+1
        NM=N-1
        DIFF=ABS(P(NP)-P(NM))

        IF(P(NP).GT.0.9*P1) THEN
        IF(DIFF.GT.DIFFX2)  THEN
        NDIFF2 = N
        DIFFX2 = DIFF
        ENDIF
        ENDIF

        IF(P(NM).LT.1.1*P1) THEN
        IF(DIFF.GT.DIFFX1)  THEN
        NDIFF1 = N
        DIFFX1 = DIFF
        ENDIF
        ENDIF

5       CONTINUE
        WRITE(*,*)'N1,N2',NDIFF1,NDIFF2
        DO 6 N=2,NCM-1

        IF(N.LT.(NDIFF2-6)) THEN
        P(N)  = P2
        RHO(N)= RHO2
        U(N)  = U2
        ENDIF

        IF(N.GT.(NDIFF2+6).AND.N.LT.(NDIFF1-6)) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        ENDIF

        IF(N.GT.(NDIFF1+6)) THEN
        P(N)  = 1.0
        RHO(N)= 1.0
        U(N)  = 0.0
        ENDIF

        IF(N.GE.(NDIFF2-6).AND.N.LE.(NDIFF2+6)) GO TO 6
        IF(N.GE.(NDIFF1-6).AND.N.LE.(NDIFF1+6)) GO TO 6

        T(N)    = P(N)/RHO(N)
        A(N)    = SQRT(1.4*T(N))
        AMACH(N)= U(N)/A(N)
        S(N)    = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)    = GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)    = GA*T(N)
        PTOT(N) =P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N) =T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
6       CONTINUE

        RETURN

8888    CONTINUE

        NMID =NA/2+2
        NDIFF=10000
        DIFFX=0.0
        DO 7 N=3,NMID
        NP=N+1
        NM=N-1
        DIFF=ABS(P(NP)-P(NM))
        IF(DIFF.GT.DIFFX) THEN
        NDIFF=N
        DIFFX=DIFF
        ENDIF
7       CONTINUE

        DO 8 N=2,NMID
        IF(N.LT.(NDIFF-6)) THEN
        P(N)  = P1
        RHO(N)= RHO1
        U(N)  = U1
        T(N)    = P(N)/RHO(N)
        A(N)    = SQRT(1.4*T(N))
        AMACH(N)= U(N)/A(N)
        S(N)    = 1.4*ALOG(T(N))-.4*ALOG(P(N))
        E(N)    = GB*P(N)+.5*RHO(N)*U(N)**2
        H(N)    = GA*T(N)
        PTOT(N) = P(N)*(1.+GD*AMACH(N)**2)**GA
        TTOT(N) = T(N)*(1.+GD*AMACH(N)**2)
        W1(N) =RHO(N)
        W2(N) =W1(N)*U(N)
        W3(N) =E(N)
        F1(N) =RHO(N)*U(N)
        F2(N) =P(N)+RHO(N)*U(N)*U(N)
        F3(N) =U(N)*(P(N)+E(N))
        FLOW(N)= RHO(N)*U(N)*AREA(N)
        FLHT(N)= U(N)*(P(N)+E(N))*AREA(N)
        ENDIF
8       CONTINUE
        RETURN

        END

C................................................................

