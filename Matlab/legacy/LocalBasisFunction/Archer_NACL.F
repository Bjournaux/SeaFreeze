#include "fintrf.h"

C This is FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2004 The MathWorks, Inc.
C $Revision: 1.9.2.1 $
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
         MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE  MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
  
          MWPOINTER M,N
          MWPOINTER PTmp,outputp
          
C-----------------------------------------------------------------------
C
          REAL*8 PTm(1000000,3),output(1000000,26)

           M = MXGETM(PRHS(1))
           N = MXGETN(PRHS(1))	

      if(M .gt. 1000000) then
        CALL MEXERRMSGTXT('Too many input points')
      endif
   
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      PTmp = MXGETPR(PRHS(1))

C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
       CALL MXCOPYPTRTOREAL8(PTmp, PTm, M*3)
	  
C CREATE A MATRIX FOR RETURN ARGUMENT
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(M,26,0)
 
C
C CREATE A MATRIX FOR RETURN ARGUMENT
C	  
       outputp = MXGETPR(PLHS(1))

C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
		
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  

C
C           write(6,400) M
 400       format(i10)
	CALL NACL_AR(M,PTm, output)
	
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
      CALL MXCOPYREAL8TOPTR(output, outputp, M*26)

      RETURN
      END


             BLOCK DATA
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            DOUBLE PRECISION MR
	    COMMON /REFS/TR,PR,MR,SH2OTR
            DATA TR,PR,MR/298.15d0,.10d0,6.d0/
            DATA SH2OTR/3.8806d0/
            END



C************************************	
	SUBROUTINE NACL_AR(np,PTm, output)
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /CPCALC/CPSAVE
      COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
      COMMON /REFS2/GEXTR,SCORR,SIONS,GH2OTR,SEXMR
      COMMON /REFS/TR,PR,MR,SH2OTR
      common /betao/bflag
        DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,MR
        Double precision PTm(np,3),output(np,26)

        CHARACTER*3 ANSWER
        DATA IAC,IAK/1,1/
        DATA B0/ 1.697521279805020e-01,2.259112226362840e-01,
     &      -2.114825336892210e-01,0,0,5.115848035497095e-02,
     &      -5.162843180495660e+01,0,6.406857762719160e-01,
     &      -1.620243674058400e+01,7.226082973204031e+01,
     &      -1.079047241164010e+00,0,0,1.015067387407250e+01,
     &      -1.909948464417210e+01,1.117650484589690e+00,
     &      -1.084847743698650e+01,0,0,2.363319007710990e+01,
     &      2.599340467700060e+00,-7.799932657592330e-01,
     &       8.134724399295809e+00/
        DATA B1/-1.794146363585010e+00,5.232546877161670e+00,0,
     &    -2.065972870059340e+01,0,0,3.923324202537330e+02,
     &    -6.459643174357840e-01,0,0,2.506258068875210e+02,
     &    0,-1.015523154403010e+02,0,0,-1.043157965021230e+02,
     &    0,0,0,0,0,0,0,0/
        DATA C0/ 1.740257414460444e-02,-1.159586306353160e-01,
     &       4.904329105911764e-02, -2.547845918784475e-02,
     &      0,0,1.213819756896780e+00,-1.143607122592849e-02,
     &     0,0,-2.220096790979990e+00,0,1.511493962926520e+00,
     &      -1.533930307449620e+01,-1.162469494510210e-01,0,
     &       1.031292777379220e-01,0,0,0,0,0,0,-1.913087079350930e-01/
        DATA C1/ 1.233396463689420e+00, -5.123875496288060e+00,
     &      1.463311468560450e+00,-2.003181141293860e+01,0,
     &      2.808779877592540e-01,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &      0,0,0,0,8.511699632329510e+00/
        DATA XV/ 1.741489217010460e+00,9.037749396513380e-01,
     &     6.789684304075950e+00,5.461964868418160e+00,
     &    -1.304039116694210e+01,1.667078379917950e-01,
     &    -9.260812944537040e+00,0,1.323605828731950e-01,0/
        DATA XC/ 1.640786577438580e+00,0,-7.791800772907960e-01,
     &     2.150272556945700e-01,0,0, -1.019009655885870e+00,
     &    -9.045299999999999e+00,-8.958900000000000e+00,
     &     1.621430000000000e-01,1.152196000000000e-01/

            SIONS=XC(11)
            LSAV=5
C           write(*,400) np
C 400       format(i10)
C   !$OMP PARALLEL DO 
C   !$OMP& DEFAULT(SHARED) 
C   !$OMP& PRIVATE(t,p,m,ap,ah,av,ak,gam,phi,vphic,comphi,cpphi,ephic,rho,si,dpr,ddt,u,h,s,cvcp,a)

        do 10 i=1,np
          t=PTm(i,2)
          p=PTm(i,1)
          m=PTm(i,3)

          bflag=0.0

          CALL DHLLN(T,P,AP,AH,AV,AC,AK,IAC,IAK)
            ALPHAW=-DRHODT/RHO
            BETAW=DRHODP/RHO
            CPW=CPSAVE
            RHOW=RHO
            DELTAW=ALPHAW*ALPHAW/RHOW/CPW



            call GAMCAL(T,P,M,GAM,PHI)
            aw = EXP(-PHI*2.d0*M*18.0152d0/1000.d0)



           call LPHI(T,P,M,HEXT)
     
          call Vphi(t,p,m,VphiC,aV)
             VS=(VPHIC*M+1.D3/RHOW)/(1.D3+M*58.443d0)
             RHOS=1.D0/VS

           call Kphi(t,p,m,COMphi,aK)
              COMSOL=(M*COMPHI+1.D3*BETAW/RHOW)/(1.D3+M*58.443d0)
              COMSOL=COMSOL*RHOS
  
          CALL CPCALC3(T,P,M,CpPHI,ac)
               CPSOL=(CPPHI*M+1.D3*CPW)/(1.D3+M*58.443d0)

          call ephi(t,p,m,ephic,av)
             ESOL=(EPHIC*M+ALPHAW*1.D3/RHOW)/(1.D3+M*58.443d0)
             ESOL=ESOL*RHOS
c    CALCULATE DIFFERENCE IN THE TWO COMPRESSIBILITIES
             DELTSL=T*ESOL*ESOL/RHOS/CPSOL
             BETAIS=COMSOL-DELTSL
             CORRIS=DELTSL*(1.D0+ALPHAW/ESOL)*EPHIC/ESOL-DELTAW*CPPHI/
     1               (CPSOL*RHOS)+(DELTAW-(DELTSL*ALPHAW/ESOL))*VPHIC


C    DO FIRST PASS THROUGH TO GET REFERENCE VALUES.
             CALL WSTEAM(TR,RHO,PR,SI,DPR,DDT,U,H,S,CV,CP,A,4)

             SCORR=SH2OTR-S
             GH2OTR=H-TR*SH2OTR      
             CALL GAMCAL(TR,PR,MR,GAMMA,PHIC)
             GEXTR=2.D0*MR*8.3144d-3*TR*(LOG(GAMMA)+1.D0-PHIC)
             CALL LPHI(TR,PR,MR,HEXMR)
             SEXMR=(HEXMR-GEXTR/MR)/TR

      
                    CALL DELG20(T,P,DELG)   
                    call s20(t,p,st)

                    delh=delg+t*st-298.15d0*sions
                    dels=st-sions           

        output(i,1)=gam
        output(i,2)=phi
        output(i,3)=aw
        output(i,4)=vphic
        output(i,5)=ephic
        output(i,6)=comphi
        output(i,7)=hext
        output(i,8)=cpphi
        output(i,9)=rhos
        output(i,10)=rhow
        output(i,11)=esol
        output(i,12)=alphaw
        output(i,13)=comsol
        output(i,14)=betaw
        output(i,15)=cpsol
        output(i,16)=cpw
        output(i,17)=delg
        output(i,18)=delh
        output(i,19)=dels*1.d3
        output(i,20)=AP
        output(i,21)=Av
        output(i,22)=Ac
        output(i,23)=AH
        output(i,24)=D
        output(i,25)=Gw
        output(i,26)=Sw

 10    CONTINUE   
C  !$OMP END PARALLEL DO 

          END        
C**********************************************************************************

      SUBROUTINE CPCALC3(T,P,M,CPPHI,ac)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XVA(10),XCA(11)
      COMMON /CPCALC/CPSAVE
        DIMENSION XB(5,24),XV(10),XC(11)
      DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,lot,lot2,lot3,lot4
     1 ,LOTA  
       DATA R,NUM,NUX,NU,ZM,ZX/8.3144D-3,1.D0,1.D0,2.D0,1.D0,-1.D0/
       DATA PS,TS,MS/.10d0,298.15d0,6.d0/
       DATA ALPHA,ALPHA2,BB/2.0D0,2.5D0,1.2D0/
          AJT=AC*R
         CP=CPSAVE

      ISS=0.5D0*MS*(NUM*ZM*ZM + NUX*ZX*ZX)
      IS=(M*(NUM*ZM*ZM+NUX*ZX*ZX))/2.D0
      SQIS=SQRT(IS)
      AI=ALPHA*SQIS
      AI2=ALPHA2*SQIS
       SQISS=SQRT(ISS)
      AIS=ALPHA*SQISS
      A2IS=ALPHA2*SQISS 
      if(m.le.0.0d0)then
         ai=0.5d0
         ai2=0.5d0
      endif
      COR=(1.D0-(1.D0+AI)*EXP(-AI))/(AI*AI)
      COR2=(6.D0-(6.D0+AI2*(6.D0+AI2*(3.D0+AI2)))*EXP(-AI2))/
     1      ((AI2*AI2)*(AI2*AI2))
      CORR=(1.D0-(1.D0+AIS)*EXP(-AIS))/(AIS*AIS)
      CORR2=(6.D0-(6.D0+A2IS*(6.D0+A2IS*(3.D0+A2IS)))*EXP(-A2IS))/
     1      ((A2IS*A2IS)*(A2IS*A2IS))

      XC(1)=1.d0
      XC(2)=300.D0/T/T
      XC(3)=T/300.d0
      XC(4)=T*T/300.d0/300.d0
      XC(6)=T*T*T/3.d2/3.d2/3.d2
      XC(7)=1.d2/T
      XC(8)=0.d0
      XC(9)=0.d0
      XC(10)=0.d0
      XC(11)=0.d0


      DO 5 I=1,10
    5 XV(I)=0.D0
      DP=P/100.D0
      DT=T/300.D0
      XV(4)=-T*6.D0*DT*(P-ps)/9.D4/1.D3
      XV(7)=-T*2.D0*(p*p-ps*ps)/9.D4/2.D0/1.D3/100.d0
      XV(10)=0.D0
      

      XB(1,1)=-2.D0*NUM*NUX*R*T*T
      XB(2,1)=XB(1,1)*(2.D0*M*COR-2.D0*MS*CORR)
      XB(3,1)=XB(1,1)*(M*M-MS*MS)*NUM*ZM
      XB(4,1)=XB(1,1)*NUM*ZM*4.D0*(M*M*COR2-MS*MS*CORR2)
      XB(1,1)=XB(1,1)*(M-MS)

        TT=T/500.D0
        LOT=t-225.D0
c NC 10/13/2010: changed from 200.D0        
c	  LOTa=t-200.D0
	  LOTa=t-215.D0
        HIT=650.d0-T
        hit2=hit*hit
        hit4=hit2*hit2
        LOT2=LOT*LOT
        lot3=lot2*lot
         lot4=lot3*lot
C NEXT SET OF LINES GIVE dBETA/dT.
      DO 1 J = 1,4
         XB(J,2)=XB(J,1)/1.d+3
         XB(J,3)=XB(J,1)*TT*2.D0/500.D0
	 XB(J,4)=-XB(J,1)/LOTA/LOTA
         XB(J,5)=-3.d0*XB(J,1)*1.D+4/LOTA/LOTA/lota/lota
         XB(J,6)=-2.d0*XB(J,1)*1.d2/LOTA/LOTA/LOTA
         XB(J,7)=-2.d0*XB(J,1)*2.D2/T/T/T
         XB(J,8)=3.D0*xb(j,1)*TT*TT/500.D0
         XB(j,9)=0.5d0*xb(j,1)/(hit**1.5d0)
         XB(J,10)=0.0D0
         XB(J,11)=-XB(J,1)*P*1.D-6*2.D2/LOT2
         XB(J,12)=3.d0*XB(J,1)*P*1.d2/HIT4
         XB(J,13)=XB(J,1)*1.D-5*P/500.D0
         XB(J,14)=XB(J,1)*P*2.D2*1.D-6/HIT2
         XB(J,15)=0.0D0
         XB(J,16)=-XB(J,1)*P*P*1.D-8*2.D2/LOT2
         XB(J,17)=3.d0*XB(J,1)*P*P/HIT4
         XB(J,18)=XB(J,1)*P*P*1.D-7/500.D0
         XB(J,19)=2.d0*XB(J,1)*p*p*TT*1.d-7/500.d0
         XB(J,20)=-2.D0*XB(J,1)*1.d-6*4.d4*P/LOT/LOT/LOT
         XB(J,21)=2.D0*XB(J,1)*1.d-5*P*T/500.d0/5.D2
         XB(J,22)=-XB(J,1)*P*P*P*1.D-10*2.D2/LOT2
         XB(J,23)=xb(j,1)*p*p*p*3.d0*1.d-2/hit4
         XB(J,24)=3.D0*XB(J,1)*2.D2/HIT4
         XB(J,24)=3.D0*XB(J,1)*2.D2/HIT4
  1   CONTINUE

C  NEXT LOOP MULTIPLIES dBETA/dT BY 2/T.
      DO 2 J = 1, 4
        DO 2 JJ = 2, 24
  2   XB(J,JJ)=XB(J,JJ)*2.d0/T

C NEXT LOOP ADDS d2BETA/dT2.
           hit5=hit4*hit
           hit3=hit2*hit
      DO 3 J = 1, 4
         XB(J,2)=XB(J,2)+0.d0 
         XB(J,3)=XB(J,3)+XB(J,1)*2.D0/500.D0/500.D0
         XB(J,4)=XB(J,4)+XB(J,1)*2.D0/LOTA/LOTA/LOTA
         XB(J,5)=xb(j,5)+12.d0*XB(J,1)*1.D+4/LOTA/LOTA/lota/lota/lota
         XB(J,6)=XB(J,6)+6.d0*XB(J,1)*1.d2/LOTA/LOTA/LOTA/LOTA
         XB(J,7)=XB(J,7)+XB(J,1)*6.d0*2.D2/T/T/T/T
         XB(J,8)=XB(J,8)+6.d0*XB(J,1)*TT/500.D0/500.D0
         XB(j,9)=xb(j,9)+0.75d0*xb(j,1)/(hit**2.5d0)
         XB(J,10)=XB(J,10)+0.0D0
         XB(J,11)=XB(J,11)+2.D0*XB(J,1)*200.D-06*P/LOT3
         XB(J,12)=XB(J,12)+12.D0*XB(J,1)*1.d2*P/HIT5
         XB(J,13)=XB(J,13)+0.0D0
         XB(J,14)=XB(J,14)+2.D0*XB(J,1)*200.D-6*P/HIT3
         XB(J,15)=0.0D0
         XB(J,16)=XB(J,16)+2.0D0*XB(J,1)*P*P*200.D-8/LOT3
         XB(J,17)=XB(J,17)+12.D0*XB(J,1)*P*P/HIT5
         XB(J,18)=XB(J,18)+0.0D0
         XB(J,19)=XB(J,19)+2.d0*XB(J,1)*p*p*1.d-7/500.d0/500.D0
         XB(J,20)=XB(J,20)+6.D0*XB(J,1)*1.d-6*4.d4*P/LOT/LOT/LOT/LOT
         XB(J,21)=XB(J,21)+2.D0*XB(J,1)*1.d-5*P/500.d0/5.D2
         XB(J,22)=XB(J,22)+2.D0*XB(J,1)*P*P*P*200.D-10/LOT3
         XB(J,23)=xb(j,23)+xb(j,1)*12.d0*p*p*p*1.d-2/hit5
         XB(J,24)=XB(J,24)+12.D0*XB(J,1)*2.D2/HIT5
         XB(J,1)=0.d0
 3    CONTINUE  


         BETA0V=0.0D0
         BETA1V=0.0D0
         CPHIV=0.0D0
         CPHIV1=0.0D0
         DCP0DP=0.0D0
         CPADD=0.0D0
         DO 20 I=1, 24
           BETA0V=BETA0V+XB(1,I)*B0(I)
 20      CONTINUE
         DO 30 I=1, 24
           BETA1V=BETA1V+XB(2,I)*B1(I)
 30      CONTINUE
         DO 40 I=1, 24
           CPHIV=CPHIV+XB(3,I)*C0(I)
 40      CONTINUE
         DO 140 I=1, 24
           CPHIV1=CPHIV1+XB(4,I)*C1(I)
 140      CONTINUE
         DO 50 I=1,10
           DCP0DP=DCP0DP+XV(I)*XVA(I)
 50      CONTINUE
         DO 60 I=1,7
           CPADD=CPADD+XC(I)*XCA(I)
 60      CONTINUE
      DHL=ABS(ZM*ZX)*(NUM+NUX)*AJT*(LOG((1.d0+BB*SQIS)/
     1     (1.d0+BB*SQISS)))/(2.d0*BB)
      DHL=DHL+BETA0V+BETA1V+CPHIV+CPHIV1+DCP0DP+CPADD

       CPPHI=-1.d0*CP/MS+DHL
       CPPHI=CPPHI*1.d3

       RETURN
      END


C************************************************************************
      subroutine ephi(t,p,m,ephic,av)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
      COMMON /CPCALC/CPSAVE
        DIMENSION XB(5,24),XVA(10),XCA(11)
        DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,lot

       DATA R,NUM,NUX,NU,ZM,ZX/8.3144D-3,1.D0,1.D0,2.D0,1.D0,-1.D0/
       DATA PS,TS,MS/.1D0,298.15D0,6.D0/
       DATA ALPHA,ALPHA2,BB/2.D0,2.5D0,1.2D0/

c  Make Ae by means of numerical differentiation.
          CALL DHLLN(T-0.01D0,P,AP,AH,AV,AC,AK,IAC,IAK)
          DRDPL=DRHODP
          AVL=AV
          CALL DHLLN(T+0.01D0,P,AP,AH,AV,AC,AK,IAC,IAK)
          DRDPDT=(DRHODP-DRDPL)/0.02D0
          AVT=(AV-AVL)/0.02D0
          CALL DHLLN(T,P,AP,AH,AV,AC,AK,IAC,IAK)
 

      ISS=0.5D0*MS*(NUM*ZM*ZM + NUX*ZX*ZX)
      IS=(M*(NUM*ZM*ZM+NUX*ZX*ZX))/2.D0
      SQIS=SQRT(IS)
      AI=ALPHA*SQIS
      AI2=ALPHA2*SQIS
       SQISS=SQRT(ISS)
      AIS=ALPHA*SQISS
      A2IS=ALPHA2*SQISS 
      IF(AI.EQ.0.d0)THEN
         ai=0.5
         ai2=0.5
      endif
      COR=(1.D0-(1.D0+AI)*EXP(-AI))/(AI*AI)
      CORR=(1.D0-(1.D0+AIS)*EXP(-AIS))/(AIS*AIS)

      DP=P/100.D0
      DT=T/300.D0
      XVA(1)=0.0D0
      XVA(2)=1.D-2/300.D0
      XVA(3)=2.D0*dt/1.d3/300.D0
      XVA(4)=3.D0*DT*DT/1.D3/3.D2
      XVA(5)=0.0D0
      XVA(6)=DP/3.D2/1.D1
      XVA(7)=2.D0*DP*DT/1.D3/3.D2
      xvA(8)=3.D0*DT*DT*DP/1.D4/3.D2
      XVA(9)=DP*DP/1.D3/3.D2
      XVA(10)=0.0D0
      xva(3)=(P+10.D0)**1.5D0/3.D2/1.D7
        DO 1000 IX=1,10
 1000           XVA(IX)=XVA(IX)*XV(IX)
C*****

C*****
      XB(1,1)=2.D0*NUM*NUX*R*T
      XB(2,1)=XB(1,1)*(2.D0*M*COR-2.D0*MS*CORR)
      XB(3,1)=XB(1,1)*(M*M-MS*MS)
      XB(4,1)=XB(1,1)*(M*M*(6.D0+(-AI2*AI2*AI2-3.D0*AI2*AI2
     1         -6.D0*AI2-6.D0)*EXP(-AI2))/(AI2**4.D0)-
     1             MS*MS*(6.D0+(-A2IS*A2IS*A2IS-3.D0*A2IS*A2IS
     1         -6.D0*A2IS-6.D0)*EXP(-A2IS))/(A2IS**4.D0))
      XB(1,1)=XB(1,1)*(M-MS)

        TT=T/500.D0
        LOT=t-225.D0
        HIT=650.d0-T
 


       DO 1 J = 1,4
         XB(J,2)=0.D0
         XB(J,3)=0.D0
         XB(J,4)=0.D0
         XB(J,5)=0.D0
	 XB(J,6)=0.D0
         XB(J,7)=0.D0
         XB(J,8)=0.D0
         XB(J,9)=0.D0
         XB(J,10)=XB(J,1)*1.d-5
         XB(J,11)=XB(J,1)*2.D2*1.d-6/LOT
         XB(J,12)=XB(J,1)*1.D2/HIT/hit/hit
         XB(J,13)=XB(J,1)*TT*1.d-5
         XB(J,14)=XB(J,1)*2.d2*1.d-6/hit
         XB(J,15)=2.D0*XB(J,1)*P*1.d-7
         XB(J,16)=2.D0*XB(J,1)*P*2.d2*1.d-8/LOT
         XB(J,17)=2.D0*XB(J,1)*P/HIT/hit/hit
         XB(J,18)=2.D0*XB(J,1)*P*TT*1.d-7
         XB(J,19)=2.d0*XB(J,1)*p*TT*TT*1.d-7
         XB(J,20)=XB(J,1)*1.d-6*4.d4/LOT/LOT
         XB(J,21)=XB(J,1)*1.d-5*T*T/500.d0/5.D2
         XB(J,22)=3.D0*XB(J,1)*P*P*1.d-10*200.D0/LOT
         XB(J,23)=3.d0*xb(j,1)*p*p*1.D-2/hit/hit/hit
         XB(J,24)=0.0D0
  1   CONTINUE
          DO 2 J=1,4
            DO 2 JJ=2,24
               XB(J,JJ)=XB(J,JJ)/T
 2     CONTINUE
       DO 3 J = 1,4
         XB(J,10)=XB(J,10)+0.0D0
         XB(J,11)=XB(J,11)-XB(J,1)*2.D2*1.d-6/LOT/lot
         XB(J,12)=XB(J,12)+3.d0*XB(J,1)*1.D2/HIT/hit/hit/hit
         XB(J,13)=XB(J,13)+XB(J,1)*1.d-5/5.d2
         XB(J,14)=XB(J,14)+XB(J,1)*2.d2*1.d-6/hit/hit
         XB(J,15)=XB(J,15)+0.0d0
         XB(J,16)=XB(J,16)-2.D0*XB(J,1)*P*2.d2*1.d-8/LOT/lot
         XB(J,17)=XB(J,17)+6.D0*XB(J,1)*P/HIT/hit/hit/hit
         XB(J,18)=XB(J,18)+2.D0*XB(J,1)*P*1.d-7/5.d2
         XB(J,19)=XB(J,19)+4.d0*XB(J,1)*P*TT*1.d-7/5.D2
         XB(J,20)=XB(J,20)-2.d0*XB(J,1)*1.d-6*4.d4/LOT/LOT/lot
         XB(J,21)=XB(J,21)+2.D0*XB(J,1)*1.d-5*T/500.d0/5.D2
         XB(J,22)=XB(J,22)-3.D0*XB(J,1)*P*P*1.d-10*200.D0/LOT/lot
         XB(J,23)=XB(J,23)+9.d0*xb(j,1)*p*p*1.D-2/hit/hit/hit/hit
         XB(J,24)=XB(J,24)+0.0D0
C         XB(J,1)=0.D0
  3   CONTINUE



         BETA0V=0.D0
         BETA1V=0.D0
         CPHIV1=0.D0
         CPHIV2=0.D0
         VMS=0.D0


         DO 10 IX=10,24
          BETA0V=BETA0V+XB(1,IX)*B0(IX)
          BETA1V=BETA1V+XB(2,IX)*B1(IX)
          CPHIV1=CPHIV1+XB(3,IX)*C0(IX)
          CPHIV2=CPHIV2+XB(4,IX)*C1(IX)
 10      CONTINUE
          DO 202 IX=1,10
            VMS=VMS+XVA(IX)
 202       CONTINUE

      DHL=ABS(ZM*ZX)*(NUM+NUX)*AVT*(DLOG(1.D0+BB*SQIS)
     1   -DLOG(1.D0+BB*SQISS))/(2.D0*BB)/1.D3

       DHL=DHL+BETA0V+BETA1V+CPHIV1+CPHIV2+VMS
       EPHIC=(1.0D0*DRHODT/(MS*RHO*RHO)+DHL)*1.D3

        RETURN
        END

C************************************************************************

      subroutine Kphi(t,p,m,COMphi,aK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
      COMMON /CPCALC/CPSAVE
        DIMENSION XB(5,24),XVA(10),XCA(11)
        DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,lot

       DATA R,NUM,NUX,NU,ZM,ZX/8.3144D-3,1.D0,1.D0,2.D0,1.D0,-1.D0/
       DATA PS,TS,MS/.1D0,298.15D0,6.D0/
       DATA ALPHA,ALPHA2,BB/2.D0,2.5D0,1.2D0/

          AVT=AK

      ISS=0.5D0*MS*(NUM*ZM*ZM + NUX*ZX*ZX)
      IS=(M*(NUM*ZM*ZM+NUX*ZX*ZX))/2.D0
      SQIS=SQRT(IS)
      AI=ALPHA*SQIS
      AI2=ALPHA2*SQIS
      SQISS=SQRT(ISS)
      AIS=ALPHA*SQISS
      A2IS=ALPHA2*SQISS 
      IF(AI.LE.0.d0)THEN
         ai=0.5
         ai2=0.5
      endif 
      COR=(1.D0-(1.D0+AI)*EXP(-AI))/(AI*AI)
      CORR=(1.D0-(1.D0+AIS)*EXP(-AIS))/(AIS*AIS)


      DP=P/1.D2
      DT=T/3.D2
      XVA(1)=0.0D0
      XVA(2)=0.0D0
      XVA(3)=0.0D0
      XVA(4)=0.0D0
      XVA(5)=-1.D0/1.D3/1.D2
      XVA(6)=-DT/1.D2/1.D1
      XVA(7)=-DT*DT/1.D3/1.D2
      xvA(8)=-DT*DT*DT/1.D4/1.D2
      XVA(9)=-2.D0*DP*DT/1.D3/1.D2
      XVA(10)=-1.5D0*SQRT(P+10.D0)/1.d6
      xva(3)=-1.5d0*sqrt(p+10.D0)*dt/1.d7

C*****
      XB(1,1)=2.D0*NUM*NUX*R*T
      XB(2,1)=XB(1,1)*(2.D0*M*COR-2.D0*MS*CORR)
      XB(3,1)=XB(1,1)*(M*M-MS*MS)
      XB(4,1)=XB(1,1)*(M*M*(6.D0+(-AI2*AI2*AI2-3.D0*AI2*AI2
     1         -6.D0*AI2-6.D0)*EXP(-AI2))/(AI2**4.D0)-
     1             MS*MS*(6.D0+(-A2IS*A2IS*A2IS-3.D0*A2IS*A2IS
     1         -6.D0*A2IS-6.D0)*EXP(-A2IS))/(A2IS**4.D0))
      XB(1,1)=XB(1,1)*(M-MS)

        TT=T/500.D0
        LOT=t-225.D0
        HIT=650.d0-T
 
       DO 1 J = 1,4
         XB(J,2)=0.D0
         XB(J,3)=0.D0
         XB(J,4)=0.D0
         XB(J,5)=0.D0
	 XB(J,6)=0.D0
         XB(J,7)=0.D0
         XB(J,8)=0.D0
         XB(J,9)=0.D0
         XB(J,10)=0.D0
         XB(J,11)=0.D0
         XB(J,12)=0.D0
         XB(J,13)=0.D0
         XB(J,14)=0.D0
         XB(J,15)=-2.D0*XB(J,1)*1.d-7
         XB(J,16)=-2.D0*XB(J,1)*2.d2*1.d-8/LOT
         XB(J,17)=-2.D0*XB(J,1)/HIT/hit/hit
         XB(J,18)=-2.D0*XB(J,1)*TT*1.d-7
         XB(J,19)=-2.d0*XB(J,1)*TT*TT*1.d-7
         XB(J,20)=0.0D0
         XB(J,21)=0.0D0
         XB(J,22)=-6.D0*XB(J,1)*P*1.d-10*200.D0/LOT
         XB(J,23)=-6.d0*xb(j,1)*1.D-2*p/hit/hit/hit
         XB(J,24)=0.0D0
         XB(J,1)=0.D0
  1   CONTINUE

         BETA0V=0.D0
         BETA1V=0.D0
         CPHIV1=0.D0
         CPHIV2=0.D0
         VMS=0.D0


         DO 10 IX=15,24
          BETA0V=BETA0V+XB(1,IX)*B0(IX)
          BETA1V=BETA1V+XB(2,IX)*B1(IX)
          CPHIV1=CPHIV1+XB(3,IX)*C0(IX)
          CPHIV2=CPHIV2+XB(4,IX)*C1(IX)
 10      CONTINUE
          DO 202 IX=1,10
            VMS=VMS+XVA(IX)*XV(IX)
 202       CONTINUE

      DHL=ABS(ZM*ZX)*(NUM+NUX)*AVT*(DLOG(1.D0+BB*SQIS)
     1 -DLOG(1.D0+BB*SQISS))/(2.D0*BB)/1.D3
       DHL=-DHL+BETA0V+BETA1V+CPHIV1+CPHIV2+VMS

      YY=-1.D0*DRHODP/(MS*RHO*RHO)+DHL
       COMPHI=YY*1000.D0

        RETURN
         END


C***************************************************************************
      subroutine Vphi(t,p,m,Vphic,avT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XVA(10),XCA(11)
      COMMON /CPCALC/CPSAVE
        DIMENSION XB(5,24),XV(10),XC(11)
        DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,lot

       DATA R,NUM,NUX,NU,ZM,ZX/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       DATA PS,TS,MS/.10d0,298.15d0,6.d0/
       DATA ALPHA,ALPHA2,BB/2.d0,2.5d0,1.2d0/

      ISS=0.5d0*MS*(NUM*ZM*ZM + NUX*ZX*ZX)
      IS=(M*(NUM*ZM*ZM+NUX*ZX*ZX))/2.d0
      SQIS=SQRT(IS)
      AI=ALPHA*SQIS
      AI2=ALPHA2*SQIS
       SQISS=SQRT(ISS)
      AIS=ALPHA*SQISS
      A2IS=ALPHA2*SQISS 
      if(ai.le.0.d0)then
         ai=0.5
         ai2=0.5
      endif
      COR=(1.d0-(1.d0+AI)*EXP(-AI))/(AI*AI)
      CORR=(1.d0-(1.d0+AIS)*EXP(-AIS))/(AIS*AIS)

      DP=P/100.d0
      DT=T/300.d0

      XV(1)=1.D-1
      XV(2)=DT/1.D2
      XV(3)=dt*dt/1.d3
      XV(4)=DT*DT*Dt/1.D3
      XV(5)=DP/1.D3
      XV(6)=DP*DT/1.d1
      XV(7)=DP*DT*DT/1.D3
      xv(8)=DT*DT*DT*DP/1.D4
      XV(9)=DP*DP*DT/1.D3
      XV(10)=(P+10.d0)**1.5D0/1.d6
      XV(3)=(P+10.d0)**1.5D0*DT/1.D7

C*****
      XB(1,1)=2.d0*NUM*NUX*R*T
      XB(2,1)=XB(1,1)*(2.d0*M*COR-2.d0*MS*CORR)
      XB(3,1)=XB(1,1)*(M*M-MS*MS)
      XB(4,1)=XB(1,1)*(M*M*(6.D0+(-AI2*AI2*AI2-3.D0*AI2*AI2
     1         -6.D0*AI2-6.D0)*EXP(-AI2))/(AI2**4.D0)-
     1             MS*MS*(6.D0+(-A2IS*A2IS*A2IS-3.D0*A2IS*A2IS
     1         -6.D0*A2IS-6.D0)*EXP(-A2IS))/(A2IS**4.D0))
      XB(1,1)=XB(1,1)*(M-MS)


        TT=T/500.d0
        LOT=T-225.d0
        HIT=650.d0-T
 
       DO 1 J = 1,4
         XB(J,2)=0.d0
         XB(J,3)=0.d0
         XB(J,4)=0.d0
         XB(J,5)=0.d0
	 XB(J,6)=0.d0
         XB(J,7)=0.d0
         XB(J,8)=0.0D0
         XB(J,9)=0.0D0
         XB(J,10)=XB(J,1)*1.d-5
         XB(J,11)=XB(J,1)*2.D2*1.d-6/LOT
         XB(J,12)=XB(J,1)*1.d2/HIT/hit/hit
         XB(J,13)=XB(J,1)*TT*1.d-5
         XB(J,14)=XB(J,1)*2.d2*1.d-6/hit
         XB(J,15)=2.D0*XB(J,1)*P*1.d-7
         XB(J,16)=2.D0*XB(J,1)*P*2.d2*1.d-8/LOT
         XB(J,17)=2.D0*XB(J,1)*P/HIT/hit/hit
         XB(J,18)=2.D0*XB(J,1)*P*TT*1.d-7
         XB(J,19)=2.d0*XB(J,1)*p*TT*TT*1.d-7
         XB(J,20)=XB(J,1)*1.d-6*4.d4/LOT/LOT
         XB(J,21)=XB(J,1)*1.d-5*T*T/500.d0/5.D2
         XB(J,22)=3.D0*XB(J,1)*P*P*1.d-10*200.D0/LOT
         XB(J,23)=3.d0*xb(j,1)*p*p*1.d-2/hit/hit/hit
         XB(J,24)=0.0D0
         XB(J,1)=0.d0
  1   CONTINUE

         BETA0V=0.d0
         BETA1V=0.d0
         CPHIV1=0.d0
         CPHIV2=0.d0
         VMS=0.d0


         DO 10 IX=10,24
          BETA0V=BETA0V+XB(1,IX)*B0(IX)
          BETA1V=BETA1V+XB(2,IX)*B1(IX)
          CPHIV1=CPHIV1+XB(3,IX)*C0(IX)
          CPHIV2=CPHIV2+XB(4,IX)*C1(IX)
 10      CONTINUE
          DO 202 IX=1,10
            VMS=VMS+XVA(IX)*XV(IX)
 202       CONTINUE

      DHL=ABS(ZM*ZX)*(NUM+NUX)*AVT*(DLOG(1.d0+BB*SQIS)
     1   -DLOG(1.d0+BB*SQISS))/(2.d0*BB)/1.d3

       DHL=DHL+BETA0V+BETA1V+CPHIV1+CPHIV2+VMS

      YY=-1.d0/(MS*RHO)+DHL
       VPHIC=YY*1000.D0

       RETURN
        END


      SUBROUTINE GAMCAL(T,P,M,GAM,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PARAMS/ B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
      COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT
     1 ,D2RDP2,D2RDT2,Gw,Sw
      COMMON /CPCALC/CPSAVE
            common /betao/bflag
        DIMENSION XB(5,24)
        DOUBLE PRECISION NUM,NUX,NU,MS,IS,M,ISS,lot,LOTA
        CHARACTER*3 ANSWER

       DATA R,NUM,NUX,NU,ZM,ZX/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       DATA PS,TS,MS/.1d0,298.15d0,6.d0/
       DATA IAC,IAK/0,0/
       DATA ALPHA,BB/2.0d0,1.2d0/

         if(m.le.0.0d0)then
               gam=1.d0
               phi=1.0d0
               return
         endif

C**************************************************************************
           CALL DHLLN(T,P,AP,AH,AV,AC,AK,IAC,IAK)
           APHIT=AP
           IS=0.5d0*M*(NUM*ZM*ZM + NUX*ZX*ZX)
           AI2=2.5D0*SQRT(M)           
           AIS=2.0D0*SQRT(M)
           AI2TRM=6.D0-(6.D0+AI2*(6.D0+AI2*(3.D0+AI2*(1.D0-0.5D0*AI2)))
     1         )*EXP(-AI2)
              AI2TRM=AI2TRM/(AI2*AI2*AI2*AI2)
      
        XB(1,1)=M*4.d0*NUM*NUX/NU
        XB(2,1)=XB(1,1)*
     1     (1.d0-(1.d0+AIS-AIS*AIS*0.5d0)*EXP(-AIS))/(AIS*AIS)
        XB(3,1)=XB(1,1)*1.5D0*NUM*ZM*M
        XB(4,1)=XB(1,1)*2.D0*NUM*ZM*M*AI2TRM

        T2=T*T/500.d0/500.d0
        LOT=T-225.d0
c       NC 10/13/2010 changed from 200 to 215 
c        LOTA=T-200.d0
        LOTA=T-215.d0
        HIT=650.d0-T
        hit3=hit*hit*hit 
       DO 1 J = 1,4
         XB(J,2)=XB(J,1)*T/1.D+3
         XB(J,3)=XB(J,1)*T2
	 XB(J,4)=XB(J,1)/LOTA
         XB(J,5)=XB(J,1)*1.D+4/LOTA/LOTA/lota
         XB(J,6)=XB(J,1)*1.d2/LOTA/LOTA
         XB(J,7)=XB(J,1)*2.d2/T/T
         XB(J,8)=XB(J,1)*T2*T/500.D0
         XB(J,9)=XB(J,1)/sqrt(HIT)
         XB(J,9)=XB(J,1)/sqrt(HIT)
         XB(J,10)=XB(J,1)*1.d-5*P
         XB(J,11)=XB(J,1)*1.d-6*P*200.D0/LOT
         XB(J,12)=XB(J,1)*P*1.d2/HIT3
         XB(J,13)=XB(J,1)*1.d-5*P*T/500.d0
         XB(J,14)=XB(J,1)*1.d-6*P*200.D0/HIT
         XB(J,15)=XB(J,1)*1.d-7*P*P
         XB(J,16)=XB(J,1)*1.d-8*P*P*200.D0/LOT
         XB(J,17)=XB(J,1)*P*P/HIT3
         XB(J,18)=XB(J,1)*1.d-7*P*P*T/500.d0
         XB(J,19)=XB(J,1)*p*p*T2*1.d-7
         XB(J,20)=XB(J,1)*1.d-6*4.d4*P/LOT/LOT
         XB(J,21)=XB(J,1)*1.d-5*P*T*T/500.d0/5.D2
         XB(J,22)=XB(J,1)*1.d-10*P*P*P*2.D2/LOT
         XB(J,23)=xb(j,1)*p*p*p*1.d-2/hit3
         XB(J,24)=XB(J,1)*2.D2/HIT/HIT/HIT
  1   CONTINUE

         BETA0V=0.d0
         BETA1V=0.d0
         CPHIV=0.d0
         CPHIV1=0.d0
         DO 2 I=1, 24
           BETA0V=BETA0V+XB(1,I)*B0(I)
 2       CONTINUE
         DO 3 I=1, 24
           BETA1V=BETA1V+XB(2,I)*B1(I)
 3       CONTINUE
         DO 4 I=1, 24
           CPHIV=CPHIV+XB(3,I)*C0(I)
 4       CONTINUE
         DO 5 I=1, 24
           CPHIV1=CPHIV1+XB(4,I)*C1(I)
 5       CONTINUE


      DHL=-ABS(ZM*ZX)*APHIT*((SQRT(IS)/(1.d0+BB*SQRT(IS))) +
     1 2.d0*LOG(1.d0+BB*SQRT(IS))/BB)
      DHL=DHL+BETA0V+BETA1V+CPHIV+CPHIV1
      GAM=EXP(DHL)


c          if((bflag.eq.1.0).AND.(M.NE.0.0D0))then
c            b0out=beta0v/xb(1,1)
c            b1out=beta1v/xb(2,1)
c            c0out=cphiv/xb(3,1)
c            c1out=cphiv1/xb(4,1)
c            print *,aphit,b0out,b1out
c             print *,c0out,c1out
c            pause
c          endif

        XB(1,1)=M*2.d0*NUM*NUX/NU
        XB(2,1)=XB(1,1)*EXP(-AIS)
        XB(3,1)=XB(1,1)*NUM*ZM*M*2.D0
        XB(4,1)=XB(3,1)*EXP(-AI2)


       DO 10 J = 1,4
         XB(J,2)=XB(J,1)*T/1.D+3
         XB(J,3)=XB(J,1)*T2
	 XB(J,4)=XB(J,1)/LOTA
         XB(J,5)=XB(J,1)*1.D+4/LOTA/LOTA/lota
         XB(J,6)=XB(J,1)*1.d2/LOTA/LOTA
         XB(J,7)=XB(J,1)*2.d2/T/T

         XB(J,8)=XB(J,1)*T2*T/500.D0
         XB(J,9)=XB(J,1)/sqrt(HIT)
         XB(J,10)=XB(J,1)*1.d-5*P
         XB(J,11)=XB(J,1)*1.d-6*P*200.D0/LOT
         XB(J,12)=XB(J,1)*P*1.d2/HIT3
         XB(J,13)=XB(J,1)*1.d-5*P*T/500.d0
         XB(J,14)=XB(J,1)*1.d-6*P*200.D0/HIT
         XB(J,15)=XB(J,1)*1.d-7*P*P
         XB(J,16)=XB(J,1)*1.d-8*P*P*200.D0/LOT
         XB(J,17)=XB(J,1)*P*P/HIT3
         XB(J,18)=XB(J,1)*1.d-7*P*P*T/500.d0
         XB(J,19)=XB(J,1)*p*p*T2*1.d-7
         XB(J,20)=XB(J,1)*1.d-6*4.d4*P/LOT/LOT
         XB(J,21)=XB(J,1)*1.d-5*P*T*T/500.d0/5.D2
         XB(J,22)=XB(J,1)*1.d-10*P*P*P*2.D2/LOT
         XB(J,23)=xb(j,1)*p*p*p*1.d-2/hit3
         XB(J,24)=XB(J,1)*2.D2/HIT/HIT/HIT
 10     CONTINUE

         BETA0V=0.d0
         BETA1V=0.d0
         CPHIV=0.d0
         CPHIV1=0.d0
         DO 20 I=1, 24
           BETA0V=BETA0V+XB(1,I)*B0(I)
 20       CONTINUE
         DO 30 I=1, 24
           BETA1V=BETA1V+XB(2,I)*B1(I)
 30       CONTINUE
         DO 40 I=1, 24
           CPHIV=CPHIV+XB(3,I)*C0(I)
 40       CONTINUE
         DO 50 I=1, 24
           CPHIV1=CPHIV1+XB(4,I)*C1(I)
 50       CONTINUE

      DHL=-ABS(ZM*ZX)*APHIT*SQRT(IS)/(1.d0+BB*SQRT(IS))
      DHL=DHL + BETA0V + BETA1V + CPHIV + CPHIV1
      PHI =  1.d0 + DHL
      RETURN
      END

C*****************************************************************************
          SUBROUTINE LPHI(T,P,M,HEXT)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
           DOUBLE PRECISION M,NUM,NUX,NU,MR,LOT,LOTA,LOT2,LOT3
            COMMON /REFS/TR,PR,MR,SH2OTR
            COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
            COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,
     1                    DRHODT,D2RDP2,D2RDT2,Gw,Sw
            DIMENSION XB(4,24)
            DATA R/.0083144D0/
            DATA NUM,NUX,NU,ZM,ZX/1.D0,1.D0,2.D0,1.D0,-1.D0/
            DATA IAC,IAK/0,0/

           if(m.le.0.0d0)then
                 hext=0.0d0
                 return
           endif

           CALL DHLLN(T,P,AP,AH,AV,AC,AK,IAC,IAK)
           AH=AH*R*T
           HEXT=AH*NU*DLOG(1.D0+1.2D0*SQRT(M))/(2.D0*1.2D0)

           AIS=2.0D0*SQRT(M)
           AI2=2.5D0*SQRT(M)           
         AI2TRM=(6.D0-(6.D0+AI2*(6.D0+AI2*(3.D0+AI2)))*EXP(-AI2)
     1       )/((AI2*AI2)*(AI2*AI2))           

      XB(1,1)=-R*T*T*M*2.D0*NUM*NUX
      XB(2,1)=XB(1,1)*2.D0*(1.D0-(1.D0+AIS)*EXP(-AIS))/(AIS*AIS)
      XB(3,1)=XB(1,1)*NUM*ZM*M
      XB(4,1)=XB(3,1)*4.D0*AI2TRM

        tt=t/500.D0
        LOT=T-225.d0
c       NC 10/13/2010 changed from 200 to 215 
c        LOTA=T-200.d0
        LOTA=T-215.d0
c        LOT=T-227.D0
        HIT=650.d0-T
        hit2=hit*hit
        hit4=hit2*hit2
        LOT2=LOT*LOT
        lot3=lot2*lot
      DO 1 J = 1,4
         XB(J,2)=XB(J,1)/1.d+3
         XB(J,3)=XB(J,1)*TT*2.D0/500.D0
	 XB(J,4)=-XB(J,1)/LOTA/LOTA
c         XB(J,5)=-XB(J,1)*4.d+4/T/T
C         XB(J,6)=-2.d0*XB(J,1)*12.d+6/T/T/T
         XB(J,5)=-3.d0*XB(J,1)*1.D+4/LOTA/LOTA/lota/lota
         XB(J,6)=-2.d0*XB(J,1)*1.D+2/LOTA/LOTA/LOTA
         XB(J,7)=-2.d0*XB(J,1)*2.D2/T/T/T
         XB(J,8)=3.D0*xb(j,1)*TT*TT/500.D0
         XB(j,9)=0.5d0*xb(j,1)/(hit**1.5d0)

         XB(J,10)=0.0D0
         XB(J,11)=-XB(J,1)*P*1.D-6*2.D2/LOT2
         XB(J,12)=3.d0*XB(J,1)*P*1.d2/HIT4
         XB(J,13)=XB(J,1)*1.D-5*P/500.D0
         XB(J,14)=XB(J,1)*P*2.D2*1.D-6/HIT2
         XB(J,15)=0.0D0
         XB(J,16)=-XB(J,1)*P*P*1.D-8*2.D2/LOT2
         XB(J,17)=3.d0*XB(J,1)*P*P/HIT4
         XB(J,18)=XB(J,1)*P*P*1.D-7/500.D0
         XB(J,19)=2.d0*XB(J,1)*p*p*TT*1.d-7/500.d0
C         XB(J,20)=XB(J,1)*0.0D0
         XB(J,20)=-2.D0*XB(J,1)*1.d-6*4.d4*P/LOT/LOT/LOT
         XB(J,21)=2.D0*XB(J,1)*1.d-5*P*T/500.d0/5.D2
C         XB(J,21)=XB(J,1)*P*P*P*1.D-10/5.D2
         XB(J,22)=-XB(J,1)*P*P*P*1.D-10*2.D2/LOT2
         XB(J,23)=xb(j,1)*p*p*p*3.d0*1.d-2/hit4
         XB(J,24)=3.D0*XB(J,1)*2.D2/HIT4
         XB(J,1)=0.0D0
  1   CONTINUE


         BETA0=0.d0
         BETA1=0.d0
         CPHI=0.d0
         CPHI1=0.d0
         DO 2 I=1, 24
           BETA0=BETA0+XB(1,I)*B0(I)
 2       CONTINUE
         DO 3 I=1, 24
           BETA1=BETA1+XB(2,I)*B1(I)
 3       CONTINUE
         DO 4 I=1, 24
           CPHI=CPHI+XB(3,I)*C0(I)
 4       CONTINUE
         DO 5 I=1, 24
           CPHI1=CPHI1+XB(4,I)*C1(I)
 5       CONTINUE

            HEXT=HEXT+BETA0+BETA1+CPHI+CPHI1
            
            
           RETURN
           END


C*****************************************************************************
            SUBROUTINE S20(T,P,S)
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            DOUBLE PRECISION MR,M
            COMMON /REFS/TR,PR,MR,SH2OTR
            COMMON /REFS2/GEXTR,SCORR,SIONS,GH2OTR,SEXMR
            COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
            COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,
     1			  DRHODT,D2RDP2,D2RDT2,Gw,Sw
		 COMMON /DEBYE/APHI,AH,AV
            COMMON /DIHYD/DGH2OT

C CALCULATE DIFFERENCE IN ENTROPY OF WATER (KG).
             CALL WSTEAM(T,RHO,P,SI,DPR,DDT,U,H,S,CV,CP,A,4)
                  SH2OT=S+scorr
                  SH2OTL=(SH2OTR-SH2OT)/MR
                  

C CALCULATE DIFFERENCES IN SEX.

                CALL GAMCAL(T,P,MR,GAMMA,PHIC)
               GEXTMR=2.d0*MR*T*.0083144D0*(LOG(GAMMA)+1.D0-PHIc)
               CALL LPHI(T,P,MR,HEXMR)
               SEXT=(HEXMR-GEXTMR/MR)/T
               SEXTL=SEXMR-SEXT
             
C  CALCULATE THE G20(T,p) - G20(Tr,pr)

C            EVALUATE INTEGRALS
               CTOT=0.0d0

C  FIRST INT (1/T**2) INT CP,PR,MR, dT dT


                TINV=1.d0/T-1.d0/TR

               CTOT=CTOT+XC(1)*DLOG(T/TR)
               CTOT=CTOT+XC(3)*(T-TR)/3.D2
               CTOT=CTOT+0.5D0*XC(4)*(T*T-TR*TR)/3.D2/3.D2
               CTOT=CTOT-XC(7)*TINV*1.D2

               DELP=P-PR
               DELP2=0.5d0*(P*P-PR*PR)
               DELP3=(P*P*P-PR*PR*PR)/3.d0
               DELPP=((P+10.D0)**2.5d0-(PR+10.D0)**2.5d0)/2.5d0
               DT=T/300.d0
               VTOT=0.0d0
            vtot=vtot+XV(2)*DELP/3.d2/1.d2
            vtot=vtot+3.d0*XV(4)*DELP*dt*dt/1.d3/3.d2
            vtot=vtot+XV(6)*DELP2/1.d2/1.d1/3.d2
            vtot=vtot+2.d0*XV(7)*DT*DELP2/1.d3/1.d2/3.d2
            vtot=vtot+XV(9)*DELP3/1.d3/1.d4/3.d2
            vtot=vtot+XV(3)*DELPP/1.D7/3.d2


               S=SIONS+SH2OTL+SEXTL-VTOT+CTOT
               RETURN
               END



C*****************************************************************************
            SUBROUTINE DELG20(T,P,DELG)
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            DOUBLE PRECISION MR,M
            COMMON /REFS/TR,PR,MR,SH2OTR
            COMMON /REFS2/GEXTR,SCORR,SIONS,GH2OTR,SEXMR
            COMMON /PARAMS/B0(24),B1(24),C0(24),C1(24),XV(10),XC(11)
            COMMON /DIELC/D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,
     1			  DRHODT,D2RDP2,D2RDT2,Gw,Sw
		 COMMON /DEBYE/APHI,AH,AV
            COMMON /DIHYD/DGH2OT
            DIMENSION XCA(11)
             CALL WSTEAM(T,RHO,P,SI,DPR,DDT,U,H,S,CV,CP,A,4)
                  GH2OT=H-T*(S+SCORR)
                  DGH2OT=(GH2OTR-GH2OT)/MR

             
C  CALCULATE THE G20(T,p) - G20(Tr,pr)

C            EVALUATE INTEGRALS
                TINV=1.d0/T-1.d0/TR

             ts=tr
            XCa(1)=T-TS
            XCA(2)=-3.D2*(1.D0/T-1.D0/TS)
            XCa(3)=0.5D0*(T*T-TS*TS)/300.D0
            XCa(4)=(T*T*T-TS*TS*TS)/27.D4
            XCa(5)=t/(680.d0-t)-ts/(680.d0-ts)-log((680.d0-ts)
     1                /(680.d0-ts))
            XCa(6)=(T*T*T*T-TS*TS*TS*TS)/27.D6/4.D0
            XCa(7)=log(t/ts)

            xca(1)=xca(1)-t*log(t/ts)
            XCA(2)=XCA(2)+T*0.5D0*3.D2*(1.D0/(T*T)-1.D0/(TS*TS))
            xca(3)=xca(3)-t*(t-ts)/3.d2
            xca(4)=xca(4)-t*0.5d0*(t*t-ts*ts)/9.d4
            xca(6)=xca(6)-t*(t*t*t-ts*ts*ts)/27.d6/3.d0
            xca(7)=xca(7)+t*(1.d0/t-1.d0/ts)
            xca(7)=xca(7)*1.d2
             ctot=0.0d0
            do 1010 i=1,7
             ctot=ctot+xc(i)*xca(i)
 1010       continue


               DELP=P-PR
               DELP2=0.5d0*(P*P-PR*PR)
               DELP3=(P*P*P-PR*PR*PR)/3.d0
               DELPP=((P+10.D0)**2.5d0-(PR+10.D0)**2.5d0)/2.5d0
               DT=T/300.d0
               VTOT=0.0d0
            vtot=vtot+XV(1)*DELP*1.d-1
            vtot=vtot+XV(2)*DELP*dt/1.d2
            vtot=vtot+XV(3)*DELP*DT*dt/1.d3
            vtot=vtot+XV(4)*DELP*dt*dt*dt/1.d3
            vtot=vtot+XV(5)*DELP2/1.d3/1.d2
            vtot=vtot+XV(6)*DT*DELP2/1.d2/1.d1
            vtot=vtot+XV(7)*DT*DT*DELP2/1.d3/1.d2
            vtot=vtot+XV(8)*DT*DT*DT*DELP2/1.d2/1.d4
            vtot=vtot+XV(9)*DT*DELP3/1.d3/1.d4
            vtot=vtot+XV(10)*delpp/1.d6
            vtot=vtot+XV(3)*DELPP*DT/1.D7



                 CALL GAMCAL(T,P,MR,GAMMA,PHIC)
                 GEXT=2.d0*MR*T*.0083144D0*(LOG(GAMMA)+1.D0-PHIc)
                 DGEXMR=(GEXTR-GEXT)/MR

               DELG=DGH2OT+DGEXMR+VTOT+CTOT
               DELG=DELG-(T-TR)*(SIONS+SH2OTR/MR+SEXMR)
               RETURN
               END


C****   SUBROUTINE TO CONVERT OLD DHLL CALLS INTO NEW CALLS.

         SUBROUTINE DHLL(TIN,PIN,IN,APHI,AH,AC,AV,iac)
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         T=TIN
         P=PIN
c         IAC=1
          iak=0
c         IAK=1
         CALL DHLLN(T,P,APHI,AH,AV,AC,AK,IAC,IAK)
c         print *,' returning from dhll  aj = ',ac
         RETURN
         END

C****   SUBROUTINE TO CONVERT HGK CALLS INTO HILL CALLS.
C***    NOT ALL VARIABLES AR PASSED.  BE CAREFUL !!!!!!!
 
        SUBROUTINE HGK(T,P,I,RHO,CP,DPD,DPDT1,DVWDT,K,D2T,DTR,D2R)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        print *, 
        LSAV=4
        CALL WSTEAM(T,RHO,P,SI,DPDR,DPDT,U,H,S,CV,CP,A,LSAV) 
          DRHODT=-DPDT/DPDR
          DVWDT=-DRHODT/RHO/RHO
          DPD=DPDR
c        WRITE(5,100)T,P,RHO
 100    FORMAT(1X,' T = ',F6.2,' P = ',E14.6,' RHO = ',F10.8)
        RETURN
        END
 
C**********************************************************************
C 	SUBROUTINE TO CALCULATE DEBYE-HUCKEL SLOPES.
C     
C                 Draft version 1.01, 10/6/89. 
C                 The user of these subroutines assumes all liability 
C                 for the accuracy of the program and the results 
C                 of any calculations performed with this code.
C
C       USAGE:
C         CALL DHLL(T,P,APHI,AH,AV,AC,AK,IAC,IAK)
C
C       ARGUMENTS:
C          T       TEMPERATURE IN KELVIN
C          P       PRESSURE IN MPa
C          APHI    LIMITING LAW SLOPE FOR OSMOTIC COEFFICIENT (UNITLESS)
C          AH      LIMITING LAW SLOPE FOR APPARENT MOLAR ENTHALPY DIVIDED
C                     BY RT (UNITLESS)
C          AV      LIMITING LAW SLOPE FOR APPARENT MOLAR VOLUME. (UNITS
C                     ARE CUBIC CENTIMETERS PER MOLE)
C          AC      LIMITING LAW SLOPE FOR APPARENT MOLAR HEAT CAPACITY
C                     DIVIDED BY R. (UNITLESS)
C          AK      LIMITING LAW SLOPE FOR APPARENT MOLAR COMPRESSIBILITY.
C                     (UNITS ARE CUBIC CENTIMETERS PER MPa PER MOLE)
C          IAC     SWITCH FOR COMPUTATION OF AC (1 = YES); INTEGER
C          IAK     SWITCH FOR COMPUTATION OF AK (1 = YES); INTEGER
C 
C       **    NOTE THAT, IN ABOVE, THERE EXISTS AN INVERSE SQUARE ROOT OF 
C       **    IONIC STRENGTH UNIT, WHICH WAS NOT EXPLICITLY STATED.
C      
C          ALL ARGUMENTS OTHER THAN IAC AND IAK ARE DOUBLE PRECISION.
C
C          THE TWO INTEGER SWITCHES CONTROL CALCULATION OF SECOND DERIVATIVE
C          FUNCTIONS.  THESE CALCULATIONS CONSUME MORE COMPUTATION TIME.
C
C          THE SUBROUTINE CONTAINS A COMMON WHICH IS NAMED DIELC.
C          THIS COMMON CONTAINS THE DIELECTRIC CONSTANT AND THE 
C          FIRST AND SECOND DERIVATIVES OF THE DIELECTRIC CONSTANT WITH
C          RESPECT TO T AND P: D, DDDT, D2DDT2, DDDP, AND D2DDP2 (BUT
C          NOT IN THIS ORDER).  THE COMMON ALSO CONTAINS THE DENSITY, RHO,
C          AND ITS DERIVATIVES: DRHODT, D2RDT2, DRHODP, AND D2RDP2 (BUT
C          NOT IN THIS ORDER).
C
C          DONALD G. ARCHER 
C          ELECTROLYTE DATA CENTER
C          NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C          GAITHERSBURG, MD 20899
C          
C
C**********************************************************************

      SUBROUTINE DHLLN(T,P,APHI,AH,AV,AC,AK,IAC,IAK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C****
C****    COMMON DIELC CONTAINS THE VALUES OF THE DERIVATIVES OF THE 
C****    DIELECTRIC CONSTANT AND WATER DENSITY AFTER THE CALL
C****    TO SUBROUTINE DIEL.
C****
      COMMON /DIELC/ D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT,
     1 D2RDP2,D2RDT2,Gw,Sw

        DOUBLE PRECISION K,MU,MU2

        DATA WM,PI,AVNO,ALPHA,K,MU/18.0153,3.14159,6.0221367E+23,
     1  1.444E-24,1.380658E-23,1.84E-18/

        DATA E,EPS0,R/1.60217733E-19,8.8541878E-12,8.314510/

C   GET THE DIELECTRIC CONSTANT, DENSITY, AND THEIR DERIVATIVES.

           CALL DIEL(T,P,IAC,IAK)

           APHI=SQRT(2.*PI*RHO*AVNO/1000.)
     1               *((E*E*100./(4.*PI*D*EPS0*K*T))**1.5)/3.

C   CALCULATE A FEW MORE DERIVATIVES.

           ALPHAV=-DRHODT/RHO
           DALPAT=ALPHAV*ALPHAV-D2RDT2/RHO
           RD2VDT=2.*DRHODT*DRHODT/(RHO*RHO)-D2RDT2/RHO
           BETAV=DRHODP/RHO
           DBETAP=-BETAV*BETAV+D2RDP2/RHO
           DLNDDT=DDDT/D
           DLNDDP=DDDP/D
           D2LNDT=-DDDT*DDDT/(D*D)+D2DDT2/D
           D2LNDP=-DDDP*DDDP/(D*D)+D2DDP2/D

C   CALCULATE SOME REQUIRED TERMS.

           BRAKET=1./T+DLNDDT+ALPHAV/3.
           BRAK2=-1./(T*T)+D2LNDT+DALPAT/3.
           VBRAK=DLNDDP-BETAV/3.

C   CALCULATE REMAINDER OF DEBYE-HUCKEL PARAMETERS.

           AH=-6.*APHI*T*BRAKET
           AV=6.*APHI*R*T*VBRAK

           IF(IAC .EQ. 1)THEN
            AC = 3.*APHI*T*(1./T+2.*DLNDDT+5.*T
     1              *DLNDDT*DLNDDT+2.*T*DLNDDT*ALPHAV+2.*ALPHAV/
     2              3.+T*ALPHAV*ALPHAV-2.*T*D2DDT2/D-2.*T
     3              *RD2VDT/3.)
            ELSE
            AC=0.0
           ENDIF
           IF(IAK .EQ. 1)THEN 
	    AK=6.*R*T*APHI*(-1.5*VBRAK*VBRAK+D2LNDP
     1             -DBETAP/3.)
           ELSE
            AK=0.0
           ENDIF

C*****    THAT'S ALL FOLKS.

        RETURN
        END

      SUBROUTINE DIEL(T,P,IAC,IAK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DIELC/ D,DDDP,D2DDP2,DDDT,D2DDT2,RHO,DRHODP,DRHODT,
     1 D2RDP2,D2RDT2,Gw,Sw
       COMMON/CPCALC/CPSAVE
        DIMENSION X(12)
        DOUBLE PRECISION K,MU,MU2

C****  Assign constants.

        DATA WM,PI,AVNO,ALPHA,K,MU/18.0153,3.1415927,6.0221367E+23,
     1 1.444E-24,1.380658E-16,1.84E-18/

C****  Assign parameters.

        data X/0.0,-3.548184,-.4044525E-1,0.0000000,
     1 -1246.311,-.6928953,.2633077E+6,103.6180,75.32165,
     2 -23.23778,-204.4473,0.0/

            TK3=3.*K*T
            MU2=MU*MU
            TR=215.
            LSAV=4

C*****   Get rho and first deriviatives for T and p.
C*****

             CALL WSTEAM(T,RHO,P,SI,DPDRHO,DPDT,U,H,S,CV,CPSAVE,
     1            A,LSAV)  
		Gw=H-T*S
		Sw=S
             DRHODP=1./DPDRHO
             DRHODT=-DPDT*DRHODP

C*****   If flagged, get second derivatives of rho with respect to T and p.
C*****
            DELP=0.01
            DELT=0.001

             IF(IAK .EQ. 1)THEN
             
                 CALL WSTEAM(T,RHO1,P-DELP,SI,DPDR1,DPDT,U,H,S,CV,CP,
     1            A,LSAV)  
                 CALL WSTEAM(T,RHO2,P+DELP,SI,DPDR2,DPDT,U,H,S,CV,CP,
     1            A,LSAV)  
           
                  DRDP1=1./DPDR1
                  DRDP2=1./DPDR2
                  D2RDP2=(DRDP2-DRDP1)/(2.*DELP)
             ELSE
                  D2RDP2=0.0
             ENDIF

             IF(IAC .EQ. 1)THEN

                 CALL WSTEAM(T-DELT,RHO1,P,SI,DPDR1,DPDT1,U,H,S,CV,CP,
     1            A,LSAV)  
                 CALL WSTEAM(T+DELT,RHO2,P,SI,DPDR2,DPDT2,U,H,S,CV,CP,
     1            A,LSAV)
  
                 DRDP1=1./DPDR1
                 DRDP2=1./DPDR2
                 DRDT1 = -DPDT1*DRDP1
                 DRDT2 = -DPDT2*DRDP2
                 D2RDT2=(DRDT2-DRDT1)/(2.*DELT)
             ELSE
                 D2RDT2=0.0
             ENDIF

C*****   Calculate g.
C*****

	 EXPTRM = EXP(X(5)/T+X(6)*P/T+X(7)/(T*T)+
     1		  X(11)*P/(T*T)+X(12)*P)
	 Y  =  X(1)+X(2)/((T-TR)**.25)+X(3)*P/T+X(8)/SQRT(T)+
     1         X(9)/((T-TR)**1.00)+X(10)/((T-TR)**0.5)+EXPTRM
	 G  =  Y*RHO+1.

C*****  Calculate first and second temperature derivatives of Y.
C*****

        DPWRDT  = (-X(5)/T-X(6)*P/T-2.*(X(7)+X(11)*P)/(T*T))/T

        DYDT    = -X(8)/(2.*(T**1.5))-X(9)/((T-TR)**2.)-X(3)*P
     1              /(T*T)-X(10)/(2.*((T-TR)**1.5))+EXPTRM*DPWRDT
     2              -X(2)/(4.*((T-TR)**1.25))

        IF(IAC .EQ. 1)THEN 

          D2PRDT  = (2.*X(5)/T+2.*X(6)*P/T+6.*(X(7)+X(11)*P)/(T*T))
     1              /(T*T)
        D2YDT2  = 3.*X(8)/(4.*(T**2.5))+2.*X(9)/((T-TR)**3.)
     1		    +1.5*X(10)/(2.*((T-TR)**2.5))+EXPTRM*DPWRDT*DPWRDT
     2              +EXPTRM*D2PRDT
     3              +2.*X(3)*P/(T*T*T)+1.25*X(2)/(4.*((T-TR)**2.25))
        ELSE
           D2PRDT=0.0
           D2YDT2=0.0
        ENDIF

C*****  Calculate first and second pressure derivatives of Y.
C*****
        DPWRDP = X(6)/T+X(12)+X(11)/(T*T)
        DYDP   = X(3)/T+EXPTRM*DPWRDP

        IF (IAK .EQ. 1)THEN
          D2YDP2 = EXPTRM*DPWRDP*DPWRDP
        ELSE
          D2YDP2=0.0
        ENDIF
        
C*****  Calculate first and second temperature derivatives of g.
C*****
        DGDT   = DYDT*RHO+Y*DRHODT
        D2GDT2 = D2YDT2*RHO+2.*DYDT*DRHODT+Y*D2RDT2

C*****  Calculate first and second pressure derivatives of g.
C*****
        DGDP   = DYDP*RHO+Y*DRHODP
        D2GDP2 = D2YDP2*RHO+2.*DYDP*DRHODP+Y*D2RDP2

C*****   Calculate A.
C*****
           A = G*MU*MU/TK3
           A = (A+ALPHA)*4.*PI*AVNO*RHO/(3.*WM)

C*****   Calculate first and second temperature derivatives of A.
C*****

         DADT   = (DRHODT*(ALPHA/MU2+G/TK3)+RHO*(DGDT-G/T)/TK3)
     1            *4.*PI*AVNO*MU2/(3.*WM)

       IF (IAC .EQ. 1)THEN  
         D2ADT2 = (D2RDT2*(ALPHA/MU2+G/TK3)+2.*DRHODT*(DGDT-G/T)
     1            /TK3+RHO*(2.*G/T-2.*DGDT+T*D2GDT2)
     2            /(TK3*T))*4.*PI*AVNO*MU2/(3.*WM)             
       ENDIF

C*****   Calculate first and second pressure derivatives of A.
C*****
	 DADP	= (DRHODP*(ALPHA/MU2+G/TK3)+RHO*(DGDP/TK3))*4.*
     1            PI*AVNO*MU2/(3.*WM)
      IF(IAK .EQ. 1)THEN
	 D2ADP2 = (D2RDP2*(ALPHA/MU2+G/TK3)+DRHODP*2.*DGDP/TK3+
     1		  RHO*D2GDP2/TK3)*4.*PI*AVNO*MU2/(3.*WM)
      ENDIF

C*****    Calculate D from A.
C***** 
          APOLY = 9.*A*A+2.*A+1.
          D     = 0.25*(1.+9.*A+3.*SQRT(APOLY))

C*****    Calculate first and second derivatives of D w.r.t. A.
C*****

          DDDA   = (9. + 1.5*(18.*A+2.)/SQRT(APOLY))/4.
       IF(IAC .EQ. 1)THEN
          D2DDA2 = (-0.75*(18.*A+2.)*(18.*A+2.)/(APOLY**1.5)+27./
     1             SQRT(APOLY))/4.
       ENDIF

C*****   Calculate first and second temperature derivatives of D.
C*****
          DDDT   = DDDA*DADT

          IF(IAC .EQ. 1) D2DDT2 = DDDA*D2ADT2 + DADT*DADT*D2DDA2

C*****   Calculate first and second pressure derivatives of D.
C*****
          DDDP   = DDDA*DADP

          IF(IAK.EQ.1) D2DDP2 = DDDA*D2ADP2 + DADP*DADP*D2DDA2

          IF(IAK .NE. 1)D2DDP2=0.0
          IF(IAC .NE. 1)D2DDT2=0.0

C****  Finished.  

	  RETURN
           END
      

C********************************************************************
C      REMAINDER OF THIS FILE IS THE EQUATION OF STATE OF HILL.
C      THIS CODE COMPILES AND EXECUTES PROPERLY UNDER MS FORTRAN AND 
C      LAHEY F77L COMPILERS.  ERRORS IN SOME CALCULATIONS HAVE BEEN
C      FOUND WHEN EXECUTING UNDER LAHEY'S EM32 COMPILER (v. 1.01).
C                                                    D.G.A. 4/26/89 
C********************************************************************



      SUBROUTINE WSTEAM(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,
     1      A,LSAV)
C
C
C********************************************************************
C
C         VARIABLES IN CALL STATEMENT
C
C********************************************************************
C         TEMK         TEMPERATURE                     K
C         RHO          DENSITY                         kJ/dm3
C         PRESS        PRESSURE                        MPa
C         SI           HELMHOLTZ FREE ENERGY           kJ/kg.K
C         DPDR         PRESSURE-DENSITY DERIVATIVE     MPa.dm3/kg
C         DPDT         PRESSURE-TEMPERATURE DERIV.     MPa/K
C         U            INTERNAL ENERGY                 kJ/kg
C         H            ENTHALPY                        kJ/kg
C         S            ENTROPY                         kJ/kg.K
C         CV           ISOCHORIC SPECIFIC HEAT         kJ/kg.K
C         CP           ISOBARIC SPECIFIC HEAT          kJ/kg.K
C         A            ISENTROPIC SONIC SPEED          m/s
C
C           ADDITIONAL VALUES PASSED IN COMMON /CARGA/
C    
C           CALL ROUTINE SWVIR FOR VIRIAL COEFF VALUES
C
C********************************************************************
C
C         OPTIONS
C
C********************************************************************
C
C         L = 1        NO OPERATION
C                           (CODE VALUE USED BY AUTOMATIC INITIALIZATION
C                            PROCEDURE ON FIRST CALL TO WSTEAM)
C
C         L = 2        GIVEN RHO AND TEMK CALCULATE ALL VARIABLES
C         L = 3        GIVEN PRESS AND TEMK ( VAPOUR SIDE )
C                            CALCULATE ALL VARIABLES
C         L = 4        GIVEN PRESS AND TEMK ( LIQUID SIDE )
C                            CALCULATE ALL VARIABLES
C         L = 5        GIVEN TEMK CALCULATE SATURATION VALUES
C
C                 [ITERATED SAT VALUES PASSED IN LABELLED COMMON]
C                 /SATLIQ/  AND  /SATVAP/
C
C
C         L = 6        GIVEN TEMK AND PRESS, CALCULATE PSAT, THEN
C                             CALCULATE ALL DERIVATIVES ON LIQ OR VAP
C                             SIDE, DEPENDING ON INPUT PRESS VALUE
C
C
C         NOTES:  ERROR MESSAGES PRINTED ON FORTRAN UNIT 6
C                APPROX SAT VALUES COULD BE OBTAINED USING FPSWAG
C                CORRELATION FOR PS WITH L=3(VAP) OR L=4 (LIQ)
C
C                FOR L=3 WITH TEMK WHERE PRESS .GT. PSAT WSTEAM MAY
C                CALCULATE METASTABLE STATES, (AND
C                 SIMILARLY FOR L=4)
C           RELEASE 2.0   (JUNE 1989)
C                COMMENTS AND SUGGESTIONS FOR IMPROVEMENT OF THESE
C                ROUTINES SHOULD BE DIRECTED TO
C                         P.G. HILL,
C                         DEPT. MECH. ENGINEERING,
C                         UNIVERSITY OF BRITISH COLUMBIA
C                         VANCOUVER, B.C., CANADA
C
C                THIS PROGRAM IS WRITTEN IN FORTRAN/77
C
C                THE APPROXIMATE RANGE OF THE INPUT DATA USED
C                TO DETERMINE THE COEFFICIENTS IS:
C
C                         273.15 <  TEMK  <  1273.15  K
C                          0.0   < PRESS  <  1000.    MPa
C                          0.0   <  RHO   <  1.       kJ/dm3
C
C              (ALSO INCLUDED IS A LIMITED AMOUNT OF DATA
C               AT 293 K, UP TO APPROX 40,000 MPa)
C
C********************************************************************
C
C         NORMALIZING CONSTANTS
C
C********************************************************************
C
C         GAS CONSTANT            RVAL = 0.46151 kJ/kg.K
C         CRITICAL TEMPERATURE      TC = 647.067 K
C         CRITICAL DENSITY        RHOC = 0.322778 kJ/dm3
C
C*******************************************************************
C
C         DIMENSIONLESS VARIABLES
C
C********************************************************************
C
C         R = RHO/RHOC                DR = R - 1
C         T = -TC/TEMK                DT = T + 1
C         P = PRESS/(RHOC*RVAL*TEMK)
C        SI = SI/(RVAL*TEMK)
C         U = U/(RVAL*TC)
C         H = H/(RVAL*TC)
C         S = S/RVAL
C        CV = CV/RVAL
C        CP = CP/RVAL
C         A = A/SQRT(RVAL*TC)
C
C********************************************************************
C
C        THE DIMENSIONLESS EQUATION OF STATE
C
C********************************************************************
C
C        SI = SIF + F*( SIN - SIF )
C        SIF = SI0 + DLOG(R) + W1 + E*W2 + G*W3 + H*W4
C                IN WHICH
C        SI0 = IDEAL GAS FUNCTION DEFINED BY SUBROUTINE SI0F
C        SIN =REVISED AND EXTENDED SCALING FUNCTION OF LEVELT-SENGERS
C             KAMGAR-PARSI, BALFOUR, AND SENGERS (SEPT 1981) WITH
C             ADJUSTMENT OF THE TWO ARBITRARY CONSTANTS
C             EVALUATED IN SUBROUTINE SINFUN
C        W1 = SUM OVER I AND J OF C1(I,J)*R**(I-2)*T**(J-1)*(1-E)
C         EXCEPT WHEN I=2 C1(I,J)*(LN(R)(1-E-R*R)+R*R/2)*T**(J-1)
C        W2 = SUM OVER I AND J OF C2(I,J)*R**(I)*T**(J-1)*E
C        W3 = SUM OVER I AND J OF C3(I,J)*R**(I+1)*T**(J+1)*G
C        W4 = SUM OVER I AND J 0F C4(I,J)*R**(I-2)*1-E)*T**(J-1)*H
C         EXCEPT WHEN I=2 C4(I,J)*(LN(R)(1-E-R*R)+R*R/2)*T**(J-1)
C        W1, W2, W3,,W4 AND THEIR DERIVATIVES EVALUATED
C        BY SUBROUTINE WFUNCS
C        E IS DEFINED BY SUBROUTINE EFUNC
C        F IS DEFINED BY SUBROUTINE FFUNC
C        G IS DEFINED BY SUBROUTINE GFUNC
C        H IS DEFINED BY SUBROUTINE HFUNC
C
C*******************************************************************
C
C        DERIVATIVES OF THE DIMENSIONLESS HELMHOLTZ FUNCTION
C
C********************************************************************
C
C         NOTATION     SIR = FIRST  PARTIAL DERIVATIVE OF SI WRT R
C                     SIRR = SECOND PARTIAL DERIVATIVE OF SI WRT R
C                     SIRT = SECOND PARTIAL DERIVATIVE OF SI WRT R AND T
C
C         PRESSURE           P = R*R*SIR
C                           PR = 2*R*SIR + R*R*SIRR
C                           PT = R*R*SIRT
C         INTERNAL ENERGY    U = -SIT
C         ENTHALPY           H = -SIT - R*SIR/T
C         ENTROPY            S = -SI + T*SIT
C         ISOCHORIC SP HT   CV = -T*T*SITT
C         ISOBARIC SP HT    CP = CV + R*R*PR*( T*PT - P )**2
C         ISENTROPIC SPEED   A = SQRT( -CP/( CV*T ) * PR )
C           OF SOUND
C
C********************************************************************
C
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DOUBLE PRECISION W(4,6),D(3),DD(3),C0(8)
          DOUBLE PRECISION QV(20)
       COMMON/CDEL/CDELTA
       COMMON/CCREAD/CREAD
        COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERRLS
          COMMON/CQV/QV
          COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
          COMMON/NORMS/TC,RHOC,RVAL
          COMMON/ITER/ITERQK
          COMMON/OSG/THET48,POSG,RFOSG,RGOSG,HFOSG,HGOSG
          COMMON/SAT/PSAT,RF,RG,HF,HG,SF,SG
        COMMON/CZERO/C0
       COMMON/CINIT/C000,C001
       COMMON/SATLIQ/CTEMKF,CRHOF,CPF,CSIF,CDPDRF,CDPDTF,CUF,CHF,
     *               CSF,CCVF,CCPF,CAF
       COMMON/SATVAP/CTEMKG,CRHOG,CPG,CSIG,CDPDRG,CDPDTG,CUG,CHG,
     *               CSG,CCVG,CCPG,CAG
C
       COMMON/CARGA/DHDPT,DTDPH,DTDPS,EXPIS,X5,X6,X7,X8,X9,X10
       COMMON/CARGF/DHDPTF,DTDPHF,DTDPSF,EXPISF,XF5,XF6,XF7,XF8,XF9,XF10
       COMMON/CARGG/DHDPTG,DTDPHG,DTDPSG,EXPISG,XG5,XG6,XG7,XG8,XG9,XG10
       COMMON/CSCNT/ICSCNT
       COMMON/CWCNT/ICWCNT
C
       DATA INIVAL/0/
C
C      CREAD SET TO 1.0D0 IF COEFFICIENTS READ FROM FORTRAN UNIT 10
C
C
C       OPEN(UNIT=10,FILE='COEF132.DAT')
       PSP=0.0D0
          L=LSAV
       IF(L.EQ.1)RETURN
       PC=22.046D0
       IF(L.NE.2) RHO= 0.0D0
       IF((L.EQ.2).OR.(L.EQ.5)) PRESS= 0.0D0
       CERRP= 0.0D0
       CERRTK=TEMK
       CERRHO=RHO
       IF((L.EQ.3).OR.(L.EQ.4).OR.(L.EQ.6)) CERRP=PRESS
       CERRL= L
C
C
       IF(INIVAL.EQ.0) ITERQK=1
       IF(INIVAL.EQ.0) THEN
       CALL WINIT(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,1)
            INIVAL=1
            ENDIF
C
C       CHECK VALUE OF CALL ARGUMENTS
C
       IF (TEMK.LE.0.0D0) THEN
                          WRITE(6,9901) TEMK
                          RETURN
                          ENDIF
9901   FORMAT('0','TEMK .LE. 0.0 ',E25.15)
C
       IF((L.EQ.2).AND.(RHO.LT.0.0D0)) THEN
                                        WRITE(6,9902)RHO
                                       RETURN
                                       ENDIF
9902   FORMAT('0','RHO .LT. 0.0 ',E25.15)
C
C
C        COMPARE INPUT PRESS ARGUMENT WITH WSTEAM SPINODAL
C                PRESSURE TO GUARD AGAINST UNSUCCESSFUL OR
C                INCORRECT ITERATION IN WRHO ROUTINE TO
C                SUPERSATURATED METASTABLE STATES
C
C        (NOTE THAT FOR L=4, WHERE THE SUPERSATURATED PRESSURE CAN BE
C                           NEGATIVE, WRHO CAN BE CALLED
C                           WITH  PRESS < 0.0)
C
       IF((L.NE.3).AND.(L.NE.4).AND.(L.NE.6)) GO TO 500
C
       IF((L.EQ.6).AND.(PRESS.GT.0.0D0))GO TO 500
       IF((L.EQ.6).AND.(PRESS.LE.0.0D0)) GO TO 599
C
       IF((L.EQ.3).AND.(PRESS.LE.0.0D0)) GO TO 599
       IF((L.EQ.3).AND.(TEMK.GE.TC)) GO TO 500
       IF((L.EQ.4).AND.(PRESS.LE.0.D0).AND.(TEMK.GE.TC))GO TO 599
       IF((L.EQ.4).AND.(TEMK.GE.TC)) GO TO 500
C
       PSW= FPSWAG(TEMK)
       IF((L.EQ.3).AND.(TEMK.LT.TC).AND.(PRESS.LE.PSW))GO TO 500
       IF(L.EQ.3)THEN
       SPING= WSPING(TEMK)
       CALL WDERIV(TEMK,SPING,PSP,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,2)
            ENDIF
       IF((L.EQ.3).AND.(TEMK.LT.TC).AND.(PRESS.GT.PSP*1.0001))THEN
                     WRITE(6,9599) TEMK,PRESS,PSP
9599      FORMAT(' ','WSTEAM: INPUT PRESS GT P(SPINODAL)',
     &              ' FOR L=3',3F15.8)
                     WRITE(6,9597)
9597      FORMAT(' IS L VALUE CORRECT?')
                     RHO=-10.
                     RETURN
                     ENDIF
       IF(L.EQ.3) GO TO 500
C
       IF((L.EQ.4).AND.(TEMK.LT.TC).AND.(PRESS.GE.PSW))GO TO 500
       IF(L.EQ.4)THEN
       SPINF= WSPINF(TEMK)
       IF((TEMK.GE.311.15).AND.(TEMK.LE.423.15))SPINF=SPINF*1.001
       CALL WDERIV(TEMK,SPINF,PSP,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,2)
                   ENDIF
C      
       IF((L.EQ.4).AND.(TEMK.LT.TC).AND.(PRESS.LT.PSP*.999))THEN
                     WRITE(6,9598) TEMK,PRESS,PSP
9598      FORMAT(' ','WSTEAM: INPUT PRESS LT P(SPINODAL)',
     &              ' FOR L=4',3F15.8)
                     WRITE(6,9597)
                     RHO=-10.
                     RETURN
                     ENDIF
       IF(L.EQ.4) GO TO 500
C
599    CONTINUE
       WRITE(6,9903)TEMK,PRESS,L
9903   FORMAT(' ','WSTEAM: PRESS VALUE OUT OF RANGE; ',2E15.6,I3)
       RHO= -10.0
       RETURN
C
500    CONTINUE
C
C
          ITERQK = 0
C      SET ITERQK = 1 TO REDUCE ITERATION TIME FOR WRHO(P,T) 
          IF( L .GE. 3 .AND. L .LE. 5 ) ITERQK = 1
          TSAVE = 0.D0
C
       IF(L.EQ.2) THEN 
      CALL WDERIV(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
           RETURN
           ENDIF
C
       IF((L.EQ.3).OR.(L.EQ.4)) THEN
      CALL WRHO(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
            RETURN
            ENDIF
C
       IF(L.EQ.5) THEN
      CALL WSAT(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
                 RETURN
                 ENDIF
C
       IF((L.EQ.6).AND.(TEMK.GT.TC)) THEN
      CALL WRHO(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,4)
                  RETURN
                  ENDIF
C
C
       IF(L.EQ.6) 
     *  CALL WSAT(TEMK,RHO,PSC,SI,DPDR,DPDT,U,H,S,CV,CP,A,5)
       IF((L.EQ.6).AND.(PRESS.GE.PSC))
     *  CALL WRHO(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,4)
C
       IF((L.EQ.6).AND.(PRESS.LT.PSC))
     *  CALL WRHO(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,3)
       IF(L.EQ.6) RETURN
C
       WRITE(6,999) L
999    FORMAT(' ','ERROR - L VALUE ACTION UNDEFINED: ',I3)
       RETURN
C
       END
       SUBROUTINE WINIT(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DOUBLE PRECISION W(4,6),D(3),DD(3),C0(8)
          DOUBLE PRECISION QV(20)
       COMMON/CCREAD/CREAD
        COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERRLS
          COMMON/CQV/QV
          COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
          COMMON/NORMS/TC,RHOC,RVAL
          COMMON/ITER/ITERQK
          COMMON/OSG/THET48,POSG,RFOSG,RGOSG,HFOSG,HGOSG
          COMMON/SAT/PSAT,RF,RG,HF,HG,SF,SG
        COMMON/CZERO/C0
       COMMON/CINIT/C000,C001
       COMMON/SATLIQ/CTEMKF,CRHOF,CPF,CSIF,CDPDRF,CDPDTF,CUF,CHF,
     *               CSF,CCVF,CCPF,CAF
       COMMON/SATVAP/CTEMKG,CRHOG,CPG,CSIG,CDPDRG,CDPDTG,CUG,CHG,
     *               CSG,CCVG,CCPG,CAG
C
       COMMON/CARGA/DHDPT,DTDPH,DTDPS,EXPIS,X5,X6,X7,X8,X9,X10
       COMMON/CARGF/DHDPTF,DTDPHF,DTDPSF,EXPISF,XF5,XF6,XF7,XF8,XF9,XF10
       COMMON/CARGG/DHDPTG,DTDPHG,DTDPSG,EXPISG,XG5,XG6,XG7,XG8,XG9,XG10
C
C********************************************************************
C
C        INITIALIZE
C
C********************************************************************
C
C      CONSTANTS FOR TRIPLE LIQ REF STATE
C
C       C0(1)= 0.707510440802D1
C       C0(2)= -.834258572635D1
C       C000= 0.168157893281D-3
C       C001=  -.390307421339D-5
C
C
 1        RVAL = 0.46151D0
          TC = 647.067D0
          RHOC = 0.322778D0
C      CONSTANT FOR E FUNCTION
          EVAL = 1.D0
C      CONSTANTS FOR F FUNCTION
          DR0 = 0.23D0
          DT0 = 0.05D0
          DELTA = 1.028667D0
C      CONSTANTS FOR G FUNCTION
          AVAL = 80.D0
          BVAL = 1.D0
          CVAL = 130.D0
          DVAL = 12.D0
          HVAL = 4.0D0
C
C      CONSTANTS FOR IDEAL GAS FUNCTION
C       IF(CREAD.EQ.1.0D0)
C     *READ(10,1010) C0(1),C0(2), (C0(I), I = 3,8 )
C 1010  FORMAT(1X,2D19.12,6D15.8)
C      WRITE(5,1010)C0
C
C
C      COEFFICIENTS FOR W1, W2, W3, AND W4 FUNCTIONS
          CALL WFUNCS( Q,Q1,W,1)

C
C      ADDITIONS TO ARBITRARY CONSTANTS OF SCALING FUNCTION
C       IF(CREAD.EQ.1.0D0)   READ(10,99) C000, C001
C 99       FORMAT(1X,10D18.10)
C
C
C
C************************  C000, C001, C0(1), C0(2) WERE ALL ENTERED
C                           INTO BLOCK DATA STATEMENT AS THEY APPEARED
C                           IN HILL'S STEAM TABLE.  THUS, NEXT LINE OF 
C                           CODE IS COMMENTED SO THAT THE CORRECTION IS
C                           MADE TO THE COEFFICIENTS.
C
C       IF(CREAD.NE.1.0D0) RETURN
C************************************************************************
C
C      THIS SECTION OF CODE CALCULATES
C      ADJUSTMENTS TO TRIPLE POINT LIQUID REFERENCE STATE
C       (REQUIRED WHEN VALUES OF COEFFICIENTS ARE CHANGED)
        TEMK = 273.16D0
       CALL WSAT(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
 100        UREF = HF - PRESS/RF
        SREF = SF
        C0(1) = C0(1) - UREF/RVAL/TC
        C0(2) = C0(2) + SREF/RVAL
        C000 = C000 + SREF/RVAL - UREF/RVAL/TC
        C001 = C001 + UREF/RVAL/TC
C        WRITE(6,101)
C        WRITE(6,102) C0(1),C0(2)
C        WRITE(6,103) C000,C001
 101        FORMAT(10X,'ARBITRARY CONSTANTS AFTER ADJUSTMENT'
     1 'TO TRIPLE POINT LIQUID STATE')
 102        FORMAT(15X,'C0(1) = ',D19.12,5X,'C0(2) = ',D19.12)
 103      FORMAT(15X,'C000  = ',D19.12,5X,'C001  = ',D19.12)
          RETURN
       END
       SUBROUTINE WDERIV(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,
     *                   A,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DOUBLE PRECISION W(4,6),D(3),DD(3),C0(8)
          DOUBLE PRECISION QV(20)
       COMMON/CCREAD/CREAD
        COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERRLS
          COMMON/CQV/QV
       COMMON/FXXX/FXVAL
          COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
          COMMON/NORMS/TC,RHOC,RVAL
          COMMON/ITER/ITERQK
          COMMON/OSG/THET48,POSG,RFOSG,RGOSG,HFOSG,HGOSG
          COMMON/SAT/PSAT,RF,RG,HF,HG,SF,SG
        COMMON/CZERO/C0
       COMMON/CINIT/C000,C001
       COMMON/SATLIQ/CTEMKF,CRHOF,CPF,CSIF,CDPDRF,CDPDTF,CUF,CHF,
     *               CSF,CCVF,CCPF,CAF
       COMMON/SATVAP/CTEMKG,CRHOG,CPG,CSIG,CDPDRG,CDPDTG,CUG,CHG,
     *               CSG,CCVG,CCPG,CAG
C
       COMMON/CARGA/DHDPT,DTDPH,DTDPS,EXPIS,X5,X6,X7,X8,X9,X10
       COMMON/CARGF/DHDPTF,DTDPHF,DTDPSF,EXPISF,XF5,XF6,XF7,XF8,XF9,XF10
       COMMON/CARGG/DHDPTG,DTDPHG,DTDPSG,EXPISG,XG5,XG6,XG7,XG8,XG9,XG10
C
C********************************************************************
C
C        GIVEN RHO AND TEMPERATURE DERIVE ALL PROPERTIES (L=2)
C
C********************************************************************
C
 20     R = RHO/RHOC
       CERRHO=RHO
        T = -TC/TEMK
        DT = T + 1.D0
C
       IF(RHO.EQ.0.0D0)THEN
            PRESS=0.0D0
            DPDT=0.0D0
            DPDR=1.0D0*RVAL*TEMK
            CALL SI0F(T,SI0,SI0T,SI0TT,2)
            CV=-T*T*SI0TT
            CP= CV+1.D0
            U= -SI0T
            H= (U-1.D0/T)
C 
            CV= CV*RVAL
            CP=CP*RVAL
            U= U *RVAL*TC
            H= H * RVAL*TC
            A= DSQRT(CP/CV*RVAL*TEMK)
            RETURN
            ENDIF
C
        CALL WFUNCS( T,R,W,2 )
C
        W1   = W(1,1)
        W1R  = W(1,2)
        W1RR = W(1,3)
        W1T  = W(1,4)
        W1TT = W(1,5)
        W1RT = W(1,6)
C
        W2   = W(2,1)
        W2R  = W(2,2)
        W2RR = W(2,3)
        W2T  = W(2,4)
        W2TT = W(2,5)
        W2RT = W(2,6)
C
        W3   = W(3,1)
        W3R  = W(3,2)
        W3RR = W(3,3)
        W3T  = W(3,4)
        W3TT = W(3,5)
        W3RT = W(3,6)
C
        W4   = W(4,1)
        W4R  = W(4,2)
        W4RR = W(4,3)
        W4T  = W(4,4)
        W4TT = W(4,5)
        W4RT = W(4,6)
C
        CALL EFUNC( R,E,ER,ERR )
        CALL FFUNC( R,T,F,FR,FRR,FT,FTT,FRT )
        CALL GFUNC( R,T,G,GR,GRR,GT,GTT,GRT )
        CALL HFUNC( R,T,H,HR,HRR,HT,HTT,HRT )
        CALL SI0F( T,SI0,SI0T,SI0TT,2 )
C
C********************************************************************
C
C       FAR FIELD FUNCTIONS
C
C********************************************************************
C
         RLOG = -170.D0
         IF( R .GT. 1.5D-37 ) RLOG = DLOG(R)
        SIF = SI0 + RLOG + W1 + E*W2 + G*W3 + H*W4
        PF = R + R*R*( W1R + ER*W2 + E*W2R + GR*W3 +G*W3R 
     1                 +HR*W4 + H*W4R )
        PFR = 1.D0 + 2.D0*R*( W1R + ER*W2 + E*W2R + GR*W3 + G*W3R 
     1        +HR*W4 + H*W4R )
     2        +R*R*( W1RR + ERR*W2 +2.D0*ER*W2R + E*W2RR
     3                    + GRR*W3 +2.D0*GR*W3R + G*W3RR 
     4                    + HRR*W4 +2.D0*HR*W4R + H*W4RR )
C
        PFT = R*R*( W1RT + ER*W2T + E*W2RT + GRT*W3 + GR*W3T
     1                                     + GT*W3R + G*W3RT 
     2                                     + HT*W4R + H*W4RT )
        UF = -SI0T - W1T - E*W2T - GT*W3 - G*W3T -HT*W4 -H*W4T
        CVF = -T*T*( SI0TT + W1TT + E*W2TT + GTT*W3 + 2.D0*GT*W3T
     1                                              + G*W3TT 
     2              +HTT*W4 + 2.D0*HT*W4T + H*W4TT )
C
       IF(F.NE.0.D0) THEN 
         CALL SINFUN(TEMK,RHO,SIN,PN,UN,PNR,PNT,CVN,2)
          GO TO 381
          ENDIF
C
        SI = SIF
        P = PF
        PR = PFR
        PT = PFT
        U = UF
        CV = CVF
       PRESS= P*RHOC*RVAL*TEMK
       SI= SI*RVAL*TEMK
       DPDR= PR*RVAL*TEMK
       IF(ITERQK.EQ.1)RETURN
        GO TO 389
C
C
C********************************************************************
C
C       COMBINED FUNCTIONS
C
C*******************************************************************
C
 381    CONTINUE
        SIN = SIN + C000 + C001*DT
        UN = UN - C001
        SI = SIF + F*( SIN - SIF)

C P IS Pbar = p/(rho(crit)*R*T) WHERE R IS THE GAS CONSTANT.
C  IN HILL'S CODE R IS THE REDUCED DENSITY.

        P = PF + F*( PN -PF ) + R*R*FR*( SIN - SIF )
        TERM1 = PFR
        TERM2 = F*( PNR - PFR )
        TERM3 = 2.D0*FR*( PN - PF )
        TERM4 = ( 2.D0*R*FR + R*R*FRR )*( SIN -SIF )

C   PR IS (dPbar/drho) AT CONSTANT T

        PR = TERM1 + TERM2 + TERM3 + TERM4
       PRESS= P*RHOC*RVAL*TEMK
       SI= SI*RVAL*TEMK
       DPDR= PR*RVAL*TEMK
       IF(ITERQK.EQ.1) RETURN
C
C
        TERM5 = PFT
        TERM6 = F*( PNT -PFT ) + FT*( PN - PF )
        TERM7 = -R*R*FR*( UN - UF )
        TERM8 = R*R*FRT*( SIN - SIF )

C  PT IS (dPbar/dT) AT CONSTANT RHO

        PT = TERM5 + TERM6 + TERM7 + TERM8
        U = UF + F*( UN -UF ) - FT*( SIN - SIF )
        TERM9 = CVF
        TERM10 = F*( CVN - CVF )
        TERM11 = 2.D0*T*T*FT*( UN -UF )
        TERM12 = -T*T*FTT*( SIN -SIF )
        CV = TERM9 + TERM10 + TERM11 + TERM12
 389    CONTINUE
        H = U - P/R/T
      CP = 1000.D10
        A = 0.D0
      DELR = R - 1.D0
      IF(DABS(DT).LE.1.D-7.AND.DABS(DELR).LE.6.D-4)GO TO 3891
C      IF(((DABS(PR).LE.1.D-20).OR.(DABS(PR).GT.1.D20)) .OR.
C     &  (DABS(T*PT-P) .GT. 1.D15)) THEN
c       PRINT *,' DABS(PR) D-20 <> D+20 ',DABS(PR)
c       PRINT *,' DABS(T*PT-P) < D+15 ',DABS(T*PT-P)

      IF(((DABS(PR).LE.1.D-30).OR.(DABS(PR).GT.1.D30)) .OR.
     &  (DABS(T*PT-P) .GT. 1.D20)) THEN
                   WRITE(6,9999)TEMK,RHO,PRESS,L
                   RETURN
                   ENDIF
        CP = CV + ( T*PT - P )**2/( R*R*PR )
         A = ( CP/CV*(-PR/T) )
       IF(A.LE.0.0D0) THEN
C          PRINT *,' A > 0.0 ',A
                      A=-10.D0
C                      WRITE(6,9999) TEMK,RHO,PRESS,L
C NEXT LINE COMMENTED 12/5/89.  A IS DIMENSIONLESS SPEED OF SOUND.
C  IN THE REMAINDER OF THIS SUBROUTINE'S CODE FURTHER CALCULATIONS OF
C   A ARE NOW MADE CONDITIONAL ON SIGN OF A.
C                      RETURN
                      ENDIF
9999   FORMAT('0','ERROR WDERIV: CANNOT CALCULATE A',3E20.10,I3/
     *    ' TK OR RHO OUT OF RANGE OR IN UNSTABLE SUPERSAT REGION')
       IF(A .GE. 0.D0)A= DSQRT(A)
 3891   CONTINUE
        S = -SI/RVAL/TEMK - U*T
C
C
C        RESTORE FROM DIMENSIONLESS VARIABLES
C
        DPDT = ( PT*RHOC*RVAL*TC + PRESS )/TEMK
        U = U*RVAL*TC
        H = H*RVAL*TC
        S = S*RVAL
        CV = CV*RVAL
        CP = CP*RVAL
        IF(A .GE. 0.D0)A = A*DSQRT( 1000.D0*RVAL*TC )
C
C
C       CALCULATION ADDITIONAL QUANTITIES
C               DHDPT:  JOULE THOMSON COEFF (ISOTHERMAL)
C               DTDPH:  JOULE THOMSON COEFF (ISENTHALPIC)
C               DTDPS:  PRESSURE TEMPERATURE COEFF (ISENTROPIC)
C               EXPIS:  ISENTROPIC EXPONENT (RHO/P*DPDRS)
C
C               QV(20): ARRAY OF VALUES USED BY PLOT PROGRAMS
C
C      DPDR= 0.0 AT CRITICAL POINT)(
C
       IF(DPDR.EQ.0.0)RETURN
       DTDPH= 1.0D0/RHO-TEMK/RHO**2*DPDT/DPDR
       DTDPH= -DTDPH/CP
       DTDPS=  TEMK/CP/RHO**2*DPDT/DPDR
       DHDPT= DTDPH*(-CP)
       EXPIS= RHO/PRESS*DPDR*CP/CV
C
       QV(1)= DTDPH
       QV(2)= DTDPS
       QV(6)= DHDPT
        RETURN
       END
C
C*******************************************************************
C
C       GIVEN PRESSURE AND TEMPERATURE FIND RHO ETC(L = 3 OR 4)
C
C********************************************************************
C
       SUBROUTINE WRHO(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,
     *                   A,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DOUBLE PRECISION W(4,6),D(3),DD(3),C0(8)
          DOUBLE PRECISION QV(20)
          COMMON/CQV/QV
          COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
          COMMON/NORMS/TC,RHOC,RVAL
          COMMON/ITER/ITERQK
          COMMON/OSG/THET48,POSG,RFOSG,RGOSG,HFOSG,HGOSG
          COMMON/SAT/PSAT,RF,RG,HF,HG,SF,SG
        COMMON/CZERO/C0
       COMMON/CINIT/C000,C001
       COMMON/SATLIQ/CTEMKF,CRHOF,CPF,CSIF,CDPDRF,CDPDTF,CUF,CHF,
     *               CSF,CCVF,CCPF,CAF
       COMMON/SATVAP/CTEMKG,CRHOG,CPG,CSIG,CDPDRG,CDPDTG,CUG,CHG,
     *               CSG,CCVG,CCPG,CAG
C
       COMMON/CARGA/DHDPT,DTDPH,DTDPS,EXPIS,X5,X6,X7,X8,X9,X10
       COMMON/CARGF/DHDPTF,DTDPHF,DTDPSF,EXPISF,XF5,XF6,XF7,XF8,XF9,XF10
       COMMON/CARGG/DHDPTG,DTDPHG,DTDPSG,EXPISG,XG5,XG6,XG7,XG8,XG9,XG10
       COMMON/CWCNT/ICWCNT
 600    CONTINUE
       RHOGSS=PRESS/TEMK/.4
       IF(RHOGSS.GT. 1.0) RHOGSS=1.0
       IF((L.EQ.4).AND.(TEMK.LT.TC))RHOGSS=FRFWAG(TEMK)
C
C      ITERATION TO RHO FROM PRESS INPUT ARG
C

       SPINF= WSPINF(TEMK)
       IF((TEMK.GE.311.15).AND.(TEMK.LE.321.15))
     &     SPINF=SPINF*.001 +SPINF
       SPING= WSPING(TEMK)
C
       DO 500 J=1,40
	
       CALL WDERIV(TEMK,RHOGSS,PCALC,SI,DPDR,Q6,Q7,Q8,Q9,Q10,Q11,Q12,L)

C	write(6,*) RHOGSS,PCALC,PRESS
       IF(PRESS.NE.0.0D0)THEN
       IF(DABS(1.0D0-PCALC/PRESS).LT..1D-8)GO TO 510
            ENDIF
      IF((L.EQ.4).AND.(DABS(PCALC-PRESS).LT..1D-8))GO TO 510
C
C      CHECK FOR DPDR NEGATIVE
C
       IF((TEMK.GE.TC).AND.(DPDR.LT.0.0D0))GO TO 599
       IF(TEMK.GE.TC) GO TO 520
       IF((L.EQ.3).AND.(RHOGSS.LT.SPING).AND.(DPDR.GE.0.0D0))
     &    GO TO 520
       IF((L.EQ.4).AND.(RHOGSS.GT.SPINF).AND.(DPDR.GT.0.0D0))
     &    GO TO 520
C
       IF((L.EQ.3).AND.(RHOGSS.GT.SPING))RHOGSS=SPING
       IF((L.EQ.4).AND.(RHOGSS.LE.SPINF))RHOGSS=SPINF
       IF(L.EQ.3) RHOGSS=RHOGSS*.98
       IF(L.EQ.4) RHOGSS=RHOGSS*1.02
       GO TO 500
520    CONTINUE
C
       DPDX= DPDR*1.1
       DPDLIM= .1
      IF(((PRESS.GE.22.046D0).AND.(PRESS.LE.23.))
     *  .AND.((TEMK.GE.647.067D0).AND.(TEMK.LT.678.15D0)))
     *     DPDLIM= .1D-8
       IF(DPDX.LT.DPDLIM) DPDX=DPDLIM
       X=(PRESS-PCALC)/DPDX
       IF(DABS(X) .GT. .1D0)X= X*.1/DABS(X)
       RHOGSS=RHOGSS+X
       IF(RHOGSS.LE..1D-12)RHOGSS=.1D-12
500    CONTINUE
C
599    CONTINUE
       WRITE(6,9999)L,TEMK,PRESS,PCALC,RHOGSS,DPDR
       RHO=-10.
9999   FORMAT('0','WRHO FAIL CONVG:L,TK,PRESS,PCALC,RGSS: ',
     *   'DPDR: '/I2,F7.2,4E20.10)
       WRITE(6,9996)
9996   FORMAT(' ','IS L SET CORRECTLY?')
       RETURN
C
C
510    CONTINUE
       RHO=RHOGSS
       ICWCNT=J
       IF(L.EQ.5) RETURN
       IF(L.EQ.1)RETURN
C
C
C      (IF L=1, OR L=5, ONLY REQUIRE CALC OF RHO,SI,DPDR FROM PRESS)
C
       ITERQK=0
       CALL WDERIV(TEMK,RHO,PCALC,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
       RETURN
       END


       SUBROUTINE WSAT(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,
     *                   A,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DOUBLE PRECISION W(4,6),D(3),DD(3),C0(8)
          DOUBLE PRECISION QV(20)
       COMMON/CDEL/CDELTA
       COMMON/CCREAD/CREAD
       COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERRLS
          COMMON/CQV/QV
          COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
          COMMON/NORMS/TC,RHOC,RVAL
          COMMON/ITER/ITERQK
          COMMON/OSG/THET48,POSG,RFOSG,RGOSG,HFOSG,HGOSG
          COMMON/SAT/PSAT,RF,RG,HF,HG,SF,SG
        COMMON/CZERO/C0
       COMMON/CINIT/C000,C001
       COMMON/SATLIQ/CTEMKF,CRHOF,CPF,CSIF,CDPDRF,CDPDTF,CUF,CHF,
     *               CSF,CCVF,CCPF,CAF
       COMMON/SATVAP/CTEMKG,CRHOG,CPG,CSIG,CDPDRG,CDPDTG,CUG,CHG,
     *               CSG,CCVG,CCPG,CAG
C
       COMMON/CARGA/DHDPT,DTDPH,DTDPS,EXPIS,X5,X6,X7,X8,X9,X10
       COMMON/CARGF/DHDPTF,DTDPHF,DTDPSF,EXPISF,XF5,XF6,XF7,XF8,XF9,XF10
       COMMON/CARGG/DHDPTG,DTDPHG,DTDPSG,EXPISG,XG5,XG6,XG7,XG8,XG9,XG10
       COMMON/CSCNT/ICSCNT
C
C********************************************************************
C
C       GIVEN TEMPERATURE FIND SATURATION STATES(L=5)
C
C********************************************************************
C
       IF(TEMK.GT.TC) THEN
       WRITE(6,705)TEMK,TC
705     FORMAT('0','WSAT: INPUT TEMK ',F10.3,' ABOVE TC ',F10.3,
     *             ' FOR L=5')
       RETURN
       ENDIF
C
C
C
        T = -TC/TEMK
       BETA= .325D0
        DT = T + 1.D0
      CF = 0.07857D0
      CG = 0.07805D0
      TCDIFF = TC -TEMK
       IF(TCDIFF.EQ.0.0D0) THEN
                           RF=RHOC
                           RG=RHOC
                           GO TO 760
                           ENDIF
       IF(TCDIFF.LE.0.0001D0) THEN
          RF = RHOC + CF*TCDIFF**BETA
          RG = RHOC - CG*TCDIFF**BETA
          GO TO 760
          ENDIF
C
C      ITERATION TO SATURATION (GDIFF= 0.0)
C
       RF= FRFWAG(TEMK)
       RG= FRGWAG(TEMK)
C
       DO 750 JSAT=1,50
       RHO=RF
       CALL WDERIV(TEMK,RHO,PF,SIF,DPDRF,Q6,Q7,Q8,Q9,Q10,Q11,Q12,2)
C
       RHO=RG
       CALL WDERIV(TEMK,RHO,PG,SIG,DPDRG,Q6,Q7,Q8,Q9,Q10,Q11,Q12,2)
C
       DP=PG-PF
       DG= SIG+PG/RG -(SIF+PF/RF)
C
       IF(DABS(DG) .LT. 0.1D-12*RVAL*TEMK) GO TO 760
C
       DP= -DP
       DG= -DG
       AA= 1.D0/RG - 1.D0/RF
       DRG= (DG-DP/RF)/(AA*DPDRG)
       DRF= (DG-DP/RG)/(AA*DPDRF)
       RG= RG+DRG
       RF= RF+DRF
C
750    CONTINUE
C
       WRITE(6,999) TEMK,PF,PG,DP,DG
999    FORMAT('0','WSAT NON-CONVG: TK,PF,PG,DP,DG; ',5e17.10)
C
760    CONTINUE
C
       ICSCNT= JSAT
        ITERQK = 0
        RHO = RG
       CALL WDERIV(TEMK,RHO,PRESS,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
      PSAT = PRESS
      RG = RHO
      HG = H
      SG = S
      AG = A
       CTEMKG=TEMK
       CRHOG=RHO
       CPG=PRESS
       CSIG=SI
       CDPDRG=DPDR
       CDPDTG=DPDT
       CUG=U
       CHG=H
       CSG=S
       CCVG=CV
       CCPG=CP
       CAG=A
       DHDPTG=DHDPT
       DTDPHG=DTDPH
       DTDPSG=DTDPS
       EXPISG=EXPIS
C
        RHO = RF
       CALL WDERIV(TEMK,RHO,PCALCF,SI,DPDR,DPDT,U,H,S,CV,CP,A,L)
      RF = RHO
      HF = H
      SF = S
      AF = A
       CTEMKF=TEMK
       CRHOF=RHO
       CPF=PRESS
       CSIF=SI
       CDPDRF= DPDR
       CDPDTF=DPDT
       CUF=U
       CHF=H
       CSF=S
       CCVF=CV
       CCPF=CP
       CAF=A
       DHDPTF=DHDPT
       DTDPHF=DTDPH
       DTDPSF=DTDPS
       EXPISF=EXPIS
        RETURN
        END
C
C********************************************************************
C********************************************************************
C
        SUBROUTINE WFUNCS( T,R,W,M )
C
C********************************************************************
C
C       OPTIONS
C
C       M = 1       INITIALIZE
C       M = 2       OBTAIN W1, W2, W3, AND ALL THEIR DERIVATIVES
C                                     ( ITERQK = 0 )
C                       OR W1, W2, W3, AND W1R, W2R, W3R
C                                     ( ITERQK = 1 )
C
C********************************************************************
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
        DOUBLE PRECISION W(4,6),C1(10,10),C2(10,12),C3(10,10),C4(10,10),
     1         R1(10),R1R(10),R1RR(10),R2(10),R2R(10),R2RR(10),
     2         T1(12),T1T(12),T1TT(12),T2(12),T2T(12),T2TT(12),
     3         R3(10),R3R(10),R3RR(10),R4(10),R4R(10),R4RR(10)
       COMMON/CCREAD/CREAD
       COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERRLS
       COMMON /CVCC/C1,C2,C3,C4,NT1,NT2,NT3,NT4
        COMMON/NORMS/TC,RHOC,RVAL
        COMMON/ITER/ITERQK
         COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
C
       DATA TSAVE/0.0D0/
C
       
        GO TO (1,200),M
 1      NR1 = 7
        NT1 = 7
        NR2 = 7
        NT2 = 12
        NR3 = 5
        NT3 = 5
        NR4 = 5
        NT4 = 10
        NRMAX = 8
        NTMAX = 12
C
        IF( CREAD.NE.1.D0)RETURN
C
C       N.B. NRMAX MUST BE .GE. NR1, NR2, AND NR3
C            NTMAX MUST BE .GE. NT1, NT2, AND NT3
C            OTHERWISE RX(I) AND TX(J) WILL BE WRONG

C
C*****************************************************************
C
C      READ NON-DIMENSIONAL COEFFICIENTS
C
C********************************************************************
           WRITE(6,101)
      DO 2 J = 1,NT1
C       IF(CREAD.EQ.1.D0)
C     *READ(10,100) (C1(I,J),I=1,NR1)
      DO 120 I=1,NR1
      IF( C1(I,J) .EQ. 0.D0 ) GO TO 120
           WRITE(6,105) I, J, C1(I,J)
 120  CONTINUE
 2    CONTINUE
           WRITE(6,102)
      DO 3 J = 1,NT2
C       IF(CREAD.EQ.1.D0)
C     *READ(10,100) (C2(I,J),I=1,NR2)
      DO 130 I=1,NR2
      IF( C2(I,J) .EQ. 0.D0 ) GO TO 130
           WRITE(6,105) I, J, C2(I,J)
 130  CONTINUE
 3    CONTINUE
           WRITE(6,103)
      DO 4 J = 1,NT3
C       IF(CREAD.EQ.1.D0)
C     *READ(10,100) (C3(I,J),I=1,NR3)
      DO 140 I=1,NR3
      IF( C3(I,J) .EQ. 0.D0 ) GO TO 140
           WRITE(6,105) I, J, C3(I,J)
 140  CONTINUE
 4    CONTINUE
           WRITE(6,103)
      DO 5 J = 1,NT4
C       IF(CREAD.EQ. 1.D0)
C     *READ(10,100) (C4(I,J),I=1,NR4)
      DO 150 I=1,NR4
      IF( C4(I,J) .EQ. 0.D0 ) GO TO 150
           WRITE(6,105) I, J, C4(I,J)
 150  CONTINUE
 5    CONTINUE
 100    FORMAT(1X,10D18.10)
 101    FORMAT(1X,'   I   J        C1(I,J)  ')
 102    FORMAT(1X,'   I   J        C2(I,J)  ')
 103    FORMAT(1X,'   I   J        C3(I,J)  ')
 105    FORMAT(1X,I4,I4,3X,D18.10)
C
        RETURN
C
C*******************************************************************
C
C       EVALUATE W FUNCTIONS
C
C********************************************************************
 200  CALL EFUNC( R,E,ER,ERR )
	NR1 = 7
        NT1 = 7
        NR2 = 7
        NT2 = 12
        NR3 = 5
        NT3 = 5
        NR4 = 5
        NT4 = 10
        NRMAX = 8
        NTMAX = 12
      ETERM = 1.D0 - E
      IF( R .LT. 0.1D-4 ) ETERM = EVAL*R*R
      DO 34 I=1,NRMAX
       ZI=I
      IF(I .EQ. 2) GO TO 31
      R1(I) = ETERM*R**(I-2)
      R1R(I) = ETERM*(ZI-2.D0)*R**(I-3) - ER*R**(I-2)
      R1RR(I) = ETERM*((ZI-2.D0)*(ZI-3.D0))*R**(I-4)
     1         -2.D0*ER*(ZI-2.D0)*R**(I-3)
     2         -ERR*R**(I-2)
      GO TO 32
 31   CONTINUE
C
       IF(R.LE.0.0D0)WRITE(6,999)CERRTK,CERRHO
999    FORMAT('0','WFUNCS R .LE. 0.0; TK,RHO: ',2E25.15)
      DLOGR= DLOG(R)
C
      R1(I) = (1.D0 - E)*DLOGR -EVAL*R*R*DLOGR + EVAL*R*R*0.5D0
      R1R(I) = (1.D0 - E)/R - ER*DLOGR - 2.D0*EVAL*R*DLOGR
      R1RR(I) = - (1.D0 - E)/R/R - 2.D0*ER/R - ERR*DLOGR
     1          - 2.D0*EVAL*DLOGR - 2.D0*EVAL
 32   CONTINUE
C
      R2(I) = R**(I)
      R2R(I) = ZI*R**(I-1)
      R2RR(I) = (ZI*(ZI-1.D0))*R**(I-2)
C
      R3(I) = R**(I+1)
      R3R(I) = (ZI+1.D0)*R**(I)
      R3RR(I) = ((ZI+1.D0)*ZI)*R**(I-1)
      R4(I) = R1(I)
      R4R(I) = R1R(I)
      R4RR(I)= R1RR(I)
 34   CONTINUE
C
       IF(T .EQ. TSAVE) GO TO 250
C
      DO 37 J = 1,NTMAX
      T1(J) = T**(J-1)
      T1T(J) = DFLOAT(J-1)*T**(J-2)
        T1TT(J) = DFLOAT((J-1)*(J-2))*T**(J-3)
      T2(J) = T**(J+1)
      T2T(J) = DFLOAT(J+1)*T**(J)
      T2TT(J) = DFLOAT((J+1)*J)*T**(J-1)
 37   CONTINUE
        TSAVE = T
 250    CONTINUE
C
        DO 301 I = 1,4
        DO 300 J = 1,6
        W(I,J) = 0.D0
 300    CONTINUE
 301    CONTINUE
C
        DO 321 I = 1,NR1
        DO 320 J = 1,NT1
        IF( C1(I,J) .EQ. 0.D0 ) GO TO 319
        W(1,1) = W(1,1) + C1(I,J)*R1(I)  *T1(J)
        W(1,2) = W(1,2) + C1(I,J)*R1R(I) *T1(J)
        W(1,3) = W(1,3) + C1(I,J)*R1RR(I)*T1(J)
        IF( ITERQK .EQ. 1 ) GO TO 319
        W(1,4) = W(1,4) + C1(I,J)*R1(I)  *T1T(J)
        W(1,5) = W(1,5) + C1(I,J)*R1(I)  *T1TT(J)
        W(1,6) = W(1,6) + C1(I,J)*R1R(I) *T1T(J)
 319    CONTINUE
 320    CONTINUE
 321    CONTINUE
C
        DO 331 I = 1,NR2
        DO 330 J = 1,NT2
        IF( C2(I,J) .EQ. 0.D0 ) GO TO 329
        W(2,1) = W(2,1) + C2(I,J)*R2(I)  *T1(J)
        W(2,2) = W(2,2) + C2(I,J)*R2R(I) *T1(J)
        W(2,3) = W(2,3) + C2(I,J)*R2RR(I)*T1(J)
        IF( ITERQK .EQ. 1 ) GO TO 329
        W(2,4) = W(2,4) + C2(I,J)*R2(I)  *T1T(J)
        W(2,5) = W(2,5) + C2(I,J)*R2(I)  *T1TT(J)
        W(2,6) = W(2,6) + C2(I,J)*R2R(I) *T1T(J)
 329    CONTINUE
 330    CONTINUE
 331    CONTINUE
C
        DO 341 I = 1,NR3
        DO 340 J = 1,NT3
        IF( C3(I,J) .EQ. 0.D0 ) GO TO 339
        W(3,1) = W(3,1) + C3(I,J)*R3(I)  *T2(J)
        W(3,2) = W(3,2) + C3(I,J)*R3R(I) *T2(J)
        W(3,3) = W(3,3) + C3(I,J)*R3RR(I)*T2(J)
        IF( ITERQK .EQ. 1 ) GO TO 339
        W(3,4) = W(3,4) + C3(I,J)*R3(I)  *T2T(J)
        W(3,5) = W(3,5) + C3(I,J)*R3(I)  *T2TT(J)
        W(3,6) = W(3,6) + C3(I,J)*R3R(I) *T2T(J)
 339    CONTINUE
 340    CONTINUE
 341    CONTINUE
C
        DO 351 I = 1,NR4
        DO 350 J = 1,NT4
        IF( C4(I,J) .EQ. 0.D0 ) GO TO 349
        W(4,1) = W(4,1) + C4(I,J)*R4(I)  *T1(J)
        W(4,2) = W(4,2) + C4(I,J)*R4R(I) *T1(J)
        W(4,3) = W(4,3) + C4(I,J)*R4RR(I)*T1(J)
        IF( ITERQK .EQ. 1 ) GO TO 349
        W(4,4) = W(4,4) + C4(I,J)*R4(I)  *T1T(J)
        W(4,5) = W(4,5) + C4(I,J)*R4(I)  *T1TT(J)
        W(4,6) = W(4,6) + C4(I,J)*R4R(I) *T1T(J)
 349    CONTINUE
 350    CONTINUE
 351    CONTINUE

        RETURN
        END
C
C********************************************************************
C********************************************************************
C
        SUBROUTINE EFUNC( R,E,ER,ERR )
C
C*******************************************************************
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
C
        ARGE = EVAL*R*R
        IF( ARGE .GT. 170.D0 ) ARGE = 170.D0
        E = DDEXP( -ARGE )
        ER = -2.D0*EVAL*R*E
        ERR = EVAL*( -2.D0 + 4.D0*EVAL*R*R )*E
C
        RETURN
        END
C
C********************************************************************
C
        SUBROUTINE FFUNC( R,T,F,FR,FRR,FT,FTT,FRT )
C
C*******************************************************************
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON/FXXX/FXVAL
        COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
       DR = R - 1.D0
       DT = T + 1.D0
       S = DSQRT( (DR/DR0)**2 + (DT/DT0)**2)
      IF( S .LT. 0.003D0 ) S = 0.003D0
C
       POWER = 4.D0
       POW1  = (POWER - 1.D0)/POWER
       POW2  = (POWER - 2.D0)/POWER
       POW3  = POWER/(POWER - 1.D0)
       POW4  = POWER*(POWER - 1.D0)
       ARGS  = (S/DELTA)**POWER
       IF(ARGS .GT. 23.5D0) ARGS = 23.5D0
       CAPS  = DDEXP( ARGS ) - 1.D0
       ARGC  = 1.D0/CAPS
       IF( ARGC .GT. 35.D0 ) ARGC = 35.D0
       EXPS  = DDEXP( -ARGC )
       FS1   = DLOG( CAPS + 1.D0 )
       FS2   = (1.D0 - CAPS - CAPS*CAPS )/CAPS/CAPS
       FS3   = (1.D0 + CAPS )/CAPS/CAPS
C
C
       F     = 1.D0 - EXPS
       IF( F .LT. 0.1D-9 ) F = 0.D0
       DFDS  = -EXPS/DELTA*FS3*POWER*FS1**POW1
       IF( DABS(DFDS) .LT. 0.1D-12 ) DFDS = 0.D0
       D2FDS2= -EXPS/DELTA/DELTA*FS3*POW4*FS1**POW2
     1         *(1.D0 + FS2*POW3*FS1)
       IF( DABS(D2FDS2) .LT. 0.1D-12 ) D2FDS2 = 0.D0
C
       DSDR = DR/S/DR0/DR0
       DSDT = DT/S/DT0/DT0
       DTFAC = (DT/DT0)**2
       DRFAC = (DR/DR0)**2
       D2SDT2 = (1.D0/S - DTFAC/(S)**3)/DT0/DT0
       D2SDR2 = (1.D0/S -DRFAC/(S)**3)/DR0/DR0
       D2SDTR = -DT*DR/(S*DT0*DR0)**2/S
C
       FT = DFDS*DSDT
       FR = DFDS*DSDR
       FTT = D2FDS2*(DSDT)**2 + DFDS*D2SDT2
       FRR = D2FDS2*(DSDR)**2 + DFDS*D2SDR2
       FRT = DFDS*D2SDTR +D2FDS2*DSDT*DSDR
       IF( DABS(FT) .LT. 0.1D-12 ) FT = 0.D0
       IF( DABS(FR) .LT. 0.1D-12 ) FR = 0.D0
       IF( DABS(FTT).LT. 0.1D-12 ) FTT = 0.D0
       IF( DABS(FRR).LT. 0.1D-12 ) FRR = 0.D0
       IF( DABS(FRT).LT. 0.1D-12 ) FRT = 0.D0
       FXVAL=F
        RETURN
        END
C
C********************************************************************
C
        SUBROUTINE GFUNC( R,T,G,GR,GRR,GT,GTT,GRT )
C
C********************************************************************
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
C
        DR = R - 1.D0
        DT = T + 1.D0
        DT2=DT*DT
        DR2=DR*DR
        A2=AVAL*AVAL
        B2=BVAL*BVAL
        S=DSQRT(A2*DT2+B2*DR2)
      IF( S .LT. 0.0001D0 ) S = 0.0001D0
        SCU=S*S*S
        TOP1 = A2*DT 
        TOP2 = B2*DR 
        DSDT=TOP1/S
        DSDR=TOP2/S
        D2SDT=A2/S - TOP1**2/SCU
        D2SDR=B2/S - TOP2**2/SCU
        DSDRT= -TOP1*TOP2/SCU
C
C
C
        GFAC = -AVAL*DT - BVAL*DR -CVAL*DT2 - DVAL*DR2
        IF(DABS(GFAC) .GT. 100.D0) GFAC=-100.D0
        G = DDEXP(GFAC)
C
C
C
        T1 = -AVAL - 2.D0*CVAL*DT
        T2 = -BVAL - 2.D0*DVAL*DR
        GT=G*T1
        GR=G*T2
C
C
C
         GTT = GT*T1 - 2.D0*CVAL*G
         GRR = GR*T2 - 2.D0*DVAL*G
         GRT = GR*T1
      IF( DABS(G)  .LE. 0.1D-12 ) G = 0.D0
      IF( DABS(GR) .LE. 0.1D-12 ) GR = 0.D0
      IF( DABS(GRR).LE. 0.1D-12 ) GRR= 0.D0
      IF( DABS(GT) .LE. 0.1D-12 ) GT = 0.D0
      IF( DABS(GTT).LE. 0.1D-12 ) GTT= 0.D0
      IF( DABS(GRT).LE. 0.1D-12 ) GRT= 0.D0
        RETURN
        END
C
       DOUBLE PRECISION FUNCTION DDEXP(X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C      THRESHOLD TO SET DEXP FUNCTION TO ZERO
C        FOR LARGE NEGATIVE ARGUMENTS  (-35)
C
         DOUBLE PRECISION X
       IF(X.LE.-35.) DDEXP= 0.0D0
       IF(X.LE.-35.) RETURN
       DDEXP= DEXP(X)
         RETURN
         END
C
C
C*********************************************************************
C
         SUBROUTINE HFUNC( R,T,H,HR,HRR,HT,HTT,HRT )
C
C*********************************************************************
C
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,BVAL,CVAL,DVAL,HVAL
      ARGH = HVAL*(T + 3.D0)
      H = DDEXP(-ARGH )
      HR = 0.D0
      HRR = 0.D0
      HT = -HVAL*H
      HTT = HVAL*HVAL*H
      HRT = 0.D0
         RETURN
         END
C**********************************************************************
C
        SUBROUTINE SI0F( T,SI0,SI0T,SI0TT,L )
C
C********************************************************************
C
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON/NORMS/TC,RHOC,RVAL
      COMMON/CZERO/C0
        DOUBLE PRECISION C0(8)
C
C
 100    SI0 = 0.D0
        SI0T = 0.D0
        SI0TT = 0.D0
        DO 110 I = 1,6
        POWER = ( -T )**(2-I)
        SI0 = SI0 + C0(I)*POWER
        SI0T = SI0T - C0(I)*DFLOAT(2-I)*POWER/(-T)
        SI0TT = SI0TT + C0(I)*DFLOAT(2-I)*DFLOAT(1-I)*POWER/T/T
 110    CONTINUE
        SI0 = SI0 + ( C0(7)*T + C0(8) )*DLOG(-T)
        SI0T = SI0T + C0(7)*( DLOG(-T) + 1.D0 ) +C0(8)/T
        SI0TT = SI0TT + C0(7)/T - C0(8)/T/T
        RETURN
        END
C
C****************************************************************
C****************************************************************
C****************************************************************
C****************************************************************
C****************************************************************
C****************************************************************
C
C     SUBROUTINES FOR THE REVISED AND EXTENDED SCALING SOLUTION
C         TO CALCULATE VALUES NEAR THE CRITICAL POINT
C
C         REF:  J.M.H. LEVELT SENGERS, B. KAMGAR-PARSI,
C               F.W. BALFOUR, J.V. SENGERS
C               "THERMODYNAMIC PROPERTIES OF STEAM IN THE 
C          CRITICAL REGION", J.PHYS.CHEM.REF.DATA VOL 12,P 1(1983)
C
C*****************************************************************
C
      SUBROUTINE SINFUN(T,DGMPCC,SI,P,U,DPDR,DPDT,CV,L)
C
C
C      Given the temperature T(K) and density D(gm/cc), this routine
C     calculates PRESSURE P(MPa), ENERGY U and ENTHALPY H(kJ/kg),
C     ENTROPY, SPECIFIC HEATS Cp and Cv(kJ/kg-K), VELOCITY of SOUND
C     Cs(m/s), and COMPRESSIBILITY COMP(1/MPa).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/COEFS/A(20),Q(20),ANAMES(20),QNAMES(20)
      COMMON/CRITS/TC,RHOC,PC,PCON,UCON,SCON,RVAL
      COMMON/SECDER/D2PDM2
      COMMON/CCMETA/CUNST
      DIMENSION S(2),XK(2),SD(2)
      EQUIVALENCE (PW1,A(5)),(PW2,A(4)),(PW3,A(2)),
     1   (AM0,A(13)),(AM1,A(14)),(AM2,A(15)),(AM3,A(16)),
     2   (P00,Q(11)),(P20,Q(12)),(P40,Q(13)),
     3   (P01,Q(18)),(P21,Q(19)),(P41,Q(20)),
     4   (AA,A(10)),(XK0,A(7)),(XK1,A(12)),(PW11,Q(9)),
     5   (ALPHA,Q(10)),(ALHI,Q(15)),(BESQ,A(9))
C
        DATA INIVAL/0/
C
       IF(L.EQ.0)RETURN
      IF  (INIVAL .NE. 0)  GO TO 50
      CALL CONGEN
      Q(3)   = A(8) * Q(7)
      TC     = Q(4) + A(8)
      RHOC   = Q(5) + A(3)
      PC     = Q(6) + Q(3)
      PCON   = PC / TC
      UCON   = 1.D3 * PC
      SCON   = UCON / TC
      DPCON  = PCON / RHOC / RHOC
      RHOCC  = 0.322778D0
      RVAL   = 0.46151D0
       INIVAL=1
C
C
 50   CONTINUE
      D  = DGMPCC * 1000.D0
      XK(1)=XK0
      XK(2)=XK1
      TEE=(T-TC)/TC
      TW=-TC/T
      DTW=1.D0+TW
      RHO=D/RHOC
C
      DELR = RHO - 1.D0
      IF(DABS(TEE).LE.1.D-7.AND.DABS(DELR).LE.0.0006D0)GO TO 55
      DTSTIN= 1.0D0 - (TC/T)
      CALL CONSIM(RHO,TEE,AMU,TH1,R1,RHO1,S,RHODI,UNSTAB,DTSTIN)
C      WRITE(9,9099) TH1,R1
9099   FORMAT(' ','THETA,R: ',2E20.8)
      DELTMU = AMU
      TT1=TH1*TH1
      TT2=TT1*TT1
      CALL SS(R1,TH1,S,SD)
      PW0=1.D0+DTW*(PW1+DTW*(PW2+DTW*PW3))
      PWMU=AMU*RHODI
C
C       pwmu = dmu * (1 + p11*dt)
C
      P0TH  = P00+P20*TT1+P40*TT2
      P1TH  = P01+P21*TT1+P41*TT2
      DPW0  = XK0*P0TH*R1**(2.D0-ALPHA)
      DPW1  = XK1*P1TH*R1**(2.D0-ALHI)
      DPW   = AA*(DPW0+DPW1)
      PW    = PW0+PWMU+DPW
      P     = PW * PCON * T
      DP0DT = PW1+DTW*(2.D0*PW2+3.D0*PW3*DTW)
      DM0DT = AM1+DTW*(2.D0*AM2+3.D0*AM3*DTW)
      UW    = DP0DT- RHO * DM0DT + PW11 * AMU+S(1)+S(2)
      HW    = PW - TW * UW
      AMW   = AMU + AM0 + DTW*(AM1+DTW*(AM2+DTW*AM3))
      SW    = HW - RHO * AMW
      U     = UW * UCON / D
      SCOND = SCON / D
      ENTROP= SW * SCOND
      SI    = U - T * ENTROP
      PDIMEN= P
      P     = P / RHOCC / RVAL / T
      U     = U / RVAL / TC
      SI    = SI / RVAL / T
      IF  (L .EQ. 1)  RETURN
C
C
      D2P0DT= 2.D0*PW2+6.D0*PW3*DTW
      D2M0DT= 2.D0*AM2+6.D0*AM3*DTW
      CALL AUX(R1,TH1,D2PDT2,D2PDMT,D2PDM2,AA,XK,SD,CVCOEX)
      DPDTCD=DP0DT+PW11*(AMU-RHO/D2PDM2)+S(1)+S(2)-D2PDMT*RHO/D2PDM2
      DPWDTW= PW-TW*DPDTCD
      CVITW2= D2P0DT-RHO*D2M0DT+D2PDT2-(PW11+D2PDMT)**2/D2PDM2
      DPDD  = DPCON * D * T / D2PDM2
      DPDR  = 1000.D0 * DPDD
C     WRITE(6,101) T,D,PCON,DPCON,D2PDM2,DPDD,DPDR
      DPDT  = PCON * DPWDTW
      CVW   = CVITW2 * TW * TW
      CV    = CVW * SCOND
      DPDR  = DPDR / RVAL / T * UNSTAB
C     WRITE(6,101) T,D,DPCON,D2PDM2,DPDD,DPDR
C101  FORMAT(1X,'T,D,Pcon,dPcon,d2Pdm2,dpdd,dpdr ',7D14.6)
      DPDT  = (DPDT * T - PDIMEN) / RHOCC / RVAL / TC
      CV    = CV / RVAL
C
C
      RETURN
 55   CONTINUE
      P = 22.0460D0/RVAL/TC/RHOCC
      U = 2014.82D0 - 468.03D0*DELR
      U = U/RVAL/TC
      ENT = 4.40498D0 - 0.8231D0*DELR
      SI = U - TC*ENT/RVAL/TC
      DPDR = 0.D0
      DPDT = 0.D0
      CV = 1000.D0
      RETURN
      END



       BLOCKDATA CLSBLK
C
C
       DOUBLE PRECISION ANAMES,QNAMES
      DOUBLE PRECISION A,Q
       COMMON/COEFS/A(20),Q(20),ANAMES(20),QNAMES(20)
C
C     This subroutine supplies the parameters used in the equation
C     of state.
C
      DATA A/-.017762D0,5.2380D0,0.D0,-25.4915D0,6.8445D0,.325D0,1.4403
     1D0,0.D0,1.3757D0,23.6666D0,4.82D0,.2942D0,-11.2317D0,
     2-22.6407D0,-17.7876D0,-4.9332D0,0.D0,0.D0,0.D0,0.D0/
      DATA Q/-.006D0,-.003D0,0.D0,647.067D0,322.778D0,22.0460D0
     1,.267D0,-1.6D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     2,0.D0,0.D0,0.D0,0.D0/
      DATA ANAMES/6HC     ,6HCOT3  ,6HDELROC,6HCOT2  ,6HDPDT  ,6HBETA  ,
     1            6HK0    ,6HDELTC ,6HB*B   ,6HA     ,6HDELTA ,6HK1    ,
     2            6HMUC   ,6HMU1   ,6HMU2   ,6HMU3   ,6HS00   ,6HS20   ,
     3            6HS01   ,6HS21   /
      DATA QNAMES/2*6H      ,6HDELPC ,6HTC    ,6HRHOC  ,6HPC    ,
     1              6HDPCDTC,6HSLOPDI,6HP11   ,6HALPHA ,6HP00   ,
     2              6HP20   ,6HP40   ,6HDELTAI,6HALPHAI,6HBETA I,
     3              6HGAMMAI,6HP01   ,6HP21   ,6HP41   /
      END



      SUBROUTINE CONGEN
C
C   This subroutine calculates all quantities not dependent on
C  R or Theta
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /COEFS/ A(20),Q(20),ANAMES(20),QNAMES(20)
C
C
      EQUIVALENCE (DELI,Q(14)),(DRHOC,A(3)),(BETA,A(6)),(DELTA,A(11))
     1,(COTEE2,A(4)),(DPDT,A(5)),(DELTC,A(8)),(BESQ,A(9))
     2,(CC,A(1)),(AA,A(10)),(TCK,Q(4)),(RHOC,Q(5)),(PC,Q(6)),(XK0,A(7))
     3,(ALHI,Q(15)),(BETI,Q(16)),(GAMI,Q(17)),(P01,Q(18)),(P21,Q(19))
     4,(ALPHA,Q(10)),(P00,Q(11)),(P20,Q(12)),(P40,Q(13)),(P41,Q(20))
     5,(S00,A(17)),(S20,A(18)),(S01,A(19)),(S21,A(20))
C
C
      DATA DELPP,DELPM,DELPPH,DELPMH/1.D-6,-1.D-6,1.000001D0,.999999D0/
      ALPHA = 2.D0 - A(6)*(A(11) + 1.D0   )
      IF  (DABS(ALPHA) .LT. 1.D-6)   ALPHA=1.D-6
      GAMMA = BETA*(DELTA - 1.D0)
      DELI = 0.50D0
      IF (DELI .LT. -BETA)   DELI=-BETA
      ALHI = ALPHA - DELI
      BETI = BETA + DELI
      GAMI = GAMMA - DELI
      IF(ALHI.LE.DELPPH.AND.ALHI.GE.DELPMH) ALHI=1.000001D0
      IF (ALHI.LE.DELPP.AND.ALHI.GE.DELPM) ALHI = .000001D0
      DELPC = Q(7)*DELTC
      ERR = 2.D0*BETA*DELTA - 1.D0
      P00 = (BETA*(DELTA-3.D0)-BESQ*ALPHA*GAMMA)/(2.D0*
     1BESQ*BESQ*(2.D0-ALPHA)*(1.D0-ALPHA)*ALPHA)
      P20=-(BETA*(DELTA-3.D0)-BESQ*ALPHA*ERR)
     1 / (2.D0 * BESQ * (1.D0 - ALPHA) * (ALPHA))
      P40=(ERR-2.D0)/2.D0/ALPHA
      S00 = (2.D0-ALPHA)*P00
      S20 = -BETA*(DELTA-3.D0)/2.D0/BESQ/ALPHA
      DA=Q(1)
      DB=Q(2)
      RA = DA/(1.D0 - BESQ)
      RB = DB/(1.D0 - BESQ)
      SW0 = S00+S20
      DR0=SW0*(RA**(1.D0-ALPHA)-RB**(1.D0-ALPHA))
      P01 = (BETA*(DELTA-3.D0)-3.D0*DELI-A(9)*ALHI*GAMI)
     1 /(2.D0*A(9)*A(9)*(2.D0-ALHI)*(1.D0-ALHI)*ALHI)
      P21 = -((BETA*(DELTA-3.D0)-3.D0*DELI-A(9)*ALHI*ERR)/
     1 (2.D0*A(9)*(1.D0-ALHI)*ALHI))
      P41 = (.5D0*ERR -1.D0)/ALHI
      S01 = (2.D0-ALHI)*P01
      S21 = -(BETA*DELTA-3.D0*BETI)/2.D0/BESQ/ALHI
      PTW = P01 + P21 + P41
      SW1 = S01+S21
      DR1=SW1*(RA**(1.D0-ALHI)-RB**(1.D0-ALHI))
      IF(XK0.LE.0.D0) XK0=0.0001D0
      D1=CC*AA*XK0*DR0
      D2=CC*AA*A(12)*DR1
      Q(9) = (Q(8)*(1.D0/(1.D0-DB)-1.D0/
     1   (1.D0 - DA)) + D1 + D2) / (DB - DA)
      RETURN
      END


      SUBROUTINE SS(R,TH,S,SD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(2),SD(2)
      COMMON /COEFS/ A(20),Q(20),ANAMES(20),QNAMES(20)
      EQUIVALENCE (ALPHA,Q(10)),(BETA,A(6)),(BESQ,A(9)),(DELTA,A(11))
     1,(DELI,Q(14)),(ALHI,Q(15)),(BETI,Q(16)),(GAMI,Q(17)),(P00,Q(11))
     2,(P01,Q(18)),(S00,A(17)),(S20,A(18)),(S01,A(19)),(S21,A(20))
      TT = TH*TH
      S(1) = S00 + S20*TT
      SD(1) = 2.D0*S20*TH
      S(2) = S01 + S21*TT
      SD(2) = 2.D0*S21*TH
      S(1)=S(1)*A(10)*A(7)*R**(1.D0-ALPHA)
      S(2)=S(2)*A(10)*A(12)*R**(1.D0-ALHI)
      RETURN
      END



      SUBROUTINE AUX(R1,TH1,D2PDT2,D2PDMT,D2PDM2,AA,XK,SD,CVCOEX)
C
C     This routine calculates second derivatives of the scaled
C     equation of state
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /COEFS/ A(20),Q(20),ANAMES(20),QNAMES(20)
      DIMENSION XK(2),S(2), SD(2),W(2),Y(2),Z(2),COEX(2)
      EQUIVALENCE (CC,A(1)),(BETA,A(6)),(BESQ,A(9)),(DELTA,A(11))
     1,(ALPHA,Q(10)),(S00,A(17)),(S20,A(18)),(S01,A(19)),(S21,A(20))
      DELI = 0.D0
      S(1)=S00+S20*TH1*TH1
      S(2)=S01+S21*TH1*TH1
      SD(1) = 2.D0*TH1*S20
      SD(2) = 2.D0*TH1*S21
      WW = 0.D0
      YY = 0.D0
      ZZ = 0.D0
      GAMMA = BETA*(DELTA -1.D0)
      TT1 = TH1*TH1
      TER = 2.D0*BETA*DELTA - 1.D0
      G = (1.D0+(BESQ*TER-3.D0)*TT1 - BESQ*(TER-2.D0)*TT1*TT1)
      CVCOEX = 0.D0
      DO 30 I = 1,2
      ALHI = ALPHA - DELI
      BETI = BETA + DELI
      GAMI = GAMMA - DELI
      W(I)=(1.D0-ALHI)*(1.D0-3.D0*TT1)*S(I)-
     1  BETA * DELTA * (1.D0 - TT1) * TH1 * SD(I)
      W(I) = (W(I) *(R1**(-ALHI)))/G
      W(I) = W(I)*XK(I)

      WW = WW + W(I)
      Y(I) = BETI*(1.D0-3.D0*TT1)*TH1 - BETA*DELTA*(1.D0-TT1)*TH1
      Y(I) = (Y(I) * (R1**(BETI - 1.D0)))*XK(I)/G
      YY = YY + Y(I)
      Z(I) = 1.D0 - BESQ*(1.D0 - (2.D0*BETI))*TT1
      Z(I) = (Z(I) *(R1**(-GAMI)))*XK(I)/G
      ZZ = ZZ + Z(I)
      A1 = (BETA*(DELTA-3.D0)-3.D0*DELI-BESQ*ALHI*GAMI)
     1/(2.D0*BESQ*BESQ*(2.D0-ALHI)*(1.D0-ALHI)*ALHI)
      A2 = 1.D0+((BETA*(DELTA-3.D0)-3.D0*DELI-BESQ*ALHI*TER)/
     1(2.D0*BESQ*(1.D0-ALHI)*ALHI))
      A2 = - A2
      A4 = 1.D0+((TER-2.D0)/(2.D0*ALHI))
      F1 = A1 + A2 + A4
      COEX(I) = ((2.D0 - ALHI)*(1.D0 - ALHI)*(R1**(-ALHI))*F1*XK(I))
      CVCOEX = CVCOEX + COEX(I)
      DELI = 0.5D0
 30   CONTINUE
      D2PDT2 = AA*WW
      D2PDMT = YY + AA*CC*WW
      D2PDM2 = ZZ/AA + 2.D0*CC*YY + (CC**2)*AA*WW
      RETURN
      END




       SUBROUTINE CONSIM(RM,T,AMU,THETA,R,RHO1S,S1,RHODI,UNSTAB,DTSTIN)
C
C      THIS ROUTINE CONVERTS CALLS IN OLD CONVER FORM TO
C                                     NEW CONVER FORM
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERLS
       COMMON/COEFS/A(20),Q(20),ANMAMES(20),QNAMES(20)
       COMMON/FXXX/FXVAL
       COMMON/CCR/CR,CTHETA
C
       DIMENSION S1(2),SD(2)
C
       AA= A(10)
       BETA=A(6)
       DELTA=A(11)
       P11=Q(9)
       TSTAR= T+1.0D0
C
       RHODI=  1.D0+P11*DTSTIN
       DELRM= RM - RHODI
C
       CALL CONVER(DTSTIN,DELRM,R,THETA,UNSTAB)
C
C      WRITE(9,9900)R,THETA
9900   FORMAT(2F10.5)
       TT= THETA*THETA
       IF(R.GT. 0.0D0)
     & AMU=  AA*R**(BETA*DELTA)*THETA*(1.0D0-TT)
C
C
C      WRITE(6,200) T,RM,R,THETA
200    FORMAT(' T,RM,R,THETA  ',4E20.8)
       CR=R
       CTHETA=THETA
C
       RETURN
       END




      SUBROUTINE CONVER(T,M,R,THETA,UNSTAB) 
*     THIS ROUTINE CONVERTS THE REDUCED TEMPERATURE AND 
*     DENSITY TO REVISED AND EXTENDED PARAMETRIC VARIABLES. 
* 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LTEST 
      DOUBLE PRECISION M,K0,K1 
      PARAMETER (ITRMX2=100,SMLNO=1.0D-08) 
       COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERLS
      COMMON /COEFFS/ A(20),Q(10) 
      COMMON /PARAMS2/ ALPHA,BETA,GAMMA,DELTA,WEG,BB,P00,P02, 
     1                P04,P10,P12,P14,S00,S20,S01,S21 
       DATA INIVAL/0/
C
      S0(TH)=S00+S20*(TH**2) 
      S1(TH)=S01+S21*(TH**2) 
C
       UNSTAB=1.0D0
       IF(INIVAL.EQ.0)THEN
                      CALL CONST
                      INIVAL=1
C                     WRITE(6,9898)S00,S20,S01,S21
9898   FORMAT(' ','S00,S20,S01,S21 ',4E20.8)
                      ENDIF
      AA=A(1) 
      K0=A(2) 
      K1=A(3) 
      C =A(4) 
      IF(K0.LE.SMLNO)THEN 
         PRINT *,  'INVALID VALUE FOR K0 --CONVER ABORTING' 
         RETURN  
      ENDIF 
      IF(K1.LT.0.0D0) K1=-K1 
***   GET INITIAL ESTIMATES FOR R & THETA *** 
      CALL RTHETA(T,M,K0,R0,TH0,UNSTAB) 
       IF(UNSTAB.EQ.-1.0)GO TO 400
*** 
***   NOW GET ACTUAL R & THETA *** 
*** 
      DO 200 K=1,ITRMX2 
C      WRITE(6,9875)T,AA,C,R0,BETA,DELTA,TH0 
9875   FORMAT(7E15.7)
         IF(R0.LE.0.0D0) GO TO 200
         P=T+AA*C*(R0**(BETA*DELTA))*TH0*(1.D0-TH0**2) 
         QQ=(M-AA*C*(R0**(1.D0-ALPHA))* 
     >      (K0*S0(TH0)+K1*S1(TH0)*R0**WEG))/ 
     >                     (1.D0+K1*(R0**WEG)/K0) 
         CALL RTHETA(P,QQ,K0,R1,TH1,UNSTAB) 
         IF(UNSTAB.EQ.-1.0)GO TO 400
         IF(DABS(R1-R0)/(1.D0+R0).LT.1.D-12.AND. 
     >      DABS(TH1-TH0)/(1.D0+DABS(TH0)).LT.1.D-12)THEN 
       R=(R1+R0)/2.D0 
             THETA=(TH1+TH0)/2.D0 
             RETURN
         ELSE 
            TH0=TH1 
      R0=R1 
         ENDIF 
200   CONTINUE 
      PRINT *,'CAN''T GET ROUTINE TO CONVERGE FOR R & THETA' 
      PRINT *,'WHEN T=',T,'   M=',M,'   RINIT=',R0,'   THETAINIT=',TH0 
      WRITE(6,150)CERRTK,CERRHO,CERRP,CERRL
150   FORMAT(' ','ORIG CALL ARGS: TK,RHO,P,L: ',3F15.8,F3.0)
      WRITE(6,151)
151   FORMAT(' ','(IN UNSTABLE SUPERSAT REGION)')
400      CONTINUE
         R=DABS(R0)
         THETA=TH0
300   RETURN  
      END 



      SUBROUTINE RTHETA(T,M,K0,R,THETA,UNSTAB) 
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON/CERR/CERRTK,CERRHO,CERRP,CERRL,CERLS
      DOUBLE PRECISION M,K0 
      PARAMETER (EPS=1.0D-12) 
      PARAMETER (BETA=0.325D0,BB=1.3757D0) 
****  THIS SUBROUTINE CONVERTS T AND M INTO THE PARAMETRIC 
***   COORDINATES R & THETA.  THE PARTICULAR MODEL USED HERE 
***   IS THE LINEAR MODEL OF SCHOFIELD; THE REFERENCE USED IS 
***   M.R. MOLDOVER, J.RES.NAT.BUR.STD. 83, 329 (JULY-AUGUST 1978). 
***   NOTE THAT THIS ROUTINE MAY BE USED ITERATIVELY TO COMPUTE 
***   MODELS WHICH DEVIATE SLIGHTLY FROM THE SYMMETRIC LINEAR CASE. 
***   T IS THE REDUCED TEMPERATURE 
***   M IS THE ORDER PARAMETER 
***   K0 IS THE NON-UNIVERSAL CONSTANT APPEARING IN THE LINEAR MODEL 
***   R IS THE PARAMETRIC VARIABLE DESCRIBING THE DISTANCE FROM 
*     THE CRITICAL POINT 
***   THETA IS THE PARAMETRIC VARIABLE SPECIFYING THE 
*     THERMODYNAMIC PATH 
***   "*" IS A STATEMENT NUMBER IN THE CALLING PROGRAM TO WHICH THE 
*     RTHETA ROUTINE SHOULD RETURN TO IF AN ERROR OCCURRS. 
      DABST=DABS(T) 
      DABSM=DABS(M) 
      B=DSQRT(BB) 
      IF(K0.LT.0.0D0)K0=-K0 
      IF(BB.LE.1.0D0)THEN 
         PRINT *,'B IS < OR = 1.0; RTHETA ROUTINE ABORTING' 
         RETURN 
      ENDIF 
***   IF THE SYSTEM IS ON THE CRITICAL ISOCHORE 
      IF(DABSM.LT.EPS)THEN 
         IF(DABST.LT.EPS)THEN 
            THETA=0.0D0 
            R=0.0D0 
            RETURN 
         ELSE IF (T.LE.-EPS) THEN 
            THETA=DSIGN(1.0D0,M) 
            R=T/(1.D0-BB) 
            RETURN 
         ELSE IF (T.GE.EPS) THEN 
            THETA=0.0D0 
            R=T 
            RETURN 
         ENDIF 
***   IF, HOWEVER, WE HAVE A NON-ZERO ORDER PARAMETER 
      ELSE 
**       CRITICAL ISOTHERM 
         IF(DABST.LT.EPS)THEN 
            THETA=DSIGN(1.D0,M)/B 
            R=(M/(K0*THETA))**(1.0D0/BETA) 
C           WRITE(6,9897)M,R,THETA
9897        FORMAT(' M,R,THETA: ',4E20.8)
            RETURN 
         ELSE IF (T.LE.-EPS)THEN 
**          CHECK FOR 2-PHASE REGION 
             ZOLD=(-T*(K0/DABSM)**(1.D0/BETA))
             ZOLD=  DSQRT(ZOLD+1.0D0)*1.00234/1.07
             IF(ZOLD.GT.B) THEN
*BYPASS THIS *       IF((-T*(K0/DABSM)**(1.D0/BETA)).GT.(BB-1.0D0))THEN 
               THETA=DSIGN(1.0D0,M) 
               R=T/(1.D0-BB) 
               UNSTAB=-1.0
               RETURN 
            ENDIF 
***         GET INITIAL ESTIMATE BY LINEAR INTERPOLATION 
***         Z=B*THETA 
            ZINIT=1.D0-(1.D0-B)*T*(K0/DABSM)**(1.D0/BETA)/(1.D0-BB) 
         ELSE IF (T.GE.EPS)THEN 
***         GET INITIAL ESTIMATE BY RATIONAL INTERPOLATION 
***      Z=B*THETA 
         ZINIT=(1.D0+T*(K0/(B*DABSM))**(1.D0/BETA))**(-BETA) 
      ENDIF 
      C=-M*B/(K0*DABST**BETA) 
      Z=DSIGN(ZINIT,M) 
***   USE NEWTON-RAPHSON METHOD TO GET EXACT Z 
C
C      WRITE(6,9898)Z,C,T,M
9898   FORMAT(3X,'Z,C,T,M: ',6E20.8)
C
      DO 10 I=1,20 
         DZ=(1.D0-Z**2)*(Z+C*DABS(1.D0-Z**2)**BETA)/ 
     >((2.D0*BETA-1.D0)*Z**2+1.D0) 
      Z=Z-DZ 
      IF(DABS(DZ)/(1.D0+DABS(Z)).LT.EPS)THEN 
      R=T/(1.0D0-Z**2) 
      THETA=Z/B 
**BYPASS**    IF(DABS(THETA).GT.(1.0D0-EPS))THETA=DSIGN(1.D0,THETA) 
      RETURN 
      ENDIF 
10    CONTINUE 
      ENDIF 
      PRINT *,' RTHETA FAILED TO CONVERGE' 
      PRINT *,' T=',T,'   M=',M,' DZ=',Z 
      WRITE(6,150)CERRTK,CERRHO,CERRP,CERRL
150   FORMAT(' ','ORIG CALL ARGS: TK,RHO,P,L: ',3F15.8,F3.0)
      WRITE(6,151)
151   FORMAT(' ','(IN UNSTABLE SUPERSAT REGION)')
       UNSTAB=-1.0D0
       R= DABS(R)
      RETURN 
      END 



      SUBROUTINE CONST 
*     THIS ROUTINE CALCULATES THE VARIOUS CONSTANTS NEEDED 
*     FOR THE NEW REVISED CONVER. 
C
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PARAMS2/ ALPHA,BETA,GAMMA,DELTA,WEG,BB,P00 
     1,P02,P04,P10,P12,P14,S00,S20,S01,S21 
C
C
       ALPHA= 0.1085D0
       BETA= 0.325D0
       GAMMA=1.2415D0
       DELTA= 4.82D0
       WEG= .5D0
       BB= 1.3757D0
C
      P00=(BETA*(DELTA-3.D0)-BB*ALPHA*GAMMA)/ 
     >(2.D0*(BB**2)*(2.D0-ALPHA)*(1.D0-ALPHA)*ALPHA) 
      P20=-(BETA*(DELTA-3.D0)-BB*ALPHA*(2.D0*BETA*DELTA-1.D0))/ 
     >(2.D0*BB*(1.D0-ALPHA)*ALPHA) 
      P40=(2.D0*BETA*DELTA-3.D0)/(2.D0*ALPHA) 
      P01=(BETA*(DELTA-3.D0)-3.D0*WEG-BB*(ALPHA-WEG)*(GAMMA-WEG))/ 
     >(2.D0*(BB**2)*(2.D0-ALPHA+WEG)*(1.D0-ALPHA+WEG)*(ALPHA-WEG)) 
      P21=-(BETA*(DELTA-3.D0)-3.D0*WEG-BB*(ALPHA-WEG)* 
     >(2.D0*BETA*DELTA-1.D0))/ 
     >(2.D0*BB*(1.D0-ALPHA+WEG)*(ALPHA-WEG)) 
      P41=(2.D0*BETA*DELTA-3.D0)/(2.D0*(ALPHA-WEG)) 
      S00=(2.D0-ALPHA)*P00 
      S20=-BETA*(DELTA-3.D0)/(2.D0*BB*ALPHA) 
      S01=(2.D0-ALPHA+WEG)*P01 
      S21=-(BETA*DELTA-3.D0*(BETA+WEG))/ 
     >(2.D0*BB*(ALPHA-WEG)) 
      RETURN 
      END 




      BLOCK DATA BCOEFF 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /COEFFS/ A(20),Q(10) 
***   THIS PROGRAM INITIALIZES THE PARAMETERS NEEDED IN THE 
***   REVISED AND EXTENDED PARAMETRIC MODEL. 
C
C      NOTE THAT THE CONSTANTS HERE ARE FOR THE NEW,REVISED
C           CONVER, AND THE ORDER IS DIFFERENT FROM THAT IN
C           BLOCKDATA CLSBLK
C     (ACTUALLY, ONLY A(1) THRU A(4) ARE REFERRED TO IN NEW CONVER)
C
*     A(1)=A 
*     A(2)=K0 
*     A(3)=K1 
*     A(4)=C 
*     A(5)=P1 
*     A(6)=P2 
*     A(7)=P3 
*     A(8)=P11 
*     A(9)=DELTA TC 
*     A(10)=DELTA RHOC 
*     A(11)=DELTA PC 
*     A(12)=MU0 
*     A(13)=MU1 
*     A(14)=MU2 
*     A(15)=MU3 
*     A(16) - A(20) ARE UNUSED 
*     Q(1)=TC 
*     Q(2)=RHOC 
*     Q(3)=PC 
*     Q(4)=DP/DT EVALUATED AT THE CRITICAL POINT IN ACTUAL UNITS 
*     Q(5)=DELTA PC CALCULATED FROM SHIFT IN TC AND DP/DT AT C.P. 
*     Q(6)=SLOPE OF THE DIAMETER IN THE CLASSICAL REGION 
*     Q(7) - Q(10) ARE UNUSED 
       DATA A/23.6666D0,1.4403D0,.2942D0,-.017762D0,6.8445D0,
     *        -25.4915D0,5.238D0,0.4918D0,0.0D0,0.0D0,
     *         0.0D0,-11.2326D0,-22.6547D0,-17.8876D0,-4.9332D0,
     *         5*0.0D0/
C
       DATA Q/647.067D0,322.778D0,22.046D0,.267D0,.0D0,
     *      -1.6D0,4*0.D0/
      END 




        SUBROUTINE SWVIR(TKD,BZ,CZ)
C
C       SUBROUTINE TO CALCULATE 
C             BZ:  SECOND VIRIAL COEFFICIENT
C             CZ:  THIRD VIRIAL  COEFFICIENT
C
C
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DOUBLE PRECISION C1(10,10),C2(10,12),C3(10,10)
        DOUBLE PRECISION C4(10,10)
        COMMON/CVCC/C1,C2,C3,C4,NT1,NT2,NT3,NT4
        COMMON/EFGVAL/EVAL,DR0,DT0,DELTA,AVAL,ABVAL,BVAL,CVAL,DVAL
        COMMON/NORMS/TC,RHOC,RVAL
C
C
        TBAR= -TC/TKD
        DTBAR= TBAR + 1.0D0
C
C       SUMMATION WHERE I=1 FOR B, I=2 FOR C;
C       COEFFICIENTS C1 CORRESPOND TO K=1
C                    C2 ~ K=2
C                    C3 ~ K=3
C
C
       BZ= 0.0
       CZ= 0.0
C
       DO 200 J=1,NT1
       BA=  EVAL * C1(1,J) * TBAR**(J-1)
       BZ= BZ+BA
       TBJ= TBAR**(J-1)
200    CONTINUE
C
       DO 220 J=1,NT2
       BB= C2(1,J) * TBAR**(J-1)
       TBJ= TBAR**(J-1)
       BZ= BZ +BB
       CZ= CZ + C2(2,J) * TBAR**(J-1)
220    CONTINUE
C
       EAVAL= (-AVAL*DTBAR-CVAL*DTBAR**2)
       DEXLIM= -35.
       DARG= EAVAL
       IF(DARG.GE.DEXLIM) EXPVAL=DEXP(DARG)
       IF(DARG.LT.DEXLIM) EXPVAL=0.0D0
       EAVAL= EXPVAL
       ECVAL= DEXP(-DVAL)
C
C      DO 230 J=1,NT3
C      BC= C3(1,J)*ECVAL*EAVAL*TBAR**(J-1)*TBAR**2
C      CZ= CZ + C3(2,J)*ECVAL*EAVAL*TBAR**(J-1)*TBAR**2
C      BZ= BZ+BC
C      TBJ= TBAR**(J-1)
C230    CONTINUE
C
       H= DEXP(-4.0D0*(TBAR+3.0D0))
C
       DO 240 J=1,NT4
       BD= (EVAL*H*C4(1,J))*TBAR**(J-1)
       CD= (EVAL*H*C4(2,J))*TBAR**(J-1)
       BZ= BZ+BD
       CZ=CZ+CD
240    CONTINUE
C
       BZ= BZ / RHOC/1000.
       CZ= CZ /(RHOC*1000.)**2
       RETURN
       END





C
C   ********************
C
C    SAUL,WAGNER (1985) H2O SATURATION FUNCTIONS
C
C     REF:  WAGNER, W, SENGERS,J
C           DRAFT RELEASE ON SATURATION PROPERTIES
C           TECHNICAL REPORT BN 1051
C
C
C                       [FPWAG,FDPWAG,RFWAG,RGWAG,FHFWAG,FHGWAG,FAWAG]
C
C   ********************************************************
C
       DOUBLE PRECISION FUNCTION FPSWAG(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Y)
       TC= 647.14D0
       IF(TK.GE.TC) THEN
                    FPSWAG=22.064D0
                    WRITE(6,999)TK
                    RETURN
999                 FORMAT(' FPSWAG: TK GE TC ',F20.10)
                    ENDIF
       TAU=  1.0D0 - TK/TC
       W= -7.858230D0 * TAU
       W= W + 1.839910D0 * TAU**1.5
       W= W + (-11.78110D0 *TAU**3)
       W= W + (22.6705D0 * TAU**3.5)
       W= W + (-15.9393D0 * TAU**4)
       W= W + (1.77516D0 * TAU**7.5)
       DEXLIM= -35.
       DARG= TC/TK*W
       IF(DARG.GE.DEXLIM) W=DEXP(DARG)
       IF(DARG.LT.DEXLIM) W=0.0D0
       FPSWAG= 22.064D0 * W
       RETURN
       END



       DOUBLE PRECISION FUNCTION FDPSW(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Y)
       PS= FPSWAG(TK)
       TC= 647.14
       IF(TK.GE.TC) THEN
                    FDPSW= 22.064D0
                    WRITE(6,999) TK
                    RETURN
999                FORMAT(' FDPSW: TK GT TC ',F20.10)
                    ENDIF
       TAU=  1.0 - TK/TC
       W= -7.858230 
       W= W + 1.839910 * TAU**0.5*1.5
       W= W + (-11.78110 *TAU**2*3.0)
       W= W + (22.6705 * TAU**2.5*3.5)
       W= W + (-15.9393 * TAU**3*4.0)
       W= W + (1.77516 * TAU**6.5*7.5)
       W= DLOG(PS/22.064D0) + W
       FDPSW= W *(-PS/TK)
       RETURN
       END



       DOUBLE PRECISION FUNCTION FRFWAG(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Y)
       TKC= 647.14
       RC= .322
       IF(TK.GE.TKC)RETURN
C
       TAU= 1.0D0 - TK/TKC
       W= 1.99206 *TAU**(1.0/3.0)
       W= W + 1.10123* TAU**(2.0/3.0)
       W= W + (-0.512506 * TAU**(5./3.))
       W= W + (-1.75263 * TAU**(16./3.))
       W= W + (-45.4485 * TAU**(43./3.))
       W= W + (-675615. * TAU**(110./3.))
       FRFWAG= (W+1.0)*RC
       RETURN
       END



       DOUBLE PRECISION FUNCTION FRGWAG(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Y)
       TC= 647.14
       IF(TK.GE.TC) RETURN
       TAU= 1.0 - TK/TC
       W=  -2.029570 * TAU**(2.0/6.0)
       W=  W + (-2.68781*TAU**(4./6.))
       W= W +  (-5.381070*TAU**(8./6.))
       W= W +  (-17.31510*TAU**(18./6.))
       W= W +  (-44.6384 *TAU**(37./6.))
       W= W +  (-64.3486 *TAU**(71./6.))
       DEXLIM= -35.
       DARG= W
       IF(DARG.GE.DEXLIM) W=DEXP(DARG)
       IF(DARG.LT.DEXLIM) W=0.0D0
       FRGWAG= .322 * W
       RETURN
       END




       DOUBLE PRECISION FUNCTION FBZLNR(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C      REF:  LEFVRE,E, NIGHTINGALE,M, ROSE, J
C            J.MECH.ENG.SCI. VOL 17,NO 5 (1975)
C
       B= 1.0D0/(1.0D0+TK/10000.D0)*.0015D0
       DARG= -1500.D0/TK
       DEXLIM= -35.
       IF(DARG.GE.DEXLIM) EXPVAL=DEXP(DARG)
       IF(DARG.LT.DEXLIM) EXPVAL=0.0D0
       C= (1.D0-EXPVAL)**(5.D0/2.D0)
       C= C*DSQRT(TK/1500.D0)*DEXP(1500.D0/TK)*(-.000942D0)
       D= 1500.D0/TK * (-.0004882D0)
       FBZLNR= B + C +D
       RETURN
       END






C
C  ********************************************
C
C      SPINODAL RHO (LIQ AND VAP) CORRELATIONS
C         OF WSTEAM TK AND RHO VALUES WHERE DPDR=0.0 
C         FOR USE IN ITERATION TO RHO FROM PRESSURE 
C         VALUE ARGUMENTS IN THE SUPERSATURATED
C         REGION, 
C
C         IN ADDITION THESE FUNCTIONS ARE USED TO 
C         CHECK INPUT PRESSURE VALUES (L=3 OR L=4)
C         TO GIVE AN ERROR MESSAGE WHEN THE RHO
C         FOR A GIVEN PRESSURE LIES IN THE UNSTABLE
C         SUPERSATURATED REGION 
C
C  ***********************************************
       DOUBLE PRECISION FUNCTION WSPINF(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       TKC= 647.067D0
       RC=  .322778D0
       WSPINF=0.0D0
      IF(TK.GE.TKC)RETURN
C
       TAU= 1.0D0 - TK/TKC
       W= -.5948463865D0 *TAU**(1.0D0/3.0D0)
       W= W + .5080947833D2* TAU**(2.0D0/3.0D0)
       W= W + (-.3797642336D3 * TAU**(3.D0/3.0D0))
       W= W + (.1249347735D4 * TAU**(4.D0/3.0D0))
       W= W + (-.1860121245D4 * TAU**(5.D0/3.0D0))
       W= W + (.1031487945D4 * TAU**(6.D0/3.0D0))
       W= W + (-.1574103153D4 * TAU**(17.D0/3.0D0))
       W= W + (.9435352115D4  * TAU**(21.D0/3.0D0))
       W= W + (-.8168784707D4 * TAU**(22.D0/3.0D0))
       WSPINF= (W+1.0D0)*RC
       RETURN
       END
       DOUBLE PRECISION FUNCTION WSPING(TK)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       TC= 647.067D0
       RC=  .322778D0
       WSPING=0.0D0
       IF(TK.GE.TC) RETURN
       TAU= 1.0 - TK/TC
       W=  -.3346401461D1 * TAU**(1.0D0/3.0D0)
       W=  W + (.1268114311D2*TAU**(2.D0/3.0D0))
       W= W +  (-.2316813580D2*TAU**(3.D0/3.0D0))
       W= W +  (-.6567012350D2*TAU**(4.D0/3.0D0))
       W= W +  (0.2166943752D3 *TAU**(5.D0/3.0D0))
       W= W +  (-.1663688422D3 *TAU**(6.D0/3.0D0))
       W= W +  (.1409844499D3 *TAU**(16.D0/3.0D0))
       W= W +  (-.1501900170D4 *TAU**(26.D0/3.0D0))
       W= W +  (.1722928004D4   *TAU**(29.D0/3.0D0))
       DEXLIM= -35.
       DARG= W
       IF(DARG.GE.DEXLIM) W=DEXP(DARG)
       IF(DARG.LT.DEXLIM) W=0.0D0
       WSPING= RC * W
       RETURN
       END




       BLOCKDATA CHIBLK
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DOUBLE PRECISION C1(10,10),C2(10,12),C3(10,10),C4(10,10)
       DOUBLE PRECISION C0(8)
        COMMON/CZERO/C0
       COMMON/CCREAD/CREAD
       COMMON/CINIT/C000,C001
       COMMON /CVCC/C1,C2,C3,C4,NT1,NT2,NT3,NT4
       DATA CREAD/0.0D0/

       DATA C0/.707501277532D1,-.834240573672D1,-.36460138D0
     *         ,-.36897043D-1,0.30338153D-2,
     *         0.39010933D-3,0.11359287D0,.24131785D1/

       DATA C000/-.3474470158D-4/
       DATA C001/-.3040231678D-4/

        DATA C1/
     &  3*0.D0,0.3384249125D+0,-.7153393406D-1,0.D0,.549368081D-3
     &  ,3*0D0,
     &  2*0.D0,.4933218501D-1,2*0.D0,-.2328491212D-1,.2402095181D-2
     &  ,3*0.D0,
     &  .7529422956D+0,0.D+0,-.2280260070D+1,7*0.D+0,
     &  0.D+0,.1142004144D+1,-.2619059624D+1,0.D+0,.4395237702D+0,
     &  -.3161046646D-1,.6814467692D-3,3*0.D+0,
     &  -.3924227294D+0,0.D+0,-.2738770648D+0,7*0.D+0,
     &  3*0.D+0,-.1943443857D-01,.3048860434D-2,5*0.D+0,
     &  2*0.D+0,.3946510403D-2,7*0.D+0
     &  ,30*0.D0/

       DATA C2/
     &  .22436101314D0,.1193250201D0,2*0.D0,.6582959348D-1,
     &   5*0.D0,
     &  4*0.D0,.1651430628D0,5*0.D0,
     &  -.2178969357D1,.2674090542D0,.8647490995D0,7*0.D0,
     &  -.1530432257D0,0.D0,.2059881454D1,2*0.D0,
     &  -.4888628703D0,.1375328753D0,3*0.D0,
     &  4*0.D0,-.9015180666D0,-.1444258609D0,.1558046279D0,3*0.D0
     &  ,3*0.D0,-.2740652563D1,0.D0,.4983771706D0,4*0.D0,
     &  3*0.D0,-.3261978564D1,.1609338784D1,5*0.D0,
     &  .3484674963D-1,-.1537646434D1,2*0.D0,.2316225257D0,5*0.D0,
     &  0.D0,-.1419249232D1,.7969984635D0,7*0.D0,
     &  4*0.D0,.7510544627D-2,5*0.D0,
     &  10*0.D0,
     &	.5364384732D-3,9*0.D0/

        DATA C3/
     & .6109381296D0,0.D0,-.1906644459D-1,0.D0,.7976092188D-2,
     & 5*0.D0,
     & .1934466766D1,9*0.D0,
     & .1921820547D1,0.D0,-.4410105919D-1,7*0.D0,
     & .6130354419D0,-.2855258689D0,2*0.D0,.2526137080D-1,5*0.D0,
     & 0.D0,-.2374074642D+0,0.D+0,.3855866402D-1,.8041672150D-2,
     & 55*0.D0/

       DATA C4/
     & 2*0.D0,-.1635439033D2,7*0.D0,
     & -.5025818675D2,9*0.D0,
     & 10*0.D0,
     & 0.D0,.1649003040D0,8*0.D0,
     & -.8499893502D0,9*0.D0,
     & 30*0.D0,
     & 0.8314382544D-2,.8781327858D-3,8*0.D0,
     & 0.D0,.1537391213D-2,-.9016873786D-3,0.D0,.3326628664D-3,
     & 5*0.D0/
       END


