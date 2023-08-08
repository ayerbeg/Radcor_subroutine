*     RADIATE
      SUBROUTINE radiate(E,TH,Tatmp,Tbtmp,Ztmp,Atmp,W,SIGNR,SIGRAD)
C     DOES NOT INCLUDE CONTRIBUTION FROM ELASTIC SCATTERING
C     Modified to use Simpson method of integration 5/14/15
C-----K. Slifer. 09/16/02
C
C     Rewrote subroutine to include external bremsstrahlung using
C     formalism from S. Stein et al. Phys. Rev. D 12 7. Equation (A82)
C     Where possible the equation number is given with each expression.
C
C     This subroutine was created from the originla nradcor.f code
C     by Ryan Zielinski and David Ruth "David  <David.Ruth@unh.edu>
C
C     For this version, small changes were done by
C     Carlos Ayerbe <gayoso@jlab.org> indicated by (CA)
C----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 SIGNR,SIGRAD,DEL,PREC,W,W0,Z,A,SPENCE
      REAL*8 xia,xib,ES_MIN,ES_MAX
      REAL*8 PI,THR,XMT,ARG
      REAL*8 QMS,D1,tr,D2,xF,Tpb,Tpa
      REAL*8 eta,xi,R,TERM1,TERM2,TERM3,TERM4
      REAL*8 X1LO,X1HI,X2LO,X2HI,ANS_Es,ANS_Ep
      REAL*8 E,E0,TH,TH0,Ta,Tb,EP_MIN,EP_MAX
      REAL*8 Tatmp,Tbtmp,Ztmp,Atmp
      REAL*8 ALPH,EMASS,xb,XM
      REAL*8 hbc2,HBARC
      REAL*8 XSSIGNR, FEP,FES
      integer NSP,KF,NUM
      EXTERNAL FEP,FES
      COMMON/PAR/E0,TH0,W0,Z,A,SPENCE,Ta,Tb
      COMMON/CONST/ALPH,EMASS,hbc2,HBARC,XM
C     the cross section sent as parameter needs to sent to two functions
C     here, so we need an auxiliary common variable to avoid confussion (CA)
      COMMON XSSIGNR 


!     Definition of constants
C     these constants were defined common and were defined in
C     the main program      

      ALPH      = 1./137.035999  
      EMASS     = 0.511          
      hbc2      = 0.389379        ! GeV^2-mbarn
      HBARC     = 197.327053      ! MeV-FM
      XM        = 931.494         ! Mass of the nucleon (u) (MeV)

      
      DEL   = 10
      NUM   = 80

      Ta=Tatmp
      Tb=Tbtmp
      Z=Ztmp
      A=Atmp
      TH0=TH
      E0=E
      W0=W

      XMT   = A*XM   ! Mass of the target
      PI    = ACOS(-1.)
      THR   = TH    ! angle in radians
      ARG   = COS(THR/2.)**2
      SPENCE= PI**2/6.-LOG(ARG)*LOG(1.-ARG)


      DO NSP=1,50
         SPENCE = SPENCE-ARG**NSP/FLOAT(NSP)**2
      END DO

        
      QMS   = 4.*E*(E-W)*SIN(THR/2.)**2
      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )    ! (A46)
      xb    = 4./3. *(1.+1./9.*((Z+1.)/(Z+eta) )/(LOG(183.*Z**(-1./3.)))  ) 
      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      tr=D1/2./xb

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-17./36. ! this term dominates D2
      D2 = D2 +0.5*(PI**2/6.-SPENCE)
      D2 = D2 -1./4.*( LOG( E/(E-W) ) )**2        ! Correct. to peak. appr.
      D2 = D2*(2.*ALPH/PI)                 
      D2 = D2+0.5772*xb*(Tb+Ta)                   ! Here D2= F-1
      xF = (1.+D2)                                ! (A44)

      Tpb = tr + Tb
      Tpa = tr + Ta
      
      XSSIGNR = SIGNR ! just defininig the value of the commom cross section
      
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) ) ! (A83)
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )              ! (A52)
      SIGRAD = SIGNR * xF
      SIGRAD = SIGRAD*( (R*DEL/E  )**(xb*Tpb) )
      SIGRAD = SIGRAD*( (DEL/(E-W))**(xb*Tpa) )
      SIGRAD = SIGRAD*(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )

      TERM1=(R*DEL/E  )**(xb*Tpb)
      TERM2=(DEL/(E-W))**(xb*Tpa)
 1    TERM3=(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )
      TERM4=xF
C
C-----Stein's 1st integral wrt dEs' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X1LO   = (E-W)*( XMT/( XMT-2.*(E-W)*(SIN(THR/2.))**2) -1.0 )
      X1HI   = W-R*DEL

      ES_MIN  = (E-W)/(1.-2.*(E-W)*SIN(THR/2.0)*SIN(THR/2.0)/XMT)    
      ES_MAX  = E -R*DEL

      ANS_Es = 0.
      IF (X1HI.GT.X1LO) THEN
          CALL SIM(ES_MIN,ES_MAX,DEL,NUM,FES,ANS_Es)
      ELSE

      ENDIF

      ANS_Es = ANS_Es
C
C-----Stein's 2nd integral wrt dEp' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X2LO   = E*( 1.0-XMT/( XMT+2.*E*(SIN(THR/2.))**2) )
      X2HI   = W-DEL

      EP_MIN  = E-W + DEL                                 
      EP_MAX  = E/(1.+2.*E*SIN(THR/2.0)*SIN(THR/2.0)/XMT) 

      ANS_Ep = 0.
      IF (X2HI.GT.X2LO) THEN
         CALL SIM(EP_MIN,EP_MAX,DEL,NUM,FEP,ANS_Ep)
      ELSE

      ENDIF
      ANS_Ep = ANS_Ep

      SIGRAD = SIGRAD + ANS_Es + ANS_Ep
 
   
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIM(FL,FU,DEL,N,F,ANS)
      IMPLICIT REAL*8(A-H,O-Z)
      NX = ((FU - FL)/DEL/2.)
      NX = 2*NX
      
      IF(NX.LT.N) NX=(N/2)*2            

      DT=(FU-FL)/FLOAT(NX)
      SI = F(FL)
      SI = SI + F(FU)
      XX = FL + DT
      SI = SI + F(XX) *4.
      AA=FL
      II=0
5     II=II+2

      IF (II.GE.NX)  GO TO 7

      AA =FLOAT (II)*DT
      AA=AA+FL
      BB=AA+DT
      SI = SI + 2.*F(AA) + 4.*F(BB)
      GO TO 5

7      ANS = DT * SI/3.0
8      FORMAT (' FL,FU,DEL,DT,N,NX=', 6(G11.4,1X), ' ANS ', G11.4)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION GAMMA(X)
      IMPLICIT REAL*8(A-H,O-Z)
       Y = X-1.

      GAMMA=1.-.57486*Y+.95123*Y*Y-.69986*Y**3
     . +.42455*Y**4-.10107*Y**5

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FES(VAR)
      IMPLICIT REAL*8(A-H,O-Z)
C     REAL *8 inelas_cxsn
      LOGICAL INTERNAL
      COMMON/PAR/E,TH,W,Z,A,SPENCE,Ta,Tb
      COMMON/CONST/ALPH,EMASS,hbc2,HBARC,XM
      COMMON XSSIGNR

      
      INTERNAL = .TRUE.
      ! Use S. Stein et al. Phys. Rev. D 12 7
        PI    = ACOS(-1.) 
        EI    = E
        EF    = E - W
        FMT   = A*XM
        THR = TH 
        ARG   = COS(THR/2.)**2
        xeta  = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
        xb    = 4./3. *(1.+1./9.*((Z+1.)/(Z+xeta) )/(LOG(183.*Z**(-1./3.)))  )
        xi    = (PI*EMASS/2./ALPH)*(Ta+Tb)
        xi    = xi/( (Z+xeta)*LOG(183.*Z**(-1./3.)) )
        R     = ( FMT+EI*(1.-COS(THR)) )/( FMT-(EF)*(1.-COS(THR)) )
        QMS1  = 4.*(VAR)*(EF)*SIN(THR/2.)**2    

        IF (INTERNAL) THEN
          TR1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)

          XF    = -14./9. + 13./12.*LOG(QMS1/EMASS**2)
          XF    = XF - 1./4.*( LOG( EI/(EF) ) )**2       ! Corr. to peak. approx.
          XF    = XF + 0.5*(PI**2/6.-SPENCE)
        ELSE
          TR1   = 0.0
          XF    = 0.0
        ENDIF

        XF      = 2.*ALPH/PI * XF + 1. +  0.5772*xb*(Tb+Ta)

        Tpb   = TR1 + Tb
        Tpa   = TR1 + Ta
      
        FES   = ( xb*Tpb/(EI-VAR) ) *phi((EI-VAR)/(EI))   
        FES   = FES + xi/(2.*(EI-VAR)**2)
        FES   = FES * XSSIGNR *XF
C     inelas_cxsn is the cross-section used in the original code
C     defined in the main code, and sent as parameter to the subroutine
C     as SIGNR, for some reason, it was called here, instead of using the
C     parameter sent to the subroutine (CA)       
C        FES   = FES * inelas_cxsn(VAR/1000.0,EF/1000.0,TH,Z,A)*XF            
        FES   = FES * ( (EI-VAR)/(EF*R) )**(xb*Tpa)    
        FES   = FES * ( (EI-VAR)/(EI)   )**(xb*Tpb)  
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CSL
      REAL*8 FUNCTION phi(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      phi=1.0-x+3./4.*x**2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FEP(VAR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL *8 inelas_cxsn
      LOGICAL INTERNAL
      COMMON/PAR/E,TH,W,Z,A,SPENCE,Ta,Tb
      COMMON/CONST/ALPH,EMASS,hbc2,HBARC,XM
      COMMON XSSIGNR
      
CSL ORIGINAL CODE
      INTERNAL = .TRUE.
      !S. Stein et al. Phys. Rev. D 12 7
        PI    = ACOS(-1.)
        EI    = E
        EF    = E - W
        FMT   = A*XM  
        THR = TH 
        ARG   = COS(THR/2.)**2

        xeta  = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
        !xb   = 4./3.  ! A45
        xb    = 4./3. *(1.+1./9.*((Z+1.)/(Z+xeta) )/(LOG(183.*Z**(-1./3.)))  )

        xi    = (PI*EMASS/2./ALPH)*(Ta+Tb)
        xi    = xi/( (Z+xeta)*LOG(183.*Z**(-1./3.)) )
        R     = ( FMT+EI*(1.-COS(THR)) )/( FMT-(EF)*(1.-COS(THR)) )

        QMS2  = 4.*EI*(VAR)*SIN(THR/2.)**2  


        IF (INTERNAL) THEN
          TR2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)

          XF    = -14./9. + 13./12.*LOG(QMS2/EMASS**2)
          XF    = XF - 1./4.*( LOG( EI/(EF) ) )**2       !KS. Correction to peak. approx.
          XF    = XF + 0.5*(PI**2/6.-SPENCE)
        ELSE
          TR2   = 0.0
          XF    = 0.0 
        ENDIF

        XF      = 2.*ALPH/PI * XF + 1.0 + 0.5772*xb*(Tb+Ta)

        Tpb   = TR2 + Tb
        Tpa   = TR2 + Ta
        FEP    = ( xb*Tpa/(VAR-EF) ) *phi((VAR-EF)/(VAR))
        FEP    = FEP + xi/(2.*(VAR-EF)**2)

C     inelas_cxsn is the cross-section used in the original code
C     defined in the main code, and sent as parameter to the subroutine
C     as SIGNR, for some reason, it was called here, instead of using the
C     parameter sent to the subroutine (CA)
C       FEP    = FEP * inelas_cxsn(EI/1000.0,VAR/1000.0,TH,Z,A)*XF
        FEP    = FEP * XSSIGNR *XF
        FEP    = FEP * ( (VAR-EF)/(VAR) )**(xb*Tpa)
        FEP    = FEP * ( (VAR-EF)*R/(EI))**(xb*Tpb)
        RETURN
      END
