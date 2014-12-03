c****************************************************************
c                 fonction qui calcule la composante de
c                     la contrainte tangentielle 
c                       sigma0xx  sur la surface
c                            de la Bulle
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Octobre 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function sigma0xxB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                   Un,dUn) 
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: c, zetap 
      double precision, intent(in) :: xi, alpha0
      double precision, dimension(Nt), intent(in) :: alpha           
      double precision, dimension(Nt+3), intent(in) :: Pn
      double precision, dimension(Nt+1), intent(in) :: Un  
      double precision, dimension(Nt+1), intent(in) :: dUn   

c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2
      double precision :: aux3
      double precision :: xn, xNt, xi2, sh, ch
      double precision :: SigmaEtaB
     
c     Declaration des variables en sorties
c     ------------------------------------      
       double precision :: sigma0xxB

c     ==================================================================       
c     on calcule la composante de la contrainte tangentielle

c     soit sigmaxxB = -P + Tau0xxB
c
c      Tau0xxB(zeta,x) = (ch-x)^(-1/2)/(2c³(1-x²)) 
c                        { Vn(1+x²-2xch) [2(ch-x)U'n- 3sh*Un ]
c                          +2(ch-x)(1-x²)V'n [sh*Un-2(ch-x)U'n] 
c                        }
c
c     ou bien suivant les Pn

c      Tau0xxB(zeta,x) = (ch-x)^(-1/2)/(2c³(1-x²)) 
c                        { [Pn-1 - Pn+1](1+x²-2xch) [2(ch-x)U'n- 3sh*Un]
c                          +2(ch-x)[nxPn-1 -nPn -(n+2)xPn+1 +(n+2)Pn+2]
c                           *[sh*Un-2(ch-x)U'n] 
c                        }
c
c     De plus, en coordonneees bispherique la pression est

c      P = (ch-x)^(-1/2)/(c³) sum{alphan*ch[(n+1/2)zeta]*Pn} 

c     avec les betan = 0 car les An=Cn=0 en raison des conditions 
c     limites sur la particule

c     Donc la contrainte totale est 
c
c      sigmaxxB = (ch-x)^(-1/2)/c³ 
c                  sum{ -alphan*ch[(n+1/2)zeta]*Pn   
c                       + (1/(2(1-x²)))
c                          *([Pn-1 - Pn+1](1+x²-2xch) [2(ch-x)U'n- 3sh*Un]
c                          +2(ch-x)[nxPn-1 -nPn -(n+2)xPn+1 +(n+2)Pn+2]
c                           *[sh*Un-2(ch-x)U'n] 
c                           )
c                     }
c                       
cc     Sachant que sur la bulle zeta = -zetap alors
c
c      sigmaxxB(-zetap,x) = (ch-x)^(-1/2)/c³ 
c                  sum{ -alphan*ch[(n+1/2)zetap]*Pn   
c                       + (1/(2(1-x²)))
c                          *([Pn-1 - Pn+1](1+x²-2xch) [2(ch-x)U'n +3sh*Un]
c                          +2(ch-x)[nxPn-1 -nPn -(n+2)xPn+1 +(n+2)Pn+2]
c                           *[-sh*Un -2(ch-x)U'n] 
c                           )
c                     }
       ! mauvaise idee de faire sh(-x)=-sh(x), on garde l'expression
       ! precedente
c                       
c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=dsinh(zetap)

c     On initialise SigmaEta_barre à O
c     --------------------------------      
      SigmaEtaB=0.d0
 
c     Calcul de la somme SigmaEta_barre
c     ----------------------------------
      aux=0.d0
      aux0=0.d0
      aux1=0.d0    
      aux2=0.d0   
      aux3=0.d0
c     on sait que Pn0=1.D0 
c     quand n= 0 U'0(-zetap) = 0 et U0(-zetap) = 0
c      Tau0xxB = 0
c     quand n= 1,  
c      Tau0zzB = (1/(2(1-x²)))*{ (1-P(2))*(1+x²-2xch)
c                                 *[2(ch-x)U'1 - 3sh*U(1)] 
c                                +2*(ch-x)*[x-P(1)-3*x*P(2)+3P(3)]
c                                *[sh*U1 -2(ch-x)U'1]  }
      xi2=xi*xi
      aux2=1.d0/(2.d0*(1.d0-xi2))      

      !Test (verif avec Mathematica OK)
c$$$        n=1
c$$$      aux1=aux2*((1.d0-Pn(2))*(1.d0+xi2-2.d0*xi*ch)
c$$$     &     *(2.d0*(ch-xi)*dUn(1)-3.d0*sh*Un(1))
c$$$     &     +2.d0*(ch-xi)*(xi-Pn(1)-3.d0*xi*Pn(2)+3.d0*Pn(3))
c$$$     &     *(-2.d0*(ch-xi)*dUn(1)+sh*Un(1)) )
c$$$       
c$$$      do 51 n=2,Nt 
c$$$       xn=dble(n)    
c$$$       aux=aux+aux2*((Pn(n-1)-Pn(n+1))*(1.d0+xi2-2.d0*xi*ch)
c$$$     &     *(2.d0*(ch-xi)*dUn(n)-3.d0*sh*Un(n))
c$$$     &     +2.d0*(ch-xi)*(xn*xi*Pn(n-1)-xn*Pn(n)
c$$$     &     -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
c$$$     &      *(-2.d0*(ch-xi)*dUn(n)+sh*Un(n)))
c$$$ 51   continue  

c     Ligne pour n=0
c     --------------
      aux0=alpha0*dcosh((1.d0/2.d0)*zetap)

c     Ligne pour n=1
c     --------------
      aux1=alpha(1)*Pn(1)*dcosh((3.d0/2.d0)*zetap)
     &     +aux2*((1.d0-Pn(2))*(1.d0+xi2-2.d0*xi*ch)
     &     *(2.d0*(ch-xi)*dUn(1)-3.d0*sh*Un(1))
     &     +2.d0*(ch-xi)*(xi-Pn(1)-3.d0*xi*Pn(2)+3.d0*Pn(3))
     &     *(-2.d0*(ch-xi)*dUn(1)+sh*Un(1)) )

c     Ligne pour n=2,Nt
c     -------------------
      do 51 n=2,Nt
       xn=dble(n)    
       aux=aux+Pn(n)*alpha(n)*dcosh((xn+(1.d0/2.d0))*zetap)
     &     +aux2*((Pn(n-1)-Pn(n+1))*(1.d0+xi2-2.d0*xi*ch)
     &     *(2.d0*(ch-xi)*dUn(n)-3.d0*sh*Un(n))
     &     +2.d0*(ch-xi)*(xn*xi*Pn(n-1)-xn*Pn(n)
     &     -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
     &      *(-2.d0*(ch-xi)*dUn(n)+sh*Un(n)))
 51   continue  
    
      SigmaEtaB=aux+aux1+aux0

c     scaling par  (ch-x)^{-1/2}/c^3 
c     ---------------------------------
      sigma0xxB=-(1.d0/(c*c*c*dsqrt(ch-xi)))          
     &           *SigmaEtaB

c      write(8,*) xi,sigma0xxB,SigmaEtaB

c     =============================
c     Fin de la fonction sigmaBerdan
c     ==============================
      end function sigma0xxB

