c****************************************************************
c                 fonction qui calcule la composante de la 
c                   derivee de la contrainte tangentielle
c                        sur la surface de la Bulle 
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Octobre 2013
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function dsigma0zzB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un)                                             
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none                       
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: xi, c, zetap     
      double precision, dimension(Nt+3), intent(in) :: Pn
      double precision, dimension(Nt+1), intent(in) :: Un
      double precision, dimension(Nt+1), intent(in) :: dUn 
      double precision, dimension(Nt+1), intent(in) :: d2Un           

c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2,aux3,aux4
      double precision :: terme1, terme2, terme3
      double precision :: terme4, terme5, aux7
      double precision :: aux5, aux6, xn2, xi2, xNt2
      double precision :: xn, xNt, sinEta, invsinEta
      double precision :: ch, sh , ch2, ch3
      double precision :: dSigmazzB
     
c     Declaration des variables en sorties
c     ------------------------------------      
       double precision :: dsigma0zzB

c     ==================================================================       
c     on calcule la de rivee de la contrainte tangentielle
c
c          dsigmazz = (-(ch-x)^(-3/2)/(8*c^3*sin(eta)))
c                     { Un [ 3(2+sin(eta))*Vn*(-3+4xch-2ch2)
c                           +sin(eta)*V'n (-6x+(5+8x²)ch-10xch+3ch3)]
c                       +4sin(eta)(ch-x)² V'n (2sh*U'n+(ch-x)U"n)
c                       +2(ch-x)*Vn (-2sh*U'n+ (2+3sin(eta))(ch-x)U"n)
c                      }

c     Soit suivant les Pn

c      dsigmazz = (-(ch-x)^(-3/2)/(8*c^3*sin(eta))) 
c                { Un [ 3(2+sin(eta))*[Pn-1 - Pn+1]*(-3+4xch-2ch2)
c                      +(sin(eta)/(1-x²))*[nxPn-1 - nPn -(n+2)xPn+1 
c                        +(n+2)Pn+2]*(-6x+(5+8x²)ch-10xch+3ch3)]
c                  +4*(ch-x)²*(sin(eta)/(1-x²))*[nxPn-1 - nPn -(n+2)xPn+1 
c                        +(n+2)Pn+2]* (2sh*U'n+(ch-x)U"n)
c                  +2(ch-x)*[Pn-1 - Pn+1]*(-2sh*U'n+ (2+3sin(eta))(ch-x)U"n)
c                 }

                   
c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=dsinh(zetap)

      ch2=dcosh(2.d0*zetap)
      ch3=dcosh(3.d0*zetap)

      

c     On initialise dSigma_barre à O
c     --------------------------------      
      dSigmazzB=0.d0
 
c     Calcul de la somme dSigma_barre5
c     -------------------------------
       aux=0.d0
       aux0=0.d0
       aux1=0.d0
       aux2=0.d0
       aux3=0.d0
       aux4=0.d0
       aux5=0.d0
       aux6=0.d0
   
     
      sinEta = dsqrt((1.d0-xi)*(1.d0+xi)) 
      invsinEta=1.d0/sinEta
      aux2=(ch-xi)
      aux3=(1.d0-xi**2.d0)
      aux4=(ch-xi)**(3.d0/2.d0)
      aux5=sinEta/aux3
      aux6=aux5*((ch-xi)**2.d0)
      xi2=xi*xi
c     On test chaque termes l'un apres l'autre

c$$$      n=3
c$$$      xn=dble(n)   
c$$$      xn2=xn*xn
c$$$      !Calcul des 3 termes dans la boucle
c$$$      terme1=(1.d0/2.d0)*Un(n)*(
c$$$     &       +3.d0*(2.d0+sinEta)*(Pn(n-1)-Pn(n+1))
c$$$     &       *(-3.d0+4.d0*xi*ch-ch2)
c$$$     &       +aux5*(xn*xi*Pn(n-1)-xn*Pn(n)
c$$$     &       -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
c$$$     &       *(-6.d0*xi+(5.d0+8.d0*xi2)*ch-10.d0*xi*ch2
c$$$     &       +3.d0*ch3))
c$$$
c$$$      terme2=+8.d0*aux6*(xn*xi*Pn(n-1)-xn*Pn(n)
c$$$     &       -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
c$$$     &       *(2.d0*sh*dUn(n)+aux2*d2Un(n))
c$$$     &       +4.d0*aux2*(Pn(n-1)-Pn(n+1))*(-2.d0*sh
c$$$     &       *dUn(n)+(2.d0+3.d0*sinEta)*aux2*d2Un(n))
c$$$
c$$$      aux=aux+terme1+terme2     

c     on sait que Pn0=1.D0 et dUn(0)=0 (U'0(0))
c     quand n= 0
c      dSigma0zxSL = 1/sinEta {sinEta P(1) +[(3/2)sinEta - 1]*P(1)}
c     quand n= 1

c     On definit les Vn, Vn' et Vn"
     
c      Vn = Pn(n-1) - Pn(n+1)
c      (1-x)Vn' = (n+1)Pn+1 + nPn-1 -(2n+1)Pn
c      V"n = - [1/(x²-1)²]*{-n(1 + x²(n+1))Pn-1 + n(2n+3)xPn
c                           + (2 -n² +6x² +5nx² +n²x²)Pn+1 
c                           - (11n +14 +2n²)xPn+2 
c                           + (6+5n+n²)Pn+3 }


c     Ligne pour n=1
c     --------------    
      aux1=(1.d0/2.d0)*Un(1)*(
     &     +3.d0*(2.d0+sinEta)*(1.d0-Pn(2))
     &     *(-3.d0+4.d0*xi*ch-ch2)
     &     +aux5*(xi-Pn(1)-3.d0*xi*Pn(2)+3.d0*Pn(3))     
     &     *(-6.d0*xi+(5.d0+8.d0*xi2)*ch-10.d0*xi*ch2
     &      +3.d0*ch3)) 
     &      +8.d0*aux6*(xi-Pn(1)-3.d0*xi*Pn(2)
     &      +3.d0*Pn(3))*(2.d0*sh*dUn(1)+aux2*d2Un(1))
     &       +4.d0*aux2*(1.d0-Pn(2))*(-2.d0*sh
     &       *dUn(1)+(2.d0+3.d0*sinEta)*aux2*d2Un(1))     

            
c     Ligne pour n=2,Nt-2
c     -------------------
      do 51 n=2,Nt+1 
       xn=dble(n)   
       xn2=xn*xn
      !Calcul des 3 termes dans la boucle
      terme1=(1.d0/2.d0)*Un(n)*(
     &       +3.d0*(2.d0+sinEta)*(Pn(n-1)-Pn(n+1))
     &       *(-3.d0+4.d0*xi*ch-ch2)
     &       +aux5*(xn*xi*Pn(n-1)-xn*Pn(n)
     &       -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
     &       *(-6.d0*xi+(5.d0+8.d0*xi2)*ch-10.d0*xi*ch2
     &       +3.d0*ch3))

      terme2=+8.d0*aux6*(xn*xi*Pn(n-1)-xn*Pn(n)
     &       -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
     &       *(2.d0*sh*dUn(n)+aux2*d2Un(n))
     &       +4.d0*aux2*(Pn(n-1)-Pn(n+1))*(-2.d0*sh
     &       *dUn(n)+(2.d0+3.d0*sinEta)*aux2*d2Un(n))


       aux=aux+terme1+terme2
 51   continue     

      dSigmazzB=aux+aux1

c     scaling par -(ch-x)^(-3/2)/(8*c^3*sin(eta))
c     -------------------------------------------
      dsigma0zzB=-(1.d0/(4.d0*c*c*c*aux4*sinEta))
     &            *dSigmazzB                                   

c      write(8,*) xi,dsigma0zzB,dSigmazzB

c     =============================
c     Fin de la fonction dsigma0zzB
c     ==============================
      end function dsigma0zzB


