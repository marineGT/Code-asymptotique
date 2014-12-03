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

      function dsigma0zxB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un,d3Un)                                             
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
      double precision, dimension(Nt+1), intent(in) :: d3Un      

c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2,aux3,aux4
      double precision :: terme1, terme2, terme3
      double precision :: terme4, terme5, aux7
      double precision :: aux5, aux6, xn2, xi2, xNt2
      double precision :: xn, xNt, sinEta, invsinEta
      double precision :: ch, sh , ch2, sh2
      double precision :: dSigmaB
     
c     Declaration des variables en sorties
c     ------------------------------------      
       double precision :: dsigma0zxB

c     ==================================================================       
c     on calcule la de rivee de la contrainte tangentielle
c
c          dsigmazx = (-(ch-x)^(-3/2)/(8*c^3*sin(eta)))
c                     { Vn*(3sh*Un(-4x²+sin(eta)+14xch-10ch²+sh²)
c                           -3(ch-x)U'n(1+2sin(eta)-4xhc+6ch²-3)
c                           +12(ch-x)²sh*U"n +8(ch-x)^3 U"'n)
c                      + 4sin(eta)(ch-x)² V"n (3sh*Un + 2(ch-x)U'n) 
c                     }

c     Soit suivant les Pn

c      dsigmazx = (-(ch-x)^(-3/2)/(8*c^3*sin(eta))) 
c                {[Pn-1 - Pn+1]*( 3sh*Un(-4x²+sin(eta)+14xch-10ch²+sh²)
c                                -3(ch-x)U'n(1+2sin(eta)-4xhc+6ch²-3)
c                                +12(ch-x)²sh*U"n +8(ch-x)^3 U"'n )
c                 +(4(ch-x)²/(1-x²)^ (3/2))
c                  *[-n(1 + x²(n+1))Pn-1 + n(2n+3)xPn
c                    +(2 -n² +6x² +5nx² +n²x²)Pn+1 -(11n +14 +2n²)xPn+2 
c                    +(6+5n+n²)Pn+3 ]*(3sh*Un + 2(ch-x)U'n)       
c                 }          
c                       
c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=dsinh(zetap)

      ch2=dcosh(zetap)*dcosh(zetap)
      sh2=dsinh(zetap)*dsinh(zetap)

c     On initialise dSigma_barre à O
c     --------------------------------      
      dSigmaB=0.d0
 
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
      aux5=(1.d0-xi**2.d0)
      aux6=(ch-xi)**(3.d0/2.d0)
      aux7=-4.d0*((ch-xi)**2.d0)/(aux5**(3.d0/2.d0))
      xi2=xi*xi
c     On test chaque termes l'un apres l'autre

c$$$      n=4
c$$$      xn=dble(n)   
c$$$      xn2=xn*xn
c$$$      !Calcul des 3 termes dans la boucle
c$$$      terme1=(Pn(n-1)-Pn(n+1))*(
c$$$     &       +3.d0*sh*Un(n)
c$$$     &       *(-4*xi2+sinEta+14*xi*ch-10*ch2+sh2) 
c$$$     &       -3.d0*(ch-xi)*dUn(n)*(1+2.d0*sinEta
c$$$     &       -4*xi*ch+6.d0*ch2-3.d0)
c$$$     &       +12.d0*((ch-xi)**2.d0)*sh*d2Un(n)
c$$$     &       +8.d0*((ch-xi)**3.d0)*d3Un(n))
c$$$
c$$$      terme2= aux7*(-(1.d0+xi2*(1.d0+xn))*xn
c$$$     +       *Pn(n-1)+xn*(2.d0*xn+3.d0)*xi*Pn(n)
c$$$     +       +(2.d0-xn2+6.d0*xi2+5.d0*xn*xi2+xi2*xn2)
c$$$     +       *Pn(n+1)-(11.d0*xn+14.d0+2.d0*xn2)*xi
c$$$     +       *Pn(n+2)+(6.d0+5.d0*xn+xn2)*Pn(n+3))
c$$$     &       *(3.d0*sh*Un(n)+2.d0*(ch-xi)*dUn(n))    
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
      aux1=(1.d0-Pn(2))*(
     &       +3.d0*sh*Un(1)
     &       *(-4*xi2+sinEta+14*xi*ch-10*ch2+sh2) 
     &       -3.d0*(ch-xi)*dUn(1)*(1+2.d0*sinEta
     &       -4*xi*ch+6.d0*ch2-3.d0)
     &       +12.d0*((ch-xi)**2.d0)*sh*d2Un(1)
     &       +8.d0*((ch-xi)**3.d0)*d3Un(1))
     &     +aux7*(-(1.d0+2.d0*xi2)+5.d0*xi*Pn(1)
     &     +(1.d0+12.d0*xi2)*Pn(2)-27.d0*xi*Pn(3)
     &      +12.d0*Pn(4))
     &      *(3.d0*sh*Un(1)+2.d0*(ch-xi)*dUn(1))
            
c     Ligne pour n=2,Nt-2
c     -------------------
      do 51 n=2,Nt+1 
       xn=dble(n)   
       xn2=xn*xn
      !Calcul des 2 termes dans la boucle
      terme1=(Pn(n-1)-Pn(n+1))*(
     &       +3.d0*sh*Un(n)
     &       *(-4*xi2+sinEta+14*xi*ch-10*ch2+sh2) 
     &       -3.d0*(ch-xi)*dUn(n)*(1+2.d0*sinEta
     &       -4*xi*ch+6.d0*ch2-3.d0)
     &       +12.d0*((ch-xi)**2.d0)*sh*d2Un(n)
     &       +8.d0*((ch-xi)**3.d0)*d3Un(n))

      terme2=aux7*(-(1.d0+xi2*(1.d0+xn))*xn
     +       *Pn(n-1)+xn*(2.d0*xn+3.d0)*xi*Pn(n)
     +       +(2.d0-xn2+6.d0*xi2+5.d0*xn*xi2+xi2*xn2)
     +       *Pn(n+1)-(11.d0*xn+14.d0+2.d0*xn2)*xi
     +       *Pn(n+2)+(6.d0+5.d0*xn+xn2)*Pn(n+3))
     &       *(3.d0*sh*Un(n)+2.d0*(ch-xi)*dUn(n))

       aux=aux+terme1+terme2
 51   continue     

      dSigmaB=aux+aux1

c     scaling par -(ch-x)^(-3/2)/(8*c^3*sin(eta))
c     -------------------------------------------
      dsigma0zxB=-(1.d0/(8.d0*c*c*c))
     &          /(aux6*sinEta)
     &          *dSigmaB                      

c      write(8,*) xi,dsigma0zxB,terme1+terme2

c     =============================
c     Fin de la fonction dsigma0zxB
c     ==============================
      end function dsigma0zxB


