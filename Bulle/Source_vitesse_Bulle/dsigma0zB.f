c****************************************************************
c                 fonction qui calcule la composante de
c                     la contrainte normale sigma0zz
c                        sur la surface libre
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Novembre 2012.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function dsigma0zx(Nt,c,xi,alpha,alpha0,Pn,dPn,dPn0,dUn)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none                       
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: xi, c, alpha0, dPn0
      double precision, dimension(Nt), intent(in) :: alpha           
      double precision, dimension(Nt+1), intent(in) :: Pn
      double precision, dimension(Nt+1), intent(in) :: dPn
      double precision, dimension(Nt+1), intent(in) :: dUn      

c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2
      double precision :: xn, xNt, sinEta, invsinEta
      double precision :: dSigma_barre
     
c     Declaration des variables en sorties
c     ------------------------------------      
       double precision :: dsigma0zx

c     ==================================================================       
c     on calcule la composante de la contrainte tangentielle

c     soit sigmazzSL = -P + Tau0zzSL
c
c       Tau0zzSL=-(1-x)^(1/2)/(2c³sin(eta))
c                 U'n(0) {(1-x)sin(eta)V'n - Vn + (3/2)sin(eta) Vn }

c     ou bien suivant les Pn

cc       Tau0zzSL=-(1-x)^(1/2)/(2c³sin(eta))
c                 U'n(0) {sin(eta)[ (n+1)Pn+1 + nPn-1 - (2n+1)Pn ]
c                          =[(3/2)sin(eta)-1] (Pn-1 - Pn+1) }      
c
c     Donc la contrainte totale est 
c
c      sigmazzSL = (1-x)^(-1/2)/c³ sum{-alphan*Pn - U'n(0)/sin(eta) 
c                       {sin(eta)[ (n+1)Pn+1 + nPn-1 - (2n+1)Pn ]
c                          =[(3/2)sin(eta)-1] (Pn-1 - Pn+1) }   
c

c     On initialise dSigma_barre à O
c     --------------------------------      
      dSigma_barre=0.d0
 
c     Calcul de la somme dSigma_barre
c     -------------------------------
       aux=0.d0
       aux0=0.d0
       aux1=0.d0
       aux2=0.d0

      sinEta = (1.d0-xi)*(1.d0+xi) 
      invsinEta=1.d0/sinEta
c     on sait que Pn0=1.D0 et U'n(0)=2n+1, 
c     donc qd n= 1 alors dUn(0)=1 (U'0(0))
c     quand n= 0
c      Tau0zzSL = 1/sinEta {sinEta P(1) +[(3/2)sinEta - 1]*P(1)}
c     quand n= 1
c      Tau0zzSL = U'1(0)/sinEta { sinEta (2P(2)+1-3P(1)) 
c                                   +[(3/2)sinEta - 1]*(1-P(2)) }
c     Ligne pour n=0
c     --------------
      aux0=alpha0-invsinEta*(sinEta*Pn(1)
     &      +((3.d0/2.d0)*sinEta-1.d0)*Pn(1))

c     Ligne pour n=1
c     --------------
      aux1=alpha(1)*Pn(1)-invsinEta*dUn(1)*(
     &     +sinEta*(2.d0*Pn(2)+1.d0 -3.d0*Pn(1))
     &     +((3.d0/2.d0)*sinEta-1.d0)*(1.d0-Pn(2)))  
                  
c     Ligne pour n=2,Nt
c     -----------------
      do 51 n=2,Nt 
       xn=dble(n)    
       aux=aux+Pn(n)*alpha(n)-invsinEta*dUn(n)*(
     &     +sinEta*((xn+1.d0)*Pn(n+1)+xn*Pn(n-1)
     &     -(2.d0*xn+1.d0)*Pn(1))
     &     +((3.d0/2.d0)*sinEta-1.d0)
     &     *(Pn(n-1)-Pn(n+1)))    
 51   continue      

      dSigma_barre=aux+aux1+aux2+aux0

c     scaling par  (-1/c^3)(1-x)^{1/2} 
c     ---------------------------------
      dsigma0zx=(1.d0/(c*c*c))
     &           *(1.d0-xi)**(1.d0/2.d0)
     &           *dSigma_barre

c     =============================
c     Fin de la fonction sigmaBerdan
c     ==============================
      end function dsigma0zx


