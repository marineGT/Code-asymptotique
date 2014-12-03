c****************************************************************
c                fonction qui calcule la composante de
c                    la vitesse uzeta sur la 
c                        surface de la Bulle
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Octobre 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function u0etaB(Nt,c,zetap,xi,Pn,Un,dUn)
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

c     declaration des variables locales
c     ---------------------------------
      integer :: m, k    
      double precision :: aux, aux0, aux1, aux2 
      double precision :: xm, sinEta, invsinEta, ch, sh
      double precision :: uEtaB
     
c     Declaration des variables en sorties
c     ------------------------------------      
      double precision :: u0etaB

c     ==================================================================      

c     Calcul de u0etaB
c     -----------------
c     on calcule la composante de la vitesse u0etaB avec zeta = -zetap
c
c      u0etaB = ((ch-x)^(-1/2)/(c²*sin(eta))) sum{3/2*sh*Un*Vn-(ch-x)*U'n*Vn}

c     Soit suivant les Pn

c      u0etaB = (-(ch-x)^(-1/2)/(c²*sin(eta)))
c                  sum{ (Pn-1 - Pn+1)*(-3/2*sh*Un + (ch-x)*U'n) }

c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=dsinh(zetap)  ! pas besoin de -zetap, la valeurs de sinh 
                       ! correspond a celle de Mathematica

c     Calcul de sin(eta) et 1/sin(eta)
c     --------------------------------   
      sinEta=dsqrt((1.d0-xi)*(1.d0+xi)) 
      invsinEta=1.d0/sinEta

c     On initialise uEtaB à O
c     ----------------------------
      uEtaB=0.d0
 
c     Calcul de la somme uEtaB
c     -----------------------------
c     La somme de uEtaB va de n=1,....,Nt+1
      aux=0.d0
      aux1=0.d0
      aux2=0.d0
c     on sait que Pn0=1.D0 
c     quand n= 0 U'0(-zetap) = 0 et U0(-zetap) = 0
c      u0etaB(n=0) = 0
c     quand n= 1,  
c      u0etaB(1) =  (-(ch-x)^(1/2)/(c²*sin(eta)))
c                *[(1-P(2))*(-3/2*sh*Un + (ch-x)*U'n)]

c     Ligne pour n=1
c     --------------
      aux1=(1.d0-Pn(2))*((-3.d0/2.d0)*sh*Un(1)
     &      +(ch-xi)*dUn(1))

c     Ligne pour n=2,Nt+1
c     -------------------
c      m=3
      do 21 m=2,Nt+1
       aux=aux+(Pn(m-1)-Pn(m+1))
     &     *(-(3.d0/2.d0)*sh*Un(m)+(ch-xi)*dUn(m))
 21   continue          
      
      uEtaB=aux1+aux

c     scaling par -(ch-x)^(1/2)/(c²*sin(eta))
c     --------------------------------------
      u0etaB=(-invsinEta/(c*c*dsqrt(ch-xi)))
     &       *uEtaB

c         write(8,*) xi, u0etaB
                  

c       write(*,*) 'c dans ueta avant sortie', c
c     =============================
c     Fin de la fonction u0etaSL
c     ==============================
      end function u0etaB


